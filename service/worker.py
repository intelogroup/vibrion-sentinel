"""Phase 1 worker: queued ingest jobs in, pipeline reports out.

Long-running process. Each iteration claims one `queued` job atomically
(POST /rest/v1/rpc/claim_next_ingest_job -> FOR UPDATE SKIP LOCKED, so
multiple workers are safe), downloads the FASTQ from R2, writes a
pipeline config.yaml using the EXISTING schema
(workflow/sentinel_lite/Snakefile), runs snakemake, and stores the
report.json verbatim on the job row + mirrored into sentinel_runs.

Failed jobs stay failed for inspection; re-queue is the explicit
POST /jobs/{id}/retry endpoint. Never retry blindly.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
import traceback
from pathlib import Path
from typing import Any, Callable

import yaml

from service.app.config import Settings
from service.app.crypto import DecryptionError, EncryptedStorage, KeyProvider
from service.app.db import JobDB, SupabaseRestDB
from service.app.storage import R2Storage, Storage

# Tail of the snakemake log kept in status_reason on failure.
LOG_TAIL_CHARS = 4000


def build_config(job: dict[str, Any], settings: Settings, local_fastq: str, outdir: str) -> dict[str, Any]:
    """Pipeline config.yaml reusing the EXISTING sentinel_lite schema.

    accession = the ingest job id, so reports join back to jobs and the
    sentinel_runs mirror trivially.
    """
    cfg: dict[str, Any] = {
        "accession": job["id"],
        "refdir": settings.pipeline_refdir,
        "outdir": outdir,
        "reads_r1": local_fastq,
        "paired": False,  # Phase 1 uploads are single-file (nanopore SE; illumina SE)
        "threads": settings.snakemake_cores,
        "platform": job["platform"],
        "tier": job.get("tier", "lite"),
    }
    if job.get("basecaller_model"):
        cfg["basecaller_model"] = job["basecaller_model"]
    return cfg


def build_snakemake_cmd(config_path: str, snakefile: str, cores: int) -> list[str]:
    return [
        "snakemake",
        "-s", snakefile,
        "--configfile", config_path,
        "--cores", str(cores),
    ]


def default_run_pipeline(workdir: str, cmd: list[str]) -> tuple[int, str]:
    """Run snakemake, capturing combined output. Returns (returncode, log)."""
    proc = subprocess.run(
        cmd,
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=6 * 3600,  # 6h wall clock; a stuck job should fail, not hang forever
    )
    log = f"$ {' '.join(cmd)}\n{proc.stdout}\n{proc.stderr}"
    return proc.returncode, log


def mirror_row(job: dict[str, Any], report: dict[str, Any]) -> dict[str, Any]:
    """Map a pipeline report.json onto the sentinel_runs schema."""
    reads = report.get("reads", {})
    classification = report.get("classification", {})
    mapping = report.get("mapping", {})
    variants = report.get("variants", {})
    qc = report.get("qc", {})
    return {
        "accession": job["id"],
        "priority": "upload",
        "pipeline_version": report.get("version"),
        "reference": report.get("reference"),
        "reads_raw": reads.get("raw"),
        "reads_qc_passed": reads.get("qc_passed"),
        "reads_vibrio_kept": reads.get("vibrio_kept"),
        "kraken_db": classification.get("kraken_db"),
        "v_cholerae_reads": classification.get("v_cholerae_reads"),
        "v_cholerae_fraction": classification.get("v_cholerae_fraction"),
        "mapped_reads": mapping.get("mapped_reads"),
        "mean_depth": mapping.get("mean_depth"),
        "breadth_pct": mapping.get("breadth_pct"),
        "snps_vs_7pet": variants.get("snps_vs_7pet"),
        "consensus_length": report.get("consensus_length"),
        "consensus_called_pct": report.get("consensus_called_pct"),
        "qc_status": qc.get("status"),
        "qc_reasons": qc.get("reasons"),
        "surveillance_loci": report.get("surveillance_loci"),
        "report": report,
    }


RunPipeline = Callable[[str, list[str]], tuple[int, str]]


def process_one_job(
    db: JobDB,
    storage: Storage,
    settings: Settings,
    run_pipeline: RunPipeline | None = None,
) -> bool:
    """Claim and process a single queued job. Returns True if one was processed."""
    run_pipeline = run_pipeline or default_run_pipeline
    job = db.claim_next_job()
    if job is None:
        return False

    job_id = job["id"]
    workdir = Path(settings.worker_workdir) / job_id
    try:
        workdir.mkdir(parents=True, exist_ok=True)
        local_fastq = str(workdir / job["filename"])
        try:
            # With EncryptedStorage this authenticates + decrypts; a missing
            # sidecar or GCM tag failure raises DecryptionError (fail closed).
            storage.download_to_file(job["r2_key"], local_fastq)
        except DecryptionError as e:
            db.transition(job, "failed", reason=f"decryption failed: {e}")
            return True

        # Optional integrity re-check at download time.
        want = (job.get("sha256") or "").strip().lower()
        if want:
            import hashlib

            h = hashlib.sha256()
            with open(local_fastq, "rb") as f:
                for chunk in iter(lambda: f.read(1 << 20), b""):
                    h.update(chunk)
            if h.hexdigest() != want:
                db.transition(job, "failed", reason="R2 download sha256 mismatch vs stored digest")
                return True

        outdir = str(workdir / "out")
        cfg = build_config(job, settings, local_fastq, outdir)
        cfg_path = str(workdir / "config.yaml")
        with open(cfg_path, "w") as f:
            yaml.safe_dump(cfg, f)

        cmd = build_snakemake_cmd(cfg_path, settings.snakefile, settings.snakemake_cores)
        rc, log = run_pipeline(str(workdir), cmd)
        if rc != 0:
            db.transition(job, "failed", reason=log[-LOG_TAIL_CHARS:])
            return True

        report_path = os.path.join(outdir, "report.json")
        with open(report_path) as f:
            report = json.load(f)
        db.update_job(job_id, {"report": report})
        db.transition({**job, "status": "running"}, "done")
        db.mirror_sentinel_run(mirror_row(job, report))
        return True
    except Exception:
        # Never let a worker crash lose the job's state: record and move on.
        reason = traceback.format_exc(limit=5)[-LOG_TAIL_CHARS:]
        try:
            db.transition({**job, "status": "running"}, "failed", reason=reason)
        except Exception:
            pass
        return True
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main() -> None:
    settings = Settings()
    settings.require_db()
    settings.require_r2()
    kek = settings.require_kek()  # fail closed: no KEK, no worker
    db = SupabaseRestDB(settings.supabase_url, settings.supabase_service_role_key)
    storage = EncryptedStorage(
        R2Storage(
            settings.r2_endpoint_url,
            settings.r2_access_key_id,
            settings.r2_secret_access_key,
            settings.r2_bucket,
        ),
        KeyProvider(db, kek),
    )
    os.makedirs(settings.worker_workdir, exist_ok=True)
    print(f"sentinel ingest worker: polling every {settings.worker_poll_interval}s", flush=True)
    while True:
        try:
            if not process_one_job(db, storage, settings):
                time.sleep(settings.worker_poll_interval)
        except Exception:
            traceback.print_exc()
            time.sleep(settings.worker_poll_interval)


if __name__ == "__main__":
    main()
