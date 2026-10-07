"""Worker tests: state machine against fakes, stubbed pipeline runner.

Covers queued->running->done with report storage + sentinel_runs mirror,
failure capture (log tail in status_reason), download integrity, and the
snakemake invocation path shape (dry — no toolchain needed).

Run: python -m unittest discover -s service/tests -t .
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from service.app.config import Settings  # noqa: E402
from service.app.storage import MemoryStorage  # noqa: E402
from service.tests.fakes import FakeDB  # noqa: E402
from service import worker  # noqa: E402


FAKE_REPORT = {
    "pipeline": "vibrion-sentinel-lite",
    "version": "0.2.0",
    "accession": "PLACEHOLDER",
    "reference": "2010EL-1786 (7PET, CP003069.1/CP003070.1)",
    "reads": {"raw": 1000, "qc_passed": 990, "vibrio_kept": 980},
    "classification": {
        "kraken_db": "kraken2_haiti_custom",
        "v_cholerae_reads": 970,
        "v_cholerae_fraction": 0.98,
    },
    "mapping": {"mapped_reads": 960, "mean_depth": 42.5, "breadth_pct": 99.1},
    "variants": {"snps_vs_7pet": 12},
    "surveillance_loci": {"ctxA": {"call": "present"}},
    "consensus_length": 4032000,
    "consensus_called_pct": 97.5,
    "qc": {"status": "pass", "reasons": []},
}


def make_payload() -> bytes:
    return b"@read1\nACGT\n+\nIIII\n" * 100


class WorkerTestCase(unittest.TestCase):
    def setUp(self):
        self.db = FakeDB()
        self.storage = MemoryStorage()
        self.tmp = tempfile.mkdtemp(prefix="sentinel-worker-test-")
        self.settings = Settings()
        self.settings.worker_workdir = os.path.join(self.tmp, "work")
        self.settings.snakemake_cores = 2
        self.settings.repo_dir = str(REPO_ROOT)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _seed_queued_job(self, **over) -> dict:
        payload = make_payload()
        job = self.db.create_job(
            {
                "org_id": "org-a",
                "filename": "run.fastq.gz",
                "size_bytes": len(payload),
                "platform": "nanopore",
                "basecaller_model": "hac",
                "tier": "lite",
                "status": "queued",
                "r2_key": "intake/org-a/JOB/run.fastq.gz",
                "sha256": hashlib.sha256(payload).hexdigest(),
                **over,
            }
        )
        # fix the r2_key to use the real job id
        job = self.db.update_job(job["id"], {"r2_key": f"intake/org-a/{job['id']}/run.fastq.gz"})
        self.storage.create(job["r2_key"])
        self.storage.append(job["r2_key"], payload)
        self.storage.finalize(job["r2_key"])
        return job

    def _ok_runner(self, report_extra=None):
        def run(workdir: str, cmd: list[str]) -> tuple[int, str]:
            # The stub plays the pipeline: writes out/report.json like the
            # real Snakefile's `report` rule would.
            cfg_path = [c for c in cmd if c.endswith("config.yaml")]
            self.assertTrue(cfg_path, f"config.yaml missing from cmd: {cmd}")
            outdir = os.path.join(workdir, "out")
            os.makedirs(outdir, exist_ok=True)
            report = dict(FAKE_REPORT)
            report.update(report_extra or {})
            with open(os.path.join(outdir, "report.json"), "w") as f:
                json.dump(report, f)
            return 0, "snakemake ok"

        return run

    # -- happy path ------------------------------------------------------
    def test_queued_to_done_stores_report_and_mirrors(self):
        job = self._seed_queued_job()
        self.assertTrue(worker.process_one_job(self.db, self.storage, self.settings,
                                               run_pipeline=self._ok_runner()))
        done = self.db.jobs[job["id"]]
        self.assertEqual(done["status"], "done")
        self.assertEqual(done["report"]["version"], "0.2.0")
        self.assertEqual(done["report"]["accession"], "PLACEHOLDER")

        self.assertEqual(len(self.db.mirrored), 1)
        mir = self.db.mirrored[0]
        self.assertEqual(mir["accession"], job["id"])
        self.assertEqual(mir["priority"], "upload")
        self.assertEqual(mir["pipeline_version"], "0.2.0")
        self.assertEqual(mir["snps_vs_7pet"], 12)
        self.assertEqual(mir["qc_status"], "pass")
        self.assertEqual(mir["report"]["version"], "0.2.0")

    def test_no_queued_job_returns_false(self):
        self.assertFalse(worker.process_one_job(self.db, self.storage, self.settings,
                                                run_pipeline=self._ok_runner()))

    def test_claims_oldest_first(self):
        j1 = self._seed_queued_job()
        j2 = self._seed_queued_job()
        worker.process_one_job(self.db, self.storage, self.settings,
                               run_pipeline=self._ok_runner())
        self.assertEqual(self.db.jobs[j1["id"]]["status"], "done")
        self.assertEqual(self.db.jobs[j2["id"]]["status"], "queued")

    # -- failures -----------------------------------------------------------
    def test_pipeline_failure_captures_log_tail(self):
        job = self._seed_queued_job()
        long_log = "".join(f"line {i}\n" for i in range(500)) + "FATAL: boom"
        def bad_runner(workdir, cmd):
            return 1, long_log
        worker.process_one_job(self.db, self.storage, self.settings, run_pipeline=bad_runner)
        done = self.db.jobs[job["id"]]
        self.assertEqual(done["status"], "failed")
        self.assertIn("FATAL: boom", done["status_reason"])
        self.assertLessEqual(len(done["status_reason"]), worker.LOG_TAIL_CHARS)
        self.assertEqual(self.db.mirrored, [])

    def test_download_sha_mismatch_fails(self):
        job = self._seed_queued_job()
        # Tamper the stored bytes after the recorded digest.
        self.storage.append(job["r2_key"], b"tamper")
        worker.process_one_job(self.db, self.storage, self.settings,
                               run_pipeline=self._ok_runner())
        done = self.db.jobs[job["id"]]
        self.assertEqual(done["status"], "failed")
        self.assertIn("sha256", done["status_reason"])

    def test_runner_exception_marks_failed_not_lost(self):
        job = self._seed_queued_job()
        def boom(workdir, cmd):
            raise RuntimeError("snakemake binary missing")
        worker.process_one_job(self.db, self.storage, self.settings, run_pipeline=boom)
        done = self.db.jobs[job["id"]]
        self.assertEqual(done["status"], "failed")
        self.assertIn("snakemake binary missing", done["status_reason"])

    # -- config / invocation shape -------------------------------------------
    def test_build_config_uses_existing_schema_keys(self):
        job = self._seed_queued_job()
        cfg = worker.build_config(job, self.settings, "/tmp/x.fastq.gz", "/tmp/out")
        # Keys the Snakefile reads: accession, refdir, outdir, reads_r1,
        # paired, threads, platform, basecaller_model, tier.
        for key in ("accession", "refdir", "outdir", "reads_r1", "paired",
                    "threads", "platform", "basecaller_model", "tier"):
            self.assertIn(key, cfg, f"missing schema key {key}")
        self.assertEqual(cfg["accession"], job["id"])
        self.assertFalse(cfg["paired"])
        self.assertEqual(cfg["platform"], "nanopore")
        self.assertEqual(cfg["basecaller_model"], "hac")

    def test_build_snakemake_cmd_shape(self):
        """Dry proof of the invocation path: binary, Snakefile, config, cores."""
        cmd = worker.build_snakemake_cmd("/w/config.yaml", "/app/workflow/sentinel_lite/Snakefile", 4)
        self.assertEqual(cmd[0], "snakemake")
        self.assertIn("-s", cmd)
        self.assertIn("/app/workflow/sentinel_lite/Snakefile", cmd)
        self.assertIn("--configfile", cmd)
        self.assertIn("/w/config.yaml", cmd)
        self.assertIn("--cores", cmd)
        self.assertIn("4", cmd)

    def test_forward_only_guard_rejects_rewind(self):
        job = self._seed_queued_job()
        worker.process_one_job(self.db, self.storage, self.settings,
                               run_pipeline=self._ok_runner())
        done = self.db.jobs[job["id"]]
        with self.assertRaises(ValueError):
            self.db.transition(done, "queued")  # done -> queued is illegal


if __name__ == "__main__":
    unittest.main()
