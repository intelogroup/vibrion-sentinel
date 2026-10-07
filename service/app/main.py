"""Phase 1 upload-ingestion API: tus 1.0.0 creation core + job status.

Endpoints (all under per-org bearer auth; other-org resources read as 404):

    POST   /uploads            tus creation (Upload-Length + Upload-Metadata)
    HEAD   /uploads/{id}       resume offsets (Upload-Offset / Upload-Length)
    PATCH  /uploads/{id}       chunk append (409 on offset mismatch)
    DELETE /uploads/{id}       cancel a partial upload
    GET    /jobs               calling org's jobs, newest first
    GET    /jobs/{id}          job detail incl. report.json when done
    POST   /jobs/{id}/retry    re-queue a failed job (explicit action only)
    GET    /health             liveness (no auth)

Use create_app(db, storage, settings) so tests can inject fakes.
"""

import os
from typing import Annotated, Any, Optional

from fastapi import Depends, FastAPI, Header, HTTPException, Request, Response

from .auth import make_org_dependency
from .config import Settings
from .db import JobDB, SupabaseError
from .storage import Storage
from .tus import parse_int_header, parse_metadata
from .validation import validate_upload_params

TUS_VERSION = "1.0.0"
FASTQ_MAGIC = b"@"  # first byte of a FASTQ header; typo catcher, not a validator


def _job_summary(job: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": job["id"],
        "org_id": job["org_id"],
        "filename": job["filename"],
        "size_bytes": job.get("size_bytes"),
        "sha256": job.get("sha256"),
        "upload_offset": job.get("upload_offset", 0),
        "platform": job["platform"],
        "basecaller_model": job.get("basecaller_model"),
        "tier": job.get("tier"),
        "status": job["status"],
        "status_reason": job.get("status_reason"),
        "created_at": job.get("created_at"),
        "updated_at": job.get("updated_at"),
    }


def _job_detail(job: dict[str, Any]) -> dict[str, Any]:
    d = _job_summary(job)
    d["report"] = job.get("report")
    return d


def create_app(
    db: JobDB,
    storage: Storage,
    settings: Optional[Settings] = None,
) -> FastAPI:
    settings = settings or Settings()
    app = FastAPI(title="Vibrion Sentinel ingest service", version="1.0.0")
    require_org = make_org_dependency(db)
    Org = Annotated[str, Depends(require_org)]

    def _not_found() -> HTTPException:
        # 404 for missing AND other-org resources: no existence oracle.
        return HTTPException(status_code=404, detail="not found")

    def _get_job_or_404(org_id: str, job_id: str) -> dict[str, Any]:
        try:
            job = db.get_job(org_id, job_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if job is None:
            raise _not_found()
        return job

    @app.get("/health")
    def health() -> dict[str, str]:
        return {"status": "ok"}

    # ------------------------------------------------------------------ tus
    @app.post("/uploads", status_code=201)
    def create_upload(
        org_id: Org,
        request: Request,
        upload_length: Annotated[Optional[str], Header(alias="Upload-Length")] = None,
        upload_metadata: Annotated[Optional[str], Header(alias="Upload-Metadata")] = None,
    ) -> Response:
        try:
            length = parse_int_header(upload_length, "Upload-Length")
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
        if length > settings.max_upload_bytes:
            raise HTTPException(
                status_code=413,
                detail=f"Upload-Length {length} exceeds limit {settings.max_upload_bytes}",
            )
        try:
            meta = parse_metadata(upload_metadata)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
        try:
            platform, basecaller_model, tier = validate_upload_params(meta)
        except ValueError as e:
            # Fail closed: nanopore without basecaller_model, bad enums, etc.
            raise HTTPException(status_code=400, detail=str(e))
        raw_filename = meta.get("filename") or "upload.fastq.gz"
        filename = os.path.basename(raw_filename) or "upload.fastq.gz"

        try:
            job = db.create_job(
                {
                    "org_id": org_id,
                    "filename": filename,
                    "size_bytes": length,
                    "platform": platform,
                    "basecaller_model": basecaller_model,
                    "tier": tier,
                    "client_checksum": meta.get("checksum"),
                    "status": "uploading",
                    # r2_key is set after we know the id (second update).
                }
            )
            r2_key = f"intake/{org_id}/{job['id']}/{filename}"
            job = db.update_job(job["id"], {"r2_key": r2_key})
            storage.create(r2_key)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")

        return Response(
            status_code=201,
            headers={
                "Location": f"/uploads/{job['id']}",
                "Tus-Resumption": TUS_VERSION,
                "Upload-Offset": "0",
            },
        )

    @app.head("/uploads/{job_id}")
    def upload_offset(org_id: Org, job_id: str) -> Response:
        job = _get_job_or_404(org_id, job_id)
        return Response(
            status_code=200,
            headers={
                "Upload-Offset": str(job.get("upload_offset", 0)),
                "Upload-Length": str(job.get("size_bytes", 0)),
                "Cache-Control": "no-store",
                "Tus-Resumption": TUS_VERSION,
            },
        )

    @app.patch("/uploads/{job_id}", status_code=204)
    async def append_chunk(
        org_id: Org,
        job_id: str,
        request: Request,
        upload_offset: Annotated[Optional[str], Header(alias="Upload-Offset")] = None,
    ) -> Response:
        job = _get_job_or_404(org_id, job_id)
        if job["status"] != "uploading":
            raise HTTPException(
                status_code=409,
                detail=f"upload is {job['status']}, not accepting chunks",
            )
        try:
            offset = parse_int_header(upload_offset, "Upload-Offset")
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
        server_offset = int(job.get("upload_offset", 0))
        if offset != server_offset:
            raise HTTPException(
                status_code=409,
                detail=f"offset mismatch: client {offset}, server {server_offset}",
            )
        chunk = await request.body()
        length = int(job.get("size_bytes", 0) or 0)
        if server_offset + len(chunk) > length:
            raise HTTPException(
                status_code=400,
                detail=f"chunk overruns declared Upload-Length {length}",
            )
        if server_offset == 0 and chunk and not chunk.startswith(FASTQ_MAGIC):
            # Typo catcher: first byte of a FASTQ record is '@'.
            db.transition(job, "failed", reason="first bytes do not look like FASTQ")
            raise HTTPException(status_code=400, detail="upload does not look like FASTQ")

        storage.append(job["r2_key"], chunk)
        new_offset = server_offset + len(chunk)

        if new_offset == length:
            # Completion: verify, then hand to the worker.
            storage.finalize(job["r2_key"])
            digest = storage.sha256(job["r2_key"])
            job = db.update_job(job_id, {"upload_offset": new_offset, "sha256": digest})
            declared = (job.get("client_checksum") or "").strip().lower()
            if declared and declared != digest:
                db.transition(job, "failed", reason="sha256 mismatch vs client-declared checksum")
                raise HTTPException(status_code=400, detail="checksum mismatch")
            job = db.transition(job, "queued")
        else:
            job = db.update_job(job_id, {"upload_offset": new_offset})

        return Response(
            status_code=204,
            headers={"Upload-Offset": str(job.get("upload_offset", new_offset))},
        )

    @app.delete("/uploads/{job_id}", status_code=204)
    def cancel_upload(org_id: Org, job_id: str) -> Response:
        job = _get_job_or_404(org_id, job_id)
        if job["status"] not in ("uploading", "queued"):
            raise HTTPException(
                status_code=409,
                detail=f"cannot cancel a {job['status']} upload",
            )
        storage.abort(job["r2_key"])
        db.transition(job, "cancelled", reason="cancelled by client")
        return Response(status_code=204)

    # ------------------------------------------------------------------ jobs
    @app.get("/jobs")
    def list_jobs(org_id: Org) -> dict[str, Any]:
        try:
            jobs = db.list_jobs(org_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        return {"jobs": [_job_summary(j) for j in jobs]}

    @app.get("/jobs/{job_id}")
    def get_job(org_id: Org, job_id: str) -> dict[str, Any]:
        return _job_detail(_get_job_or_404(org_id, job_id))

    @app.post("/jobs/{job_id}/retry")
    def retry_job(org_id: Org, job_id: str) -> dict[str, Any]:
        job = _get_job_or_404(org_id, job_id)
        if job["status"] != "failed":
            raise HTTPException(
                status_code=409,
                detail=f"only failed jobs can be retried (this one is {job['status']})",
            )
        try:
            job = db.transition(job, "queued", reason=None)
        except (SupabaseError, ValueError) as e:
            raise HTTPException(status_code=502, detail=str(e))
        return _job_summary(job)

    return app


def build_app() -> FastAPI:
    """Production entrypoint: real Supabase + R2 from env."""
    from .storage import R2Storage

    settings = Settings()
    settings.require_db()
    settings.require_r2()
    db = SupabaseRestDB(settings.supabase_url, settings.supabase_service_role_key)
    storage = R2Storage(
        settings.r2_endpoint_url,
        settings.r2_access_key_id,
        settings.r2_secret_access_key,
        settings.r2_bucket,
    )
    return create_app(db, storage, settings)


def asgi_app() -> FastAPI:
    """Uvicorn factory: `uvicorn --factory service.app.main:asgi_app`."""
    return build_app()
