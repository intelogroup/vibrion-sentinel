"""Phase 2 multi-tenancy API: roles on keys, key management, envelope
encryption, deletion workflows, audit log.

Roles (on API keys; the service is machine-to-machine):

    admin  — everything
    member — upload + read + delete own org's jobs (no key management, no DEK rotation)
    viewer — read-only

Cross-org: 404, never 403 (no existence oracle). Within an org, role
violations are 403.

Endpoints (all under per-org bearer auth unless noted):

    POST   /uploads            tus creation                    (admin, member)
    HEAD   /uploads/{id}       resume offsets                  (admin, member)
    PATCH  /uploads/{id}       chunk append                    (admin, member)
    DELETE /uploads/{id}       cancel a partial upload         (admin, member)
    GET    /jobs               calling org's jobs              (all roles)
    GET    /jobs/{id}          job detail incl. report         (all roles)
    POST   /jobs/{id}/retry    re-queue a failed job           (admin, member)
    DELETE /jobs/{id}          delete job + its bytes          (admin, member)
    POST   /org/keys           create key (plaintext once)     (admin)
    GET    /org/keys           list keys                       (admin)
    DELETE /org/keys/{key_id}  revoke key                      (admin)
    POST   /org/keys/rotate-dek  new DEK version               (admin)
    DELETE /org/data           full org purge (confirm=org_id)  (admin)
    GET    /org/audit          audit log, paginated            (admin)
    GET    /health             liveness (no auth)

Use create_app(db, storage, settings, key_provider=None) so tests can
inject fakes. In production the storage MUST be an EncryptedStorage —
build_app() wires that and fails closed without SENTINEL_KEK.
"""

import os
import secrets
from typing import Annotated, Any, Optional

from fastapi import Depends, FastAPI, Header, HTTPException, Request, Response
from fastapi.responses import HTMLResponse, RedirectResponse
from jinja2 import Environment, FileSystemLoader, select_autoescape

from .auth import VALID_ROLES, Actor, hash_key, make_org_dependency
from .config import Settings
from .crypto import EncryptedStorage, KeyProvider
from .db import JobDB, SupabaseError
from .human_auth import SESSION_COOKIE, create_human_auth_router, make_user_dependency
from .reports import render_html, render_pdf
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
    key_provider: Optional[KeyProvider] = None,
) -> FastAPI:
    settings = settings or Settings()
    app = FastAPI(title="Vibrion Sentinel ingest service", version="3.0.0")
    require_actor = make_org_dependency(db)
    ActorDep = Annotated[Actor, Depends(require_actor)]
    get_portal_user = make_user_dependency(db)

    def _forbidden(detail: str = "insufficient role for this action") -> HTTPException:
        # Within an org, role violations are 403. (Between orgs the 404
        # oracle rule still applies: _get_job_or_404 etc. filter by org.)
        return HTTPException(status_code=403, detail=detail)

    def _require(actor: Actor, *roles: str) -> None:
        if actor.role not in roles:
            raise _forbidden(
                f"role {actor.role!r} cannot perform this action "
                f"(requires one of {roles})"
            )

    def _audit(actor: Actor, action: str, target: str, detail: Optional[dict] = None) -> None:
        try:
            db.audit(actor.org_id, actor.key_id, action, target, detail or {})
        except SupabaseError:
            # Audit failure must not mask the primary action's result, but
            # it must not pass silently either: surface it in the response?
            # Decision: log-and-continue would hide a broken audit trail, so
            # we let it raise as a 502 — the action already happened, and
            # the operator sees the audit gap explicitly.
            raise HTTPException(status_code=502, detail="audit log write failed")

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
        actor: ActorDep,
        request: Request,
        upload_length: Annotated[Optional[str], Header(alias="Upload-Length")] = None,
        upload_metadata: Annotated[Optional[str], Header(alias="Upload-Metadata")] = None,
    ) -> Response:
        _require(actor, "admin", "member")
        org_id = actor.org_id
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
    def upload_offset(actor: ActorDep, job_id: str) -> Response:
        _require(actor, "admin", "member")
        job = _get_job_or_404(actor.org_id, job_id)
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
        actor: ActorDep,
        job_id: str,
        request: Request,
        upload_offset: Annotated[Optional[str], Header(alias="Upload-Offset")] = None,
    ) -> Response:
        _require(actor, "admin", "member")
        org_id = actor.org_id
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
    def cancel_upload(actor: ActorDep, job_id: str) -> Response:
        _require(actor, "admin", "member")
        job = _get_job_or_404(actor.org_id, job_id)
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
    def list_jobs(actor: ActorDep) -> dict[str, Any]:
        try:
            jobs = db.list_jobs(actor.org_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        return {"jobs": [_job_summary(j) for j in jobs]}

    @app.get("/jobs/{job_id}")
    def get_job(actor: ActorDep, job_id: str) -> dict[str, Any]:
        return _job_detail(_get_job_or_404(actor.org_id, job_id))

    @app.post("/jobs/{job_id}/retry")
    def retry_job(actor: ActorDep, job_id: str) -> dict[str, Any]:
        _require(actor, "admin", "member")
        job = _get_job_or_404(actor.org_id, job_id)
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

    @app.delete("/jobs/{job_id}", status_code=204)
    def delete_job(actor: ActorDep, job_id: str) -> Response:
        """Delete a job: its bytes (object + sidecar), its job row.

        The sentinel_runs mirror row is KEPT and tombstoned with deleted_at:
        population-level aggregates must not silently rewrite history when a
        source sample is deleted. See README Phase 2.
        """
        _require(actor, "admin", "member")
        job = _get_job_or_404(actor.org_id, job_id)
        if job.get("r2_key"):
            storage.delete(job["r2_key"])  # EncryptedStorage also drops the sidecar
        try:
            db.delete_job(job_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        try:
            db.mark_sentinel_run_deleted(job_id)
        except SupabaseError:
            pass  # no mirror row (job never ran) — not an error
        _audit(actor, "job.deleted", job_id, {"filename": job.get("filename")})
        return Response(status_code=204)

    # ------------------------------------------------------------------ org
    def _key_summary(row: dict[str, Any]) -> dict[str, Any]:
        return {
            "id": row["id"],
            "name": row.get("name"),
            "role": row.get("role"),
            "key_hash_prefix": (row.get("key_hash") or "")[:10],
            "created_at": row.get("created_at"),
            "revoked_at": row.get("revoked_at"),
        }

    @app.post("/org/keys", status_code=201)
    def create_api_key(actor: ActorDep, body: dict) -> dict[str, Any]:
        _require(actor, "admin")
        name = (body or {}).get("name") or ""
        role = (body or {}).get("role") or "member"
        if role not in VALID_ROLES:
            raise HTTPException(
                status_code=400,
                detail=f"unknown role {role!r}; expected one of {VALID_ROLES}",
            )
        plaintext = "sk_" + secrets.token_urlsafe(32)
        try:
            row = db.insert_api_key(
                actor.org_id, hash_key(plaintext), name, role
            )
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        _audit(actor, "key.created", row["id"], {"name": name, "role": role})
        # The plaintext is returned ONCE and never stored.
        return {"id": row["id"], "name": name, "role": role, "api_key": plaintext}

    @app.get("/org/keys")
    def list_api_keys(actor: ActorDep) -> dict[str, Any]:
        _require(actor, "admin")
        try:
            rows = db.list_api_keys(actor.org_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        return {"keys": [_key_summary(r) for r in rows]}

    @app.delete("/org/keys/{key_id}", status_code=204)
    def revoke_api_key(actor: ActorDep, key_id: str) -> Response:
        _require(actor, "admin")
        try:
            row = db.get_api_key(actor.org_id, key_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if row is None:
            # 404 for other orgs' keys too: no existence oracle.
            raise _not_found()
        if row.get("revoked_at"):
            raise HTTPException(status_code=409, detail="key already revoked")
        try:
            db.revoke_api_key(actor.org_id, key_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        _audit(actor, "key.revoked", key_id, {"name": row.get("name")})
        return Response(status_code=204)

    @app.post("/org/keys/rotate-dek")
    def rotate_dek(actor: ActorDep) -> dict[str, Any]:
        _require(actor, "admin")
        if key_provider is None:
            raise HTTPException(
                status_code=503, detail="encryption not configured on this deployment"
            )
        try:
            new_version = key_provider.rotate(actor.org_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        _audit(actor, "dek.rotated", actor.org_id, {"dek_version": new_version})
        return {"org_id": actor.org_id, "dek_version": new_version}

    @app.delete("/org/data")
    async def purge_org_data(actor: ActorDep, request: Request) -> dict[str, Any]:
        """Full org purge: the "we're leaving" button.

        Requires {"confirm": "<org_id>"} in the body (else 400). Deletes all
        R2 objects under the org prefix, all job rows, and revokes all keys.
        DEK versions are retained (decrypting the audit trail / backups is a
        separate, deliberate process). Loud in the audit log.
        """
        _require(actor, "admin")
        try:
            body = await request.json()
        except Exception:
            body = {}
        if not isinstance(body, dict) or body.get("confirm") != actor.org_id:
            raise HTTPException(
                status_code=400,
                detail=f'body must be {{"confirm": "{actor.org_id}"}}',
            )
        prefix = f"intake/{actor.org_id}/"
        try:
            keys = storage.list_prefix(prefix)
            for k in keys:
                storage.delete(k)
            n_jobs = db.delete_jobs_for_org(actor.org_id)
            n_keys = db.revoke_all_keys(actor.org_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        _audit(
            actor,
            "org.purged",
            actor.org_id,
            {
                "objects_deleted": len(keys),
                "jobs_deleted": n_jobs,
                "keys_revoked": n_keys,
            },
        )
        return {
            "org_id": actor.org_id,
            "objects_deleted": len(keys),
            "jobs_deleted": n_jobs,
            "keys_revoked": n_keys,
        }

    @app.get("/org/audit")
    def get_audit(
        actor: ActorDep, limit: int = 50, offset: int = 0
    ) -> dict[str, Any]:
        _require(actor, "admin")
        limit = min(max(int(limit), 1), 200)
        offset = max(int(offset), 0)
        try:
            rows = db.list_audit(actor.org_id, limit, offset)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        return {"audit": rows, "limit": limit, "offset": offset}

    # ------------------------------------------------------- Phase 3: reports
    def _resolve_report_access(
        request: Request, job_id: str
    ) -> tuple[dict[str, Any], str]:
        """(job, role) for report viewing: API key OR portal session.

        API keys keep working unchanged. Humans opening a dashboard link
        authenticate with the session cookie; they must be members of the
        job's org. Cross-org: 404, never an oracle."""
        auth = request.headers.get("authorization", "")
        if auth.lower().startswith("bearer "):
            try:
                row = db.lookup_key(hash_key(auth[7:].strip()))
            except SupabaseError as e:
                raise HTTPException(status_code=502, detail=f"database error: {e}")
            if row is None or row.get("revoked_at"):
                raise HTTPException(status_code=401, detail="invalid API key")
            job = _get_job_or_404(row["org_id"], job_id)
            return job, row.get("role") or "member"
        user = get_portal_user(request)
        if user is None:
            raise HTTPException(status_code=401, detail="login required")
        try:
            job = db.get_job_any(job_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if job is None:
            raise _not_found()
        try:
            mem = db.get_org_member(user.id, job["org_id"])
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if mem is None:
            raise _not_found()
        return job, mem.get("role") or "member"

    def _report_or_404(job: dict[str, Any]) -> dict[str, Any]:
        report = job.get("report")
        if not isinstance(report, dict):
            raise HTTPException(status_code=404, detail="report not ready yet")
        return report

    @app.get("/jobs/{job_id}/report", response_class=HTMLResponse)
    def get_report_html(request: Request, job_id: str, lang: str = "en") -> str:
        job, _role = _resolve_report_access(request, job_id)
        report = _report_or_404(job)
        return render_html(
            report, job_id, job["org_id"], lang=lang,
            uploaded_at=job.get("created_at"), completed_at=job.get("updated_at"),
        )

    @app.get("/jobs/{job_id}/report.pdf")
    def get_report_pdf(request: Request, job_id: str, lang: str = "en") -> Response:
        job, _role = _resolve_report_access(request, job_id)
        report = _report_or_404(job)
        try:
            pdf = render_pdf(
                report, job_id, job["org_id"], lang=lang,
                uploaded_at=job.get("created_at"),
                completed_at=job.get("updated_at"),
            )
        except RuntimeError as e:
            raise HTTPException(status_code=503, detail=str(e))
        return Response(
            content=pdf,
            media_type="application/pdf",
            headers={
                "Content-Disposition":
                    f'inline; filename="sentinel-report-{job_id[:8]}.pdf"'
            },
        )

    @app.get("/jobs/{job_id}/tree")
    def get_tree(request: Request, job_id: str) -> Any:
        """Phylogeny link-out (Auspice). Embedding is a follow-up; until the
        Auspice deployment exists this documents the pending integration."""
        _resolve_report_access(request, job_id)  # auth first, then integrate
        base = (settings.auspice_base_url or "").strip().rstrip("/")
        if base:
            return RedirectResponse(f"{base}/tree?job={job_id}", status_code=302)
        return {
            "tree": None,
            "followup": "Auspice integration pending — see README Phase 3; "
                        "set AUSPICE_BASE_URL to enable the link-out.",
        }

    # ------------------------------------------------------- Phase 3: portal
    app.include_router(
        create_human_auth_router(db, cookie_secure=settings.session_cookie_secure)
    )

    _dash_env = Environment(
        loader=FileSystemLoader(
            str(__import__("pathlib").Path(__file__).parent / "templates")
        ),
        autoescape=select_autoescape(default_for_string=True, default=False),
    )

    @app.get("/login", response_class=HTMLResponse)
    def login_page() -> str:
        return _dash_env.get_template("login.html").render()

    @app.get("/signup", response_class=HTMLResponse)
    def signup_page() -> str:
        return _dash_env.get_template("signup.html").render()

    @app.get("/dashboard", response_class=HTMLResponse)
    def dashboard(request: Request) -> str:
        user = get_portal_user(request)
        if user is None:
            return RedirectResponse("/login", status_code=302)
        try:
            memberships = db.get_org_memberships(user.id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        orgs = []
        for m in memberships:
            try:
                jobs = db.list_jobs(m["org_id"], limit=100)
            except SupabaseError:
                jobs = []
            for j in jobs:
                rep = j.get("report") or {}
                qc = rep.get("qc") or {}
                j["qc_status"] = qc.get("status")
            orgs.append({"org_id": m["org_id"], "role": m["role"], "jobs": jobs})
        return _dash_env.get_template("dashboard.html").render(
            user_email=user.email, orgs=orgs
        )

    @app.post("/dashboard/jobs/{job_id}/retry")
    def dashboard_retry(request: Request, job_id: str) -> Any:
        """Dashboard retry (session auth). Read-only dashboard except this."""
        user = get_portal_user(request)
        if user is None:
            return RedirectResponse("/login", status_code=302)
        try:
            job = db.get_job_any(job_id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if job is None:
            raise _not_found()
        try:
            mem = db.get_org_member(user.id, job["org_id"])
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if mem is None:
            raise _not_found()
        if mem.get("role") not in ("admin", "member"):
            raise _forbidden()
        if job["status"] != "failed":
            raise HTTPException(status_code=409, detail="only failed jobs can be retried")
        try:
            db.transition(job, "queued", reason=None)
        except (SupabaseError, ValueError) as e:
            raise HTTPException(status_code=502, detail=str(e))
        return RedirectResponse("/dashboard", status_code=303)

    return app


def build_app() -> FastAPI:
    """Production entrypoint: real Supabase + R2 from env.

    Fail-closed: SENTINEL_KEK must be set (the storage layer encrypts every
    upload with the org's DEK; without a KEK the service refuses to start
    rather than writing plaintext).
    """
    from .storage import R2Storage

    settings = Settings()
    settings.require_db()
    settings.require_r2()
    kek = settings.require_kek()  # raises RuntimeError if missing/malformed
    db = SupabaseRestDB(settings.supabase_url, settings.supabase_service_role_key)
    keys = KeyProvider(db, kek)
    storage = EncryptedStorage(
        R2Storage(
            settings.r2_endpoint_url,
            settings.r2_access_key_id,
            settings.r2_secret_access_key,
            settings.r2_bucket,
        ),
        keys,
    )
    return create_app(db, storage, settings, keys)


def asgi_app() -> FastAPI:
    """Uvicorn factory: `uvicorn --factory service.app.main:asgi_app`."""
    return build_app()
