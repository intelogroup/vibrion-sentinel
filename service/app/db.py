"""PostgREST client for Supabase, following the repo's existing convention.

The sentinel-watch workflow already writes to Supabase via the REST API
(`$SUPABASE_URL/rest/v1/...` with the service-role key); this client does
the same from Python so the service needs no direct Postgres connection
(which does not traverse this project's egress proxy anyway).

All methods return plain dicts parsed from JSON. Errors raise
SupabaseError with the status code attached.
"""

from __future__ import annotations

import json
from typing import Any, Optional

import httpx


class SupabaseError(RuntimeError):
    def __init__(self, message: str, status: int = 0):
        super().__init__(message)
        self.status = status


class JobDB:
    """Persistence interface for ingest jobs + API keys.

    The FastAPI app and the worker code against this interface; tests use
    service.tests.fakes.FakeDB. The forward-only status discipline is
    enforced here in code (TRANSITIONS) AND by the database trigger in
    sql/ingest_jobs.sql (defense in depth).
    """

    # Legal transitions. `failed -> queued` exists only for the explicit
    # retry endpoint; nothing rewinds a job otherwise.
    TRANSITIONS: dict[str, set[str]] = {
        "uploading": {"queued", "failed", "cancelled"},
        "queued": {"running", "failed", "cancelled"},
        "running": {"done", "failed"},
        "failed": {"queued"},
        "done": set(),
        "cancelled": set(),
    }

    def create_job(self, row: dict[str, Any]) -> dict[str, Any]:
        raise NotImplementedError

    def get_job(self, org_id: str, job_id: str) -> Optional[dict[str, Any]]:
        raise NotImplementedError

    def list_jobs(self, org_id: str, limit: int = 50) -> list[dict[str, Any]]:
        raise NotImplementedError

    def update_job(self, job_id: str, fields: dict[str, Any]) -> dict[str, Any]:
        raise NotImplementedError

    def transition(self, job: dict[str, Any], to: str, reason: Optional[str] = None) -> dict[str, Any]:
        """Move a job forward one legal step, enforcing the discipline."""
        frm = job["status"]
        if to != frm and to not in self.TRANSITIONS.get(frm, set()):
            raise ValueError(f"illegal ingest job transition {frm} -> {to}")
        fields: dict[str, Any] = {"status": to}
        if reason is not None:
            fields["status_reason"] = reason
        return self.update_job(job["id"], fields)

    def claim_next_job(self) -> Optional[dict[str, Any]]:
        """Atomically claim one queued job (sets it running)."""
        raise NotImplementedError

    def lookup_key(self, key_hash: str) -> Optional[dict[str, Any]]:
        raise NotImplementedError

    def insert_api_key(self, org_id: str, key_hash: str, name: str) -> dict[str, Any]:
        raise NotImplementedError

    def mirror_sentinel_run(self, row: dict[str, Any]) -> None:
        """Upsert a row into sentinel_runs (accession = job id)."""
        raise NotImplementedError


class SupabaseRestDB(JobDB):
    def __init__(self, base_url: str, service_role_key: str, timeout: float = 30.0):
        self._base = base_url.rstrip("/")
        self._client = httpx.Client(
            base_url=self._base,
            timeout=timeout,
            headers={
                "apikey": service_role_key,
                "Authorization": f"Bearer {service_role_key}",
                "Content-Type": "application/json",
            },
        )

    def _raise(self, resp: httpx.Response, what: str) -> None:
        try:
            detail = resp.json()
        except Exception:
            detail = resp.text[:500]
        raise SupabaseError(f"{what}: HTTP {resp.status_code}: {detail}", status=resp.status_code)

    def create_job(self, row: dict[str, Any]) -> dict[str, Any]:
        r = self._client.post(
            "/rest/v1/ingest_jobs",
            content=json.dumps(row),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code not in (200, 201):
            self._raise(r, "create_job")
        return r.json()[0]

    def get_job(self, org_id: str, job_id: str) -> Optional[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/ingest_jobs",
            params={"id": f"eq.{job_id}", "org_id": f"eq.{org_id}", "limit": 1},
        )
        if r.status_code != 200:
            self._raise(r, "get_job")
        rows = r.json()
        return rows[0] if rows else None

    def list_jobs(self, org_id: str, limit: int = 50) -> list[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/ingest_jobs",
            params={
                "org_id": f"eq.{org_id}",
                "order": "created_at.desc",
                "limit": limit,
            },
        )
        if r.status_code != 200:
            self._raise(r, "list_jobs")
        return r.json()

    def update_job(self, job_id: str, fields: dict[str, Any]) -> dict[str, Any]:
        r = self._client.patch(
            "/rest/v1/ingest_jobs",
            params={"id": f"eq.{job_id}"},
            content=json.dumps(fields),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code != 200:
            self._raise(r, "update_job")
        rows = r.json()
        if not rows:
            raise SupabaseError(f"update_job: no such job {job_id}", status=404)
        return rows[0]

    def claim_next_job(self) -> Optional[dict[str, Any]]:
        r = self._client.post("/rest/v1/rpc/claim_next_ingest_job")
        if r.status_code != 200:
            self._raise(r, "claim_next_ingest_job")
        rows = r.json()
        return rows[0] if rows else None

    def lookup_key(self, key_hash: str) -> Optional[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/org_api_keys",
            params={"key_hash": f"eq.{key_hash}", "limit": 1},
        )
        if r.status_code != 200:
            self._raise(r, "lookup_key")
        rows = r.json()
        return rows[0] if rows else None

    def insert_api_key(self, org_id: str, key_hash: str, name: str) -> dict[str, Any]:
        r = self._client.post(
            "/rest/v1/org_api_keys",
            content=json.dumps({"org_id": org_id, "key_hash": key_hash, "name": name}),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code not in (200, 201):
            self._raise(r, "insert_api_key")
        return r.json()[0]

    def mirror_sentinel_run(self, row: dict[str, Any]) -> None:
        r = self._client.post(
            "/rest/v1/sentinel_runs",
            content=json.dumps(row),
            headers={"Prefer": "resolution=merge-duplicates,return=minimal"},
            params={"on_conflict": "accession"},
        )
        if r.status_code not in (200, 201, 204):
            self._raise(r, "mirror_sentinel_run")
