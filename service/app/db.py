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

    def insert_api_key(
        self, org_id: str, key_hash: str, name: str, role: str = "member"
    ) -> dict[str, Any]:
        raise NotImplementedError

    def list_api_keys(self, org_id: str) -> list[dict[str, Any]]:
        """All keys for the org (for admin key management)."""
        raise NotImplementedError

    def get_api_key(self, org_id: str, key_id: str) -> Optional[dict[str, Any]]:
        raise NotImplementedError

    def revoke_api_key(self, org_id: str, key_id: str) -> Optional[dict[str, Any]]:
        """Soft revoke (sets revoked_at). Returns the row, or None."""
        raise NotImplementedError

    def revoke_all_keys(self, org_id: str) -> int:
        """Revoke every key for the org. Returns the count revoked."""
        raise NotImplementedError

    # -- data keys (Phase 2: envelope encryption) ---------------------------
    def get_data_key(self, org_id: str) -> Optional[dict[str, Any]]:
        """The org's CURRENT DEK row (org_data_keys)."""
        raise NotImplementedError

    def get_data_key_version(
        self, org_id: str, version: int
    ) -> Optional[dict[str, Any]]:
        """A specific DEK version (org_data_key_versions)."""
        raise NotImplementedError

    def store_data_key(
        self, org_id: str, version: int, dek_wrapped_b64: str, kek_id: str
    ) -> dict[str, Any]:
        """Persist a DEK version and make it current."""
        raise NotImplementedError

    # -- deletion (Phase 2) --------------------------------------------------
    def delete_job(self, job_id: str) -> None:
        raise NotImplementedError

    def delete_jobs_for_org(self, org_id: str) -> int:
        """Delete every ingest job row for the org. Returns the count."""
        raise NotImplementedError

    def mark_sentinel_run_deleted(self, accession: str) -> None:
        """Tombstone the sentinel_runs mirror row (kept for aggregates)."""
        raise NotImplementedError

    # -- audit log (Phase 2, append-only) ------------------------------------
    def audit(
        self,
        org_id: str,
        actor_key_id: str,
        action: str,
        target: str,
        detail: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        raise NotImplementedError

    def list_audit(
        self, org_id: str, limit: int = 50, offset: int = 0
    ) -> list[dict[str, Any]]:
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

    def insert_api_key(
        self, org_id: str, key_hash: str, name: str, role: str = "member"
    ) -> dict[str, Any]:
        r = self._client.post(
            "/rest/v1/org_api_keys",
            content=json.dumps(
                {"org_id": org_id, "key_hash": key_hash, "name": name, "role": role}
            ),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code not in (200, 201):
            self._raise(r, "insert_api_key")
        return r.json()[0]

    def list_api_keys(self, org_id: str) -> list[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/org_api_keys",
            params={"org_id": f"eq.{org_id}", "order": "created_at.desc"},
        )
        if r.status_code != 200:
            self._raise(r, "list_api_keys")
        return r.json()

    def get_api_key(self, org_id: str, key_id: str) -> Optional[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/org_api_keys",
            params={"id": f"eq.{key_id}", "org_id": f"eq.{org_id}", "limit": 1},
        )
        if r.status_code != 200:
            self._raise(r, "get_api_key")
        rows = r.json()
        return rows[0] if rows else None

    def revoke_api_key(self, org_id: str, key_id: str) -> Optional[dict[str, Any]]:
        from datetime import datetime, timezone

        r = self._client.patch(
            "/rest/v1/org_api_keys",
            params={"id": f"eq.{key_id}", "org_id": f"eq.{org_id}"},
            content=json.dumps(
                {"revoked_at": datetime.now(timezone.utc).isoformat()}
            ),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code != 200:
            self._raise(r, "revoke_api_key")
        rows = r.json()
        return rows[0] if rows else None

    def revoke_all_keys(self, org_id: str) -> int:
        from datetime import datetime, timezone

        # Count first (PostgREST bulk PATCH can return the rows instead).
        keys = self.list_api_keys(org_id)
        live = [k for k in keys if not k.get("revoked_at")]
        if not live:
            return 0
        r = self._client.patch(
            "/rest/v1/org_api_keys",
            params={"org_id": f"eq.{org_id}", "revoked_at": "is.null"},
            content=json.dumps(
                {"revoked_at": datetime.now(timezone.utc).isoformat()}
            ),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code != 200:
            self._raise(r, "revoke_all_keys")
        return len(r.json())

    # -- data keys -----------------------------------------------------------
    # NOTE on bytea over PostgREST: bytea columns are sent as base64 in JSON
    # payloads (the postgrest-py convention) and read back the same way. If
    # Supabase ever rejects this, switch the columns to text holding base64.
    def get_data_key(self, org_id: str) -> Optional[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/org_data_keys",
            params={"org_id": f"eq.{org_id}", "limit": 1},
        )
        if r.status_code != 200:
            self._raise(r, "get_data_key")
        rows = r.json()
        return rows[0] if rows else None

    def get_data_key_version(
        self, org_id: str, version: int
    ) -> Optional[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/org_data_key_versions",
            params={
                "org_id": f"eq.{org_id}",
                "dek_version": f"eq.{version}",
                "limit": 1,
            },
        )
        if r.status_code != 200:
            self._raise(r, "get_data_key_version")
        rows = r.json()
        return rows[0] if rows else None

    def store_data_key(
        self, org_id: str, version: int, dek_wrapped_b64: str, kek_id: str
    ) -> dict[str, Any]:
        row = {
            "org_id": org_id,
            "dek_version": version,
            "dek_wrapped": dek_wrapped_b64,
            "kek_id": kek_id,
        }
        r = self._client.post(
            "/rest/v1/org_data_key_versions",
            content=json.dumps(row),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code not in (200, 201):
            self._raise(r, "store_data_key (versions)")
        # Make it current (upsert on org_id).
        r = self._client.post(
            "/rest/v1/org_data_keys",
            content=json.dumps(row),
            headers={"Prefer": "resolution=merge-duplicates,return=representation"},
            params={"on_conflict": "org_id"},
        )
        if r.status_code not in (200, 201):
            self._raise(r, "store_data_key (current)")
        return r.json()[0]

    # -- deletion --------------------------------------------------------------
    def delete_job(self, job_id: str) -> None:
        r = self._client.delete(
            "/rest/v1/ingest_jobs", params={"id": f"eq.{job_id}"}
        )
        if r.status_code not in (200, 204):
            self._raise(r, "delete_job")

    def delete_jobs_for_org(self, org_id: str) -> int:
        jobs = self.list_jobs(org_id, limit=10000)
        if not jobs:
            return 0
        r = self._client.delete(
            "/rest/v1/ingest_jobs", params={"org_id": f"eq.{org_id}"}
        )
        if r.status_code not in (200, 204):
            self._raise(r, "delete_jobs_for_org")
        return len(jobs)

    def mark_sentinel_run_deleted(self, accession: str) -> None:
        from datetime import datetime, timezone

        r = self._client.patch(
            "/rest/v1/sentinel_runs",
            params={"accession": f"eq.{accession}"},
            content=json.dumps(
                {"deleted_at": datetime.now(timezone.utc).isoformat()}
            ),
        )
        if r.status_code not in (200, 204):
            self._raise(r, "mark_sentinel_run_deleted")

    # -- audit -------------------------------------------------------------------
    def audit(
        self,
        org_id: str,
        actor_key_id: str,
        action: str,
        target: str,
        detail: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        r = self._client.post(
            "/rest/v1/audit_log",
            content=json.dumps(
                {
                    "org_id": org_id,
                    "actor_key_id": actor_key_id,
                    "action": action,
                    "target": target,
                    "detail": detail or {},
                }
            ),
            headers={"Prefer": "return=representation"},
        )
        if r.status_code not in (200, 201):
            self._raise(r, "audit")
        return r.json()[0]

    def list_audit(
        self, org_id: str, limit: int = 50, offset: int = 0
    ) -> list[dict[str, Any]]:
        r = self._client.get(
            "/rest/v1/audit_log",
            params={
                "org_id": f"eq.{org_id}",
                "order": "at.desc",
                "limit": limit,
                "offset": offset,
            },
        )
        if r.status_code != 200:
            self._raise(r, "list_audit")
        return r.json()

    def mirror_sentinel_run(self, row: dict[str, Any]) -> None:
        r = self._client.post(
            "/rest/v1/sentinel_runs",
            content=json.dumps(row),
            headers={"Prefer": "resolution=merge-duplicates,return=minimal"},
            params={"on_conflict": "accession"},
        )
        if r.status_code not in (200, 201, 204):
            self._raise(r, "mirror_sentinel_run")
