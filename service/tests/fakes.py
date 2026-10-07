"""In-memory fakes for the Phase 1 service tests.

FakeDB implements service.app.db.JobDB without any network. It enforces
the same forward-only transition discipline as production (via the shared
JobDB.transition), so tests prove the state machine, not just the mocks.
"""

from __future__ import annotations

import copy
import uuid
from datetime import datetime, timezone
from typing import Any, Optional

from service.app.db import JobDB


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class FakeDB(JobDB):
    def __init__(self) -> None:
        self.jobs: dict[str, dict[str, Any]] = {}
        self.keys: dict[str, dict[str, Any]] = {}  # key_hash -> row
        self.mirrored: list[dict[str, Any]] = []  # sentinel_runs rows

    # -- jobs -----------------------------------------------------------
    def create_job(self, row: dict[str, Any]) -> dict[str, Any]:
        job = {
            "id": str(uuid.uuid4()),
            "upload_offset": 0,
            "status": "uploading",
            "tier": "lite",
            "basecaller_model": None,
            "size_bytes": None,
            "sha256": None,
            "client_checksum": None,
            "status_reason": None,
            "r2_key": None,
            "report": None,
            "created_at": _now(),
            "updated_at": _now(),
            **copy.deepcopy(row),
        }
        self.jobs[job["id"]] = job
        return copy.deepcopy(job)

    def get_job(self, org_id: str, job_id: str) -> Optional[dict[str, Any]]:
        job = self.jobs.get(job_id)
        if job is None or job["org_id"] != org_id:
            return None
        return copy.deepcopy(job)

    def list_jobs(self, org_id: str, limit: int = 50) -> list[dict[str, Any]]:
        rows = [j for j in self.jobs.values() if j["org_id"] == org_id]
        rows.sort(key=lambda j: j["created_at"], reverse=True)
        return copy.deepcopy(rows[:limit])

    def update_job(self, job_id: str, fields: dict[str, Any]) -> dict[str, Any]:
        job = self.jobs[job_id]
        job.update(copy.deepcopy(fields))
        job["updated_at"] = _now()
        return copy.deepcopy(job)

    def claim_next_job(self) -> Optional[dict[str, Any]]:
        queued = [j for j in self.jobs.values() if j["status"] == "queued"]
        if not queued:
            return None
        queued.sort(key=lambda j: j["created_at"])
        job = queued[0]
        # Same legal move the SQL claim function makes: queued -> running.
        self.transition(job, "running")
        return copy.deepcopy(self.jobs[job["id"]])

    # -- api keys --------------------------------------------------------
    def lookup_key(self, key_hash: str) -> Optional[dict[str, Any]]:
        row = self.keys.get(key_hash)
        return copy.deepcopy(row) if row else None

    def insert_api_key(self, org_id: str, key_hash: str, name: str) -> dict[str, Any]:
        row = {
            "id": str(uuid.uuid4()),
            "org_id": org_id,
            "key_hash": key_hash,
            "name": name,
            "created_at": _now(),
            "revoked_at": None,
        }
        self.keys[key_hash] = row
        return copy.deepcopy(row)

    def revoke_key(self, key_hash: str) -> None:
        self.keys[key_hash]["revoked_at"] = _now()

    # -- sentinel_runs mirror ---------------------------------------------
    def mirror_sentinel_run(self, row: dict[str, Any]) -> None:
        self.mirrored = [r for r in self.mirrored if r["accession"] != row["accession"]]
        self.mirrored.append(copy.deepcopy(row))
