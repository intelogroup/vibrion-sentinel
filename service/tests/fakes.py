"""In-memory fakes for the Phase 1 service tests.

FakeDB implements service.app.db.JobDB without any network. It enforces
the same forward-only transition discipline as production (via the shared
JobDB.transition), so tests prove the state machine, not just the mocks.
"""

from __future__ import annotations

import copy
import threading
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
        self.keys_by_id: dict[str, dict[str, Any]] = {}  # id -> row
        self.mirrored: list[dict[str, Any]] = []  # sentinel_runs rows
        self.data_keys: dict[str, dict[str, Any]] = {}  # org_id -> current row
        self.data_key_versions: dict[tuple[str, int], dict[str, Any]] = {}
        self.audit_rows: list[dict[str, Any]] = []
        self.users: dict[str, dict[str, Any]] = {}          # id -> row
        self.users_by_email: dict[str, dict[str, Any]] = {}  # email -> row
        self.sessions: dict[str, dict[str, Any]] = {}       # token_hash -> row
        self.members: dict[tuple[str, str], dict[str, Any]] = {}  # (user_id, org_id) -> row
        self.invites: dict[str, dict[str, Any]] = {}       # token_hash -> row
        self.invites_by_id: dict[str, dict[str, Any]] = {}
        self._lock = threading.Lock()

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
        # Lock: models the atomic FOR UPDATE SKIP LOCKED claim so the
        # double-claim test is meaningful.
        with self._lock:
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

    def insert_api_key(
        self, org_id: str, key_hash: str, name: str, role: str = "member"
    ) -> dict[str, Any]:
        row = {
            "id": str(uuid.uuid4()),
            "org_id": org_id,
            "key_hash": key_hash,
            "name": name,
            "role": role,
            "created_at": _now(),
            "revoked_at": None,
        }
        self.keys[key_hash] = row
        self.keys_by_id[row["id"]] = row
        return copy.deepcopy(row)

    def list_api_keys(self, org_id: str) -> list[dict[str, Any]]:
        rows = [k for k in self.keys_by_id.values() if k["org_id"] == org_id]
        rows.sort(key=lambda k: k["created_at"], reverse=True)
        return copy.deepcopy(rows)

    def get_api_key(self, org_id: str, key_id: str) -> Optional[dict[str, Any]]:
        row = self.keys_by_id.get(key_id)
        if row is None or row["org_id"] != org_id:
            return None
        return copy.deepcopy(row)

    def revoke_api_key(self, org_id: str, key_id: str) -> Optional[dict[str, Any]]:
        row = self.keys_by_id.get(key_id)
        if row is None or row["org_id"] != org_id:
            return None
        row["revoked_at"] = _now()
        return copy.deepcopy(row)

    def revoke_all_keys(self, org_id: str) -> int:
        n = 0
        for row in self.keys_by_id.values():
            if row["org_id"] == org_id and not row.get("revoked_at"):
                row["revoked_at"] = _now()
                n += 1
        return n

    def revoke_key(self, key_hash: str) -> None:
        # Legacy helper kept for older tests.
        self.keys[key_hash]["revoked_at"] = _now()

    # -- data keys -------------------------------------------------------
    def get_data_key(self, org_id: str) -> Optional[dict[str, Any]]:
        row = self.data_keys.get(org_id)
        return copy.deepcopy(row) if row else None

    def get_data_key_version(
        self, org_id: str, version: int
    ) -> Optional[dict[str, Any]]:
        row = self.data_key_versions.get((org_id, version))
        return copy.deepcopy(row) if row else None

    def store_data_key(
        self, org_id: str, version: int, dek_wrapped_b64: str, kek_id: str
    ) -> dict[str, Any]:
        row = {
            "org_id": org_id,
            "dek_version": version,
            "dek_wrapped": dek_wrapped_b64,
            "kek_id": kek_id,
            "created_at": _now(),
            "rotated_at": None,
        }
        prev = self.data_keys.get(org_id)
        if prev is not None:
            prev["rotated_at"] = _now()
        self.data_key_versions[(org_id, version)] = row
        self.data_keys[org_id] = row
        return copy.deepcopy(row)

    # -- deletion ----------------------------------------------------------
    def delete_job(self, job_id: str) -> None:
        self.jobs.pop(job_id, None)

    def delete_jobs_for_org(self, org_id: str) -> int:
        doomed = [jid for jid, j in self.jobs.items() if j["org_id"] == org_id]
        for jid in doomed:
            del self.jobs[jid]
        return len(doomed)

    def mark_sentinel_run_deleted(self, accession: str) -> None:
        for r in self.mirrored:
            if r["accession"] == accession:
                r["deleted_at"] = _now()

    # -- audit -----------------------------------------------------------------
    def audit(
        self,
        org_id: str,
        actor_key_id: str,
        action: str,
        target: str,
        detail: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        row = {
            "id": str(uuid.uuid4()),
            "org_id": org_id,
            "actor_key_id": actor_key_id,
            "action": action,
            "target": target,
            "detail": copy.deepcopy(detail or {}),
            "at": _now(),
        }
        self.audit_rows.append(row)
        return copy.deepcopy(row)

    def list_audit(
        self, org_id: str, limit: int = 50, offset: int = 0
    ) -> list[dict[str, Any]]:
        rows = [r for r in self.audit_rows if r["org_id"] == org_id]
        rows.sort(key=lambda r: r["at"], reverse=True)
        return copy.deepcopy(rows[offset : offset + limit])

    # -- human auth (Phase 3) ------------------------------------------------
    def create_user(self, email: str, password_hash: str) -> dict[str, Any]:
        if email in self.users_by_email:
            raise ValueError("email already registered")
        row = {
            "id": str(uuid.uuid4()),
            "email": email,
            "password_hash": password_hash,
            "created_at": _now(),
        }
        self.users[row["id"]] = row
        self.users_by_email[email] = row
        return copy.deepcopy(row)

    def get_user(self, user_id: str) -> Optional[dict[str, Any]]:
        row = self.users.get(user_id)
        return copy.deepcopy(row) if row else None

    def get_user_by_email(self, email: str) -> Optional[dict[str, Any]]:
        row = self.users_by_email.get(email)
        return copy.deepcopy(row) if row else None

    def update_user_password(self, user_id: str, password_hash: str) -> None:
        self.users[user_id]["password_hash"] = password_hash

    def update_user_last_login(self, user_id: str) -> None:
        self.users[user_id]["last_login_at"] = _now()

    def create_session(
        self, user_id: str, token_hash: str, expires_at: str
    ) -> dict[str, Any]:
        row = {
            "id": str(uuid.uuid4()),
            "user_id": user_id,
            "token_hash": token_hash,
            "created_at": _now(),
            "expires_at": expires_at,
            "revoked_at": None,
        }
        self.sessions[token_hash] = row
        return copy.deepcopy(row)

    def get_session(self, token_hash: str) -> Optional[dict[str, Any]]:
        row = self.sessions.get(token_hash)
        return copy.deepcopy(row) if row else None

    def revoke_session(self, token_hash: str) -> None:
        row = self.sessions.get(token_hash)
        if row:
            row["revoked_at"] = _now()

    def revoke_all_sessions(self, user_id: str) -> int:
        n = 0
        for row in self.sessions.values():
            if row["user_id"] == user_id and not row.get("revoked_at"):
                row["revoked_at"] = _now()
                n += 1
        return n

    def add_org_member(
        self, user_id: str, org_id: str, role: str
    ) -> dict[str, Any]:
        row = {
            "user_id": user_id,
            "org_id": org_id,
            "role": role,
            "created_at": _now(),
        }
        self.members[(user_id, org_id)] = row
        return copy.deepcopy(row)

    def get_org_member(
        self, user_id: str, org_id: str
    ) -> Optional[dict[str, Any]]:
        row = self.members.get((user_id, org_id))
        return copy.deepcopy(row) if row else None

    def get_org_memberships(self, user_id: str) -> list[dict[str, Any]]:
        return copy.deepcopy(
            [m for (uid, _), m in self.members.items() if uid == user_id]
        )

    def create_invite(
        self,
        org_id: str,
        email: str,
        role: str,
        token_hash: str,
        expires_at: str,
    ) -> dict[str, Any]:
        row = {
            "id": str(uuid.uuid4()),
            "org_id": org_id,
            "email": email,
            "role": role,
            "token_hash": token_hash,
            "created_at": _now(),
            "expires_at": expires_at,
            "used_at": None,
        }
        self.invites[token_hash] = row
        self.invites_by_id[row["id"]] = row
        return copy.deepcopy(row)

    def get_invite(self, token_hash: str) -> Optional[dict[str, Any]]:
        row = self.invites.get(token_hash)
        return copy.deepcopy(row) if row else None

    def use_invite(self, invite_id: str) -> None:
        row = self.invites_by_id.get(invite_id)
        if row:
            row["used_at"] = _now()

    def get_job_any(self, job_id: str) -> Optional[dict[str, Any]]:
        job = self.jobs.get(job_id)
        return copy.deepcopy(job) if job else None

    # -- sentinel_runs mirror ---------------------------------------------
    def mirror_sentinel_run(self, row: dict[str, Any]) -> None:
        self.mirrored = [r for r in self.mirrored if r["accession"] != row["accession"]]
        self.mirrored.append(copy.deepcopy(row))
