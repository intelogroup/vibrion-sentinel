"""Bearer API-key auth, per-org scoping.

Keys are SHA-256 hashed at rest (sql/org_api_keys.sql); the plaintext is
shown once at creation. Every request resolves to an org_id, and all data
access is filtered by it. A key for org A can never observe org B's
resources: missing/other-org resources return 404, never 401/403, so the
API is not an existence oracle.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Annotated

from fastapi import Depends, HTTPException
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

from .db import JobDB

VALID_ROLES = ("admin", "member", "viewer")


@dataclass(frozen=True)
class Actor:
    """The authenticated caller: org identity + key identity + role.

    Roles live on API keys (Phase 2 decision: the service is
    machine-to-machine; human accounts/SSO are a Phase 3 portal concern).
    """

    org_id: str
    key_id: str
    role: str

_bearer = HTTPBearer(auto_error=False)


def hash_key(plaintext: str) -> str:
    return hashlib.sha256(plaintext.encode()).hexdigest()


def _unauthorized() -> HTTPException:
    return HTTPException(
        status_code=401,
        detail="invalid or missing API key",
        headers={"WWW-Authenticate": "Bearer"},
    )


def make_org_dependency(db: JobDB):
    async def _get_actor(
        credentials: Annotated[HTTPAuthorizationCredentials | None, Depends(_bearer)],
    ) -> Actor:
        if credentials is None or credentials.scheme.lower() != "bearer":
            raise _unauthorized()
        digest = hash_key(credentials.credentials)
        row = db.lookup_key(digest)
        if row is None or row.get("revoked_at"):
            # Deliberately no distinction between unknown and revoked.
            # (A timing oracle on hash compare is not a concern here: the
            # lookup is a DB equality check, not a secret comparison.)
            raise _unauthorized()
        return Actor(
            org_id=row["org_id"],
            key_id=row["id"],
            role=row.get("role") or "member",
        )

    return _get_actor
