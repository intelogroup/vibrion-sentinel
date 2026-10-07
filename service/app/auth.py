"""Bearer API-key auth, per-org scoping.

Keys are SHA-256 hashed at rest (sql/org_api_keys.sql); the plaintext is
shown once at creation. Every request resolves to an org_id, and all data
access is filtered by it. A key for org A can never observe org B's
resources: missing/other-org resources return 404, never 401/403, so the
API is not an existence oracle.
"""

from __future__ import annotations

import hashlib
from typing import Annotated

from fastapi import Depends, HTTPException
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

from .db import JobDB

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
    async def _get_org(
        credentials: Annotated[HTTPAuthorizationCredentials | None, Depends(_bearer)],
    ) -> str:
        if credentials is None or credentials.scheme.lower() != "bearer":
            raise _unauthorized()
        digest = hash_key(credentials.credentials)
        row = db.lookup_key(digest)
        if row is None or row.get("revoked_at"):
            # Deliberately no distinction between unknown and revoked.
            # (A timing oracle on hash compare is not a concern here: the
            # lookup is a DB equality check, not a secret comparison.)
            raise _unauthorized()
        return row["org_id"]

    return _get_org
