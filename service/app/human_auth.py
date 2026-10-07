"""Human auth for the Phase 3 portal: password accounts, sessions, invites.

Design (spec §3):
- Password auth (simpler ops than magic-link email delivery; documented tradeoff).
- argon2 password hashing. Never MD5/SHA.
- Sessions: opaque tokens (sha256 at rest), 30-day expiry, HttpOnly cookie
  `sentinel_session`. Revoked on password change and on logout.
- org_members(user_id, org_id, role) — roles reuse the key role set
  (admin/member/viewer). API keys keep working unchanged.
- Signup creates a personal org + admin membership + one admin API key
  (returned once) so the new human can use the API immediately.
- Invites: admin creates a token (single-use, 7-day expiry); a logged-in user
  accepts it to join the org.

This module holds the service logic + FastAPI router. Persistence goes
through the JobDB interface (methods added in Phase 3); tests use FakeDB.
"""


import hashlib
import re
import secrets
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Annotated, Any, Optional

from argon2 import PasswordHasher
from argon2.exceptions import VerifyMismatchError
from fastapi import APIRouter, Depends, HTTPException, Request, Response
from fastapi.responses import JSONResponse

from .auth import VALID_ROLES, Actor, hash_key
from .db import JobDB, SupabaseError

SESSION_COOKIE = "sentinel_session"
SESSION_DAYS = 30
INVITE_DAYS = 7
MIN_PASSWORD_LEN = 12

_ph = PasswordHasher()


def hash_password(plaintext: str) -> str:
    return _ph.hash(plaintext)


def verify_password(password_hash: str, plaintext: str) -> bool:
    try:
        return _ph.verify(password_hash, plaintext)
    except VerifyMismatchError:
        return False
    except Exception:
        return False


def new_session_token() -> tuple[str, str]:
    """(plaintext token, sha256 hex for storage)."""
    token = secrets.token_urlsafe(32)
    return token, hashlib.sha256(token.encode()).hexdigest()


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _org_slug(email: str) -> str:
    local = email.split("@")[0].lower()
    slug = re.sub(r"[^a-z0-9]+", "-", local).strip("-") or "org"
    return f"{slug}-{secrets.token_hex(3)}"


def _valid_email(email: str) -> bool:
    return bool(re.match(r"^[^@\s]+@[^@\s]+\.[^@\s]+$", email or ""))


@dataclass(frozen=True)
class PortalUser:
    id: str
    email: str


def make_user_dependency(db: JobDB):
    """Resolve the logged-in user from the session cookie (None if absent)."""

    def _get_user(request: Request) -> Optional[PortalUser]:
        token = request.cookies.get(SESSION_COOKIE)
        if not token:
            return None
        try:
            row = db.get_session(hashlib.sha256(token.encode()).hexdigest())
        except SupabaseError:
            return None
        if row is None or row.get("revoked_at"):
            return None
        exp = row.get("expires_at")
        try:
            exp_dt = datetime.fromisoformat(str(exp).replace("Z", "+00:00"))
        except Exception:
            return None
        if exp_dt.tzinfo is None:
            exp_dt = exp_dt.replace(tzinfo=timezone.utc)
        if exp_dt < _now():
            return None
        try:
            user = db.get_user(row["user_id"])
        except SupabaseError:
            return None
        if user is None:
            return None
        return PortalUser(id=user["id"], email=user["email"])

    return _get_user


def _set_session_cookie(resp: Response, token: str, secure: bool) -> None:
    resp.set_cookie(
        SESSION_COOKIE,
        token,
        max_age=SESSION_DAYS * 86400,
        httponly=True,
        samesite="lax",
        secure=secure,
    )


def create_human_auth_router(
    db: JobDB, cookie_secure: bool = False
) -> APIRouter:
    router = APIRouter()
    get_user = make_user_dependency(db)
    UserDep = Annotated[Optional[PortalUser], Depends(get_user)]

    def _require_user(user: Optional[PortalUser]) -> PortalUser:
        if user is None:
            raise HTTPException(status_code=401, detail="login required")
        return user

    @router.post("/auth/signup", status_code=201)
    def signup(body: dict) -> Response:
        body = body or {}
        email = (body.get("email") or "").strip().lower()
        password = body.get("password") or ""
        if not _valid_email(email):
            raise HTTPException(status_code=400, detail="invalid email")
        if len(password) < MIN_PASSWORD_LEN:
            raise HTTPException(
                status_code=400,
                detail=f"password must be at least {MIN_PASSWORD_LEN} characters",
            )
        try:
            if db.get_user_by_email(email) is not None:
                raise HTTPException(status_code=409, detail="email already registered")
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")

        org_id = _org_slug(email)
        try:
            user = db.create_user(email, hash_password(password))
            db.add_org_member(user["id"], org_id, "admin")
            # One admin API key so the human can use the API immediately.
            plaintext = "sk_" + secrets.token_urlsafe(32)
            key_row = db.insert_api_key(
                org_id, hash_key(plaintext), "personal key", "admin"
            )
            token, token_hash = new_session_token()
            db.create_session(
                user["id"],
                token_hash,
                (_now() + timedelta(days=SESSION_DAYS)).isoformat(),
            )
            try:
                db.audit(org_id, f"user:{user['id']}", "user.signed_up", org_id, {})
            except SupabaseError:
                pass
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")

        resp = JSONResponse(
            {
                "id": user["id"],
                "email": email,
                "org_id": org_id,
                "api_key": plaintext,
                "session_token": token,
            },
            status_code=201,
        )
        _set_session_cookie(resp, token, cookie_secure)
        return resp

    @router.post("/auth/login")
    def login(body: dict) -> Response:
        body = body or {}
        email = (body.get("email") or "").strip().lower()
        password = body.get("password") or ""
        try:
            user = db.get_user_by_email(email)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if user is None or not verify_password(user["password_hash"], password):
            # No distinction between unknown email and wrong password.
            raise HTTPException(status_code=401, detail="invalid email or password")
        token, token_hash = new_session_token()
        try:
            db.update_user_last_login(user["id"])
            db.create_session(
                user["id"],
                token_hash,
                (_now() + timedelta(days=SESSION_DAYS)).isoformat(),
            )
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        resp = JSONResponse(
            {"id": user["id"], "email": user["email"], "session_token": token}
        )
        _set_session_cookie(resp, token, cookie_secure)
        return resp

    @router.post("/auth/logout", status_code=204)
    def logout(request: Request) -> Response:
        token = request.cookies.get(SESSION_COOKIE)
        if token:
            try:
                db.revoke_session(hashlib.sha256(token.encode()).hexdigest())
            except SupabaseError:
                pass
        resp = Response(status_code=204)
        resp.delete_cookie(SESSION_COOKIE)
        return resp

    @router.get("/auth/me")
    def me(user: UserDep) -> dict[str, Any]:
        u = _require_user(user)
        try:
            memberships = db.get_org_memberships(u.id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        return {
            "id": u.id,
            "email": u.email,
            "orgs": [
                {"org_id": m["org_id"], "role": m["role"]} for m in memberships
            ],
        }

    @router.post("/auth/change-password")
    def change_password(user: UserDep, body: dict) -> Response:
        u = _require_user(user)
        body = body or {}
        current = body.get("current_password") or ""
        new = body.get("new_password") or ""
        if len(new) < MIN_PASSWORD_LEN:
            raise HTTPException(
                status_code=400,
                detail=f"password must be at least {MIN_PASSWORD_LEN} characters",
            )
        try:
            row = db.get_user(u.id)
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if row is None or not verify_password(row["password_hash"], current):
            raise HTTPException(status_code=401, detail="current password is wrong")
        try:
            db.update_user_password(u.id, hash_password(new))
            n = db.revoke_all_sessions(u.id)  # all sessions die, incl. this one
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        resp = JSONResponse({"status": "ok", "sessions_revoked": n})
        resp.delete_cookie(SESSION_COOKIE)
        return resp

    # -- invites ---------------------------------------------------------
    @router.post("/org/invites", status_code=201)
    def create_invite(
        user: UserDep, request: Request, body: dict
    ) -> dict[str, Any]:
        """Admin invites someone to one of their orgs. Token returned once."""
        # Accept either a session user (admin of the org) or an admin API key.
        org_id, role = _resolve_admin(request, user, db)
        body = body or {}
        email = (body.get("email") or "").strip().lower()
        inv_role = body.get("role") or "member"
        if not _valid_email(email):
            raise HTTPException(status_code=400, detail="invalid email")
        if inv_role not in VALID_ROLES:
            raise HTTPException(
                status_code=400, detail=f"unknown role {inv_role!r}"
            )
        token = secrets.token_urlsafe(24)
        try:
            db.create_invite(
                org_id,
                email,
                inv_role,
                hashlib.sha256(token.encode()).hexdigest(),
                (_now() + timedelta(days=INVITE_DAYS)).isoformat(),
            )
            db.audit(org_id, "portal", "invite.created", email, {"role": inv_role})
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        return {"org_id": org_id, "email": email, "role": inv_role,
                "invite_token": token}

    @router.post("/auth/accept-invite")
    def accept_invite(user: UserDep, body: dict) -> dict[str, Any]:
        u = _require_user(user)
        token = (body or {}).get("token") or ""
        try:
            inv = db.get_invite(hashlib.sha256(token.encode()).hexdigest())
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if inv is None or inv.get("used_at"):
            raise HTTPException(status_code=404, detail="invite not found")
        try:
            exp = datetime.fromisoformat(
                str(inv["expires_at"]).replace("Z", "+00:00")
            )
        except Exception:
            raise HTTPException(status_code=400, detail="invite is invalid")
        if exp.tzinfo is None:
            exp = exp.replace(tzinfo=timezone.utc)
        if exp < _now():
            raise HTTPException(status_code=410, detail="invite expired")
        try:
            if db.get_org_member(u.id, inv["org_id"]) is None:
                db.add_org_member(u.id, inv["org_id"], inv["role"])
            db.use_invite(inv["id"])
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        return {"org_id": inv["org_id"], "role": inv["role"]}

    return router


def _resolve_admin(
    request: Request, user: Optional[PortalUser], db: JobDB
) -> tuple[str, str]:
    """Return (org_id, role) for an admin caller, via session or API key.

    The invite endpoint needs an org context: for session users it comes
    from the X-Org-Id header (must be an org they admin); for API keys the
    key's own org is used (key must be admin).
    """
    # API key path.
    auth = request.headers.get("authorization", "")
    if auth.lower().startswith("bearer "):
        try:
            row = db.lookup_key(hash_key(auth[7:].strip()))
        except SupabaseError as e:
            raise HTTPException(status_code=502, detail=f"database error: {e}")
        if row is None or row.get("revoked_at"):
            raise HTTPException(status_code=401, detail="invalid API key")
        if (row.get("role") or "member") != "admin":
            raise HTTPException(status_code=403, detail="admin required")
        return row["org_id"], "admin"
    # Session path.
    if user is None:
        raise HTTPException(status_code=401, detail="login required")
    org_id = request.headers.get("x-org-id") or ""
    if not org_id:
        raise HTTPException(
            status_code=400, detail="X-Org-Id header required for session callers"
        )
    try:
        mem = db.get_org_member(user.id, org_id)
    except SupabaseError as e:
        raise HTTPException(status_code=502, detail=f"database error: {e}")
    if mem is None:
        raise HTTPException(status_code=404, detail="not found")
    if mem.get("role") != "admin":
        raise HTTPException(status_code=403, detail="admin required")
    return org_id, "admin"
