"""Phase 2 envelope encryption: per-org DEKs, KEK-wrapped, AES-256-GCM.

Key hierarchy:

    SENTINEL_KEK (env, 32 bytes base64)  -- key-encryption key, "env-v1"
        wraps
    DEK per org (AES-256, random per org, versioned)
        encrypts
    upload bytes (AES-256-GCM, fresh 96-bit nonce per file)

Documented loudly: an env-held KEK is the stepping stone, not the end
state. A real KMS (per-deploy KEK in a managed HSM/KMS, with rotation and
audit) is the follow-up. What envelope encryption buys even with an
env KEK: per-org cryptographic isolation of bytes at rest (a leaked R2
credential or bucket listing yields ciphertext, and org A's DEK never
decrypts org B's objects), plus versioned rotation without re-encrypting
history.

DEK storage: `org_data_keys` (current pointer) + `org_data_key_versions`
(all versions, old ones retained for decrypting history). The wrapped DEK
is stored base64-encoded; see sql/org_data_keys.sql.

Per-file metadata lives in a sidecar object (`{key}.meta.json`), NOT R2
object metadata. Rationale: the sidecar works uniformly across backends
(MemoryStorage, R2Storage, anything implementing the Storage interface),
is visible/debuggable with plain object listing, and has no provider
metadata size/charset constraints. Content:

    {"encryption": "aes256-gcm-v1", "dek_version": 2,
     "nonce": "<base64>", "sha256_plaintext": "<hex>",
     "plaintext_size": 12345}

The `encryption` marker is explicit so "unencrypted" is never ambiguous:
a missing sidecar, or a sidecar naming another scheme, is a hard
decryption failure (fail closed), never a silent plaintext fallback.
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from typing import Any, Optional

from cryptography.exceptions import InvalidTag
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

from .db import JobDB
from .storage import Storage

ENCRYPTION_SCHEME = "aes256-gcm-v1"
KEK_ID = "env-v1"  # bump when the KEK itself rotates (a KMS-era concern)
NONCE_BYTES = 12  # 96-bit nonce, fresh per file
GCM_TAG_BYTES = 16

_CHUNK = 1 << 20


class DecryptionError(RuntimeError):
    """Raised when ciphertext cannot be authenticated/decrypted.

    Callers must treat this as a failed job, never as corrupt-but-usable
    data. In particular: never fall back to interpreting the bytes as
    plaintext.
    """


# ------------------------------------------------------------------ KEK/DEK

def load_kek() -> bytes:
    """Load the KEK from SENTINEL_KEK (base64, 32 bytes). Fail closed."""
    raw = os.environ.get("SENTINEL_KEK", "")
    if not raw:
        raise RuntimeError(
            "SENTINEL_KEK is not set: refusing to start without a key-encryption key. "
            "Generate one with: python -c \"import os,base64; "
            "print(base64.b64encode(os.urandom(32)).decode())\""
        )
    try:
        kek = base64.b64decode(raw, validate=True)
    except Exception as e:
        raise RuntimeError(f"SENTINEL_KEK is not valid base64: {e}")
    if len(kek) != 32:
        raise RuntimeError(
            f"SENTINEL_KEK must decode to 32 bytes (AES-256); got {len(kek)}"
        )
    return kek


def generate_dek() -> bytes:
    return os.urandom(32)


def wrap_dek(dek: bytes, kek: bytes) -> bytes:
    """Wrap a DEK with the KEK. Output: nonce(12) || ciphertext || tag(16).

    NOTE: for GCM, Encryptor.finalize() returns b"" — the tag is only
    available as encryptor.tag AFTER finalize(). Forgetting this silently
    drops authentication.
    """
    nonce = os.urandom(NONCE_BYTES)
    enc = Cipher(algorithms.AES(kek), modes.GCM(nonce)).encryptor()
    ct = enc.update(dek) + enc.finalize()
    return nonce + ct + enc.tag


def unwrap_dek(wrapped: bytes, kek: bytes) -> bytes:
    """Unwrap. Any tampering raises InvalidTag (an authentication failure)."""
    if len(wrapped) < NONCE_BYTES + GCM_TAG_BYTES:
        raise DecryptionError("wrapped DEK too short")
    nonce, blob = wrapped[:NONCE_BYTES], wrapped[NONCE_BYTES:]
    tag, ct = blob[-GCM_TAG_BYTES:], blob[:-GCM_TAG_BYTES]
    decryptor = Cipher(algorithms.AES(kek), modes.GCM(nonce, tag)).decryptor()
    try:
        return decryptor.update(ct) + decryptor.finalize()
    except InvalidTag as e:
        raise DecryptionError("DEK unwrap authentication failed") from e


class KeyProvider:
    """Per-org DEK lifecycle backed by the JobDB data-key tables."""

    def __init__(self, db: JobDB, kek: bytes, kek_id: str = KEK_ID):
        self._db = db
        self._kek = kek
        self._kek_id = kek_id

    def current(self, org_id: str) -> tuple[int, bytes]:
        """(version, dek) for the org's current DEK; creates v1 if none."""
        row = self._db.get_data_key(org_id)
        if row is None:
            return self._new_version(org_id, 1)
        return row["dek_version"], unwrap_dek(
            base64.b64decode(row["dek_wrapped"]), self._kek
        )

    def version(self, org_id: str, dek_version: int) -> bytes:
        row = self._db.get_data_key_version(org_id, dek_version)
        if row is None:
            raise DecryptionError(
                f"no DEK version {dek_version} for org {org_id!r}"
            )
        return unwrap_dek(base64.b64decode(row["dek_wrapped"]), self._kek)

    def rotate(self, org_id: str) -> int:
        """Generate a new DEK version; old versions stay for history decrypt.

        No eager re-encryption: it is O(all bytes) and buys nothing — old
        objects keep decrypting with their recorded version, new uploads use
        the new DEK. Rotation bounds *future* exposure after a suspected DEK
        compromise; it does not retroactively protect already-written bytes
        (nothing can, short of re-encryption).
        """
        row = self._db.get_data_key(org_id)
        nxt = (row["dek_version"] + 1) if row else 1
        self._new_version(org_id, nxt)
        return nxt

    def _new_version(self, org_id: str, version: int) -> tuple[int, bytes]:
        dek = generate_dek()
        wrapped = base64.b64encode(wrap_dek(dek, self._kek)).decode()
        self._db.store_data_key(org_id, version, wrapped, self._kek_id)
        return version, dek


# ------------------------------------------------------- encrypted storage

@dataclass
class _Session:
    encryptor: Any
    nonce: bytes
    version: int
    sha: Any  # hashlib._Hash
    size: int


class EncryptedStorage(Storage):
    """AES-256-GCM envelope encryption over any Storage backend.

    Write path streams: chunks are encrypted as they arrive (no whole-file
    buffering), the 16-byte GCM tag is appended at finalize(), and the
    sidecar records {dek_version, nonce, plaintext sha256/size}.
    Read path requires the sidecar; a missing/foreign sidecar is a
    DecryptionError, never a plaintext fallback.

    The org is parsed from the key layout `intake/{org_id}/...`, which the
    API assigns at upload creation.
    """

    def __init__(self, inner: Storage, keys: KeyProvider):
        self._inner = inner
        self._keys = keys
        self._sessions: dict[str, _Session] = {}

    # -- key layout helpers ------------------------------------------------
    @staticmethod
    def _org_of(key: str) -> str:
        parts = key.split("/")
        if len(parts) < 3 or parts[0] != "intake" or not parts[1]:
            raise ValueError(f"key {key!r} is not an intake object")
        return parts[1]

    @staticmethod
    def _meta_key(key: str) -> str:
        return key + ".meta.json"

    # -- Storage interface ---------------------------------------------------
    def create(self, key: str) -> None:
        version, dek = self._keys.current(self._org_of(key))
        nonce = os.urandom(NONCE_BYTES)
        encryptor = Cipher(algorithms.AES(dek), modes.GCM(nonce)).encryptor()
        self._sessions[key] = _Session(encryptor, nonce, version, hashlib.sha256(), 0)
        self._inner.create(key)

    def append(self, key: str, data: bytes) -> None:
        s = self._sessions[key]
        self._inner.append(key, s.encryptor.update(data))
        s.sha.update(data)
        s.size += len(data)

    def finalize(self, key: str) -> None:
        s = self._sessions.pop(key)
        # GCM: finalize() returns b""; the tag is encryptor.tag afterwards.
        tail = s.encryptor.finalize()
        tag = s.encryptor.tag
        if tail:
            self._inner.append(key, tail)
        self._inner.append(key, tag)
        self._inner.finalize(key)
        meta = {
            "encryption": ENCRYPTION_SCHEME,
            "dek_version": s.version,
            "nonce": base64.b64encode(s.nonce).decode(),
            "sha256_plaintext": s.sha.hexdigest(),
            "plaintext_size": s.size,
        }
        mk = self._meta_key(key)
        blob = json.dumps(meta).encode()
        self._inner.create(mk)
        self._inner.append(mk, blob)
        self._inner.finalize(mk)

    def abort(self, key: str) -> None:
        self._sessions.pop(key, None)
        self._inner.abort(key)

    def delete(self, key: str) -> None:
        self._sessions.pop(key, None)
        self._inner.delete(key)
        # Sidecar delete is idempotent on both backends (pop/delete_object).
        self._inner.delete(self._meta_key(key))

    def list_prefix(self, prefix: str) -> list[str]:
        return self._inner.list_prefix(prefix)

    def size(self, key: str) -> int:
        s = self._sessions.get(key)
        if s is not None:
            return s.size
        return int(self._read_meta(key)["plaintext_size"])

    def sha256(self, key: str) -> str:
        """Plaintext digest: from the live session, else from the sidecar."""
        s = self._sessions.get(key)
        if s is not None:
            return s.sha.hexdigest()
        return str(self._read_meta(key)["sha256_plaintext"])

    def download_to_file(self, key: str, dest_path: str) -> None:
        meta = self._read_meta(key)
        if meta.get("encryption") != ENCRYPTION_SCHEME:
            raise DecryptionError(
                f"object {key!r} uses unknown scheme {meta.get('encryption')!r}"
            )
        dek = self._keys.version(self._org_of(key), int(meta["dek_version"]))
        nonce = base64.b64decode(meta["nonce"])
        tmp = dest_path + ".enc"
        try:
            self._inner.download_to_file(key, tmp)
            self._stream_decrypt(tmp, dest_path, dek, nonce)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)

    # -- internals -----------------------------------------------------------
    def _read_meta(self, key: str) -> dict[str, Any]:
        mk = self._meta_key(key)
        fd, path = tempfile.mkstemp(prefix="sentinel-meta-")
        os.close(fd)
        try:
            try:
                self._inner.download_to_file(mk, path)
            except Exception as e:
                raise DecryptionError(
                    f"missing encryption metadata for {key!r}: refusing to "
                    f"interpret bytes as plaintext ({e})"
                )
            with open(path, "rb") as f:
                meta = json.loads(f.read().decode())
            for field in ("dek_version", "nonce", "sha256_plaintext", "plaintext_size"):
                if field not in meta:
                    raise DecryptionError(
                        f"encryption metadata for {key!r} lacks {field!r}"
                    )
            return meta
        finally:
            if os.path.exists(path):
                os.remove(path)

    def _stream_decrypt(self, src: str, dest: str, dek: bytes, nonce: bytes) -> None:
        n = os.path.getsize(src)
        if n < GCM_TAG_BYTES:
            raise DecryptionError("ciphertext shorter than the GCM tag")
        try:
            with open(src, "rb") as f:
                f.seek(n - GCM_TAG_BYTES)
                tag = f.read(GCM_TAG_BYTES)
                f.seek(0)
                decryptor = Cipher(
                    algorithms.AES(dek), modes.GCM(nonce, tag)
                ).decryptor()
                remaining = n - GCM_TAG_BYTES
                with open(dest, "wb") as out:
                    while remaining > 0:
                        chunk = f.read(min(_CHUNK, remaining))
                        if not chunk:
                            raise DecryptionError("truncated ciphertext")
                        out.write(decryptor.update(chunk))
                        remaining -= len(chunk)
                    decryptor.finalize()  # raises InvalidTag on tampering
        except InvalidTag as e:
            if os.path.exists(dest):
                os.remove(dest)
            raise DecryptionError("GCM authentication failed: bytes were tampered with") from e
        except DecryptionError:
            if os.path.exists(dest):
                os.remove(dest)
            raise
