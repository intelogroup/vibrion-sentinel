"""Encrypt-in-place migration for Phase 1 (unencrypted) intake objects.

Phase 1 wrote plaintext FASTQ to R2. Phase 2 encrypts everything with the
org's DEK. This script closes the gap: for each org, it ensures a DEK
exists, then encrypts every object under intake/{org_id}/ that lacks an
encryption sidecar, marking it encryption=aes256-gcm-v1.

Safety rules:
- Objects that already have a .meta.json sidecar are skipped (idempotent;
  re-running is safe).
- Objects belonging to jobs still in 'uploading' status are skipped (partial
  uploads must never be "completed" by migration).
- Sidecar files themselves (*.meta.json) are never treated as data.
- --dry-run reports the plan and changes nothing.

Usage:
    export SUPABASE_URL=... SUPABASE_SERVICE_ROLE_KEY=... \\
           R2_ENDPOINT_URL=... R2_ACCESS_KEY_ID=... R2_SECRET_ACCESS_KEY=... \\
           R2_BUCKET=... SENTINEL_KEK=...
    python -m service.migrate_encrypt --dry-run
    python -m service.migrate_encrypt --org-id mirebalais
"""

from __future__ import annotations

import argparse
import os
import tempfile
from typing import Any, Optional

from service.app.config import Settings
from service.app.crypto import EncryptedStorage, KeyProvider, load_kek
from service.app.db import JobDB, SupabaseRestDB
from service.app.storage import R2Storage, Storage

_CHUNK = 1 << 20


def discover_orgs(storage: Storage) -> list[str]:
    """Org ids from the intake/ prefix layout."""
    orgs: set[str] = set()
    for key in storage.list_prefix("intake/"):
        parts = key.split("/")
        if len(parts) >= 3 and parts[0] == "intake" and parts[1]:
            orgs.add(parts[1])
    return sorted(orgs)


def plan_org(
    db: JobDB, storage: Storage, org_id: str
) -> tuple[list[str], list[str], list[str]]:
    """(to_encrypt, already_done, skipped_partial) for one org.

    Reads through the RAW storage (these objects are plaintext by
    construction — that is what this migration fixes).
    """
    prefix = f"intake/{org_id}/"
    keys = storage.list_prefix(prefix)
    listed = set(keys)
    uploading = {
        j["r2_key"]
        for j in db.list_jobs(org_id, limit=10000)
        if j.get("r2_key") and j.get("status") == "uploading"
    }
    to_encrypt, already_done, skipped = [], [], []
    for k in keys:
        if k.endswith(".meta.json"):
            continue
        if k + ".meta.json" in listed:
            already_done.append(k)
        elif k in uploading:
            skipped.append(k)
        else:
            to_encrypt.append(k)
    return sorted(to_encrypt), sorted(already_done), sorted(skipped)


def encrypt_one(raw: Storage, enc: EncryptedStorage, key: str) -> None:
    """Download plaintext via raw, stream-encrypt via enc (same key)."""
    fd, tmp = tempfile.mkstemp(prefix="sentinel-migrate-")
    os.close(fd)
    try:
        raw.download_to_file(key, tmp)
        # enc.create truncates the destination via inner.create — the
        # plaintext is already safe in tmp.
        enc.create(key)
        with open(tmp, "rb") as f:
            while True:
                chunk = f.read(_CHUNK)
                if not chunk:
                    break
                enc.append(key, chunk)
        enc.finalize(key)  # writes the sidecar marking aes256-gcm-v1
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def migrate(
    db: JobDB,
    raw_storage: Storage,
    enc_storage: EncryptedStorage,
    org_ids: Optional[list[str]] = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    report: dict[str, Any] = {"dry_run": dry_run, "orgs": {}}
    for org_id in org_ids or discover_orgs(raw_storage):
        to_encrypt, already_done, skipped = plan_org(db, raw_storage, org_id)
        done: list[str] = []
        if not dry_run:
            for key in to_encrypt:
                encrypt_one(raw_storage, enc_storage, key)
                done.append(key)
        report["orgs"][org_id] = {
            "encrypted": done if not dry_run else [],
            "would_encrypt": to_encrypt if dry_run else [],
            "already_encrypted": already_done,
            "skipped_partial": skipped,
        }
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description="Encrypt Phase 1 plaintext objects in place")
    ap.add_argument("--org-id", action="append", default=None,
                    help="limit to an org (repeatable); default: all orgs found")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the plan without changing anything")
    args = ap.parse_args()

    settings = Settings()
    settings.require_db()
    settings.require_r2()
    kek = load_kek()  # fail closed

    db = SupabaseRestDB(settings.supabase_url, settings.supabase_service_role_key)
    raw = R2Storage(
        settings.r2_endpoint_url,
        settings.r2_access_key_id,
        settings.r2_secret_access_key,
        settings.r2_bucket,
    )
    enc = EncryptedStorage(raw, KeyProvider(db, kek))

    import json

    print(json.dumps(migrate(db, raw, enc, args.org_id, args.dry_run), indent=2))


if __name__ == "__main__":
    main()
