"""Create a per-org API key for the ingest service.

The plaintext key is printed ONCE and never stored; only its SHA-256 hash
goes into org_api_keys. Store the printed key somewhere safe (a password
manager, the lab's secrets store) — it cannot be recovered.

Usage:
    python -m service.make_key --org-id mirebalais --name mirebalais-lab-uploader
"""

from __future__ import annotations

import argparse
import secrets

from service.app.auth import hash_key
from service.app.config import Settings
from service.app.db import SupabaseRestDB


def main() -> None:
    ap = argparse.ArgumentParser(description="Create an ingest-service API key")
    ap.add_argument("--org-id", required=True, help="org identifier, e.g. mirebalais")
    ap.add_argument("--name", default="", help="human label for the key")
    args = ap.parse_args()

    settings = Settings()
    settings.require_db()
    db = SupabaseRestDB(settings.supabase_url, settings.supabase_service_role_key)

    plaintext = "sk_" + secrets.token_urlsafe(32)
    row = db.insert_api_key(args.org_id, hash_key(plaintext), args.name)
    print(f"org_id:  {args.org_id}")
    print(f"key id:  {row['id']}")
    print(f"API key (shown once — store it now):\n{plaintext}")


if __name__ == "__main__":
    main()
