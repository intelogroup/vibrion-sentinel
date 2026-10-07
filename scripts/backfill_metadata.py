#!/usr/bin/env python3
"""Backfill collection_date / country on sentinel_runs from a cohort TSV.

The watcher rules need temporal + geographic metadata that older pipeline runs
never recorded. This script PATCHes existing sentinel_runs rows from a manifest
TSV with run_accession / collection_date / country columns (e.g.
cohorts/haiti-baseline.tsv). Rows whose accession is not in sentinel_runs are
reported and skipped; existing non-null values are never overwritten unless
--overwrite is given.

Usage:
    python3 scripts/backfill_metadata.py --tsv cohorts/haiti-baseline.tsv --dry-run
    python3 scripts/backfill_metadata.py --tsv cohorts/haiti-baseline.tsv   # needs SUPABASE_URL + SUPABASE_SERVICE_ROLE_KEY
"""

import argparse
import csv
import json
import os
import sys
import urllib.parse
import urllib.request


def load_mapping(path):
    """accession -> (collection_date|None, country|None) from a cohort TSV."""
    mapping = {}
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#") and l.strip()]
    for r in csv.DictReader(lines, delimiter="\t"):
        acc = (r.get("run_accession") or "").strip()
        if not acc:
            continue
        date = (r.get("collection_date") or "").strip() or None
        country = (r.get("country") or "").strip() or None
        mapping[acc] = (date, country)
    return mapping


def patch_row(base, key, accession, date, country, overwrite):
    """PATCH one sentinel_runs row. Returns (status, detail)."""
    body = {}
    if date is not None:
        body["collection_date"] = date
    if country is not None:
        body["country"] = country
    if not body:
        return "skip", "no metadata in TSV"
    if not overwrite:
        # Only fill nulls: fetch current values first.
        q = urllib.parse.urlencode(
            {"accession": f"eq.{accession}", "select": "collection_date,country"})
        req = urllib.request.Request(
            f"{base}/rest/v1/sentinel_runs?{q}",
            headers={"apikey": key, "Authorization": f"Bearer {key}"})
        with urllib.request.urlopen(req, timeout=30) as r:
            rows = json.loads(r.read().decode())
        if not rows:
            return "skip", "accession not in sentinel_runs"
        cur = rows[0]
        if date is not None and cur.get("collection_date"):
            body.pop("collection_date", None)
        if country is not None and cur.get("country"):
            body.pop("country", None)
        if not body:
            return "skip", "already populated"
    q = urllib.parse.urlencode({"accession": f"eq.{accession}"})
    req = urllib.request.Request(
        f"{base}/rest/v1/sentinel_runs?{q}",
        data=json.dumps(body).encode(), method="PATCH",
        headers={"apikey": key, "Authorization": f"Bearer {key}",
                 "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        r.read()
    return "updated", json.dumps(body)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--overwrite", action="store_true",
                    help="overwrite existing non-null values (default: fill nulls only)")
    args = ap.parse_args()

    mapping = load_mapping(args.tsv)
    print(f"mapping: {len(mapping)} accessions from {args.tsv}")
    if args.dry_run:
        n_date = sum(1 for d, _ in mapping.values() if d)
        n_cty = sum(1 for _, c in mapping.values() if c)
        print(f"dry-run: {n_date} with collection_date, {n_cty} with country")
        return

    url = os.environ.get("SUPABASE_URL", "").rstrip("/")
    key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY", "")
    if not url or not key:
        sys.exit("SUPABASE_URL / SUPABASE_SERVICE_ROLE_KEY not set")

    from collections import Counter
    stats = Counter()
    for acc, (date, country) in mapping.items():
        try:
            st, detail = patch_row(url, key, acc, date, country, args.overwrite)
        except Exception as e:  # noqa: BLE001 - report per-row, keep going
            st, detail = "error", str(e)[:120]
        stats[st] += 1
        if st in ("error",):
            print(f"  {acc}: {st} ({detail})")
    print("result:", dict(stats))


if __name__ == "__main__":
    main()
