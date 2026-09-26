#!/usr/bin/env python3
"""Vibrion Sentinel watcher — find new V. cholerae SRA accessions.

Queries NCBI E-utilities for recent Vibrio cholerae sequencing runs,
diffs against state/seen_accessions.txt, and emits a GitHub Actions
matrix of new accessions to process.

Priority tiers:
  - haiti: accession linked to Haiti/Caribbean (any date, not yet seen)
  - global: new in the last --days days

Usage:
    python3 scripts/watch_sra.py --state state/seen_accessions.txt \
        --days 7 --max-new 10 --out-matrix /tmp/matrix.json
"""

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "vibrion-sentinel"
EMAIL = "vibrion-sentinel@intelogroup.example"  # NCBI etiquette; replace with real contact


def eu(path, params):
    params = dict(params, tool=TOOL, email=EMAIL, retmode="json")
    url = f"{EUTILS}/{path}?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def eu_xml(path, params):
    params = dict(params, tool=TOOL, email=EMAIL)
    url = f"{EUTILS}/{path}?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=60) as r:
        return ET.fromstring(r.read())


def esearch_ids(term, reldate=None):
    """Return SRA ids matching term (paged)."""
    ids = []
    retstart, retmax = 0, 500
    while True:
        p = {"db": "sra", "term": term, "retstart": retstart, "retmax": retmax}
        if reldate:
            p.update(datetype="pdat", reldate=str(reldate))
        d = eu("esearch.fcgi", p)
        r = d["esearchresult"]
        ids.extend(r.get("idlist", []))
        total = int(r["count"])
        retstart += retmax
        if retstart >= total or not r.get("idlist"):
            break
        time.sleep(0.4)  # stay under 3 req/s without API key
    return ids


def extract_run_accessions(esummary_root):
    """Run accessions (SRR/ERR/DRR...) from a parsed esummary XML root
    (containing one or more DocSum elements). Pulled out of
    summarize_run_ids so the escaped-XML parsing can be tested with a
    synthetic tree, no network call needed.

    NCBI nests the actual run accessions inside a "Runs" Item as escaped
    XML text: <Run acc="SRR..." .../> (one experiment may carry several
    runs). An earlier version of this function looked for an Item named
    "Run" (singular) instead, which never matches -- esummary_root.iter()
    is find-based, so a wrong name fails silently rather than raising,
    and the watcher returned zero new accessions on every scheduled run."""
    accs = []
    for doc in esummary_root.iter("DocSum"):
        for item in doc.iter("Item"):
            if item.get("Name") == "Runs":
                accs.extend(re.findall(r'<Run acc="([^"]+)"', item.text or ""))
    return accs


def summarize_run_ids(ids):
    """Map SRA internal ids -> run accessions (SRR...)."""
    accs = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        root = eu_xml("esummary.fcgi", {"db": "sra", "id": ",".join(chunk)})
        accs.extend(extract_run_accessions(root))
        time.sleep(0.4)
    return accs


def load_seen(path):
    """Set of already-processed accessions from state/seen_accessions.txt.
    A missing file means nothing has been seen yet, not an error."""
    try:
        with open(path) as f:
            return {l.strip() for l in f if l.strip() and not l.startswith("#")}
    except FileNotFoundError:
        return set()


def build_new_items(geo_accs, recent_accs, seen, max_new):
    """Haiti-tier accessions first, then global-tier, each deduplicated
    against `seen`, against each other (an accession already added as
    "haiti" is not re-added as "global"), and against repeats within the
    same tier (esummary chunking could plausibly return one accession
    twice) -- capped at max_new total."""
    new_items = []
    already = set()
    for acc in geo_accs:
        if acc not in seen and acc not in already:
            new_items.append({"accession": acc, "priority": "haiti"})
            already.add(acc)
    for acc in recent_accs:
        if acc not in seen and acc not in already:
            new_items.append({"accession": acc, "priority": "global"})
            already.add(acc)
    return new_items[:max_new]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--days", type=int, default=7)
    ap.add_argument("--max-new", type=int, default=10)
    ap.add_argument("--out-matrix", required=True)
    args = ap.parse_args()

    seen = set()
    try:
        with open(args.state) as f:
            seen = {l.strip() for l in f if l.strip() and not l.startswith("#")}
    except FileNotFoundError:
        pass
    print(f"seen accessions: {len(seen)}", file=sys.stderr)

    # Tier 1: Haiti/Caribbean-linked (any date)
    geo_term = '("Vibrio cholerae"[Organism]) AND (Haiti[All Fields] OR Caribbean[All Fields])'
    geo_ids = esearch_ids(geo_term)
    print(f"Haiti/Caribbean SRA records: {len(geo_ids)}", file=sys.stderr)

    # Tier 2: global, recent
    recent_ids = esearch_ids('"Vibrio cholerae"[Organism]', reldate=args.days)
    print(f"global SRA records (last {args.days}d): {len(recent_ids)}", file=sys.stderr)

    new_items = build_new_items(
        summarize_run_ids(geo_ids), summarize_run_ids(recent_ids), seen, args.max_new)
    print(f"new accessions to process: {len(new_items)}", file=sys.stderr)

    matrix = {"include": new_items}
    with open(args.out_matrix, "w") as f:
        json.dump(matrix, f)
    # also emit accessions list for the state-update job
    print("MATRIX_JSON=" + json.dumps(matrix))


if __name__ == "__main__":
    main()
