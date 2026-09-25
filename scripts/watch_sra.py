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


def summarize_run_ids(ids):
    """Map SRA internal ids -> run accessions (SRR...)."""
    accs = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        root = eu_xml("esummary.fcgi", {"db": "sra", "id": ",".join(chunk)})
        for doc in root.iter("DocSum"):
            # Runs are in a "Runs" item as escaped XML: <Run acc="SRR..." .../>
            # (one experiment may carry several runs)
            for item in doc.iter("Item"):
                if item.get("Name") == "Runs":
                    accs.extend(re.findall(r'<Run acc="([^"]+)"', item.text or ""))
        time.sleep(0.4)
    return accs


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

    new_items = []
    for acc in summarize_run_ids(geo_ids):
        if acc not in seen:
            new_items.append({"accession": acc, "priority": "haiti"})
    for acc in summarize_run_ids(recent_ids):
        if acc not in seen and all(i["accession"] != acc for i in new_items):
            new_items.append({"accession": acc, "priority": "global"})

    new_items = new_items[:args.max_new]
    print(f"new accessions to process: {len(new_items)}", file=sys.stderr)

    matrix = {"include": new_items}
    with open(args.out_matrix, "w") as f:
        json.dump(matrix, f)
    # also emit accessions list for the state-update job
    print("MATRIX_JSON=" + json.dumps(matrix))


if __name__ == "__main__":
    main()
