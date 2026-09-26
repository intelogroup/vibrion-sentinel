#!/usr/bin/env python3
"""Build a reproducible cohort manifest from ENA.

Cohorts are committed TSVs so a back-test runs against a fixed sample list
rather than whatever ENA returns on the day. Uses ENA's structured fields
(country, collection_date, instrument_platform, library_strategy) rather
than the free-text matching in scripts/watch_sra.py -- 'Haiti[All Fields]'
also hits papers and submitter addresses, ENA's country field does not.

Usage:
    python3 scripts/build_cohort.py --cohort haiti-baseline --out cohorts/haiti-baseline.tsv
    python3 scripts/build_cohort.py --cohort backtest-2022 --limit 25 \
        --out cohorts/backtest-2022.tsv
"""

import argparse
import sys
import urllib.parse
import urllib.request

import build_cohort_lib

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"

FIELDS = [
    "run_accession",
    "sample_accession",
    "study_accession",
    "country",
    "collection_date",
    "first_public",
    "instrument_platform",
    "library_strategy",
    "read_count",
]

# V. cholerae taxid 666. Illumina WGS only: the lite pipeline maps short
# reads with bwa mem, so long-read or amplicon runs would map badly and
# produce junk rather than failing loudly.
SHORT_READ_WGS = 'instrument_platform="ILLUMINA" AND library_strategy="WGS"'

COHORTS = {
    # Everything ENA has for Haiti -- the 2010-to-date baseline.
    "haiti-baseline": f'tax_eq(666) AND country="Haiti*" AND {SHORT_READ_WGS}',
    # Back-test target: the published Haiti 2022 resurgence project.
    # Verified via ENA: 96 runs, all Haiti, collections 2017-06..2022-11.
    "backtest-2022": (
        f'study_accession="PRJNA900623" AND {SHORT_READ_WGS} '
        'AND collection_date>=2022-01-01'
    ),
}


def ena_search(query, limit=0):
    params = {
        "result": "read_run",
        "query": query,
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "limit": str(limit),
    }
    url = ENA_SEARCH + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=120) as r:
        text = r.read().decode()
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[1:]]
    return header, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=sorted(COHORTS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0,
                    help="cap the cohort size (0 = all); applied after sorting "
                         "by collection_date so a subset is deterministic")
    args = ap.parse_args()

    query = COHORTS[args.cohort]
    print(f"cohort: {args.cohort}", file=sys.stderr)
    print(f"query:  {query}", file=sys.stderr)
    header, rows = ena_search(query)
    print(f"ENA returned {len(rows)} runs", file=sys.stderr)
    if not rows:
        sys.exit("no runs returned -- refusing to write an empty manifest")

    rows = build_cohort_lib.sort_and_limit(header, rows, args.limit)
    if args.limit:
        print(f"limited to {len(rows)} runs", file=sys.stderr)

    with open(args.out, "w") as f:
        f.write(build_cohort_lib.format_manifest(args.cohort, query, header, rows))
    print(f"wrote {args.out} ({len(rows)} runs)", file=sys.stderr)

    dr = build_cohort_lib.date_range(header, rows)
    if dr:
        print(f"collection dates: {dr[0]} .. {dr[1]}", file=sys.stderr)


if __name__ == "__main__":
    main()
