#!/usr/bin/env python3
"""Recurrence / cluster watcher rules over sentinel_runs rows.

Two rules, evaluated over the Haiti tier:

  R1 recurrence: a QC-pass toxigenic O1 detection in a region with no toxigenic
                 O1 seen for >= N quiet years (default 2).
  R2 cluster:    >= 2 QC-pass toxigenic O1 samples in the same region whose SNP
                 distances differ by <= 5 and whose collection dates fall
                 within 14 days.

All functions here are pure (no network, no Supabase): rows in, findings out.
I/O (fetching rows, writing sentinel_alerts, ntfy delivery) lives in main()
only, so the rules are unit-testable and validatable against the back-test
fixture without credentials.

Row shape (dict), mirroring the sentinel_runs columns the rules need:
    accession, collection_date ("YYYY-MM-DD" | None), country (str | None),
    priority ("haiti" | "global" | None), qc_status ("pass" | "fail"),
    snps_vs_7pet (int | None), v_cholerae_fraction (float | None),
    surveillance_loci: {locus: {"call": "present" | ...}}

Dates are parsed to datetime.date before any comparison -- never string-compared
('2018-06-23' > '2018' lexicographically, the classic trap). Rows without a
parseable full collection_date are excluded from temporal logic and counted in
the returned `skipped_undated` tally, never silently dropped.
"""

import argparse
import datetime
import json
import os
import sys
import urllib.parse
import urllib.request

DEFAULT_QUIET_YEARS = 2
DEFAULT_SNP_TOL = 5
DEFAULT_DAY_TOL = 14
MIN_PURITY = 0.9


# --------------------------------------------------------------------------
# primitives
# --------------------------------------------------------------------------

def parse_date(raw):
    """datetime.date from an ISO 'YYYY-MM-DD' string, else None."""
    if not raw:
        return None
    try:
        return datetime.date.fromisoformat(str(raw).strip()[:10])
    except ValueError:
        return None


def locus_present(row, locus):
    loci = row.get("surveillance_loci") or {}
    return (loci.get(locus) or {}).get("call") == "present"


def is_toxigenic_o1(row):
    """QC-pass, V. cholerae-pure, cholera-toxin positive."""
    if row.get("qc_status") != "pass":
        return False
    frac = row.get("v_cholerae_fraction")
    if frac is None or float(frac) < MIN_PURITY:
        return False
    return locus_present(row, "ctxA") and locus_present(row, "ctxB")


def is_haiti(row):
    """Haiti tier: explicit priority tag, or a Haiti country value."""
    if (row.get("priority") or "").lower() == "haiti":
        return True
    return (row.get("country") or "").lower().startswith("haiti")


def dated_detections(rows):
    """Haiti-tier toxigenic O1 detections with parseable dates, oldest first.

    Returns (detections, skipped_undated). Each detection is
    (date, accession, snps_vs_7pet).
    """
    dets, skipped = [], 0
    for r in rows:
        if not (is_haiti(r) and is_toxigenic_o1(r)):
            continue
        d = parse_date(r.get("collection_date"))
        if d is None:
            skipped += 1
            continue
        dets.append((d, r["accession"], r.get("snps_vs_7pet")))
    dets.sort(key=lambda t: (t[0], t[1]))
    return dets, skipped


# --------------------------------------------------------------------------
# R1: recurrence after a quiet period
# --------------------------------------------------------------------------

def detect_recurrence(rows, n_years=DEFAULT_QUIET_YEARS):
    """First sample ending a quiet period of >= n_years with no toxigenic O1.

    Returns (finding, skipped_undated). finding is None when nothing qualifies.
    The very first detection in history never fires: with no previous detection
    there is no measurable gap (this also makes cold-start safe).
    """
    dets, skipped = dated_detections(rows)
    if len(dets) < 2:
        return None, skipped
    quiet_days = int(n_years * 365)
    prev_date = dets[0][0]
    for d, acc, snps in dets[1:]:
        gap = (d - prev_date).days
        if gap >= quiet_days:
            return {
                "rule": "R1",
                "accession": acc,
                "collection_date": d.isoformat(),
                "gap_days": gap,
                "gap_years": round(gap / 365.25, 1),
                "previous_detection": prev_date.isoformat(),
                "snps_vs_7pet": snps,
                "dedupe_key": f"R1:haiti:{acc}",
            }, skipped
        prev_date = d
    return None, skipped


# --------------------------------------------------------------------------
# R2: tight spatiotemporal cluster
# --------------------------------------------------------------------------

def detect_clusters(rows, snp_tol=DEFAULT_SNP_TOL, day_tol=DEFAULT_DAY_TOL):
    """Groups of >= 2 detections within day_tol days and snp_tol SNP distance.

    Relatedness uses |snps_vs_7pet(a) - snps_vs_7pet(b)| <= snp_tol: a
    triangle-inequality proxy for pairwise distance, not a true pairwise
    measurement. Upgrade path: pairwise consensus distances from the
    R2-archived QC-pass consensuses.
    Returns (clusters, skipped_undated); each cluster is a dict with sorted
    accessions, date range, and SNP range.
    """
    dets, skipped = dated_detections(rows)
    dets = [t for t in dets if t[2] is not None]
    parent = list(range(len(dets)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    for i in range(len(dets)):
        for j in range(i + 1, len(dets)):
            if abs((dets[j][0] - dets[i][0]).days) > day_tol:
                continue
            if abs(dets[j][2] - dets[i][2]) <= snp_tol:
                union(i, j)

    groups = {}
    for i, (d, acc, snps) in enumerate(dets):
        groups.setdefault(find(i), []).append((d, acc, snps))
    clusters = []
    for members in groups.values():
        if len(members) < 2:
            continue
        members.sort()
        dates = [m[0] for m in members]
        snps = [m[2] for m in members]
        accs = sorted(m[1] for m in members)
        clusters.append({
            "rule": "R2",
            "accessions": accs,
            "n": len(accs),
            "date_from": min(dates).isoformat(),
            "date_to": max(dates).isoformat(),
            "snps_min": min(snps),
            "snps_max": max(snps),
            "dedupe_key": "R2:haiti:" + "-".join(accs),
        })
    clusters.sort(key=lambda c: c["date_from"])
    return clusters, skipped


# --------------------------------------------------------------------------
# I/O: Supabase fetch, sentinel_alerts upsert, ntfy delivery (main only)
# --------------------------------------------------------------------------

def fetch_rows():
    url = os.environ.get("SUPABASE_URL", "").rstrip("/")
    key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY", "")
    if not url or not key:
        sys.exit("SUPABASE_URL / SUPABASE_SERVICE_ROLE_KEY not set")
    cols = ("accession,collection_date,country,priority,qc_status,snps_vs_7pet,"
            "v_cholerae_fraction,surveillance_loci")
    q = urllib.parse.urlencode({"select": cols, "order": "accession"})
    req = urllib.request.Request(
        f"{url}/rest/v1/sentinel_runs?{q}",
        headers={"apikey": key, "Authorization": f"Bearer {key}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def upsert_alert(finding):
    url = os.environ.get("SUPABASE_URL", "").rstrip("/")
    key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY", "")
    row = {
        "rule_id": finding["rule"],
        "region": "haiti",
        "dedupe_key": finding["dedupe_key"],
        "accessions": ([finding["accession"]] if "accession" in finding
                       else finding["accessions"]),
        "payload": finding,
    }
    data = json.dumps(row).encode()
    req = urllib.request.Request(
        f"{url}/rest/v1/sentinel_alerts",
        data=data, method="POST",
        headers={"apikey": key, "Authorization": f"Bearer {key}",
                 "Content-Type": "application/json",
                 "Prefer": "resolution=merge-duplicates"})
    with urllib.request.urlopen(req, timeout=30) as r:
        r.read()
    print(f"alert recorded: {finding['rule']} {finding['dedupe_key']}")


def notify_ntfy(text):
    server = os.environ.get("NTFY_SERVER", "https://ntfy-alerts.clixen.app").rstrip("/")
    topic = os.environ.get("NTFY_TOPIC", "")
    if not topic:
        print("NTFY_TOPIC not set -- skipping push notification")
        return
    req = urllib.request.Request(
        f"{server}/{urllib.parse.quote(topic)}",
        data=text.encode(), method="POST",
        headers={"Title": "Vibrion Sentinel alert"})
    with urllib.request.urlopen(req, timeout=30) as r:
        r.read()
    print(f"ntfy sent to {topic}")


def format_finding(f):
    if f["rule"] == "R1":
        return (f"RECURRENCE: toxigenic O1 {f['accession']} collected "
                f"{f['collection_date']} after {f['gap_years']} quiet years "
                f"(previous: {f['previous_detection']}; {f['snps_vs_7pet']} SNPs vs 2010EL-1786).")
    c = f
    return (f"CLUSTER: {c['n']} related toxigenic O1 samples "
            f"{c['date_from']}..{c['date_to']} "
            f"({c['snps_min']}-{c['snps_max']} SNPs vs 2010EL-1786): "
            f"{', '.join(c['accessions'])}.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet-years", type=float, default=DEFAULT_QUIET_YEARS)
    ap.add_argument("--snp-tol", type=int, default=DEFAULT_SNP_TOL)
    ap.add_argument("--day-tol", type=int, default=DEFAULT_DAY_TOL)
    ap.add_argument("--dry-run", action="store_true",
                    help="evaluate and print, do not write alerts or notify")
    args = ap.parse_args()

    rows = fetch_rows()
    print(f"rows fetched: {len(rows)}")
    rec, skipped_r = detect_recurrence(rows, n_years=args.quiet_years)
    clusters, skipped_c = detect_clusters(rows, snp_tol=args.snp_tol,
                                         day_tol=args.day_tol)
    print(f"skipped undated: {max(skipped_r, skipped_c)}")
    findings = ([rec] if rec else []) + clusters
    if not findings:
        print("no rule firings")
        return
    for f in findings:
        print(format_finding(f))
        if args.dry_run:
            continue
        upsert_alert(f)
        notify_ntfy(format_finding(f))


if __name__ == "__main__":
    main()
