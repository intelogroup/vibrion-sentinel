#!/usr/bin/env python3
"""Compare pipeline output for a cohort against published expectations.

Reads per-run reports either from a local results tree (report.json per
accession) or from the Supabase sentinel_runs table, and checks them
against a cohort expectations file (see cohorts/*.expect.yaml).

A failed check is a finding about the pipeline. It is reported as a
failure; it is not a cue to loosen the threshold until it passes.

Usage:
    python3 scripts/backtest_report.py --expect cohorts/backtest-2022.expect.yaml \
        --results-dir results
    python3 scripts/backtest_report.py --expect cohorts/backtest-2022.expect.yaml \
        --supabase   # needs SUPABASE_URL + SUPABASE_SERVICE_ROLE_KEY
"""

import argparse
import json
import os
import statistics
import sys
import urllib.parse
import urllib.request


def load_expect(path):
    """Minimal YAML subset reader (avoids a pyyaml dep for CI)."""
    import re
    data = {}
    stack = [(0, data)]
    pending_block = None
    for raw in open(path):
        line = raw.rstrip("\n")
        if pending_block is not None:
            if line.startswith("  ") or not line.strip():
                continue
            pending_block = None
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        key, _, val = line.strip().partition(":")
        val = val.strip()
        while stack and indent < stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if val in ("", ">", "|"):
            if val in (">", "|"):
                pending_block = key
                parent[key] = ""
                continue
            parent[key] = {}
            stack.append((indent + 2, parent[key]))
        else:
            if re.fullmatch(r"-?\d+", val):
                parent[key] = int(val)
            elif re.fullmatch(r"-?\d*\.\d+", val):
                parent[key] = float(val)
            else:
                parent[key] = val.strip("'\"")
    return data


def cohort_accessions(manifest):
    accs = []
    for line in open(manifest):
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if parts[0] == "run_accession":
            continue
        accs.append(parts[0])
    return accs


def from_results_dir(accs, results_dir):
    rows = []
    for a in accs:
        p = os.path.join(results_dir, a, "report.json")
        if os.path.exists(p):
            rows.append(json.load(open(p)))
    return rows


def from_supabase(accs):
    url = os.environ.get("SUPABASE_URL", "").rstrip("/")
    key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY", "")
    if not url or not key:
        sys.exit("SUPABASE_URL / SUPABASE_SERVICE_ROLE_KEY not set")
    q = urllib.parse.urlencode(
        {"accession": f"in.({','.join(accs)})", "select": "report"})
    req = urllib.request.Request(
        f"{url}/rest/v1/sentinel_runs?{q}",
        headers={"apikey": key, "Authorization": f"Bearer {key}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return [row["report"] for row in json.loads(r.read().decode())
                if row.get("report")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expect", required=True)
    ap.add_argument("--results-dir")
    ap.add_argument("--supabase", action="store_true")
    ap.add_argument("--json-out")
    args = ap.parse_args()

    exp = load_expect(args.expect)
    manifest = exp["manifest"]
    accs = cohort_accessions(manifest)

    if args.supabase:
        reports = from_supabase(accs)
    elif args.results_dir:
        reports = from_results_dir(accs, args.results_dir)
    else:
        sys.exit("pass --results-dir or --supabase")

    print(f"cohort: {exp.get('cohort')}")
    print(f"manifest: {manifest} ({len(accs)} accessions)")
    print(f"reports found: {len(reports)}")
    if not reports:
        sys.exit("no reports found -- has the cohort been processed?")

    checks = []

    # --- QC pass rate ---
    passed = [r for r in reports if r.get("qc", {}).get("status") == "pass"]
    rate = len(passed) / len(reports)
    floor = exp.get("qc", {}).get("min_pass_rate", 0)
    checks.append((
        "qc pass rate", f"{rate:.2f} ({len(passed)}/{len(reports)})",
        f">= {floor}", rate >= floor))

    if not passed:
        print("no QC-pass samples; cannot evaluate biological expectations")
    else:
        # --- locus presence across QC-pass samples ---
        for locus, min_frac in (exp.get("loci_present") or {}).items():
            present = sum(
                1 for r in passed
                if (r.get("surveillance_loci", {}).get(locus, {}) or {}).get("call")
                == "present")
            frac = present / len(passed)
            checks.append((
                f"locus {locus} present", f"{frac:.2f} ({present}/{len(passed)})",
                f">= {min_frac}", frac >= min_frac))

        # --- SNP spread (the testable form of "homogeneous") ---
        snps = sorted(r.get("variants", {}).get("snps_vs_7pet", 0) for r in passed)
        se = exp.get("snps_vs_7pet") or {}
        if snps:
            med = statistics.median(snps)
            if len(snps) >= 4:
                q = statistics.quantiles(snps, n=4)
                iqr = q[2] - q[0]
            else:
                iqr = max(snps) - min(snps)
            if "max_median" in se:
                checks.append(("snps median", f"{med:g}",
                               f"<= {se['max_median']}", med <= se["max_median"]))
            if "max_iqr" in se:
                checks.append((f"snps spread (IQR, n={len(snps)})", f"{iqr:g}",
                               f"<= {se['max_iqr']}", iqr <= se["max_iqr"]))
            print(f"snps_vs_7pet across QC-pass samples: {snps}")

    width = max(len(c[0]) for c in checks)
    print()
    print(f"{'check'.ljust(width)}  {'observed':>22}  {'expected':>12}  result")
    n_fail = 0
    for name, obs, want, ok in checks:
        n_fail += 0 if ok else 1
        print(f"{name.ljust(width)}  {obs:>22}  {want:>12}  {'PASS' if ok else 'FAIL'}")
    print()
    print(f"{len(checks) - n_fail}/{len(checks)} checks passed")

    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump({
                "cohort": exp.get("cohort"),
                "reports": len(reports),
                "qc_pass": len(passed),
                "checks": [
                    {"check": n, "observed": o, "expected": w, "pass": ok}
                    for n, o, w, ok in checks
                ],
                "failed": n_fail,
            }, f, indent=2)

    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
