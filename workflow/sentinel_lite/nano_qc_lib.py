"""Nanopore QC helpers for the sentinel_lite pipeline (Phase 0).

NanoPlot writes NanoStats.txt (a "key:\\tvalue" text file); chopper does the
filtering. The report rule expects a fastp-shaped JSON summary, so
write_qc_summary() adapts: {"summary": {"before_filtering": {"total_reads": N},
"after_filtering": {"total_reads": M}}}. Pure functions, unit-tested.
"""

import gzip
import json


def parse_nanostats(path):
    """{key: value} from a NanoPlot NanoStats.txt file. Values stay strings;
    callers convert. Missing file -> {} (never raises)."""
    stats = {}
    try:
        with open(path) as f:
            for line in f:
                if ":\t" not in line:
                    continue
                key, _, value = line.partition(":\t")
                stats[key.strip()] = value.strip()
    except OSError:
        pass
    return stats


def count_fastq_reads_gz(path):
    """Read count of a gzipped FASTQ (4 lines per read). 0 on any error --
    a missing cleaned file is a QC signal, not a crash."""
    n = 0
    try:
        with gzip.open(path, "rt") as f:
            for _ in f:
                n += 1
    except OSError:
        return 0
    return n // 4


def write_qc_summary(nanostats_path, cleaned_fastq_path, out_path):
    """Write the fastp-shaped QC summary JSON.

    before_filtering.total_reads comes from NanoStats 'number_of_reads'
    (raw input); after_filtering.total_reads is counted from the cleaned
    FASTQ. Returns the dict written.
    """
    stats = parse_nanostats(nanostats_path)
    try:
        raw_reads = int(float(stats.get("number_of_reads", 0)))
    except (ValueError, TypeError):
        raw_reads = 0
    kept_reads = count_fastq_reads_gz(cleaned_fastq_path)
    summary = {
        "summary": {
            "before_filtering": {"total_reads": raw_reads},
            "after_filtering": {"total_reads": kept_reads},
        }
    }
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    return summary
