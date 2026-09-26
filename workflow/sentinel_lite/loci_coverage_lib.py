"""Pure logic for the sentinel_lite `loci_coverage` rule: converting BED
intervals to samtools region strings. Extracted because this is exactly
the class of off-by-one bug this pipeline has hit before (BED is 0-based
half-open; samtools coverage/region syntax is 1-based inclusive).
"""


def iter_bed_regions(bed_path):
    """Yield (bed_fields, region_string) for each data line in a BED
    file, skipping blank lines and comments. bed_fields is the raw
    tab-split line (so callers can still get the locus name at index 3);
    region_string is what `samtools coverage -r` expects.
    """
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            region = f"{p[0]}:{int(p[1]) + 1}-{p[2]}"
            yield p, region
