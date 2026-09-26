# Merge consecutive low-depth positions from `samtools depth -aa` output
# (columns: chrom, pos [1-based], depth) into a BED of masked intervals
# (0-based, end-exclusive). Invoked with -v m=<MIN_SITE_DEPTH>.
#
# Extracted from workflow/sentinel_lite/Snakefile's depth_mask rule so it
# can be tested independent of Snakemake; see
# tests/sentinel_lite/test_depth_mask_awk.py.
BEGIN { OFS = "\t" }
$3 < m {
    if ($1 == c && $2 == e + 1) { e = $2 }
    else { if (c != "") print c, s - 1, e; c = $1; s = $2; e = $2 }
    next
}
{ if (c != "") print c, s - 1, e; c = "" }
END { if (c != "") print c, s - 1, e }
