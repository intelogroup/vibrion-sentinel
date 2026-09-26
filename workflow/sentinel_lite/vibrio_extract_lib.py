"""Pure logic for the sentinel_lite `extract_vibrio` rule, extracted out
of the Snakefile's `run:` block so it can be unit tested without needing
real (gzipped, multi-line) FASTQ fixtures wired through Snakemake.
"""

import gzip


def build_keep_set(kraken_out_path, cholerae_taxid="666"):
    """Read ids (without a trailing /1 or /2) to keep for mapping, from a
    Kraken2 `--output` file (columns: C/U status, read id, taxid, ...).

    Kraken is the coarse filter; BWA mapping to the cholerae reference
    downstream is the confirmatory step, so unclassified ("U") reads are
    kept too, to rescue divergent strains rather than discard them here.

    Returns (keep_set, n_classified_lines, n_kept).
    """
    keep = set()
    n_class = n_keep = 0
    with open(kraken_out_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            status, rid, taxid = parts[0], parts[1], parts[2]
            n_class += 1
            base = rid.rsplit("/", 1)[0] if "/" in rid else rid
            if status == "C" and taxid == cholerae_taxid:
                keep.add(base)
                n_keep += 1
            elif status == "U":
                keep.add(base)
                n_keep += 1
    return keep, n_class, n_keep


def filter_fastq(in_path, out_path, keep):
    """Write only the 4-line FASTQ records whose read id (header before
    any whitespace, with a trailing /1 or /2 stripped) is in `keep`.

    Returns (n_in, n_out). A truncated/malformed input (a header with no
    full following 3 lines) is not explicitly guarded here -- readline()
    returning "" for seq/plus/qual would write a short record. Fixed-width
    FASTQ readers commonly make this same assumption; see
    test_extract_vibrio_lib.py for what happens on a truncated file."""
    n_in = n_out = 0
    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            seq, plus, qual = fin.readline(), fin.readline(), fin.readline()
            n_in += 1
            base = h[1:].split()[0].rsplit("/", 1)[0]
            if base in keep:
                fout.write(h + seq + plus + qual)
                n_out += 1
    return n_in, n_out
