"""Pure logic for the sentinel_lite `report` rule, extracted out of the
Snakefile's `run:` block so it can be unit tested. No file I/O paths are
hardcoded here; every function takes an open-able path or string and
returns plain data. The Snakefile calls these and only does the
Snakemake-specific wiring (input/output paths, config lookups).
"""

import re


def parse_kraken_report(path):
    """Return (cholerae_reads, unclassified_reads, root_reads, total_kraken).

    Kraken2 report columns: pct, clade_reads, direct_reads, rank_code,
    taxid, name. total_kraken = unclassified (taxid "0") + root (taxid "1")
    clade counts -- every read Kraken processed. A prior version used
    max(clade_reads) over all rows, which is not the total and made the
    resulting fraction meaningless regardless of database size.
    """
    cholerae_reads = unclassified_reads = root_reads = 0
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            n_reads, taxid = parts[1], parts[4]
            try:
                n_reads = int(n_reads)
            except ValueError:
                continue
            if taxid == "666":
                cholerae_reads = n_reads
            elif taxid == "0":
                unclassified_reads = n_reads
            elif taxid == "1":
                root_reads = n_reads
    total_kraken = unclassified_reads + root_reads
    return cholerae_reads, unclassified_reads, root_reads, total_kraken


def species_purity(cholerae_reads, total_kraken):
    """Fraction (0-1) and percent (0-100) of Kraken-processed reads
    classified as V. cholerae. 0.0 when total_kraken is 0 (nothing to
    divide by), not an error -- callers decide whether that itself is a
    QC failure."""
    fraction = cholerae_reads / total_kraken if total_kraken else 0.0
    return fraction, 100.0 * fraction


def genome_wide_coverage(path):
    """Length-weighted mean depth and breadth (%) across all contigs from
    a `samtools coverage` table (columns: rname startpos endpos numreads
    covbases coverage meandepth ...). Weighting by contig length is what
    makes this correct across chr1+chr2 rather than reporting chr1 alone."""
    tot_len = tot_cov = 0
    depth_sum = 0.0
    with open(path) as f:
        for line in f:
            if line.startswith("#rname"):
                continue
            p = line.strip().split("\t")
            if len(p) < 7:
                continue
            try:
                length = int(p[2]) - int(p[1]) + 1
                tot_len += length
                tot_cov += int(p[4])
                depth_sum += float(p[6]) * length
            except ValueError:
                continue
    mean_depth = depth_sum / tot_len if tot_len else 0.0
    breadth = 100.0 * tot_cov / tot_len if tot_len else 0.0
    return mean_depth, breadth


def parse_flagstat_mapped(text):
    """Mapped-read count from `samtools flagstat` text, e.g.
    '12345 + 0 mapped (99.00% : N/A)'. 0 if the line isn't found."""
    m = re.search(r"(\d+) \+ \d+ mapped", text)
    return int(m.group(1)) if m else 0


def parse_loci_bed(path):
    """{(chrom, start_1based, end) : locus_name} from a BED file. BED is
    0-based half-open; samtools coverage -r's startpos is 1-based, so the
    key uses start+1 to match coverage rows directly."""
    bed_labels = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.strip().split("\t")
            bed_labels[(p[0], int(p[1]) + 1, int(p[2]))] = p[3]
    return bed_labels


def call_loci(path, bed_labels, present_breadth, absent_breadth):
    """{locus_name: {call, present, mean_depth, breadth_pct}} from a
    per-locus `samtools coverage -r` table, joined to bed_labels by
    coordinates (samtools coverage prints rname/startpos/endpos, not the
    BED name).

    Called on breadth (% of the locus with >=1 read), not mean depth:
    AT-rich loci lose Illumina depth disproportionately (observed: toxT,
    ~28% GC vs ~47.5% genome, drops to ~4x against ~50-78x flanks with no
    soft-clips or split reads -- consistent with GC bias, not a deletion),
    so a depth threshold called a present gene absent in a real sample.
    """
    loci_calls = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#rname") or not line.strip() or "failed" in line:
                continue
            p = line.strip().split("\t")
            if len(p) < 7:
                continue
            try:
                key = (p[0], int(p[1]), int(p[2]))
                depth, cov = float(p[6]), float(p[5])
            except ValueError:
                continue
            name = bed_labels.get(key, f"{p[0]}:{p[1]}-{p[2]}")
            if cov >= present_breadth:
                call = "present"
            elif cov < absent_breadth:
                call = "absent"
            else:
                call = "partial"
            loci_calls[name] = {
                "call": call,
                "present": call == "present",
                "mean_depth": round(depth, 2),
                "breadth_pct": round(cov, 2),
            }
    return loci_calls


def consensus_stats(path):
    """(length, n_count, called_pct) for a (possibly multi-record) FASTA.
    called_pct is 0.0 for a zero-length file, not a ZeroDivisionError."""
    cons_len = cons_n = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip()
                cons_len += len(seq)
                cons_n += seq.upper().count("N")
    called_pct = 100.0 * (cons_len - cons_n) / cons_len if cons_len else 0.0
    return cons_len, cons_n, called_pct


def read_snp_count(path):
    """int from a snp_count.txt file; 0 (not an exception) if the file
    holds something non-numeric -- call_variants' shell rule already
    falls back to writing 0 on its own failure, so this mirrors that
    rather than crashing the report rule on top of an upstream failure."""
    try:
        with open(path) as f:
            return int(f.read().strip())
    except ValueError:
        return 0


def evaluate_qc(mean_depth, called_pct, species_purity_pct, thresholds):
    """(status, reasons) where status is 'pass' or 'fail'. thresholds is
    {min_mean_depth, min_called_pct, min_species_purity}. Each failing
    check appends one human-readable reason; qc.reasons in report.json is
    this list."""
    reasons = []
    if mean_depth < thresholds["min_mean_depth"]:
        reasons.append(f"mean_depth {mean_depth:.2f} < {thresholds['min_mean_depth']}")
    if called_pct < thresholds["min_called_pct"]:
        reasons.append(f"called_pct {called_pct:.2f} < {thresholds['min_called_pct']}")
    if species_purity_pct < thresholds["min_species_purity"]:
        reasons.append(
            f"species_purity {species_purity_pct:.2f} < {thresholds['min_species_purity']}"
        )
    return ("fail" if reasons else "pass"), reasons
