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


# ---------------------------------------------------------------------------
# Phase 0 (Nanopore branch): platform/tier config validation, QC defaults,
# assembly-block parsing, tool-version collection.
# ---------------------------------------------------------------------------

VALID_PLATFORMS = ("illumina", "nanopore")
VALID_TIERS = ("lite", "assembly")
VALID_BASECALLERS = ("fast", "hac", "sup")

PIPELINE_VERSION = "0.2.0"

# Nanopore: homopolymer noise means stricter depth, looser per-base
# expectations (ends of long reads drop out). Illumina values are the
# pre-Phase-0 defaults, unchanged.
PLATFORM_QC_DEFAULTS = {
    "illumina": {"min_site_depth": 10, "qc_min_mean_depth": 20, "qc_min_called_pct": 90},
    "nanopore": {"min_site_depth": 15, "qc_min_mean_depth": 30, "qc_min_called_pct": 85},
}

# Assembly QC gate (fail loudly): a fragmented assembly must never produce
# lineage or AMR calls.
ASM_MIN_LENGTH_BP = 3800000
ASM_MAX_LENGTH_BP = 4400000
ASM_MAX_CONTIGS = 50
ASM_MIN_N50_BP = 200000

# Load-bearing RUO disclaimer: must appear in every report carrying AMR calls
# and in every downstream rendering of them.
AMR_NOTE = "Gene detected/not detected only. Not a susceptibility prediction."

# Medaka consensus models for R10.4.1 400bps chemistry, keyed by basecaller.
# Confirm against `medaka tools list` for the installed Medaka version; a
# mismatched model is a silent accuracy killer, which is why basecaller_model
# is required (fail-closed) rather than guessed.
MEDAKA_MODELS = {
    "fast": "r1041_e82_400bps_fast_g615",
    "hac": "r1041_e82_400bps_hac_g615",
    "sup": "r1041_e82_400bps_sup_g615",
}


def validate_platform_config(config):
    """(platform, basecaller_model, tier) from a config dict, or raise
    ValueError. Fail-closed: platform=nanopore without basecaller_model
    refuses to run, because Medaka models are basecaller-matched and a
    wrong model silently degrades accuracy."""
    platform = config.get("platform", "illumina")
    if platform not in VALID_PLATFORMS:
        raise ValueError(
            f"unknown platform {platform!r}; expected one of {VALID_PLATFORMS}"
        )
    tier = config.get("tier", "lite")
    if tier not in VALID_TIERS:
        raise ValueError(f"unknown tier {tier!r}; expected one of {VALID_TIERS}")
    basecaller = config.get("basecaller_model")
    if platform == "nanopore":
        if not basecaller:
            raise ValueError(
                "platform=nanopore requires basecaller_model "
                f"(one of {VALID_BASECALLERS}); refusing to run"
            )
        if basecaller not in VALID_BASECALLERS:
            raise ValueError(
                f"unknown basecaller_model {basecaller!r}; "
                f"expected one of {VALID_BASECALLERS}"
            )
    if tier == "assembly" and platform != "nanopore":
        raise ValueError("tier=assembly is nanopore-only in Phase 0")
    return platform, basecaller, tier


def platform_qc_thresholds(config, platform):
    """QC thresholds for a platform: platform defaults, with explicit
    config values winning. Keys: min_site_depth, qc_min_mean_depth,
    qc_min_called_pct."""
    defaults = PLATFORM_QC_DEFAULTS[platform]
    return {k: config.get(k, v) for k, v in defaults.items()}


def medaka_model_for(basecaller):
    """Medaka consensus model name for a basecaller; KeyError on unknown."""
    return MEDAKA_MODELS[basecaller]


def flye_preset_for(basecaller):
    """Flye --nano-* preset: hac/sup reads are high-quality, fast are not."""
    if basecaller in ("hac", "sup"):
        return "--nano-hq"
    return "--nano-raw"


def parse_assembly_qc(path):
    """{"length_bp", "contigs", "n50_bp"} from an assembly_qc.json file
    written by the assembly_qc rule. Missing file -> zeros (gate fails)."""
    import json as _json

    try:
        with open(path) as f:
            d = _json.load(f)
        return {
            "length_bp": int(d.get("length_bp", 0)),
            "contigs": int(d.get("contigs", 0)),
            "n50_bp": int(d.get("n50_bp", 0)),
        }
    except (OSError, ValueError):
        return {"length_bp": 0, "contigs": 0, "n50_bp": 0}


def assembly_qc_gate(stats):
    """(status, reasons) for the assembly QC gate. A failing assembly must
    never produce lineage or AMR calls -- callers must check status first."""
    reasons = []
    n = stats["length_bp"]
    if not (ASM_MIN_LENGTH_BP <= n <= ASM_MAX_LENGTH_BP):
        reasons.append(
            f"assembly length {n} outside {ASM_MIN_LENGTH_BP}-{ASM_MAX_LENGTH_BP}"
        )
    if stats["contigs"] > ASM_MAX_CONTIGS:
        reasons.append(f"contig count {stats['contigs']} > {ASM_MAX_CONTIGS}")
    if stats["n50_bp"] < ASM_MIN_N50_BP:
        reasons.append(f"N50 {stats['n50_bp']} < {ASM_MIN_N50_BP}")
    return ("fail" if reasons else "pass"), reasons


def parse_mlst(path):
    """ST string from tseemann/mlst TSV output (file, scheme, ST, alleles...).
    '-' or missing -> None. Never raises on malformed input."""
    try:
        with open(path) as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[2] not in ("-", ""):
                    return parts[2]
                return None
    except OSError:
        return None
    return None


def parse_vibecheck(path):
    """Lineage string (e.g. 'T12') from vibecheck output. Tolerant: tries a
    JSON {"lineage": ...} first, then falls back to the first T\\d+ token in
    text. None when nothing is found. The exact vibecheck CLI/output format
    must be confirmed against https://github.com/cholgen/vibecheck."""
    import json as _json

    try:
        with open(path) as f:
            text = f.read()
    except OSError:
        return None
    try:
        d = _json.loads(text)
        lin = d.get("lineage")
        if lin:
            return str(lin)
    except ValueError:
        pass
    m = re.search(r"\b(T\d{1,2})\b", text)
    return m.group(1) if m else None


def parse_amrfinder(path):
    """Sorted, deduplicated AMR gene symbols from AMRFinderPlus TSV output
    ('Gene symbol' column). Empty list on missing/malformed input. Gene
    presence only -- never a susceptibility prediction."""
    genes = set()
    try:
        with open(path) as f:
            header = None
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = parts
                    try:
                        idx = header.index("Gene symbol")
                    except ValueError:
                        return []
                    continue
                if idx < len(parts) and parts[idx].strip():
                    genes.add(parts[idx].strip())
    except OSError:
        return []
    return sorted(genes)


def build_assembly_block(stats, mlst_st, lineage, amr_genes, gate_status):
    """The report.json 'assembly' block. lineage/amr_genes are None/[]
    unless gate_status is 'pass' -- a failing assembly never yields calls."""
    block = {
        "length_bp": stats["length_bp"],
        "contigs": stats["contigs"],
        "n50_bp": stats["n50_bp"],
        "qc_gate": gate_status,
        "mlst_st": mlst_st if gate_status == "pass" else None,
        "vibecheck_lineage": lineage if gate_status == "pass" else None,
        "amr_genes": sorted(set(amr_genes)) if gate_status == "pass" else [],
        "amr_note": AMR_NOTE,
    }
    return block


def collect_tool_versions(specs, runner=None):
    """{tool: version_string} for [(name, argv), ...]. runner defaults to
    subprocess.run; inject a fake in tests. A tool that errors or is absent
    records 'unknown' instead of raising -- version collection must never
    fail a run."""
    import subprocess as _sp

    if runner is None:
        def runner(argv):
            return _sp.run(argv, capture_output=True, text=True, timeout=60)

    versions = {}
    for name, argv in specs:
        try:
            r = runner(argv)
            out = (r.stdout or "") + (r.stderr or "")
            first = out.strip().splitlines()[0] if out.strip() else ""
            versions[name] = first[:200] if first else "unknown"
        except Exception:
            versions[name] = "unknown"
    return versions


def amrfinder_db_version(raw_version_output):
    """Database version from `amrfinder --version` output, which prints both
    software and database versions. 'unknown' when not found."""
    m = re.search(r"[Dd]atabase version:?\s*(\S+)", raw_version_output or "")
    return m.group(1) if m else "unknown"
