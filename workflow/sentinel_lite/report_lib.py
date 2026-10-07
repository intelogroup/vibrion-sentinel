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
    """(status, reasons) where status is 'pass', 'provisional', or 'fail'.

    Tiered gates, each traceable to a published standard (see the FP/FN
    reduction plan):
    - hard 'fail' (excluded from all outputs): extreme contamination
      (species_purity_pct < hard_min_species_purity, default 95 -- >5%
      non-Vibrio reads, FWD-AMR-RefLabCap 2022) or <50% breadth
      (hard_min_called_pct, default 50 -- Kenya 2022-23 outbreak paper).
    - 'provisional' (usable for early alerting, flagged; excluded from
      public-facing summaries): below the strict cutoffs -- mean_depth <
      min_mean_depth (default 30, FWD-AMR-RefLabCap), called_pct <
      min_called_pct (default 90), or species_purity_pct <
      min_species_purity (default 99, PulseNet's 1% secondary-abundance bar).
    - 'pass' otherwise.
    thresholds keys: min_mean_depth, min_called_pct, min_species_purity,
    hard_min_called_pct, hard_min_species_purity. Reasons list every
    breached gate with its threshold, so the report shows its work.
    """
    hard = []
    if species_purity_pct < thresholds["hard_min_species_purity"]:
        hard.append(f"species_purity {species_purity_pct:.2f} < "
                    f"{thresholds['hard_min_species_purity']} (hard fail)")
    if called_pct < thresholds["hard_min_called_pct"]:
        hard.append(f"called_pct {called_pct:.2f} < "
                    f"{thresholds['hard_min_called_pct']} (hard fail)")
    if hard:
        return "fail", hard
    reasons = []
    if mean_depth < thresholds["min_mean_depth"]:
        reasons.append(f"mean_depth {mean_depth:.2f} < {thresholds['min_mean_depth']}")
    if called_pct < thresholds["min_called_pct"]:
        reasons.append(f"called_pct {called_pct:.2f} < {thresholds['min_called_pct']}")
    if species_purity_pct < thresholds["min_species_purity"]:
        reasons.append(f"species_purity {species_purity_pct:.2f} < "
                       f"{thresholds['min_species_purity']}")
    if reasons:
        return "provisional", reasons
    return "pass", []


# ---------------------------------------------------------------------------
# Phase 0 (Nanopore branch): platform/tier config validation, QC defaults,
# assembly-block parsing, tool-version collection.
# ---------------------------------------------------------------------------

VALID_PLATFORMS = ("illumina", "nanopore")
VALID_TIERS = ("lite", "assembly")
VALID_BASECALLERS = ("fast", "hac", "sup")

PIPELINE_VERSION = "0.3.0"  # 5B: alleles + ICE profile + mobile screens

# Tiered QC gates (see evaluate_qc): provisional cutoffs, plus hard-fail floors.
# Provisional defaults: FWD-AMR-RefLabCap (>=30x; >5% off-species contaminated),
# PulseNet PT SOP (<=1% secondary species). Hard floors: Kenya 2022-23 outbreak
# paper (<50% breadth excluded).
PLATFORM_QC_DEFAULTS = {
    "illumina": {"min_site_depth": 10, "qc_min_mean_depth": 30, "qc_min_called_pct": 90,
                 "qc_min_species_purity": 99,
                 "qc_hard_min_called_pct": 50, "qc_hard_min_species_purity": 95},
    "nanopore": {"min_site_depth": 15, "qc_min_mean_depth": 30, "qc_min_called_pct": 85,
                 "qc_min_species_purity": 99,
                 "qc_hard_min_called_pct": 50, "qc_hard_min_species_purity": 95},
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
    qc_min_called_pct, qc_min_species_purity, qc_hard_min_called_pct,
    qc_hard_min_species_purity."""
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


# ---------------------------------------------------------------------------
# Phase 5B: allele calling, ICE profiling, plasmid screening, mobile elements.
# All functions are pure (paths in, plain data out) and unit-testable.
# Global RUO rule: reports state sequence observations with citations, never
# virulence, transmissibility, or clinical predictions.
# ---------------------------------------------------------------------------

import gzip as _gzip

RUO_ALLELE_TEMPLATE = (
    "Sequence-level observation: this sample's {feature} matches the published "
    "{variant} ({citation}). This is an observation, not a prediction of "
    "virulence, transmissibility, or clinical outcome."
)

_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def load_alleles(path):
    """Parse alleles.yaml. Returns the table dict; raises on missing keys so
    a malformed table fails loudly instead of silently dropping alleles."""
    import yaml
    with open(path) as f:
        table = yaml.safe_load(f)
    for key in ("codons", "cds_scans"):
        if key not in table or not table[key]:
            raise ValueError(f"alleles table {path} missing '{key}'")
    return table


def load_mobile_elements(path):
    """Parse mobile_elements.yaml. Returns the element list (possibly empty).
    Raises on a malformed entry so a bad manifest fails loudly instead of
    silently screening nothing."""
    import yaml
    with open(path) as f:
        manifest = yaml.safe_load(f) or {}
    elements = manifest.get("elements") or []
    for el in elements:
        for key in ("id", "ref_fasta", "citation"):
            if not el.get(key):
                raise ValueError(
                    f"mobile_elements manifest {path}: element missing '{key}'")
    return elements


def _parse_vcf_alt_calls(path):
    """{(chrom, pos): (ref, alt)} for haploid ALT calls in a bcftools VCF.

    The lite `call_variants` rule already filters QUAL>=30 and DP>=min site
    depth, so records here passed those floors. GT must be haploid-alt
    ("1", "1/1", "1|1"); anything else is ignored, never guessed."""
    calls = {}
    opener = _gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            chrom, pos, _id, ref, alts = p[0], int(p[1]), p[2], p[3], p[4]
            alt = alts.split(",")[0]
            fmt = p[8].split(":")
            vals = p[9].split(":")
            try:
                gt = vals[fmt.index("GT")]
            except (ValueError, IndexError):
                continue
            if gt in ("1", "1/1", "1|1"):
                calls[(chrom, pos)] = (ref, alt)
    return calls


def _parse_fasta_records(path):
    """{record_name: sequence} (upper-cased). First word of the header."""
    records, name, chunks = {}, None, []
    opener = _gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def _translate_codon(codon):
    """Single-letter AA, or None when the codon is incomplete/ambiguous."""
    codon = (codon or "").upper()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return None
    return _CODON_TABLE[codon]


def call_alleles(vcf_path, consensus_path, allele_table):
    """Allele calls from the filtered VCF + N-masked consensus.

    For each codon entry: per-position called bases (VCF ALT where a haploid
    variant record exists, else the consensus base — which is N below the
    depth floor), assembled in mRNA order (complemented on minus strand),
    translated, and matched against the entry's interpretations.

    Returns {"codons": {id: {...}}, "genes": {gene: {...}}, "cds_scans": [...]}.
    Any site that cannot be established -> "indeterminate", never guessed.
    """
    alt_calls = _parse_vcf_alt_calls(vcf_path)
    consensus = _parse_fasta_records(consensus_path)

    codon_results = {}
    for entry in allele_table["codons"]:
        cid = entry["id"]
        chrom = entry["chrom"]
        strand = entry.get("strand", "+")
        bases = []
        ok = True
        for pos in entry["positions"]:
            b = None
            if (chrom, pos) in alt_calls:
                _ref, alt = alt_calls[(chrom, pos)]
                # Only single-base ALTs resolve a codon site; anything
                # else (indel/MNP) makes the codon indeterminate.
                b = alt if alt and len(alt) == 1 else None
            else:
                seq = consensus.get(chrom)
                if seq is not None and 1 <= pos <= len(seq) and seq[pos - 1] != "N":
                    b = seq[pos - 1]
            if b is None:
                ok = False
                break
            bases.append(b)
        if not ok:
            codon_results[cid] = {
                "call": "indeterminate",
                "observed_aa": None,
                "reason": "site below depth/QUAL floor or unresolvable",
            }
            continue
        codon = "".join(_COMP[b] if strand == "-" else b for b in bases)
        aa = _translate_codon(codon)
        if aa is None:
            codon_results[cid] = {
                "call": "indeterminate",
                "observed_aa": None,
                "reason": "ambiguous codon",
            }
            continue
        interp = (entry.get("interpretations") or {}).get(aa)
        codon_results[cid] = {
            "call": interp["allele"] if interp else "other",
            "observed_aa": aa,
            "observed_codon": codon,
            "ref_aa": entry.get("ref_aa"),
            "citation": interp["citation"] if interp else None,
            "ruo": RUO_ALLELE_TEMPLATE.format(
                feature=f"{entry['gene']} codon {cid.split('_')[-1]}",
                variant=interp["allele"] if interp else f"non-canonical {aa}",
                citation=interp["citation"] if interp else entry.get("derivation", ""),
            ),
        }

    genes = {}
    for gene, rollup in (allele_table.get("gene_rollups") or {}).items():
        parts = [codon_results.get(c, {}).get("call") for c in rollup["codons"]]
        if any(p == "indeterminate" for p in parts):
            call = "indeterminate"
        elif all(p == rollup["all_ctxB7"] for p in parts):
            call = rollup["all_ctxB7"]
        elif all(p == rollup["all_classical-like"] for p in parts):
            call = rollup["all_classical-like"]
        else:
            call = "other"
        genes[gene] = {
            "call": call,
            "codon_calls": {c: codon_results.get(c, {}).get("call") for c in rollup["codons"]},
            "citation": rollup.get("citation"),
        }

    scans = []
    for scan in allele_table.get("cds_scans", []):
        chrom = scan["chrom"]
        start, end = scan["span_1based"]
        variants = []
        for (c, pos), (ref, alt) in sorted(alt_calls.items()):
            if c == chrom and start <= pos <= end:
                kind = "snp" if len(ref) == 1 and len(alt) == 1 else "indel"
                variants.append({
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "kind": kind,
                    "frameshift_candidate": kind == "indel" and abs(len(ref) - len(alt)) % 3 != 0,
                })
        scans.append({
            "id": scan["id"],
            "span": scan["span_1based"],
            "n_variants": len(variants),
            "variants": variants,
            "note": scan.get("note"),
            "citation": scan.get("citation"),
        })

    return {"codons": codon_results, "genes": genes, "cds_scans": scans}


# --- Phase 5B.6: SXT-ICE segment profiling -----------------------------------

def profile_ice_segments(segment_breadths, present_breadth=80.0, absent_breadth=20.0):
    """ICE structural pattern from per-segment coverage breadth.

    segment_breadths: {segment_name: breadth_pct}; expected segments
    SXT_cargo1 (floR..sul2), SXT_backbone1 (traI/traD), SXT_backbone2 (traC),
    SXT_cargo2 (dfrA1). A segment is 'low' below absent_breadth, 'high' at or
    above present_breadth, else 'partial'.

    The XDR AFR13 genotype (Nat Commun Lebanon 2024) layers a ~10 kb deletion
    in ICEVchInd5 that removes cargo block 1 while keeping the backbone and
    dfrA1. This reports the coverage PATTERN consistent with that deletion --
    short-read breadth cannot prove a deletion, so the call is 'consistent
    with', never 'deletion confirmed'.
    """
    segs = {}
    for name, breadth in (segment_breadths or {}).items():
        if breadth is None:
            segs[name] = "indeterminate"
        elif breadth >= present_breadth:
            segs[name] = "high"
        elif breadth < absent_breadth:
            segs[name] = "low"
        else:
            segs[name] = "partial"

    def is_state(name, *states):
        return segs.get(name) in states

    if any(v == "indeterminate" for v in segs.values()):
        pattern = "indeterminate"
    elif (is_state("SXT_cargo1", "high") and is_state("SXT_backbone1", "high")
          and is_state("SXT_backbone2", "high") and is_state("SXT_cargo2", "high")):
        pattern = "intact"
    elif (is_state("SXT_cargo1", "low") and is_state("SXT_backbone1", "high")
          and is_state("SXT_backbone2", "high") and is_state("SXT_cargo2", "high")):
        pattern = "ICEVchInd5-like"
    else:
        pattern = "other"

    notes = {
        "intact": "ICE backbone and cargo blocks all high-breadth.",
        "ICEVchInd5-like": (
            "Coverage pattern consistent with the published ICEVchInd5 "
            "~10 kb deletion (cargo block 1 low, backbone and dfrA1 retained; "
            "Rouard et al. via Nat Commun Lebanon 2024). Short-read breadth "
            "cannot confirm a deletion."
        ),
        "other": "Segment pattern matches neither intact nor ICEVchInd5-like.",
        "indeterminate": "A segment had no usable coverage data.",
    }
    return {
        "pattern": pattern,
        "segments": segs,
        "citation": "Nat Commun 2024 (Lebanon XDR AFR13 ICEVchInd5 deletion)",
        "note": notes[pattern],
    }


# --- Phase 5B.7: IncC YemVchMDRI plasmid gene-set screen ----------------------

YEMVCHMDRI_GENES = ["blaPER-7", "mph(A)", "mph(E)", "msr(E)", "aadA2", "qac", "sul1"]
YEMVCHMDRI_ALIASES = {
    "blaPER7": "blaPER-7", "blaPER_7": "blaPER-7",
    "mphA": "mph(A)", "mphE": "mph(E)", "msrE": "msr(E)",
    "qacE": "qac", "qacEΔ1": "qac", "qacEdelta1": "qac", "qacEΔ1": "qac",
}


def plasmid_screen(amr_genes):
    """YemVchMDRI gene-set co-occurrence on AMRFinderPlus gene symbols.

    The XDR AFR13 plasmid pCNRVC190243 carries blaPER-7, mph(A), mph(E),
    msr(E), aadA2, qac, sul1 (Rouard et al., NEJM Dec 2024). AMRFinderPlus
    detects genes, not plasmids: the call is 'consistent with', never
    'plasmid confirmed' -- short reads cannot prove plasmid carriage.
    """
    normalized = set()
    for g in amr_genes or []:
        g = (g or "").strip()
        normalized.add(YEMVCHMDRI_ALIASES.get(g, g))
    detected = [g for g in YEMVCHMDRI_GENES if g in normalized]
    missing = [g for g in YEMVCHMDRI_GENES if g not in normalized]
    if not detected:
        call = "absent"
    elif not missing:
        call = "consistent"
    else:
        call = "partial"
    return {
        "call": call,
        "genes_detected": detected,
        "genes_missing": missing,
        "citation": "Rouard et al., NEJM Dec 2024 (XDR AFR13, IncC pCNRVC190243)",
        "note": (
            "Gene-set co-detection only. Consistent with the published IncC "
            "plasmid gene set; plasmid carriage cannot be proven from "
            "short-read gene detection. " + AMR_NOTE
        ),
    }


# --- Phase 5B.3/5B.4: mobile-element screens (config-gated) --------------------

def summarize_mobile_elements(element_rows, present_breadth=80.0, absent_breadth=20.0):
    """{element: {call, breadth_pct, mean_depth}} from per-element coverage
    rows: [{element, breadth_pct, mean_depth}]. Same breadth semantics as
    call_loci: presence of reference-like sequence only. Elements are mobile
    by nature -- call presence per element, never genome completeness."""
    out = {}
    for row in element_rows or []:
        name = row.get("element")
        breadth = row.get("breadth_pct")
        if not name or breadth is None:
            continue
        if breadth >= present_breadth:
            call = "present"
        elif breadth < absent_breadth:
            call = "absent"
        else:
            call = "partial"
        out[name] = {
            "call": call,
            "breadth_pct": round(float(breadth), 2),
            "mean_depth": round(float(row.get("mean_depth") or 0), 2),
            "citation": row.get("citation"),
        }
    return out
