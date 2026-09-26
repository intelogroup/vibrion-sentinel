"""Pure logic for workflow/cohort_phylo/Snakefile's stack_alignment and
snp_distances rules, extracted so it can be unit tested without real
consensus FASTA fixtures wired through Snakemake and R2.
"""

import gzip
import os


def load_mask(mask_bed_path):
    """{contig: [(start0, end), ...]} from a BED file, or {} if the path
    doesn't exist (an unset/missing mask is not an error -- it just masks
    nothing)."""
    masked = {}
    if not os.path.exists(mask_bed_path):
        return masked
    with open(mask_bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            masked.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    return masked


def read_consensus(consensus_dir, acc):
    """{contig_name: sequence} for one accession's consensus FASTA,
    trying .fasta.gz then .fasta in consensus_dir. None if neither
    exists (caller decides whether that's an exclusion reason)."""
    for pat, is_gz in ((f"{consensus_dir}/{acc}.fasta.gz", True),
                       (f"{consensus_dir}/{acc}.fasta", False)):
        if os.path.exists(pat):
            op = gzip.open if is_gz else open
            seqs, name = {}, None
            with op(pat, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        name = line[1:].split()[0]
                        seqs[name] = []
                    elif name:
                        seqs[name].append(line.strip())
            return {k: "".join(v) for k, v in seqs.items()}
    return None


def load_and_filter_samples(accessions, consensus_dir, min_sample_called):
    """(loaded, excluded) where loaded is {acc: (contig_dict, called_frac)}
    for samples that both have a consensus file and clear
    min_sample_called, and excluded is [(acc, reason), ...] for the rest.

    A sample missing entirely from consensus_dir (not yet processed, or
    QC-fail and therefore never archived) is excluded with a distinct
    reason from a low-completeness one, so alignment_stats.json /
    excluded_samples.tsv can tell the two apart.
    """
    loaded, excluded = {}, []
    for acc in accessions:
        c = read_consensus(consensus_dir, acc)
        if not c:
            excluded.append((acc, "no consensus available"))
            continue
        seq = "".join(c[k] for k in sorted(c))
        called = 1.0 - (seq.upper().count("N") / len(seq)) if seq else 0.0
        if called < min_sample_called:
            excluded.append((acc, f"called {called:.3f} < {min_sample_called}"))
            continue
        loaded[acc] = (c, called)
    return loaded, excluded


def validate_coordinates(loaded):
    """(contigs, lengths) if every sample's contigs are the same set at
    the same lengths -- the invariant that makes stacking sequences
    positionally valid without an MSA step. Raises RuntimeError, naming
    the offending accession, the moment that invariant breaks (this is
    what caught a real bug this session: bcftools consensus applying
    indels made chr1 length differ by 1-4bp between samples)."""
    if not loaded:
        raise RuntimeError("no samples to validate coordinates against")
    contigs = sorted(next(iter(loaded.values()))[0])
    lengths = {ct: len(next(iter(loaded.values()))[0][ct]) for ct in contigs}
    for acc, (c, _called) in loaded.items():
        for ct in contigs:
            if len(c.get(ct, "")) != lengths[ct]:
                raise RuntimeError(
                    f"{acc} contig {ct} length {len(c.get(ct, ''))} != "
                    f"{lengths[ct]}: consensus sequences are not on identical "
                    "reference coordinates, so stacking them is invalid")
    return contigs, lengths


def build_skip_mask(masked, contigs, lengths):
    """({contig: bytearray marking masked positions}, n_masked_positions)
    for the given exclusion-BED intervals, contigs and per-contig
    lengths. Intervals on a contig not in `contigs`, or extending past
    `lengths`, are silently clipped rather than raising -- the mask file
    is allowed to describe a superset of what this cohort's reference
    actually has."""
    skip = {ct: bytearray(lengths[ct]) for ct in contigs}
    n_masked = 0
    for ct, ivs in masked.items():
        if ct not in skip:
            continue
        for s, e in ivs:
            for i in range(max(0, s), min(e, lengths[ct])):
                if not skip[ct][i]:
                    skip[ct][i] = 1
                    n_masked += 1
    return skip, n_masked


def build_alignment_columns(loaded, contigs, lengths, skip, max_site_missing):
    """(accs, cols, total_sites, n_missing_dropped, n_invariant).

    accs: sorted accession list (column order in each returned column).
    cols: one list[base] per kept site, in genome order across contigs --
        variable (>=2 distinct real bases among samples) and with no more
        than max_site_missing fraction of non-ACGT calls.
    The three counts partition every non-masked reference site into
    exactly one bucket (considered = dropped_missing + invariant + kept),
    which is what lets alignment_stats.json be checked for internal
    consistency, not just plausibility.
    """
    accs = sorted(loaded)
    cols = []
    total_sites = n_missing_dropped = n_invariant = 0
    for ct in contigs:
        seqs = [loaded[a][0][ct].upper() for a in accs]
        for i in range(lengths[ct]):
            if skip[ct][i]:
                continue
            total_sites += 1
            col = [s[i] for s in seqs]
            miss = sum(1 for b in col if b not in "ACGT")
            if miss / len(col) > max_site_missing:
                n_missing_dropped += 1
                continue
            if len({b for b in col if b in "ACGT"}) < 2:
                n_invariant += 1
                continue
            cols.append(col)
    return accs, cols, total_sites, n_missing_dropped, n_invariant


def format_alignment_fasta(accs, cols):
    """FASTA text: one record per accession, its sequence being the
    kept columns read across in order (column j's base at row i is
    accs[i]'s allele at that site)."""
    lines = []
    for j, a in enumerate(accs):
        lines.append(f">{a}")
        lines.append("".join(c[j] for c in cols))
    return "\n".join(lines) + ("\n" if lines else "")


def format_excluded_tsv(excluded):
    lines = ["accession\treason"]
    lines += [f"{acc}\t{why}" for acc, why in excluded]
    return "\n".join(lines) + "\n"


def pairwise_snp_distances(seqs):
    """(sorted_accessions, {(a, b): distance}) -- the count of aligned
    positions where both sequences have a real base (A/C/G/T) and they
    differ. Positions where either side is N/gap/anything else are
    excluded from both the denominator and the count, not treated as a
    mismatch."""
    accs = sorted(seqs)
    dist = {}
    for a in accs:
        for b in accs:
            d = sum(
                1 for x, y in zip(seqs[a], seqs[b])
                if x in "ACGT" and y in "ACGT" and x != y
            )
            dist[(a, b)] = d
    return accs, dist


def format_distance_matrix(accs, dist):
    lines = ["\t" + "\t".join(accs)]
    for a in accs:
        lines.append("\t".join([a] + [str(dist[(a, b)]) for b in accs]))
    return "\n".join(lines) + "\n"


def parse_fasta(path):
    """{record_name: sequence} from a (small, in-memory-sized) FASTA
    file -- used to read core_alignment.fasta back in for snp_distances,
    since it's the alignment's own output rather than a per-sample
    consensus (so it doesn't need read_consensus's gzip handling)."""
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = ""
            elif name is not None:
                seqs[name] += line.strip()
    return seqs
