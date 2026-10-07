#!/usr/bin/env python3
"""True pairwise SNP distances between reference-mapped consensus genomes.

Phase 4 of the recurrence/cluster watcher. The lite pipeline's consensus is
emitted in 2010EL-1786 reference coordinates, so two consensuses align
position-for-position with no aligner: the true pairwise SNP distance is the
count of sites where both samples have a called base (ACGT) and the bases
differ. Sites where either sample is N (below min_site_depth) carry no
information and are excluded -- pairwise-complete deletion, the standard
choice for this kind of comparison.

Relationship to the R2 proxy: |snps_vs_7pet(a) - snps_vs_7pet(b)| <= 5 is a
triangle-inequality screen -- necessary but not sufficient for "a and b are
<= 5 true SNPs apart". The proxy never misses a true cluster, but it
over-includes (e.g. the first two 2022 resurgence samples: proxy 3, true 6).
This module confirms screen hits with the true distance; it does not replace
the screen, and confirmation annotates the alert rather than gating it, so
the validated R2 behavior is preserved while the true measurement is visible.
"""

CALLED = frozenset("ACGT")


def load_consensus(path):
    """Read a consensus FASTA (.fasta or .fasta.gz), return the sequence."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    seq = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def pairwise_snp_distance(seq_a, seq_b):
    """(distance, n_compared) between two equal-length consensus strings.

    Counts positions where both bases are called (ACGT) and differ.
    Raises ValueError on length mismatch -- consensuses from the same
    reference must align; a mismatch is data corruption, not a distance.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError(
            f"consensus length mismatch: {len(seq_a)} vs {len(seq_b)}")
    dist = 0
    compared = 0
    for x, y in zip(seq_a, seq_b):
        if x in CALLED and y in CALLED:
            compared += 1
            if x != y:
                dist += 1
    return dist, compared


def distance_matrix(seqs):
    """{(acc_a, acc_b): (dist, n_compared)} for sorted accession pairs."""
    accs = sorted(seqs)
    out = {}
    for i in range(len(accs)):
        for j in range(i + 1, len(accs)):
            a, b = accs[i], accs[j]
            out[(a, b)] = pairwise_snp_distance(seqs[a], seqs[b])
    return out


def confirm_cluster(accessions, seqs, snp_tol=5):
    """Annotate a proxy-screened cluster with true pairwise distances.

    accessions: the cluster's accession list (from detect_clusters).
    seqs: {accession: consensus string} for whichever members have a
        consensus available (QC-pass only in production).
    Returns a dict with:
      pairwise_snps: {(a, b): dist} for pairs with both consensuses
      pairwise_compared: {(a, b): n_compared}
      pairwise_diameter: max true distance over confirmed pairs (None if none)
      pairwise_confirmed: True if every confirmed pair is within snp_tol,
        False if any exceeds it, None if no pair could be confirmed.
      missing: accessions with no consensus available.
    """
    members = sorted(accessions)
    have = {a: seqs[a] for a in members if a in seqs}
    missing = [a for a in members if a not in seqs]
    mat = distance_matrix(have)
    pairwise_snps = {k: v[0] for k, v in mat.items()}
    pairwise_compared = {k: v[1] for k, v in mat.items()}
    if not mat:
        confirmed = None
        diameter = None
    else:
        confirmed = all(d <= snp_tol for d in pairwise_snps.values())
        diameter = max(pairwise_snps.values())
    return {
        "pairwise_snps": pairwise_snps,
        "pairwise_compared": pairwise_compared,
        "pairwise_diameter": diameter,
        "pairwise_confirmed": confirmed,
        "missing_consensus": missing,
    }
