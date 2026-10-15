"""Regression: HL2601, the novel-ST (ST2040) non-O1/non-O139 strain.

HL2601 (Xiamen bloodstream isolate, Front. Microbiol. 2026; assembly
GCF_059628145.1) is the case that exposed a classification gap: its MLST
profile matched nothing in PubMLST, and Sentinel reduced that to
mlst_st=null -- indistinguishable from a failed typing run, i.e. the
strain was silently dropped as a negative.

This test feeds the REAL deposited genome (tests/fixtures/hl2601/,
content-verified at fixture entry -- see the fixture README for the
verification record) through the report layer and pins the outcome:

  V. cholerae, unassigned ST; non-toxigenic; ctxAB/tcpA absent,
  rtxA/hlyA present, T6SS present.

Loci calls here are derived from the assembly sequence itself (distinct
31-mer containment of each 2010EL-1786 panel gene: >=0.5 present,
<=0.05 absent, else partial), standing in for the pipeline's mapping-
breadth calls, which need the full kraken/bwa stack. The functions
under test -- assembly_qc_gate, parse_mlst_result,
build_assembly_block, virulence_markers -- are the production ones.
"""
import gzip
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import report_lib  # noqa: E402

FIXTURES = Path(__file__).parent.parent / "fixtures" / "hl2601"
K = 31
_BASE2 = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMP2 = {"A": 3, "C": 2, "G": 1, "T": 0}


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    records, name, chunks = {}, None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def canon_kmers(seq):
    """Distinct canonical (strand-free) K-mers as 2-bit-encoded ints."""
    mask = (1 << (2 * K)) - 1
    top = 2 * (K - 1)
    fwd = rev = run = 0
    out = set()
    for ch in seq:
        b = _BASE2.get(ch)
        if b is None:
            fwd = rev = run = 0
            continue
        fwd = ((fwd << 2) | b) & mask
        rev = (rev >> 2) | (_COMP2[ch] << top)
        run += 1
        if run >= K:
            out.add(fwd if fwd < rev else rev)
    return out


def n50(lengths):
    total = sum(lengths)
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc >= total / 2:
            return length
    return 0


class TestHL2601(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = read_fasta(FIXTURES / "genomic.fasta.gz")
        cls.gene_refs = read_fasta(FIXTURES / "marker_refs.fasta")
        genome_kmers = set()
        for seq in cls.genome.values():
            genome_kmers |= canon_kmers(seq)
        cls.containment = {}
        for gene, seq in cls.gene_refs.items():
            gk = canon_kmers(seq)
            cls.containment[gene] = (
                len(gk & genome_kmers) / len(gk) if gk else 0.0)

    def _stats(self):
        lengths = [len(s) for s in self.genome.values()]
        return {"length_bp": sum(lengths), "contigs": len(lengths),
                "n50_bp": n50(lengths)}

    def _loci_calls(self):
        calls = {}
        for gene, frac in self.containment.items():
            call = ("present" if frac >= 0.5
                    else "absent" if frac <= 0.05 else "partial")
            calls[gene] = {"call": call, "present": call == "present",
                           "mean_depth": 0.0, "breadth_pct": 0.0}
        return calls

    def test_assembly_passes_qc_gate(self):
        # The deposited assembly itself: V. cholerae-scale and clean
        # enough that the assembly tier would type it, not suppress it.
        stats = self._stats()
        self.assertEqual(stats["length_bp"], 4297354)
        status, reasons = report_lib.assembly_qc_gate(stats)
        self.assertEqual(status, "pass", reasons)

    def test_classification_is_unassigned_st(self):
        mlst = report_lib.parse_mlst_result(FIXTURES / "mlst.tsv")
        self.assertIsNone(mlst["st"])
        self.assertEqual(mlst["st_status"], "unassigned")
        block = report_lib.build_assembly_block(
            self._stats(), mlst["st"], None, [], "pass",
            mlst_st_status=mlst["st_status"])
        self.assertEqual(block["classification"],
                         "V. cholerae, unassigned ST")
        self.assertEqual(block["st_status"], "unassigned")

    def test_marker_line_is_non_toxigenic_novc(self):
        markers = report_lib.virulence_markers(self._loci_calls())
        self.assertEqual(markers["ctxAB"]["call"], "absent")
        self.assertIs(markers["toxigenic"], False)
        self.assertEqual(markers["tcpA"]["call"], "absent")
        self.assertEqual(markers["rtxA"]["call"], "present")
        self.assertEqual(markers["hlyA"]["call"], "present")
        # T6SS: HL2601's hcp is a diverged paralog of the reference
        # hcp-1 (exact-31-mer containment 0.227 -- far above the 0.000
        # of the truly absent genes, below the 0.5 'present' bar), and
        # PGAP annotates hcp-2 + the vas cluster in this assembly
        # (fixture README). The marker must not call it absent.
        self.assertGreater(self.containment["hcp"], 0.1)
        self.assertNotEqual(markers["t6ss"]["call"], "absent")


if __name__ == "__main__":
    unittest.main()
