"""Phase 5B unit tests: allele calling, ICE profiling, plasmid screening,
mobile-element summaries, manifest loading.

Conventions: synthetic VCFs are bgzipped (as the pipeline produces);
the synthetic consensus is N-masked except at the curated codon sites.
"""
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..",
                                "workflow", "sentinel_lite"))
import report_lib

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
ALLELES_YAML = os.path.join(REPO, "workflow", "sentinel_lite", "alleles.yaml")

# Plus-strand reference bases at the 15 curated codon sites (derived
# 2026-10-07 from CP003069.1; see alleles.yaml header).
REF_SITES = {
    1041555: "T", 1041554: "T", 1041553: "A",   # ctxB codon 20 (minus)
    1041498: "G", 1041497: "T", 1041496: "A",   # ctxB codon 39 (minus)
    1041411: "T", 1041410: "G", 1041409: "A",   # ctxB codon 68 (minus)
    368214: "A", 368215: "G", 368216: "T",       # tcpA codon 89 (plus)
    1040742: "T", 1040743: "G", 1040744: "A",   # rtxA codon 4534 (plus)
}
CHROM = "CP003069.1"
CONS_LEN = 1042000

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def make_consensus(extra_n=()):
    seq = ["N"] * CONS_LEN
    for pos, base in REF_SITES.items():
        seq[pos - 1] = base
    for pos in extra_n:
        seq[pos - 1] = "N"
    return {CHROM: "".join(seq)}


def write_vcf_gz(records):
    """records: [(pos, ref, alt)] -> path of a temp bgzipped VCF."""
    fd, path = tempfile.mkstemp(suffix=".vcf.gz")
    os.close(fd)
    with gzip.open(path, "wt") as f:
        f.write(VCF_HEADER)
        for pos, ref, alt in records:
            f.write(f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t"
                    f"DP=50\tGT:DP\t1:50\n")
    return path


def write_consensus(seqs):
    fd, path = tempfile.mkstemp(suffix=".fasta")
    os.close(fd)
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")
    return path


class TestLoadAlleles(unittest.TestCase):
    def test_loads_curated_table(self):
        table = report_lib.load_alleles(ALLELES_YAML)
        self.assertEqual(len(table["codons"]), 5)
        self.assertEqual(len(table["cds_scans"]), 2)
        self.assertIn("ctxB", table["gene_rollups"])
        ids = {c["id"] for c in table["codons"]}
        self.assertEqual(ids, {"ctxB_20", "ctxB_39", "ctxB_68",
                               "tcpA_89", "rtxA_4534"})

    def test_ctxB20_coordinates(self):
        table = report_lib.load_alleles(ALLELES_YAML)
        c20 = next(c for c in table["codons"] if c["id"] == "ctxB_20")
        self.assertEqual(c20["positions"], [1041555, 1041554, 1041553])
        self.assertEqual(c20["strand"], "-")
        self.assertEqual(c20["ref_aa"], "N")

    def test_malformed_table_raises(self):
        fd, path = tempfile.mkstemp(suffix=".yaml")
        os.close(fd)
        with open(path, "w") as f:
            f.write("codons: []\n")
        with self.assertRaises(ValueError):
            report_lib.load_alleles(path)
        os.unlink(path)


class TestCallAlleles(unittest.TestCase):
    def setUp(self):
        self.table = report_lib.load_alleles(ALLELES_YAML)
        self.consensus = write_consensus(make_consensus())
        self.empty_vcf = write_vcf_gz([])

    def tearDown(self):
        for p in (self.consensus, self.empty_vcf):
            if os.path.exists(p):
                os.unlink(p)

    def test_reference_like_is_ctxB7(self):
        res = report_lib.call_alleles(self.empty_vcf, self.consensus,
                                      self.table)
        self.assertEqual(res["genes"]["ctxB"]["call"], "ctxB7")
        self.assertEqual(res["codons"]["ctxB_20"]["observed_aa"], "N")
        self.assertEqual(res["codons"]["ctxB_39"]["observed_aa"], "H")
        self.assertEqual(res["codons"]["ctxB_68"]["observed_aa"], "T")
        self.assertEqual(res["codons"]["tcpA_89"]["call"], "Wave-3-like")
        self.assertEqual(res["codons"]["tcpA_89"]["observed_aa"], "S")
        self.assertEqual(res["codons"]["rtxA_4534"]["observed_aa"], "*")
        self.assertIn("Weng", res["codons"]["ctxB_20"]["citation"])

    def test_classical_variant_at_ctxB20(self):
        # Plus-strand 1041555 T->G flips minus-strand codon base 1 A->C:
        # AAT(N) -> CAT(H) = classical-like.
        vcf = write_vcf_gz([(1041555, "T", "G")])
        try:
            res = report_lib.call_alleles(vcf, self.consensus, self.table)
            self.assertEqual(res["codons"]["ctxB_20"]["observed_aa"], "H")
            self.assertEqual(res["codons"]["ctxB_20"]["call"], "classical-like")
            # Mixed with the two ctxB7 codons -> gene call "other".
            self.assertEqual(res["genes"]["ctxB"]["call"], "other")
        finally:
            os.unlink(vcf)

    def test_n_masked_site_is_indeterminate(self):
        cons = write_consensus(make_consensus(extra_n=(1041555,)))
        try:
            res = report_lib.call_alleles(self.empty_vcf, cons, self.table)
            self.assertEqual(res["codons"]["ctxB_20"]["call"], "indeterminate")
            self.assertEqual(res["genes"]["ctxB"]["call"], "indeterminate")
            # Untouched codons still resolve.
            self.assertEqual(res["codons"]["tcpA_89"]["call"], "Wave-3-like")
        finally:
            os.unlink(cons)

    def test_indel_at_codon_site_is_indeterminate(self):
        vcf = write_vcf_gz([(368214, "A", "AG")])
        try:
            res = report_lib.call_alleles(vcf, self.consensus, self.table)
            self.assertEqual(res["codons"]["tcpA_89"]["call"], "indeterminate")
        finally:
            os.unlink(vcf)

    def test_hapR_scan_reports_variants(self):
        # SNP inside the hapR span (1-based 3023604-3024215).
        vcf = write_vcf_gz([(3024000, "A", "G")])
        try:
            res = report_lib.call_alleles(vcf, self.consensus, self.table)
            hapr = next(s for s in res["cds_scans"] if s["id"] == "hapR")
            self.assertEqual(hapr["n_variants"], 1)
            self.assertEqual(hapr["variants"][0]["pos"], 3024000)
            self.assertFalse(hapr["variants"][0]["frameshift_candidate"])
        finally:
            os.unlink(vcf)

    def test_frameshift_candidate_flagged(self):
        vcf = write_vcf_gz([(2586000, "AT", "A")])  # 1-bp del in vexB span
        try:
            res = report_lib.call_alleles(vcf, self.consensus, self.table)
            vexb = next(s for s in res["cds_scans"] if s["id"] == "vexB")
            self.assertEqual(vexb["n_variants"], 1)
            self.assertTrue(vexb["variants"][0]["frameshift_candidate"])
            self.assertEqual(vexb["variants"][0]["kind"], "indel")
        finally:
            os.unlink(vcf)

    def test_ruo_wording_present(self):
        res = report_lib.call_alleles(self.empty_vcf, self.consensus,
                                      self.table)
        ruo = res["codons"]["ctxB_20"]["ruo"]
        self.assertIn("observation", ruo)
        self.assertIn("not a prediction", ruo)


class TestProfileIceSegments(unittest.TestCase):
    def test_intact(self):
        r = report_lib.profile_ice_segments({
            "SXT_cargo1": 95.0, "SXT_backbone1": 99.0,
            "SXT_backbone2": 98.0, "SXT_cargo2": 96.0})
        self.assertEqual(r["pattern"], "intact")

    def test_icevchind5_like(self):
        r = report_lib.profile_ice_segments({
            "SXT_cargo1": 5.0, "SXT_backbone1": 99.0,
            "SXT_backbone2": 98.0, "SXT_cargo2": 96.0})
        self.assertEqual(r["pattern"], "ICEVchInd5-like")
        self.assertIn("consistent with", r["note"])
        self.assertNotIn("confirmed", r["note"])

    def test_other(self):
        r = report_lib.profile_ice_segments({
            "SXT_cargo1": 5.0, "SXT_backbone1": 5.0,
            "SXT_backbone2": 98.0, "SXT_cargo2": 96.0})
        self.assertEqual(r["pattern"], "other")

    def test_indeterminate(self):
        r = report_lib.profile_ice_segments({
            "SXT_cargo1": 95.0, "SXT_backbone1": None,
            "SXT_backbone2": 98.0, "SXT_cargo2": 96.0})
        self.assertEqual(r["pattern"], "indeterminate")


class TestPlasmidScreen(unittest.TestCase):
    FULL = ["blaPER-7", "mph(A)", "mph(E)", "msr(E)", "aadA2", "qac", "sul1"]

    def test_consistent(self):
        r = report_lib.plasmid_screen(self.FULL + ["strA"])
        self.assertEqual(r["call"], "consistent")
        self.assertEqual(r["genes_missing"], [])

    def test_partial(self):
        r = report_lib.plasmid_screen(["blaPER-7", "sul1"])
        self.assertEqual(r["call"], "partial")
        self.assertEqual(len(r["genes_missing"]), 5)

    def test_absent(self):
        r = report_lib.plasmid_screen(["strA", "floR"])
        self.assertEqual(r["call"], "absent")

    def test_alias(self):
        r = report_lib.plasmid_screen(
            ["blaPER-7", "mph(A)", "mph(E)", "msr(E)", "aadA2",
             "qacE", "sul1"])
        self.assertEqual(r["call"], "consistent")

    def test_never_claims_plasmid(self):
        r = report_lib.plasmid_screen(self.FULL)
        self.assertIn("cannot be proven", r["note"])


class TestSummarizeMobileElements(unittest.TestCase):
    def test_calls(self):
        rows = [
            {"element": "PLE1", "breadth_pct": 95.0, "mean_depth": 40.0,
             "citation": "x"},
            {"element": "PLE11", "breadth_pct": 5.0, "mean_depth": 1.0,
             "citation": "y"},
            {"element": "WonAB", "breadth_pct": 50.0, "mean_depth": 10.0,
             "citation": "z"},
            {"element": "Ghost", "breadth_pct": None, "mean_depth": None,
             "citation": ""},
        ]
        out = report_lib.summarize_mobile_elements(rows)
        self.assertEqual(out["PLE1"]["call"], "present")
        self.assertEqual(out["PLE11"]["call"], "absent")
        self.assertEqual(out["WonAB"]["call"], "partial")
        self.assertNotIn("Ghost", out)


class TestLoadMobileElements(unittest.TestCase):
    def test_empty_manifest(self):
        fd, path = tempfile.mkstemp(suffix=".yaml")
        os.close(fd)
        with open(path, "w") as f:
            f.write("elements: []\n")
        try:
            els = report_lib.load_mobile_elements(path)
            self.assertEqual(els, [])
        finally:
            os.unlink(path)

    def test_malformed_entry_raises(self):
        fd, path = tempfile.mkstemp(suffix=".yaml")
        os.close(fd)
        with open(path, "w") as f:
            f.write("elements:\n  - id: PLE1\n")
        try:
            with self.assertRaises(ValueError):
                report_lib.load_mobile_elements(path)
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
