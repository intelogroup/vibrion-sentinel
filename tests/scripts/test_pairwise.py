import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))
import pairwise  # noqa: E402


class TestPairwiseSnpDistance(unittest.TestCase):
    def test_identical_is_zero(self):
        d, n = pairwise.pairwise_snp_distance("ACGTACGT", "ACGTACGT")
        self.assertEqual((d, n), (0, 8))

    def test_counts_differences(self):
        d, n = pairwise.pairwise_snp_distance("ACGT", "ACGA")
        self.assertEqual((d, n), (1, 4))

    def test_n_sites_excluded(self):
        # N vs called, called vs N, N vs N: none of them count.
        d, n = pairwise.pairwise_snp_distance("ANGT", "NCNT")
        self.assertEqual(d, 0)
        self.assertEqual(n, 1)  # only position 2 (G vs G)

    def test_all_n_compares_nothing(self):
        d, n = pairwise.pairwise_snp_distance("NNNN", "ACGT")
        self.assertEqual((d, n), (0, 0))

    def test_length_mismatch_raises(self):
        with self.assertRaises(ValueError):
            pairwise.pairwise_snp_distance("ACGT", "ACG")

    def test_load_consensus_skips_headers_and_uppercases(self):
        with tempfile.NamedTemporaryFile("w", suffix=".fasta",
                                         delete=False) as f:
            f.write(">SRR1 consensus\nacgt\nTGCA\n")
            p = f.name
        self.assertEqual(pairwise.load_consensus(p), "ACGTTGCA")


class TestDistanceMatrix(unittest.TestCase):
    def test_pairs_and_sorting(self):
        seqs = {"b": "ACGT", "a": "ACGT", "c": "ACGA"}
        m = pairwise.distance_matrix(seqs)
        self.assertEqual(set(m), {("a", "b"), ("a", "c"), ("b", "c")})
        self.assertEqual(m[("a", "b")], (0, 4))
        self.assertEqual(m[("a", "c")], (1, 4))


class TestConfirmCluster(unittest.TestCase):
    def test_confirmed_when_all_within_tol(self):
        r = pairwise.confirm_cluster(["a", "b"], {"a": "ACGT", "b": "ACGA"})
        self.assertTrue(r["pairwise_confirmed"])
        self.assertEqual(r["pairwise_diameter"], 1)
        self.assertEqual(r["missing_consensus"], [])

    def test_not_confirmed_when_any_pair_exceeds(self):
        # proxy screen can over-include: true distance above tol.
        r = pairwise.confirm_cluster(
            ["a", "b", "c"],
            {"a": "AAAAAAAA", "b": "AAAAAAAT", "c": "TTTTTTTT"},
            snp_tol=5)
        self.assertFalse(r["pairwise_confirmed"])
        self.assertEqual(r["pairwise_diameter"], 8)

    def test_none_when_no_consensus_available(self):
        r = pairwise.confirm_cluster(["a", "b"], {})
        self.assertIsNone(r["pairwise_confirmed"])
        self.assertIsNone(r["pairwise_diameter"])
        self.assertEqual(r["missing_consensus"], ["a", "b"])

    def test_partial_availability(self):
        r = pairwise.confirm_cluster(["a", "b", "c"], {"a": "ACGT", "b": "ACGT"})
        self.assertTrue(r["pairwise_confirmed"])  # only pair (a,b) confirmed
        self.assertEqual(r["missing_consensus"], ["c"])
        self.assertEqual(set(r["pairwise_snps"]), {("a", "b")})


if __name__ == "__main__":
    unittest.main()
