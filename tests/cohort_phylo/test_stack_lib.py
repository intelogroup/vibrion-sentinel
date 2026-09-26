import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "cohort_phylo"))
import stack_lib  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)


class TestLoadMask(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_file_returns_empty_dict_not_exception(self):
        self.assertEqual(stack_lib.load_mask(Path(self.tmp.name) / "nope.bed"), {})

    def test_parses_intervals_grouped_by_contig(self):
        p = Path(self.tmp.name) / "mask.bed"
        write(p, "chr1\t10\t20\trRNA\nchr1\t30\t40\tmobile\nchr2\t5\t15\trRNA\n")
        masked = stack_lib.load_mask(p)
        self.assertEqual(masked, {"chr1": [(10, 20), (30, 40)], "chr2": [(5, 15)]})

    def test_comments_and_blank_lines_skipped(self):
        p = Path(self.tmp.name) / "mask.bed"
        write(p, "# comment\n\nchr1\t10\t20\trRNA\n")
        self.assertEqual(stack_lib.load_mask(p), {"chr1": [(10, 20)]})


class TestReadConsensus(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_prefers_gz_over_plain(self):
        write_gz(Path(self.tmp.name) / "SRR1.fasta.gz", ">chr1\nACGT\n")
        write(Path(self.tmp.name) / "SRR1.fasta", ">chr1\nTTTT\n")
        c = stack_lib.read_consensus(self.tmp.name, "SRR1")
        self.assertEqual(c, {"chr1": "ACGT"})

    def test_plain_fasta_when_no_gz(self):
        write(Path(self.tmp.name) / "SRR2.fasta", ">chr1\nACGT\n")
        c = stack_lib.read_consensus(self.tmp.name, "SRR2")
        self.assertEqual(c, {"chr1": "ACGT"})

    def test_multi_contig_and_multi_line_sequence(self):
        write_gz(Path(self.tmp.name) / "SRR3.fasta.gz",
                 ">chr1\nACGT\nACGT\n>chr2\nTTTT\n")
        c = stack_lib.read_consensus(self.tmp.name, "SRR3")
        self.assertEqual(c, {"chr1": "ACGTACGT", "chr2": "TTTT"})

    def test_missing_accession_returns_none(self):
        self.assertIsNone(stack_lib.read_consensus(self.tmp.name, "SRR404"))


class TestLoadAndFilterSamples(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_consensus_excluded_with_distinct_reason(self):
        loaded, excluded = stack_lib.load_and_filter_samples(
            ["SRR404"], self.tmp.name, min_sample_called=0.8)
        self.assertEqual(loaded, {})
        self.assertEqual(excluded, [("SRR404", "no consensus available")])

    def test_low_completeness_excluded_with_distinct_reason(self):
        write_gz(Path(self.tmp.name) / "SRR1.fasta.gz", ">chr1\nACNN\n")  # 50% called
        loaded, excluded = stack_lib.load_and_filter_samples(
            ["SRR1"], self.tmp.name, min_sample_called=0.8)
        self.assertEqual(loaded, {})
        self.assertEqual(len(excluded), 1)
        self.assertEqual(excluded[0][0], "SRR1")
        self.assertIn("called", excluded[0][1])

    def test_sample_clearing_threshold_is_loaded(self):
        write_gz(Path(self.tmp.name) / "SRR1.fasta.gz", ">chr1\nACGT\n")  # 100% called
        loaded, excluded = stack_lib.load_and_filter_samples(
            ["SRR1"], self.tmp.name, min_sample_called=0.8)
        self.assertEqual(excluded, [])
        self.assertIn("SRR1", loaded)
        seqs, called = loaded["SRR1"]
        self.assertEqual(seqs, {"chr1": "ACGT"})
        self.assertEqual(called, 1.0)

    def test_mixed_cohort_partitions_correctly(self):
        write_gz(Path(self.tmp.name) / "SRR1.fasta.gz", ">chr1\nACGT\n")
        write_gz(Path(self.tmp.name) / "SRR2.fasta.gz", ">chr1\nNNNN\n")  # 0% called
        loaded, excluded = stack_lib.load_and_filter_samples(
            ["SRR1", "SRR2", "SRR3"], self.tmp.name, min_sample_called=0.8)
        self.assertEqual(set(loaded), {"SRR1"})
        excluded_accs = {a for a, _ in excluded}
        self.assertEqual(excluded_accs, {"SRR2", "SRR3"})


class TestValidateCoordinates(unittest.TestCase):
    def test_matching_contigs_and_lengths_pass(self):
        loaded = {
            "SRR1": ({"chr1": "ACGT", "chr2": "TT"}, 1.0),
            "SRR2": ({"chr1": "AAGT", "chr2": "GG"}, 1.0),
        }
        contigs, lengths = stack_lib.validate_coordinates(loaded)
        self.assertEqual(contigs, ["chr1", "chr2"])
        self.assertEqual(lengths, {"chr1": 4, "chr2": 2})

    def test_length_mismatch_raises_naming_the_accession(self):
        # this is the real bug from earlier this session: bcftools
        # consensus applying indels made contig lengths differ by sample
        loaded = {
            "SRR1": ({"chr1": "ACGT"}, 1.0),
            "SRR2": ({"chr1": "ACGTA"}, 1.0),  # one base longer
        }
        with self.assertRaises(RuntimeError) as ctx:
            stack_lib.validate_coordinates(loaded)
        self.assertIn("SRR2", str(ctx.exception))
        self.assertIn("chr1", str(ctx.exception))

    def test_missing_contig_in_one_sample_raises(self):
        loaded = {
            "SRR1": ({"chr1": "ACGT", "chr2": "TT"}, 1.0),
            "SRR2": ({"chr1": "ACGT"}, 1.0),  # chr2 entirely absent
        }
        with self.assertRaises(RuntimeError):
            stack_lib.validate_coordinates(loaded)

    def test_empty_loaded_raises_rather_than_indexerror(self):
        with self.assertRaises(RuntimeError):
            stack_lib.validate_coordinates({})


class TestBuildSkipMask(unittest.TestCase):
    def test_masked_positions_marked(self):
        masked = {"chr1": [(2, 4)]}
        skip, n_masked = stack_lib.build_skip_mask(masked, ["chr1"], {"chr1": 10})
        self.assertEqual(list(skip["chr1"]), [0, 0, 1, 1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(n_masked, 2)

    def test_overlapping_intervals_counted_once(self):
        masked = {"chr1": [(2, 5), (3, 6)]}
        skip, n_masked = stack_lib.build_skip_mask(masked, ["chr1"], {"chr1": 10})
        self.assertEqual(n_masked, 4)  # positions 2,3,4,5 -- not 6 (3+3)

    def test_contig_not_in_alignment_is_ignored(self):
        masked = {"chr_unrelated": [(0, 5)]}
        skip, n_masked = stack_lib.build_skip_mask(masked, ["chr1"], {"chr1": 10})
        self.assertEqual(n_masked, 0)
        self.assertEqual(list(skip["chr1"]), [0] * 10)

    def test_interval_extending_past_contig_length_is_clipped(self):
        masked = {"chr1": [(8, 20)]}
        skip, n_masked = stack_lib.build_skip_mask(masked, ["chr1"], {"chr1": 10})
        self.assertEqual(n_masked, 2)  # positions 8,9 only

    def test_no_mask_at_all(self):
        skip, n_masked = stack_lib.build_skip_mask({}, ["chr1"], {"chr1": 5})
        self.assertEqual(n_masked, 0)
        self.assertEqual(list(skip["chr1"]), [0] * 5)


class TestBuildAlignmentColumns(unittest.TestCase):
    def test_variable_site_kept(self):
        loaded = {"SRR1": ({"chr1": "AAAA"}, 1.0), "SRR2": ({"chr1": "AAGA"}, 1.0)}
        skip = {"chr1": bytearray(4)}
        accs, cols, total, dropped_missing, invariant = stack_lib.build_alignment_columns(
            loaded, ["chr1"], {"chr1": 4}, skip, max_site_missing=0.1)
        self.assertEqual(accs, ["SRR1", "SRR2"])
        self.assertEqual(total, 4)
        self.assertEqual(invariant, 3)  # positions 0,1,3
        self.assertEqual(cols, [["A", "G"]])  # position 2

    def test_invariant_site_dropped(self):
        loaded = {"SRR1": ({"chr1": "AAAA"}, 1.0), "SRR2": ({"chr1": "AAAA"}, 1.0)}
        skip = {"chr1": bytearray(4)}
        _, cols, _, _, invariant = stack_lib.build_alignment_columns(
            loaded, ["chr1"], {"chr1": 4}, skip, max_site_missing=0.1)
        self.assertEqual(cols, [])
        self.assertEqual(invariant, 4)

    def test_high_missingness_site_dropped_not_counted_as_invariant(self):
        # one N among two samples at a variable-looking site: with
        # max_site_missing=0.1 (needs <=10% N), 1/2 = 50% N drops it
        loaded = {"SRR1": ({"chr1": "A"}, 1.0), "SRR2": ({"chr1": "N"}, 1.0)}
        skip = {"chr1": bytearray(1)}
        _, cols, total, dropped_missing, invariant = stack_lib.build_alignment_columns(
            loaded, ["chr1"], {"chr1": 1}, skip, max_site_missing=0.1)
        self.assertEqual(cols, [])
        self.assertEqual(dropped_missing, 1)
        self.assertEqual(invariant, 0)

    def test_masked_site_excluded_from_every_bucket(self):
        # a masked position must not be counted in total_sites,
        # dropped_missing, invariant, or cols -- it's outside the
        # partition entirely
        loaded = {"SRR1": ({"chr1": "AG"}, 1.0), "SRR2": ({"chr1": "AC"}, 1.0)}
        skip = {"chr1": bytearray([1, 0])}  # position 0 masked
        _, cols, total, dropped_missing, invariant = stack_lib.build_alignment_columns(
            loaded, ["chr1"], {"chr1": 2}, skip, max_site_missing=0.1)
        self.assertEqual(total, 1)  # only position 1 considered
        self.assertEqual(cols, [["G", "C"]])

    def test_column_order_follows_sorted_accession_order(self):
        loaded = {"SRR2": ({"chr1": "G"}, 1.0), "SRR1": ({"chr1": "A"}, 1.0)}
        skip = {"chr1": bytearray(1)}
        accs, cols, *_ = stack_lib.build_alignment_columns(
            loaded, ["chr1"], {"chr1": 1}, skip, max_site_missing=0.1)
        self.assertEqual(accs, ["SRR1", "SRR2"])
        self.assertEqual(cols, [["A", "G"]])  # SRR1's base first


class TestFormatAlignmentFasta(unittest.TestCase):
    def test_basic_format(self):
        accs = ["SRR1", "SRR2"]
        cols = [["A", "G"], ["C", "C"]]
        text = stack_lib.format_alignment_fasta(accs, cols)
        self.assertEqual(text, ">SRR1\nAC\n>SRR2\nGC\n")

    def test_no_columns_still_writes_headers_with_empty_sequence(self):
        text = stack_lib.format_alignment_fasta(["SRR1"], [])
        self.assertEqual(text, ">SRR1\n\n")

    def test_no_accessions_yields_empty_string(self):
        self.assertEqual(stack_lib.format_alignment_fasta([], []), "")


class TestFormatExcludedTsv(unittest.TestCase):
    def test_header_and_rows(self):
        text = stack_lib.format_excluded_tsv([("SRR1", "no consensus available")])
        self.assertEqual(text, "accession\treason\nSRR1\tno consensus available\n")

    def test_empty_list_still_has_header(self):
        self.assertEqual(stack_lib.format_excluded_tsv([]), "accession\treason\n")


class TestPairwiseSnpDistances(unittest.TestCase):
    def test_identical_sequences_zero_distance(self):
        accs, dist = stack_lib.pairwise_snp_distances({"SRR1": "ACGT", "SRR2": "ACGT"})
        self.assertEqual(dist[("SRR1", "SRR2")], 0)
        self.assertEqual(dist[("SRR1", "SRR1")], 0)

    def test_counts_only_real_base_mismatches(self):
        accs, dist = stack_lib.pairwise_snp_distances({"SRR1": "ACGT", "SRR2": "ACGA"})
        self.assertEqual(dist[("SRR1", "SRR2")], 1)

    def test_n_positions_excluded_from_denominator_not_counted_as_mismatch(self):
        # position 3: SRR1 has N, SRR2 has A -- must not count as a diff
        accs, dist = stack_lib.pairwise_snp_distances({"SRR1": "ACGN", "SRR2": "ACGA"})
        self.assertEqual(dist[("SRR1", "SRR2")], 0)

    def test_matrix_is_symmetric(self):
        accs, dist = stack_lib.pairwise_snp_distances({"SRR1": "ACGT", "SRR2": "ACGA", "SRR3": "TCGA"})
        for a in accs:
            for b in accs:
                self.assertEqual(dist[(a, b)], dist[(b, a)])

    def test_accessions_returned_sorted(self):
        accs, _ = stack_lib.pairwise_snp_distances({"SRR2": "A", "SRR1": "A"})
        self.assertEqual(accs, ["SRR1", "SRR2"])


class TestFormatDistanceMatrix(unittest.TestCase):
    def test_basic_matrix(self):
        accs = ["SRR1", "SRR2"]
        dist = {("SRR1", "SRR1"): 0, ("SRR1", "SRR2"): 5,
                ("SRR2", "SRR1"): 5, ("SRR2", "SRR2"): 0}
        text = stack_lib.format_distance_matrix(accs, dist)
        self.assertEqual(text, "\tSRR1\tSRR2\nSRR1\t0\t5\nSRR2\t5\t0\n")


class TestParseFasta(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_multi_record_multi_line(self):
        p = Path(self.tmp.name) / "aln.fasta"
        write(p, ">SRR1\nACGT\nACGT\n>SRR2\nTTTT\n")
        seqs = stack_lib.parse_fasta(p)
        self.assertEqual(seqs, {"SRR1": "ACGTACGT", "SRR2": "TTTT"})

    def test_empty_file(self):
        p = Path(self.tmp.name) / "aln.fasta"
        write(p, "")
        self.assertEqual(stack_lib.parse_fasta(p), {})


if __name__ == "__main__":
    unittest.main()
