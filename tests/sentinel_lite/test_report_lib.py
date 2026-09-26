import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import report_lib  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


class TestParseKrakenReport(unittest.TestCase):
    """This parser replaced a real bug: the old code used
    max(clade_reads) over every row as its denominator, which is not the
    total and made v_cholerae_fraction meaningless regardless of database
    size."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _report(self, rows):
        p = Path(self.tmp.name) / "kraken.report"
        write(p, "\n".join(rows) + "\n")
        return p

    def test_pure_cholerae_sample(self):
        p = self._report([
            "0.50\t500\t500\tU\t0\tunclassified",
            "99.50\t99500\t0\tR\t1\troot",
            "99.40\t99400\t99400\tS\t666\tVibrio cholerae",
        ])
        chol, unclass, root, total = report_lib.parse_kraken_report(p)
        self.assertEqual(chol, 99400)
        self.assertEqual(unclass, 500)
        self.assertEqual(root, 99500)
        self.assertEqual(total, 100000)  # this is the number that was wrong before

    def test_mixed_sample_low_purity(self):
        # what the fixed denominator was for: species_purity must reflect
        # cholerae's share of ALL reads, not just the classified ones.
        p = self._report([
            "60.00\t60000\t60000\tU\t0\tunclassified",
            "40.00\t40000\t0\tR\t1\troot",
            "0.06\t56\t56\tS\t666\tVibrio cholerae",
            "39.94\t39944\t39944\tS\t644\tAeromonas hydrophila",
        ])
        chol, unclass, root, total = report_lib.parse_kraken_report(p)
        self.assertEqual(total, 100000)
        frac, pct = report_lib.species_purity(chol, total)
        self.assertAlmostEqual(pct, 0.056, places=2)

    def test_no_cholerae_row_at_all(self):
        p = self._report([
            "100.00\t1000\t1000\tU\t0\tunclassified",
            "0.00\t0\t0\tR\t1\troot",
        ])
        chol, unclass, root, total = report_lib.parse_kraken_report(p)
        self.assertEqual(chol, 0)
        self.assertEqual(total, 1000)

    def test_malformed_lines_are_skipped_not_fatal(self):
        p = self._report([
            "not\tenough\tcolumns",
            "0.50\t500\t500\tU\t0\tunclassified",
            "99.50\tNaN\t0\tR\t1\troot",  # non-integer clade count
        ])
        chol, unclass, root, total = report_lib.parse_kraken_report(p)
        self.assertEqual(unclass, 500)
        self.assertEqual(root, 0)  # malformed row skipped, not crashed on

    def test_empty_report_gives_zero_total(self):
        p = self._report([])
        chol, unclass, root, total = report_lib.parse_kraken_report(p)
        self.assertEqual(total, 0)


class TestSpeciesPurity(unittest.TestCase):
    def test_zero_total_is_zero_not_zerodivisionerror(self):
        frac, pct = report_lib.species_purity(0, 0)
        self.assertEqual((frac, pct), (0.0, 0.0))

    def test_full_purity(self):
        frac, pct = report_lib.species_purity(1000, 1000)
        self.assertEqual((frac, pct), (1.0, 100.0))


class TestGenomeWideCoverage(unittest.TestCase):
    """Regression coverage for the chr1-only bug: an earlier version broke
    out of the loop after the first data row, silently ignoring chr2."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _cov(self, rows):
        p = Path(self.tmp.name) / "coverage.txt"
        header = "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
        write(p, header + "\n" + "\n".join(rows) + "\n")
        return p

    def test_two_contigs_are_length_weighted(self):
        # chr1: 100bp, 100% covered, depth 10. chr2: 100bp, 100% covered, depth 50.
        # A chr1-only read would report depth 10, not the weighted 30.
        p = self._cov([
            "CP003069.1\t1\t100\t500\t100\t100.0\t10.0\t30\t60",
            "CP003070.1\t1\t100\t500\t100\t100.0\t50.0\t30\t60",
        ])
        mean_depth, breadth = report_lib.genome_wide_coverage(p)
        self.assertAlmostEqual(mean_depth, 30.0)
        self.assertAlmostEqual(breadth, 100.0)

    def test_partial_breadth_on_one_contig(self):
        p = self._cov([
            "CP003069.1\t1\t100\t500\t50\t50.0\t10.0\t30\t60",   # half covered
            "CP003070.1\t1\t100\t500\t100\t100.0\t10.0\t30\t60",  # fully covered
        ])
        mean_depth, breadth = report_lib.genome_wide_coverage(p)
        self.assertAlmostEqual(breadth, 75.0)  # (50+100)/200

    def test_empty_table_gives_zero_not_crash(self):
        p = self._cov([])
        mean_depth, breadth = report_lib.genome_wide_coverage(p)
        self.assertEqual((mean_depth, breadth), (0.0, 0.0))


class TestParseFlagstatMapped(unittest.TestCase):
    def test_typical_flagstat_text(self):
        text = (
            "700172 + 0 in total (QC-passed reads + QC-failed reads)\n"
            "700172 + 0 mapped (100.00% : N/A)\n"
        )
        self.assertEqual(report_lib.parse_flagstat_mapped(text), 700172)

    def test_no_matching_line_returns_zero(self):
        self.assertEqual(report_lib.parse_flagstat_mapped("nothing here"), 0)


class TestLociCalling(unittest.TestCase):
    """This is the fix for the AT-rich-locus bug: toxT (~28% GC) was
    called absent under a mean-depth>=5x rule despite being present, with
    no soft-clips or split reads at its boundaries. Breadth doesn't have
    that GC-bias failure mode."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _bed(self, rows):
        p = Path(self.tmp.name) / "loci.bed"
        write(p, "\n".join(rows) + "\n")
        return p

    def _loci(self, rows):
        p = Path(self.tmp.name) / "loci_coverage.txt"
        header = "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
        write(p, header + "\n" + "\n".join(rows) + "\n")
        return p

    def test_present_by_breadth_despite_low_depth(self):
        bed = self._bed(["CP003069.1\t100\t200\ttoxT\t.\t+"])
        # startpos is 1-based = bed start + 1 = 101
        loci = self._loci(["CP003069.1\t101\t200\t50\t100\t100.0\t4.3\t30\t60"])
        bed_labels = report_lib.parse_loci_bed(bed)
        calls = report_lib.call_loci(loci, bed_labels, present_breadth=80, absent_breadth=20)
        self.assertEqual(calls["toxT"]["call"], "present")
        self.assertTrue(calls["toxT"]["present"])
        self.assertEqual(calls["toxT"]["mean_depth"], 4.3)  # depth still reported

    def test_absent_below_threshold(self):
        bed = self._bed(["CP003069.1\t100\t200\tctxA\t.\t+"])
        loci = self._loci(["CP003069.1\t101\t200\t0\t0\t0.0\t0.0\t0\t0"])
        bed_labels = report_lib.parse_loci_bed(bed)
        calls = report_lib.call_loci(loci, bed_labels, present_breadth=80, absent_breadth=20)
        self.assertEqual(calls["ctxA"]["call"], "absent")
        self.assertFalse(calls["ctxA"]["present"])

    def test_partial_between_thresholds(self):
        bed = self._bed(["CP003069.1\t100\t200\twbeT\t.\t+"])
        loci = self._loci(["CP003069.1\t101\t200\t20\t50\t50.0\t3.0\t30\t60"])
        bed_labels = report_lib.parse_loci_bed(bed)
        calls = report_lib.call_loci(loci, bed_labels, present_breadth=80, absent_breadth=20)
        self.assertEqual(calls["wbeT"]["call"], "partial")
        self.assertFalse(calls["wbeT"]["present"])

    def test_unmatched_coordinates_fall_back_to_coordinate_name(self):
        bed = self._bed(["CP003069.1\t100\t200\tctxA\t.\t+"])
        # loci row coordinates don't match any BED entry
        loci = self._loci(["CP003069.1\t500\t600\t10\t50\t50.0\t3.0\t30\t60"])
        bed_labels = report_lib.parse_loci_bed(bed)
        calls = report_lib.call_loci(loci, bed_labels, present_breadth=80, absent_breadth=20)
        self.assertIn("CP003069.1:500-600", calls)

    def test_failed_marker_line_is_skipped(self):
        bed = self._bed(["CP003069.1\t100\t200\tctxA\t.\t+"])
        loci_path = Path(self.tmp.name) / "loci_coverage.txt"
        write(loci_path, "loci coverage failed (check BED)\n")
        bed_labels = report_lib.parse_loci_bed(bed)
        calls = report_lib.call_loci(loci_path, bed_labels, present_breadth=80, absent_breadth=20)
        self.assertEqual(calls, {})


class TestConsensusStats(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _fa(self, text):
        p = Path(self.tmp.name) / "consensus.fasta"
        write(p, text)
        return p

    def test_fully_called_no_n(self):
        p = self._fa(">chr1\nACGTACGT\n")
        length, n, pct = report_lib.consensus_stats(p)
        self.assertEqual((length, n), (8, 0))
        self.assertEqual(pct, 100.0)

    def test_masked_positions_reduce_called_pct(self):
        p = self._fa(">chr1\nACGTNNNN\n")
        length, n, pct = report_lib.consensus_stats(p)
        self.assertEqual((length, n), (8, 4))
        self.assertEqual(pct, 50.0)

    def test_lowercase_n_is_counted(self):
        p = self._fa(">chr1\nacgtnnnn\n")
        _, n, pct = report_lib.consensus_stats(p)
        self.assertEqual(n, 4)

    def test_multi_record_fasta_sums_across_contigs(self):
        p = self._fa(">chr1\nACGT\n>chr2\nNNNN\n")
        length, n, pct = report_lib.consensus_stats(p)
        self.assertEqual((length, n), (8, 4))
        self.assertEqual(pct, 50.0)

    def test_empty_file_is_zero_not_crash(self):
        p = self._fa("")
        length, n, pct = report_lib.consensus_stats(p)
        self.assertEqual((length, n, pct), (0, 0, 0.0))


class TestReadSnpCount(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_normal_count(self):
        p = Path(self.tmp.name) / "snp_count.txt"
        write(p, "42\n")
        self.assertEqual(report_lib.read_snp_count(p), 42)

    def test_non_numeric_content_falls_back_to_zero(self):
        p = Path(self.tmp.name) / "snp_count.txt"
        write(p, "loci coverage failed (check BED)\n")
        self.assertEqual(report_lib.read_snp_count(p), 0)


class TestEvaluateQc(unittest.TestCase):
    THRESHOLDS = {"min_mean_depth": 20, "min_called_pct": 90, "min_species_purity": 90}

    def test_all_pass(self):
        status, reasons = report_lib.evaluate_qc(50.0, 98.0, 99.0, self.THRESHOLDS)
        self.assertEqual(status, "pass")
        self.assertEqual(reasons, [])

    def test_fails_on_depth_only(self):
        status, reasons = report_lib.evaluate_qc(0.03, 98.0, 99.0, self.THRESHOLDS)
        self.assertEqual(status, "fail")
        self.assertEqual(len(reasons), 1)
        self.assertIn("mean_depth", reasons[0])

    def test_fails_on_all_three_independently(self):
        # this is the real ERR11684929 case: every one of these fails at once
        status, reasons = report_lib.evaluate_qc(0.01, 0.0, 0.08, self.THRESHOLDS)
        self.assertEqual(status, "fail")
        self.assertEqual(len(reasons), 3)

    def test_exactly_at_threshold_passes(self):
        # boundary: >= threshold passes, not > threshold
        status, reasons = report_lib.evaluate_qc(20.0, 90.0, 90.0, self.THRESHOLDS)
        self.assertEqual(status, "pass")

    def test_just_below_threshold_fails(self):
        status, reasons = report_lib.evaluate_qc(19.99, 90.0, 90.0, self.THRESHOLDS)
        self.assertEqual(status, "fail")
        self.assertEqual(len(reasons), 1)


if __name__ == "__main__":
    unittest.main()
