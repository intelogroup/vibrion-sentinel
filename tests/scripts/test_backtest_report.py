import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))
import backtest_report  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


class TestLoadExpect(unittest.TestCase):
    """The minimal YAML subset reader. This already had one real bug this
    session (an UnboundLocalError from `data, stack = {}, [(0, data)]`
    referencing `data` on the same line it's assigned) -- these tests
    exercise every construct cohorts/*.expect.yaml actually uses."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _yaml(self, text):
        p = Path(self.tmp.name) / "test.yaml"
        write(p, text)
        return p

    def test_scalar_types(self):
        p = self._yaml("name: backtest-2022\ncount: 42\nfrac: 0.9\n")
        d = backtest_report.load_expect(p)
        self.assertEqual(d, {"name": "backtest-2022", "count": 42, "frac": 0.9})
        self.assertIsInstance(d["count"], int)
        self.assertIsInstance(d["frac"], float)

    def test_nested_mapping(self):
        p = self._yaml("qc:\n  min_pass_rate: 0.8\n")
        d = backtest_report.load_expect(p)
        self.assertEqual(d, {"qc": {"min_pass_rate": 0.8}})

    def test_multiple_keys_in_nested_mapping(self):
        p = self._yaml("snps_vs_7pet:\n  max_median: 400\n  max_iqr: 150\n")
        d = backtest_report.load_expect(p)
        self.assertEqual(d, {"snps_vs_7pet": {"max_median": 400, "max_iqr": 150}})

    def test_dedent_ends_nested_mapping(self):
        p = self._yaml("qc:\n  min_pass_rate: 0.8\ncohort: backtest-2022\n")
        d = backtest_report.load_expect(p)
        self.assertEqual(d, {"qc": {"min_pass_rate": 0.8}, "cohort": "backtest-2022"})

    def test_comments_and_blank_lines_ignored(self):
        p = self._yaml("# a comment\ncohort: backtest-2022\n\n# another\n")
        d = backtest_report.load_expect(p)
        self.assertEqual(d, {"cohort": "backtest-2022"})

    def test_quoted_string_value_unquoted(self):
        p = self._yaml("manifest: 'cohorts/backtest-2022.tsv'\n")
        d = backtest_report.load_expect(p)
        self.assertEqual(d["manifest"], "cohorts/backtest-2022.tsv")

    def test_block_scalar_with_greater_than_is_skipped_not_crashed(self):
        # cohorts/backtest-2022.expect.yaml has a `notes: >` block; this
        # reader doesn't need to reconstruct its text, just not choke on it
        p = self._yaml("cohort: backtest-2022\nnotes: >\n  line one\n  line two\n")
        d = backtest_report.load_expect(p)
        self.assertEqual(d["cohort"], "backtest-2022")
        self.assertEqual(d["notes"], "")

    def test_real_expect_file_parses(self):
        # exercises the actual committed file end to end
        real = Path(__file__).parent.parent.parent / "cohorts" / "backtest-2022.expect.yaml"
        d = backtest_report.load_expect(real)
        self.assertEqual(d["cohort"], "backtest-2022")
        self.assertIn("loci_present", d)
        self.assertIn("ctxA", d["loci_present"])
        self.assertEqual(d["qc"]["min_pass_rate"], 0.8)


class TestCohortAccessions(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_header_and_comments_skipped(self):
        p = Path(self.tmp.name) / "cohort.tsv"
        write(p, "# cohort: x\n# runs: 2\nrun_accession\tother\nSRR1\tfoo\nSRR2\tbar\n")
        self.assertEqual(backtest_report.cohort_accessions(p), ["SRR1", "SRR2"])

    def test_empty_manifest(self):
        p = Path(self.tmp.name) / "cohort.tsv"
        write(p, "# cohort: x\nrun_accession\n")
        self.assertEqual(backtest_report.cohort_accessions(p), [])


class TestEvaluateChecks(unittest.TestCase):
    """The check math itself: pass rate, per-locus presence fraction,
    SNP median/IQR -- the numbers a back-test result actually turns on."""

    def report(self, qc_status="pass", loci=None, snps=0):
        return {
            "qc": {"status": qc_status},
            "surveillance_loci": loci or {},
            "variants": {"snps_vs_7pet": snps},
        }

    def get_check(self, checks, name):
        for n, obs, want, ok in checks:
            if n == name:
                return obs, want, ok
        self.fail(f"no check named {name!r} in {[c[0] for c in checks]}")

    def test_qc_pass_rate_below_floor_fails(self):
        reports = [self.report("pass"), self.report("fail"), self.report("fail")]
        exp = {"qc": {"min_pass_rate": 0.8}}
        checks = backtest_report.evaluate_checks(reports, exp)
        _, _, ok = self.get_check(checks, "qc pass rate")
        self.assertFalse(ok)

    def test_no_qc_pass_samples_returns_only_the_qc_check(self):
        reports = [self.report("fail"), self.report("fail")]
        exp = {"qc": {"min_pass_rate": 0.8}, "loci_present": {"ctxA": 0.9}}
        checks = backtest_report.evaluate_checks(reports, exp)
        self.assertEqual(len(checks), 1)
        self.assertEqual(checks[0][0], "qc pass rate")

    def test_locus_present_fraction_computed_over_qc_pass_only(self):
        reports = [
            self.report("pass", loci={"ctxA": {"call": "present"}}),
            self.report("pass", loci={"ctxA": {"call": "absent"}}),
            # a QC-fail sample's loci calls must not count either way
            self.report("fail", loci={"ctxA": {"call": "absent"}}),
        ]
        exp = {"qc": {"min_pass_rate": 0.0}, "loci_present": {"ctxA": 0.9}}
        checks = backtest_report.evaluate_checks(reports, exp)
        obs, want, ok = self.get_check(checks, "locus ctxA present")
        self.assertIn("0.50", obs)
        self.assertIn("(1/2)", obs)
        self.assertFalse(ok)

    def test_missing_locus_call_counts_as_not_present(self):
        # a locus entirely absent from surveillance_loci (e.g. not in the
        # BED at all) must not crash and must not count as present
        reports = [self.report("pass", loci={})]
        exp = {"qc": {"min_pass_rate": 0.0}, "loci_present": {"ctxA": 0.9}}
        checks = backtest_report.evaluate_checks(reports, exp)
        obs, want, ok = self.get_check(checks, "locus ctxA present")
        self.assertIn("(0/1)", obs)
        self.assertFalse(ok)

    def test_snps_median_and_iqr_with_fewer_than_four_samples(self):
        # statistics.quantiles needs >=4 points; below that this falls
        # back to max-min as the spread, per the code's own comment
        reports = [self.report("pass", snps=s) for s in (40, 42, 44)]
        exp = {"qc": {"min_pass_rate": 0.0}, "snps_vs_7pet": {"max_median": 400, "max_iqr": 150}}
        checks = backtest_report.evaluate_checks(reports, exp)
        obs, want, ok = self.get_check(checks, "snps median")
        self.assertEqual(obs, "42")
        obs, want, ok = self.get_check(checks, "snps spread (IQR, n=3)")
        self.assertEqual(obs, "4")  # max(44)-min(40)

    def test_snps_median_over_floor_fails(self):
        reports = [self.report("pass", snps=s) for s in (500, 600, 700, 800)]
        exp = {"qc": {"min_pass_rate": 0.0}, "snps_vs_7pet": {"max_median": 400}}
        checks = backtest_report.evaluate_checks(reports, exp)
        _, _, ok = self.get_check(checks, "snps median")
        self.assertFalse(ok)

    def test_absent_thresholds_in_expect_omit_that_check(self):
        # snps_vs_7pet with only max_median (no max_iqr) must not emit an
        # IQR check at all
        reports = [self.report("pass", snps=s) for s in (1, 2, 3, 4)]
        exp = {"qc": {"min_pass_rate": 0.0}, "snps_vs_7pet": {"max_median": 400}}
        checks = backtest_report.evaluate_checks(reports, exp)
        names = [c[0] for c in checks]
        self.assertIn("snps median", names)
        self.assertNotIn("snps spread (IQR, n=4)", names)

    def test_no_reports_at_all_gives_zero_pass_rate_not_crash(self):
        checks = backtest_report.evaluate_checks([], {"qc": {"min_pass_rate": 0.0}})
        self.assertEqual(len(checks), 1)
        obs, want, ok = checks[0][1], checks[0][2], checks[0][3]
        self.assertIn("(0/0)", obs)
        self.assertTrue(ok)  # 0.0 >= 0.0


if __name__ == "__main__":
    unittest.main()
