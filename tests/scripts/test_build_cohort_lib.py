import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))
import build_cohort_lib  # noqa: E402


class TestGetField(unittest.TestCase):
    def test_normal_lookup(self):
        header = ["run_accession", "collection_date"]
        row = ["SRR1", "2022-10-03"]
        self.assertEqual(build_cohort_lib.get_field(header, row, "collection_date"), "2022-10-03")

    def test_missing_column_name_returns_empty_string(self):
        header = ["run_accession"]
        row = ["SRR1"]
        self.assertEqual(build_cohort_lib.get_field(header, row, "collection_date"), "")

    def test_short_row_returns_empty_string_not_indexerror(self):
        # ENA has, in practice, returned rows shorter than the header
        header = ["run_accession", "collection_date", "country"]
        row = ["SRR1"]
        self.assertEqual(build_cohort_lib.get_field(header, row, "country"), "")


class TestSortAndLimit(unittest.TestCase):
    HEADER = ["run_accession", "collection_date"]

    def test_sorted_by_date_then_accession(self):
        rows = [["SRR3", "2022-10-01"], ["SRR1", "2022-09-01"], ["SRR2", "2022-09-01"]]
        result = build_cohort_lib.sort_and_limit(self.HEADER, rows)
        self.assertEqual([r[0] for r in result], ["SRR1", "SRR2", "SRR3"])

    def test_limit_zero_keeps_everything(self):
        rows = [["SRR1", "2022-09-01"], ["SRR2", "2022-09-02"]]
        result = build_cohort_lib.sort_and_limit(self.HEADER, rows, limit=0)
        self.assertEqual(len(result), 2)

    def test_limit_caps_after_sorting_not_before(self):
        # this is what makes a --limit N subset deterministic: it must
        # take the earliest N by date, not an arbitrary N from ENA's
        # unspecified result order
        rows = [["SRR3", "2022-10-01"], ["SRR1", "2022-08-01"], ["SRR2", "2022-09-01"]]
        result = build_cohort_lib.sort_and_limit(self.HEADER, rows, limit=2)
        self.assertEqual([r[0] for r in result], ["SRR1", "SRR2"])

    def test_empty_rows(self):
        self.assertEqual(build_cohort_lib.sort_and_limit(self.HEADER, [], limit=10), [])

    def test_missing_dates_sort_first_not_crash(self):
        # empty string sorts before any real date, which is an acceptable
        # (documented by test, not just assumed) tie-break rather than a
        # crash on a row ENA returned with no collection_date at all
        rows = [["SRR2", "2022-01-01"], ["SRR1", ""]]
        result = build_cohort_lib.sort_and_limit(self.HEADER, rows)
        self.assertEqual([r[0] for r in result], ["SRR1", "SRR2"])


class TestDateRange(unittest.TestCase):
    HEADER = ["run_accession", "collection_date"]

    def test_normal_range(self):
        rows = [["SRR1", "2022-10-03"], ["SRR2", "2022-11-16"], ["SRR3", "2022-10-20"]]
        self.assertEqual(build_cohort_lib.date_range(self.HEADER, rows), ("2022-10-03", "2022-11-16"))

    def test_no_dates_returns_none(self):
        rows = [["SRR1", ""], ["SRR2", ""]]
        self.assertIsNone(build_cohort_lib.date_range(self.HEADER, rows))

    def test_empty_rows_returns_none(self):
        self.assertIsNone(build_cohort_lib.date_range(self.HEADER, []))

    def test_mixed_missing_and_present_dates(self):
        rows = [["SRR1", ""], ["SRR2", "2022-06-01"]]
        self.assertEqual(build_cohort_lib.date_range(self.HEADER, rows), ("2022-06-01", "2022-06-01"))


class TestFormatManifest(unittest.TestCase):
    def test_header_and_metadata_lines(self):
        header = ["run_accession", "collection_date"]
        rows = [["SRR1", "2022-10-03"]]
        text = build_cohort_lib.format_manifest("backtest-2022", "tax_eq(666)", header, rows)
        lines = text.splitlines()
        self.assertEqual(lines[0], "# cohort: backtest-2022")
        self.assertEqual(lines[1], "# ena_query: tax_eq(666)")
        self.assertEqual(lines[2], "# runs: 1")
        self.assertEqual(lines[3], "run_accession\tcollection_date")
        self.assertEqual(lines[4], "SRR1\t2022-10-03")

    def test_run_count_matches_actual_row_count(self):
        header = ["run_accession"]
        rows = [["SRR1"], ["SRR2"], ["SRR3"]]
        text = build_cohort_lib.format_manifest("x", "q", header, rows)
        self.assertIn("# runs: 3", text)

    def test_empty_rows_still_produces_header(self):
        header = ["run_accession"]
        text = build_cohort_lib.format_manifest("x", "q", header, [])
        self.assertIn("# runs: 0", text)
        self.assertTrue(text.rstrip("\n").endswith("run_accession"))

    def test_ends_with_single_trailing_newline(self):
        header = ["run_accession"]
        rows = [["SRR1"]]
        text = build_cohort_lib.format_manifest("x", "q", header, rows)
        self.assertTrue(text.endswith("\n"))
        self.assertFalse(text.endswith("\n\n"))


if __name__ == "__main__":
    unittest.main()
