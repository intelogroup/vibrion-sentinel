import datetime
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))
import watch_rules as wr  # noqa: E402

FIXTURE = Path(__file__).parent.parent / "fixtures" / "backtest_rows.json"


def row(acc, date, snps=40, qc="pass", priority="haiti", ctx=True, frac=0.99,
        o1=True):
    call = "present" if ctx else "absent"
    o1_call = "present" if o1 else "absent"
    return {
        "accession": acc,
        "collection_date": date,
        "country": "Haiti",
        "priority": priority,
        "qc_status": qc,
        "snps_vs_7pet": snps,
        "v_cholerae_fraction": frac,
        "surveillance_loci": {
            "ctxA": {"call": call},
            "ctxB": {"call": call},
            "wbeT": {"call": o1_call},
            "rfbV": {"call": o1_call},
        },
    }


class TestParseDate(unittest.TestCase):
    def test_full_iso(self):
        self.assertEqual(wr.parse_date("2022-10-03"), datetime.date(2022, 10, 3))

    def test_none_and_empty(self):
        self.assertIsNone(wr.parse_date(None))
        self.assertIsNone(wr.parse_date(""))

    def test_malformed_never_raises(self):
        # Regression guard: dates are parsed, never string-compared.
        # '2018-06-23' > '2018' lexicographically, which once silently
        # excluded full 2018 dates from a <= '2018' filter.
        self.assertIsNone(wr.parse_date("not-a-date"))
        self.assertIsNone(wr.parse_date("2022"))  # year-only -> excluded


class TestIsToxigenicO1(unittest.TestCase):
    def test_positive(self):
        self.assertTrue(wr.is_toxigenic_o1(row("A", "2022-01-01")))

    def test_qc_fail_excluded(self):
        self.assertFalse(wr.is_toxigenic_o1(row("A", "2022-01-01", qc="fail")))

    def test_provisional_qc_eligible_for_alerting(self):
        # Tiered gates: provisional is usable for early alerting (flagged).
        self.assertTrue(wr.is_toxigenic_o1(row("A", "2022-01-01", qc="provisional")))

    def test_missing_toxin_excluded(self):
        self.assertFalse(wr.is_toxigenic_o1(row("A", "2022-01-01", ctx=False)))

    def test_non_o1_toxigenic_excluded(self):
        # O139-like: toxin genes present, O1 antigen genes absent --
        # must not be reported as O1.
        self.assertFalse(wr.is_toxigenic_o1(row("A", "2022-01-01", o1=False)))

    def test_missing_o1_loci_excluded(self):
        # Panel without wbeT/rfbV cannot evidence O1: conservative exclusion.
        r = row("A", "2022-01-01")
        del r["surveillance_loci"]["wbeT"]
        del r["surveillance_loci"]["rfbV"]
        self.assertFalse(wr.is_toxigenic_o1(r))

    def test_purity_boundary(self):
        self.assertTrue(wr.is_toxigenic_o1(row("A", "2022-01-01", frac=0.9)))
        self.assertFalse(wr.is_toxigenic_o1(row("A", "2022-01-01", frac=0.899)))


class TestHaitiScope(unittest.TestCase):
    def test_priority_haiti(self):
        self.assertTrue(wr.is_haiti(row("A", "2022-01-01", priority="haiti")))

    def test_country_haiti(self):
        r = row("A", "2022-01-01", priority="global")
        r["country"] = "Haiti: Port-au-Prince"
        self.assertTrue(wr.is_haiti(r))

    def test_global_excluded(self):
        r = row("A", "2022-01-01", priority="global")
        r["country"] = "Bangladesh"
        self.assertFalse(wr.is_haiti(r))


class TestRecurrence(unittest.TestCase):
    def test_fires_after_quiet_gap(self):
        rows = [row("OLD", "2018-06-30", snps=41), row("NEW", "2022-10-03", snps=42)]
        f, skipped = wr.detect_recurrence(rows, n_years=2)
        self.assertIsNotNone(f)
        self.assertEqual(f["accession"], "NEW")
        self.assertGreaterEqual(f["gap_days"], 2 * 365)

    def test_no_fire_when_continuous(self):
        rows = [row("A", "2018-01-01"), row("B", "2019-06-01"), row("C", "2020-12-01")]
        f, _ = wr.detect_recurrence(rows, n_years=2)
        self.assertIsNone(f)

    def test_no_fire_on_first_ever_detection(self):
        # Cold start: no history means no measurable gap.
        f, _ = wr.detect_recurrence([row("A", "2022-10-03")], n_years=2)
        self.assertIsNone(f)

    def test_boundary_exactly_two_years(self):
        rows = [row("A", "2020-01-01"), row("B", "2021-12-31")]
        f, _ = wr.detect_recurrence(rows, n_years=2)
        self.assertIsNotNone(f)  # 730 days >= 730
        rows = [row("A", "2020-01-01"), row("B", "2021-12-30")]
        f, _ = wr.detect_recurrence(rows, n_years=2)
        self.assertIsNone(f)  # 729 days < 730

    def test_undated_skipped_not_crashed(self):
        rows = [row("OLD", "2018-06-30"), row("NODATE", None), row("NEW", "2022-10-03")]
        f, skipped = wr.detect_recurrence(rows, n_years=2)
        self.assertIsNotNone(f)
        self.assertEqual(skipped, 1)

    def test_non_toxigenic_does_not_reset_gap(self):
        rows = [row("OLD", "2018-06-30"), row("MID", "2020-01-01", ctx=False),
                row("NEW", "2022-10-03")]
        f, _ = wr.detect_recurrence(rows, n_years=2)
        self.assertIsNotNone(f)
        self.assertEqual(f["previous_detection"], "2018-06-30")

    def test_fires_once_per_event(self):
        rows = [row("OLD", "2018-06-30"), row("N1", "2022-10-03"), row("N2", "2022-10-04")]
        f, _ = wr.detect_recurrence(rows, n_years=2)
        self.assertEqual(f["accession"], "N1")


class TestClusters(unittest.TestCase):
    def test_pair_within_tolerances(self):
        rows = [row("A", "2022-10-03", snps=42), row("B", "2022-10-04", snps=45)]
        clusters, _ = wr.detect_clusters(rows)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["accessions"], ["A", "B"])

    def test_snp_diff_boundary(self):
        rows = [row("A", "2022-10-03", snps=40), row("B", "2022-10-04", snps=45)]
        clusters, _ = wr.detect_clusters(rows, snp_tol=5)
        self.assertEqual(len(clusters), 1)
        rows = [row("A", "2022-10-03", snps=40), row("B", "2022-10-04", snps=46)]
        clusters, _ = wr.detect_clusters(rows, snp_tol=5)
        self.assertEqual(len(clusters), 0)

    def test_day_boundary(self):
        rows = [row("A", "2022-10-01", snps=42), row("B", "2022-10-15", snps=43)]
        clusters, _ = wr.detect_clusters(rows, day_tol=14)
        self.assertEqual(len(clusters), 1)
        rows = [row("A", "2022-10-01", snps=42), row("B", "2022-10-16", snps=43)]
        clusters, _ = wr.detect_clusters(rows, day_tol=14)
        self.assertEqual(len(clusters), 0)

    def test_chain_links_via_union_find(self):
        rows = [row("A", "2022-10-01", snps=40), row("B", "2022-10-05", snps=44),
                row("C", "2022-10-09", snps=48)]
        clusters, _ = wr.detect_clusters(rows)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["accessions"], ["A", "B", "C"])

    def test_qc_fail_excluded(self):
        rows = [row("A", "2022-10-03", snps=42), row("B", "2022-10-04", snps=43, qc="fail")]
        clusters, _ = wr.detect_clusters(rows)
        self.assertEqual(len(clusters), 0)

    def test_singleton_no_cluster(self):
        clusters, _ = wr.detect_clusters([row("A", "2022-10-03")])
        self.assertEqual(len(clusters), 0)


class TestBacktestAcceptance(unittest.TestCase):
    """The plan's acceptance criteria, run against the real 69-row fixture."""

    @classmethod
    def setUpClass(cls):
        cls.rows = json.loads(FIXTURE.read_text())

    def test_recurrence_fires_once_on_first_2022_sample(self):
        f, skipped = wr.detect_recurrence(self.rows, n_years=2)
        self.assertIsNotNone(f)
        self.assertEqual(f["accession"], "SRR22265444")
        self.assertEqual(f["collection_date"], "2022-10-03")
        self.assertGreater(f["gap_years"], 4.0)
        self.assertEqual(f["previous_detection"], "2018-06-30")

    def test_no_recurrence_false_fire_in_anchor(self):
        anchor = [r for r in self.rows
                  if (r.get("collection_date") or "") <= "2018-12-31"
                  and r.get("collection_date")]
        f, _ = wr.detect_recurrence(anchor, n_years=2)
        self.assertIsNone(f)

    def test_cluster_fires_on_second_2022_sample(self):
        clusters, skipped = wr.detect_clusters(self.rows)
        oct4 = [c for c in clusters if "SRR22265443" in c["accessions"]]
        self.assertTrue(oct4)
        self.assertIn("SRR22265444", oct4[0]["accessions"])
        self.assertEqual(skipped, 2)  # the two undated samples

    def test_2022_cluster_is_tight(self):
        clusters, _ = wr.detect_clusters(self.rows)
        big = max(clusters, key=lambda c: c["n"])
        self.assertGreaterEqual(big["n"], 20)
        self.assertLessEqual(big["snps_max"] - big["snps_min"], 5)


if __name__ == "__main__":
    unittest.main()
