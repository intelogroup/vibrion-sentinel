import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import report_lib  # noqa: E402


def loci(**calls):
    """surveillance_loci-shaped dict from {locus: call} pairs."""
    return {name: {"call": call, "present": call == "present",
                   "mean_depth": 50.0, "breadth_pct": 100.0}
            for name, call in calls.items()}


class TestReferenceMarkerPanel(unittest.TestCase):
    """The panel is a VIEW: it groups four loci calls that already exist
    in surveillance_loci (ctxA, tcpA, rtxA, hlyA) under reference-
    similarity semantics. It computes nothing new and screens nothing."""

    def test_groups_exactly_the_four_loci(self):
        panel = report_lib.reference_marker_panel(loci(
            ctxA="present", tcpA="present", rtxA="absent", hlyA="partial",
            ctxB="present", wbeT="absent", tcpP="present"))
        self.assertEqual(panel["markers"], {
            "ctxA": "present", "tcpA": "present",
            "rtxA": "absent", "hlyA": "partial",
        })

    def test_calls_pass_through_unchanged(self):
        panel = report_lib.reference_marker_panel(loci(
            ctxA="absent", tcpA="absent", rtxA="present", hlyA="present"))
        self.assertEqual(panel["markers"]["ctxA"], "absent")
        self.assertEqual(panel["markers"]["rtxA"], "present")

    def test_missing_locus_is_not_assessed(self):
        panel = report_lib.reference_marker_panel(loci(ctxA="present"))
        self.assertEqual(panel["markers"]["tcpA"], "not_assessed")
        self.assertEqual(panel["markers"]["hlyA"], "not_assessed")

    def test_empty_and_none_panel(self):
        for empty in ({}, None):
            panel = report_lib.reference_marker_panel(empty)
            self.assertEqual(set(panel["markers"].values()), {"not_assessed"})

    def test_semantics_labelled_reference_similarity(self):
        # The BED's semantics are mapping breadth vs the O1 reference
        # allele; the block must carry that label and caveat, so a
        # divergent non-O1 allele scoring 'absent' is never read as
        # gene absence.
        panel = report_lib.reference_marker_panel(loci())
        self.assertEqual(panel["basis"], "reference_similarity:2010EL-1786")
        self.assertIn("mapping breadth", panel["note"])
        self.assertIn("divergent non-O1 allele", panel["note"])

    def test_no_rollup_verdicts(self):
        # No combined ctxAB call, no toxigenic flag, no T6SS entry:
        # grouping only.
        panel = report_lib.reference_marker_panel(loci(
            ctxA="present", ctxB="present", tcpA="present",
            rtxA="present", hlyA="present"))
        self.assertEqual(set(panel), {"markers", "basis", "note"})
        self.assertNotIn("ctxAB", panel["markers"])


if __name__ == "__main__":
    unittest.main()
