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


class TestVirulenceMarkers(unittest.TestCase):
    """The five-marker line is a rollup of loci calls already produced:
    ctxAB (both genes required), tcpA, rtxA, hlyA, T6SS via hcp."""

    def test_toxigenic_o1_like_profile(self):
        m = report_lib.virulence_markers(loci(
            ctxA="present", ctxB="present", tcpA="present",
            rtxA="present", hlyA="present", hcp="present"))
        self.assertEqual(m["ctxAB"]["call"], "present")
        self.assertEqual(m["ctxAB"]["loci"], ["ctxA", "ctxB"])
        self.assertEqual(m["tcpA"]["call"], "present")
        self.assertEqual(m["rtxA"]["call"], "present")
        self.assertEqual(m["hlyA"]["call"], "present")
        self.assertEqual(m["t6ss"]["call"], "present")
        self.assertEqual(m["t6ss"]["locus"], "hcp")
        self.assertIs(m["toxigenic"], True)

    def test_non_toxigenic_novc_profile(self):
        # HL2601-shaped: no ctxAB, no tcpA, but rtxA/hlyA/T6SS present.
        m = report_lib.virulence_markers(loci(
            ctxA="absent", ctxB="absent", tcpA="absent",
            rtxA="present", hlyA="present", hcp="present"))
        self.assertEqual(m["ctxAB"]["call"], "absent")
        self.assertIs(m["toxigenic"], False)
        self.assertEqual(m["tcpA"]["call"], "absent")
        self.assertEqual(m["rtxA"]["call"], "present")
        self.assertEqual(m["hlyA"]["call"], "present")
        self.assertEqual(m["t6ss"]["call"], "present")

    def test_ctxab_needs_both_genes(self):
        m = report_lib.virulence_markers(loci(ctxA="present", ctxB="absent"))
        self.assertEqual(m["ctxAB"]["call"], "partial")
        self.assertIsNone(m["toxigenic"])

    def test_ctxab_partial_when_either_partial(self):
        m = report_lib.virulence_markers(loci(ctxA="partial", ctxB="present"))
        self.assertEqual(m["ctxAB"]["call"], "partial")
        self.assertIsNone(m["toxigenic"])

    def test_missing_loci_are_not_assessed_not_absent(self):
        m = report_lib.virulence_markers(loci(ctxA="present"))
        self.assertEqual(m["ctxAB"]["call"], "not_assessed")
        self.assertIsNone(m["toxigenic"])
        self.assertEqual(m["tcpA"]["call"], "not_assessed")
        self.assertEqual(m["t6ss"]["call"], "not_assessed")

    def test_empty_panel(self):
        m = report_lib.virulence_markers({})
        for key in ("ctxAB", "tcpA", "rtxA", "hlyA", "t6ss"):
            self.assertEqual(m[key]["call"], "not_assessed")
        self.assertIsNone(m["toxigenic"])

    def test_none_panel(self):
        m = report_lib.virulence_markers(None)
        self.assertEqual(m["ctxAB"]["call"], "not_assessed")

    def test_note_carries_ruo_posture(self):
        m = report_lib.virulence_markers(loci())
        self.assertIn("not a clinical or virulence prediction", m["note"])


if __name__ == "__main__":
    unittest.main()
