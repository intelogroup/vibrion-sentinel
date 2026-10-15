import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import report_lib  # noqa: E402


def loci(**calls):
    return {name: {"call": call, "present": call == "present",
                   "mean_depth": 50.0, "breadth_pct": 100.0}
            for name, call in calls.items()}


class TestStrainPattern(unittest.TestCase):
    """One rule, one label: QC pass (species gate included) + ctxA/tcpA
    absent + SNP distance to the 7PET reference at/above the operational
    floor => the non-pandemic V. cholerae pattern. Every leg is a signal
    the lite report already produces."""

    def test_hl2601_shaped_signals_match(self):
        p = report_lib.strain_pattern(
            "pass", loci(ctxA="absent", tcpA="absent", rtxA="present"),
            25000)
        self.assertTrue(p["matched"])
        self.assertEqual(p["name"], "non_pandemic_v_cholerae")
        self.assertEqual(p["label"], "Non-pandemic V. cholerae pattern")
        self.assertTrue(all(p["signals"][k] for k in (
            "qc_pass", "ctxA_absent", "tcpA_absent", "snps_far_from_7pet")))

    def test_7pet_shaped_signals_do_not_match(self):
        # Toxigenic O1 near the reference: three legs fail at once.
        p = report_lib.strain_pattern(
            "pass", loci(ctxA="present", tcpA="present"), 120)
        self.assertFalse(p["matched"])
        self.assertIsNone(p["name"])
        self.assertIsNone(p["label"])
        self.assertFalse(p["signals"]["ctxA_absent"])
        self.assertFalse(p["signals"]["snps_far_from_7pet"])

    def test_each_leg_is_required(self):
        base = loci(ctxA="absent", tcpA="absent")
        self.assertFalse(report_lib.strain_pattern(
            "fail", base, 25000)["matched"])
        self.assertFalse(report_lib.strain_pattern(
            "provisional", base, 25000)["matched"])
        self.assertFalse(report_lib.strain_pattern(
            "pass", loci(ctxA="partial", tcpA="absent"), 25000)["matched"])
        self.assertFalse(report_lib.strain_pattern(
            "pass", loci(ctxA="absent", tcpA="present"), 25000)["matched"])
        self.assertFalse(report_lib.strain_pattern(
            "pass", base, 999)["matched"])

    def test_snp_floor_boundary(self):
        base = loci(ctxA="absent", tcpA="absent")
        at = report_lib.strain_pattern("pass", base, 1000)
        self.assertTrue(at["signals"]["snps_far_from_7pet"])
        self.assertTrue(at["matched"])
        below = report_lib.strain_pattern("pass", base, 999)
        self.assertFalse(below["signals"]["snps_far_from_7pet"])

    def test_missing_snp_count_never_matches(self):
        p = report_lib.strain_pattern(
            "pass", loci(ctxA="absent", tcpA="absent"), None)
        self.assertFalse(p["signals"]["snps_far_from_7pet"])
        self.assertFalse(p["matched"])
        self.assertIsNone(p["signals"]["snps_vs_7pet"])

    def test_missing_loci_never_match(self):
        p = report_lib.strain_pattern("pass", {}, 25000)
        self.assertFalse(p["matched"])
        self.assertFalse(p["signals"]["ctxA_absent"])

    def test_note_is_ruo_safe_and_scoped(self):
        p = report_lib.strain_pattern("pass", loci(), None)
        self.assertIn("not a serogroup", p["note"])
        self.assertNotIn("virulence", p["note"].lower())


if __name__ == "__main__":
    unittest.main()
