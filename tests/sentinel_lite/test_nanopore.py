"""Phase 0 (Nanopore branch) unit tests.

Covers the pure-Python additions only: platform/tier config validation
(including fail-closed behavior), platform QC thresholds, assembly-block
parsing/rendering, and the NanoPlot QC adapter. Nothing here needs the
bioinformatics toolchain -- CI runs `python -m unittest discover tests`.
"""

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import nano_qc_lib  # noqa: E402
import report_lib  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


class TestValidatePlatformConfig(unittest.TestCase):
    def test_defaults_are_illumina_lite(self):
        platform, basecaller, tier = report_lib.validate_platform_config({})
        self.assertEqual((platform, basecaller, tier), ("illumina", None, "lite"))

    def test_nanopore_requires_basecaller_model(self):
        # Fail-closed: no basecaller_model -> refuse, don't guess.
        with self.assertRaises(ValueError):
            report_lib.validate_platform_config({"platform": "nanopore"})
        with self.assertRaises(ValueError):
            report_lib.validate_platform_config(
                {"platform": "nanopore", "tier": "assembly"}
            )

    def test_nanopore_with_basecaller_ok(self):
        platform, basecaller, tier = report_lib.validate_platform_config(
            {"platform": "nanopore", "basecaller_model": "hac", "tier": "assembly"}
        )
        self.assertEqual((platform, basecaller, tier), ("nanopore", "hac", "assembly"))

    def test_bad_basecaller_rejected(self):
        with self.assertRaises(ValueError):
            report_lib.validate_platform_config(
                {"platform": "nanopore", "basecaller_model": "ultra"}
            )

    def test_bad_platform_rejected(self):
        with self.assertRaises(ValueError):
            report_lib.validate_platform_config({"platform": "pacbio"})

    def test_bad_tier_rejected(self):
        with self.assertRaises(ValueError):
            report_lib.validate_platform_config({"tier": "ultra"})

    def test_assembly_tier_is_nanopore_only(self):
        with self.assertRaises(ValueError):
            report_lib.validate_platform_config(
                {"platform": "illumina", "tier": "assembly"}
            )


class TestPlatformQcThresholds(unittest.TestCase):
    def test_nanopore_defaults(self):
        t = report_lib.platform_qc_thresholds({}, "nanopore")
        self.assertEqual(t["min_site_depth"], 15)
        self.assertEqual(t["qc_min_mean_depth"], 30)
        self.assertEqual(t["qc_min_called_pct"], 85)

    def test_illumina_tiered_defaults(self):
        # Tiered QC gates (2026-10-07): provisional cutoffs grounded in
        # FWD-AMR-RefLabCap / PulseNet; hard floors from Kenya 2022-23.
        t = report_lib.platform_qc_thresholds({}, "illumina")
        self.assertEqual(t["min_site_depth"], 10)
        self.assertEqual(t["qc_min_mean_depth"], 30)
        self.assertEqual(t["qc_min_called_pct"], 90)
        self.assertEqual(t["qc_min_species_purity"], 99)
        self.assertEqual(t["qc_hard_min_called_pct"], 50)
        self.assertEqual(t["qc_hard_min_species_purity"], 95)

    def test_explicit_config_wins(self):
        t = report_lib.platform_qc_thresholds({"qc_min_mean_depth": 50}, "nanopore")
        self.assertEqual(t["qc_min_mean_depth"], 50)
        self.assertEqual(t["min_site_depth"], 15)  # default still applies


class TestModelMaps(unittest.TestCase):
    def test_medaka_models(self):
        self.assertEqual(
            report_lib.medaka_model_for("hac"), "r1041_e82_400bps_hac_g615"
        )
        self.assertEqual(
            report_lib.medaka_model_for("sup"), "r1041_e82_400bps_sup_g615"
        )
        self.assertEqual(
            report_lib.medaka_model_for("fast"), "r1041_e82_400bps_fast_g615"
        )

    def test_flye_presets(self):
        self.assertEqual(report_lib.flye_preset_for("hac"), "--nano-hq")
        self.assertEqual(report_lib.flye_preset_for("sup"), "--nano-hq")
        self.assertEqual(report_lib.flye_preset_for("fast"), "--nano-raw")


class TestAssemblyQc(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _qc(self, obj):
        p = Path(self.tmp.name) / "assembly_qc.json"
        import json

        write(p, json.dumps(obj))
        return report_lib.parse_assembly_qc(p)

    def test_gate_pass(self):
        stats = self._qc({"length_bp": 4038120, "contigs": 2, "n50_bp": 2961102})
        status, reasons = report_lib.assembly_qc_gate(stats)
        self.assertEqual(status, "pass")
        self.assertEqual(reasons, [])

    def test_gate_fails_short(self):
        stats = self._qc({"length_bp": 2000000, "contigs": 2, "n50_bp": 2961102})
        status, reasons = report_lib.assembly_qc_gate(stats)
        self.assertEqual(status, "fail")
        self.assertTrue(any("length" in r for r in reasons))

    def test_gate_fails_fragmented(self):
        stats = self._qc({"length_bp": 4038120, "contigs": 120, "n50_bp": 50000})
        status, reasons = report_lib.assembly_qc_gate(stats)
        self.assertEqual(status, "fail")
        self.assertEqual(len(reasons), 2)  # contigs + N50

    def test_missing_file_fails_gate(self):
        stats = report_lib.parse_assembly_qc(Path(self.tmp.name) / "nope.json")
        status, _ = report_lib.assembly_qc_gate(stats)
        self.assertEqual(status, "fail")


class TestAssemblyBlock(unittest.TestCase):
    def test_pass_block_carries_calls_and_disclaimer(self):
        stats = {"length_bp": 4038120, "contigs": 2, "n50_bp": 2961102}
        block = report_lib.build_assembly_block(
            stats, "69", "T12", ["sul2", "floR", "sul2"], "pass"
        )
        self.assertEqual(block["mlst_st"], "69")
        self.assertEqual(block["vibecheck_lineage"], "T12")
        self.assertEqual(block["amr_genes"], ["floR", "sul2"])  # deduped + sorted
        self.assertIn("Not a susceptibility prediction", block["amr_note"])

    # -- mlst_call tri-state: raw TSV state x gate status ----------------
    # mlst_st null is overloaded (unassigned profile / failed run / gate
    # suppression); mlst_call must never be read off the ST value alone.

    def _result(self, st, alleles, row_present=True, scheme="vcholerae"):
        return {"st": st, "scheme": scheme, "alleles": alleles,
                "row_present": row_present}

    def test_mlst_call_assigned(self):
        stats = {"length_bp": 4038120, "contigs": 2, "n50_bp": 2961102}
        block = report_lib.build_assembly_block(
            stats, "69", "T12", [], "pass",
            mlst_result=self._result("69", ["7", "5", "3"]))
        self.assertEqual(block["mlst_call"], "assigned")
        self.assertEqual(block["mlst_st"], "69")
        self.assertEqual(block["mlst_scheme"], "vcholerae")
        self.assertEqual(block["mlst_alleles"], ["7", "5", "3"])

    def test_mlst_call_unassigned_on_dash_with_gate_pass(self):
        # A real tseemann row with ST '-': loci typed, profile matched
        # no known ST (the HL2601 shape). The raw profile is preserved.
        stats = {"length_bp": 4038120, "contigs": 2, "n50_bp": 2961102}
        profile = ["51", "41", "4", "79", "72", "1", "83"]
        block = report_lib.build_assembly_block(
            stats, None, None, [], "pass",
            mlst_result=self._result(None, profile))
        self.assertEqual(block["mlst_call"], "unassigned")
        self.assertIsNone(block["mlst_st"])
        self.assertEqual(block["mlst_alleles"], profile)

    def test_mlst_call_gate_fail_is_not_interpretable_even_with_st(self):
        stats = {"length_bp": 2000000, "contigs": 120, "n50_bp": 50000}
        block = report_lib.build_assembly_block(
            stats, "69", "T12", [], "fail",
            mlst_result=self._result("69", ["7", "5", "3"]))
        self.assertEqual(block["mlst_call"], "not_interpretable")
        self.assertIsNone(block["mlst_st"])
        # The raw profile survives as evidence; it is not a call.
        self.assertEqual(block["mlst_alleles"], ["7", "5", "3"])

    def test_mlst_call_missing_tsv_is_not_unassigned(self):
        # No usable TSV row (missing file / failed-run fallback) must
        # read as not_interpretable, distinctly from unassigned.
        stats = {"length_bp": 4038120, "contigs": 2, "n50_bp": 2961102}
        block = report_lib.build_assembly_block(
            stats, None, None, [], "pass",
            mlst_result=self._result(None, [], row_present=False,
                                     scheme=None))
        self.assertEqual(block["mlst_call"], "not_interpretable")
        self.assertNotEqual(block["mlst_call"], "unassigned")
        self.assertEqual(block["mlst_alleles"], [])

    def test_mlst_call_legacy_caller_without_result(self):
        stats = {"length_bp": 4038120, "contigs": 2, "n50_bp": 2961102}
        block = report_lib.build_assembly_block(stats, "69", "T12", [], "pass")
        self.assertEqual(block["mlst_call"], "assigned")
        block = report_lib.build_assembly_block(stats, None, None, [], "pass")
        self.assertEqual(block["mlst_call"], "not_interpretable")

    def test_fail_gate_suppresses_calls(self):
        # A fragmented assembly must never yield lineage or AMR calls.
        stats = {"length_bp": 2000000, "contigs": 120, "n50_bp": 50000}
        block = report_lib.build_assembly_block(
            stats, "69", "T12", ["sul2"], "fail"
        )
        self.assertIsNone(block["mlst_st"])
        self.assertEqual(block["mlst_call"], "not_interpretable")
        self.assertIsNone(block["vibecheck_lineage"])
        self.assertEqual(block["amr_genes"], [])
        # the disclaimer is still present even when calls are suppressed
        self.assertIn("Not a susceptibility prediction", block["amr_note"])


class TestParsers(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_parse_mlst(self):
        p = Path(self.tmp.name) / "mlst.tsv"
        write(p, "polished.fasta\tvcholerae\t69\ttnaA_1\ttcpA_2\n")
        self.assertEqual(report_lib.parse_mlst(p), "69")

    def test_parse_mlst_dash_is_none(self):
        p = Path(self.tmp.name) / "mlst.tsv"
        write(p, "polished.fasta\t-\t-\n")
        self.assertIsNone(report_lib.parse_mlst(p))

    def test_parse_mlst_result_assigned_row(self):
        p = Path(self.tmp.name) / "mlst.tsv"
        write(p, "polished.fasta\tvcholerae\t69\t7\t5\t3\n")
        r = report_lib.parse_mlst_result(p)
        self.assertEqual(r, {"st": "69", "scheme": "vcholerae",
                             "alleles": ["7", "5", "3"], "row_present": True})

    def test_parse_mlst_result_dash_row_preserves_profile(self):
        # Raw state only: interpretation (unassigned vs not) happens in
        # build_assembly_block, crossed with the gate. The full allele
        # profile survives verbatim.
        p = Path(self.tmp.name) / "mlst.tsv"
        write(p, "HL2601.fasta\tvcholerae\t-\t51\t41\t4\t79\t72\t1\t83\n")
        r = report_lib.parse_mlst_result(p)
        self.assertIsNone(r["st"])
        self.assertEqual(r["scheme"], "vcholerae")
        self.assertEqual(r["alleles"], ["51", "41", "4", "79", "72", "1", "83"])
        self.assertTrue(r["row_present"])

    def test_parse_mlst_result_partial_profile_kept_raw(self):
        # '-' allele entries are data (incomplete recovery), not noise:
        # they stay in the preserved profile.
        p = Path(self.tmp.name) / "mlst.tsv"
        write(p, "x.fasta\tvcholerae\t-\t51\t-\t4\n")
        r = report_lib.parse_mlst_result(p)
        self.assertEqual(r["alleles"], ["51", "-", "4"])
        self.assertTrue(r["row_present"])

    def test_parse_mlst_result_fallback_row_is_no_row(self):
        # The Snakefile's failure fallback ("<file>\t-\t-\n", written when
        # the mlst run produced nothing) has no scheme: not a typing row.
        p = Path(self.tmp.name) / "mlst.tsv"
        write(p, "polished.fasta\t-\t-\n")
        r = report_lib.parse_mlst_result(p)
        self.assertFalse(r["row_present"])
        self.assertIsNone(r["st"])
        self.assertEqual(r["alleles"], [])

    def test_parse_mlst_result_missing_and_empty_files(self):
        r = report_lib.parse_mlst_result(Path(self.tmp.name) / "nope.tsv")
        self.assertFalse(r["row_present"])
        p = Path(self.tmp.name) / "empty.tsv"
        write(p, "")
        r = report_lib.parse_mlst_result(p)
        self.assertFalse(r["row_present"])

    def test_parse_vibecheck_json(self):
        p = Path(self.tmp.name) / "vibecheck.txt"
        write(p, '{"lineage": "T12", "confidence": 0.99}')
        self.assertEqual(report_lib.parse_vibecheck(p), "T12")

    def test_parse_vibecheck_text_fallback(self):
        p = Path(self.tmp.name) / "vibecheck.txt"
        write(p, "assigned lineage: T5 (bootstrap 87)")
        self.assertEqual(report_lib.parse_vibecheck(p), "T5")

    def test_parse_vibecheck_none(self):
        p = Path(self.tmp.name) / "vibecheck.txt"
        write(p, "vibecheck: no lineage assigned")
        self.assertIsNone(report_lib.parse_vibecheck(p))

    def test_parse_amrfinder(self):
        p = Path(self.tmp.name) / "amrfinder.tsv"
        write(
            p,
            "Gene symbol\tSequence name\tScope\n"
            "sul2\tsulfonamide\tcore\n"
            "floR\tphenicol\tcore\n"
            "sul2\tsulfonamide\tcore\n",
        )
        self.assertEqual(report_lib.parse_amrfinder(p), ["floR", "sul2"])

    def test_parse_amrfinder_header_only(self):
        p = Path(self.tmp.name) / "amrfinder.tsv"
        write(p, "Gene symbol\tSequence name\n")
        self.assertEqual(report_lib.parse_amrfinder(p), [])

    def test_amrfinder_db_version(self):
        out = "AMRFinderPlus version 3.12.8\nDatabase version: 2024-10-22.1\n"
        self.assertEqual(report_lib.amrfinder_db_version(out), "2024-10-22.1")
        self.assertEqual(report_lib.amrfinder_db_version("garbage"), "unknown")


class TestCollectToolVersions(unittest.TestCase):
    def test_fake_runner(self):
        def fake(argv):
            class R:
                stdout = "flye 2.9.3\n"
                stderr = ""

            return R()

        versions = report_lib.collect_tool_versions([("flye", ["flye", "--version"])], runner=fake)
        self.assertEqual(versions["flye"], "flye 2.9.3")

    def test_failing_tool_is_unknown_not_fatal(self):
        def boom(argv):
            raise FileNotFoundError("nope")

        versions = report_lib.collect_tool_versions([("nope", ["nope"])], runner=boom)
        self.assertEqual(versions["nope"], "unknown")


class TestNanoQcLib(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_parse_nanostats(self):
        p = Path(self.tmp.name) / "NanoStats.txt"
        write(p, "number_of_reads:\t12345\nn50:\t6789\n")
        stats = nano_qc_lib.parse_nanostats(p)
        self.assertEqual(stats["number_of_reads"], "12345")

    def test_write_qc_summary_shape(self):
        ns = Path(self.tmp.name) / "NanoStats.txt"
        write(ns, "number_of_reads:\t100\n")
        fq = Path(self.tmp.name) / "cleaned.fastq.gz"
        with gzip.open(fq, "wt") as f:
            for i in range(80):
                f.write(f"@r{i}\nACGT\n+\nIIII\n")
        out = Path(self.tmp.name) / "qc_summary.json"
        summary = nano_qc_lib.write_qc_summary(ns, fq, out)
        # fastp-shaped: the report rule reads these exact keys
        self.assertEqual(summary["summary"]["before_filtering"]["total_reads"], 100)
        self.assertEqual(summary["summary"]["after_filtering"]["total_reads"], 80)
        import json

        self.assertEqual(json.loads(out.read_text()), summary)


if __name__ == "__main__":
    unittest.main()
