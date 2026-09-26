"""Tests the depth-masking awk script directly by piping synthetic
`samtools depth -aa` output through it and checking the merged BED,
without going through Snakemake or samtools at all.
"""

import shutil
import subprocess
import unittest
from pathlib import Path

AWK_SCRIPT = Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite" / "scripts" / "depth_mask.awk"


@unittest.skipUnless(shutil.which("awk"), "awk not available")
class TestDepthMaskAwk(unittest.TestCase):
    def run_mask(self, depth_lines, min_depth=10):
        """depth_lines: list of (chrom, pos, depth) 1-based tuples."""
        stdin = "\n".join(f"{c}\t{p}\t{d}" for c, p, d in depth_lines) + "\n"
        result = subprocess.run(
            ["awk", "-v", f"m={min_depth}", "-f", str(AWK_SCRIPT)],
            input=stdin, capture_output=True, text=True, check=True,
        )
        rows = [r.split("\t") for r in result.stdout.splitlines() if r]
        return [(c, int(s), int(e)) for c, s, e in rows]

    def test_all_positions_above_threshold_masks_nothing(self):
        depths = [("chr1", p, 20) for p in range(1, 11)]
        self.assertEqual(self.run_mask(depths), [])

    def test_single_low_depth_run_becomes_one_interval(self):
        depths = [("chr1", p, 20) for p in range(1, 6)]
        depths += [("chr1", p, 3) for p in range(6, 11)]   # positions 6-10 low
        depths += [("chr1", p, 20) for p in range(11, 16)]
        # 0-based, end-exclusive: positions 6..10 (1-based) -> [5, 10)
        self.assertEqual(self.run_mask(depths), [("chr1", 5, 10)])

    def test_two_separate_low_runs_stay_separate(self):
        depths = (
            [("chr1", p, 20) for p in (1, 2)]
            + [("chr1", 3, 2)]
            + [("chr1", p, 20) for p in (4, 5)]
            + [("chr1", 6, 1)]
            + [("chr1", p, 20) for p in (7, 8)]
        )
        self.assertEqual(self.run_mask(depths), [("chr1", 2, 3), ("chr1", 5, 6)])

    def test_low_depth_run_extending_to_end_of_contig(self):
        # this exercises the END block: a run with no higher-depth
        # position after it to trigger the flush inside the main loop
        depths = [("chr1", p, 20) for p in (1, 2)] + [("chr1", p, 3) for p in (3, 4, 5)]
        self.assertEqual(self.run_mask(depths), [("chr1", 2, 5)])

    def test_entire_contig_below_threshold(self):
        # e.g. a contig with zero reads mapped at all -- samtools depth -aa
        # still emits every position with depth 0
        depths = [("chr2", p, 0) for p in range(1, 6)]
        self.assertEqual(self.run_mask(depths), [("chr2", 0, 5)])

    def test_low_run_does_not_merge_across_contigs(self):
        # adjacent low-depth positions on different contigs must not be
        # merged into one interval just because they're consecutive lines
        depths = [("chr1", 100, 2)] + [("chr2", 1, 2)]
        self.assertEqual(self.run_mask(depths), [("chr1", 99, 100), ("chr2", 0, 1)])

    def test_threshold_is_strict_less_than(self):
        # depth exactly at the threshold is NOT masked (mirrors
        # report_lib.evaluate_qc's >= convention: >= threshold is fine)
        depths = [("chr1", 1, 10)]
        self.assertEqual(self.run_mask(depths, min_depth=10), [])
        depths = [("chr1", 1, 9)]
        self.assertEqual(self.run_mask(depths, min_depth=10), [("chr1", 0, 1)])

    def test_empty_input_produces_empty_output(self):
        self.assertEqual(self.run_mask([]), [])


if __name__ == "__main__":
    unittest.main()
