import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import loci_coverage_lib  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


class TestIterBedRegions(unittest.TestCase):
    """This is the off-by-one that has bitten this pipeline before: BED is
    0-based half-open, samtools region syntax is 1-based inclusive."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _bed(self, text):
        p = Path(self.tmp.name) / "loci.bed"
        write(p, text)
        return p

    def test_basic_conversion(self):
        # BED [100, 200) -> samtools region 101-200
        p = self._bed("CP003069.1\t100\t200\tctxA\t.\t+\n")
        results = list(loci_coverage_lib.iter_bed_regions(p))
        self.assertEqual(len(results), 1)
        fields, region = results[0]
        self.assertEqual(region, "CP003069.1:101-200")
        self.assertEqual(fields[3], "ctxA")

    def test_zero_based_start_at_origin(self):
        # a locus starting at position 0 -> region starts at 1, not 0
        p = self._bed("CP003069.1\t0\t50\tfirst_gene\t.\t+\n")
        _, region = next(loci_coverage_lib.iter_bed_regions(p))
        self.assertEqual(region, "CP003069.1:1-50")

    def test_single_base_interval(self):
        p = self._bed("CP003069.1\t99\t100\tsnp_site\t.\t+\n")
        _, region = next(loci_coverage_lib.iter_bed_regions(p))
        self.assertEqual(region, "CP003069.1:100-100")

    def test_comments_and_blank_lines_skipped(self):
        p = self._bed(
            "# header comment\n"
            "\n"
            "CP003069.1\t100\t200\tctxA\t.\t+\n"
            "# another comment\n"
            "CP003070.1\t50\t150\thlyA\t.\t+\n"
        )
        results = list(loci_coverage_lib.iter_bed_regions(p))
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0][1], "CP003069.1:101-200")
        self.assertEqual(results[1][1], "CP003070.1:51-150")

    def test_second_chromosome_preserved_in_region(self):
        # regression guard for the earlier bug where only chr1 was
        # considered anywhere in this pipeline's coverage logic
        p = self._bed("CP003070.1\t1010364\t1012590\thlyA\t.\t+\n")
        _, region = next(loci_coverage_lib.iter_bed_regions(p))
        self.assertTrue(region.startswith("CP003070.1:"))

    def test_empty_bed_yields_nothing(self):
        p = self._bed("")
        self.assertEqual(list(loci_coverage_lib.iter_bed_regions(p)), [])


if __name__ == "__main__":
    unittest.main()
