import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "sentinel_lite"))
import vibrio_extract_lib  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)


def read_gz(path):
    with gzip.open(path, "rt") as f:
        return f.read()


class TestBuildKeepSet(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _kraken_out(self, rows):
        p = Path(self.tmp.name) / "kraken.out"
        write(p, "\n".join(rows) + "\n")
        return p

    def test_cholerae_and_unclassified_are_kept(self):
        p = self._kraken_out([
            "C\tread1/1\t666\t150\t0:116",
            "U\tread2/1\t0\t150\t0:116",
            "C\tread3/1\t644\t150\t0:116",  # not cholerae -> dropped
        ])
        keep, n_class, n_keep = vibrio_extract_lib.build_keep_set(p)
        self.assertEqual(keep, {"read1", "read2"})
        self.assertEqual(n_class, 3)
        self.assertEqual(n_keep, 2)

    def test_read_id_without_pair_suffix(self):
        # single-end reads have no /1 or /2 at all
        p = self._kraken_out(["C\treadA\t666\t150\t0:116"])
        keep, _, _ = vibrio_extract_lib.build_keep_set(p)
        self.assertEqual(keep, {"readA"})

    def test_custom_taxid(self):
        p = self._kraken_out(["C\tread1/1\t644\t150\t0:116"])
        keep, _, n_keep = vibrio_extract_lib.build_keep_set(p, cholerae_taxid="644")
        self.assertEqual(keep, {"read1"})
        self.assertEqual(n_keep, 1)

    def test_short_lines_skipped(self):
        p = self._kraken_out(["C\tonly_two_fields"])
        keep, n_class, n_keep = vibrio_extract_lib.build_keep_set(p)
        self.assertEqual((keep, n_class, n_keep), (set(), 0, 0))

    def test_empty_file(self):
        p = self._kraken_out([])
        keep, n_class, n_keep = vibrio_extract_lib.build_keep_set(p)
        self.assertEqual((keep, n_class, n_keep), (set(), 0, 0))


class TestFilterFastq(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_only_kept_reads_written(self):
        inp = Path(self.tmp.name) / "in.fastq.gz"
        out = Path(self.tmp.name) / "out.fastq.gz"
        write_gz(inp,
            "@read1/1\nACGT\n+\nIIII\n"
            "@read2/1\nTTTT\n+\nIIII\n"
            "@read3/1\nGGGG\n+\nIIII\n")
        n_in, n_out = vibrio_extract_lib.filter_fastq(inp, out, {"read1", "read3"})
        self.assertEqual((n_in, n_out), (3, 2))
        content = read_gz(out)
        self.assertIn("@read1/1", content)
        self.assertIn("@read3/1", content)
        self.assertNotIn("@read2/1", content)

    def test_header_whitespace_and_pair_suffix_both_stripped(self):
        inp = Path(self.tmp.name) / "in.fastq.gz"
        out = Path(self.tmp.name) / "out.fastq.gz"
        write_gz(inp, "@read1/2 some comment\nACGT\n+\nIIII\n")
        n_in, n_out = vibrio_extract_lib.filter_fastq(inp, out, {"read1"})
        self.assertEqual((n_in, n_out), (1, 1))

    def test_nothing_kept_produces_valid_empty_gzip(self):
        inp = Path(self.tmp.name) / "in.fastq.gz"
        out = Path(self.tmp.name) / "out.fastq.gz"
        write_gz(inp, "@read1/1\nACGT\n+\nIIII\n")
        n_in, n_out = vibrio_extract_lib.filter_fastq(inp, out, set())
        self.assertEqual((n_in, n_out), (1, 0))
        self.assertEqual(read_gz(out), "")

    def test_empty_input(self):
        inp = Path(self.tmp.name) / "in.fastq.gz"
        out = Path(self.tmp.name) / "out.fastq.gz"
        write_gz(inp, "")
        n_in, n_out = vibrio_extract_lib.filter_fastq(inp, out, {"whatever"})
        self.assertEqual((n_in, n_out), (0, 0))

    def test_keep_set_from_build_keep_set_round_trips_with_filter(self):
        # integration of the two functions in this module: a taxid that
        # should be dropped really is dropped end to end
        kraken_path = Path(self.tmp.name) / "kraken.out"
        write(kraken_path,
              "C\treadA/1\t666\t150\t0:116\n"
              "C\treadB/1\t644\t150\t0:116\n")
        keep, _, _ = vibrio_extract_lib.build_keep_set(kraken_path)

        inp = Path(self.tmp.name) / "in.fastq.gz"
        out = Path(self.tmp.name) / "out.fastq.gz"
        write_gz(inp, "@readA/1\nACGT\n+\nIIII\n@readB/1\nTTTT\n+\nIIII\n")
        n_in, n_out = vibrio_extract_lib.filter_fastq(inp, out, keep)
        self.assertEqual((n_in, n_out), (2, 1))
        self.assertIn("@readA/1", read_gz(out))


if __name__ == "__main__":
    unittest.main()
