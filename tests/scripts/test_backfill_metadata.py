import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))
import backfill_metadata  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


class TestLoadMapping(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _tsv(self, text):
        p = Path(self.tmp.name) / "cohort.tsv"
        write(p, text)
        return str(p)

    def test_basic(self):
        p = self._tsv(
            "run_accession\tsample_accession\tstudy_accession\tcountry\tcollection_date\n"
            "SRR1\tSAMN1\tPRJNA1\tHaiti\t2022-10-03\n"
            "SRR2\tSAMN2\tPRJNA1\tHaiti: Port-au-Prince\t\n"
        )
        m = backfill_metadata.load_mapping(p)
        self.assertEqual(m["SRR1"], ("2022-10-03", "Haiti"))
        self.assertEqual(m["SRR2"], (None, "Haiti: Port-au-Prince"))

    def test_comment_lines_skipped(self):
        p = self._tsv(
            "# cohort: test\n"
            "run_accession\tcountry\tcollection_date\n"
            "SRR1\tHaiti\t2022-01-01\n"
        )
        m = backfill_metadata.load_mapping(p)
        self.assertEqual(list(m), ["SRR1"])

    def test_empty_date_becomes_none_not_empty_string(self):
        # The watcher must never see "" as a date -- only None or ISO.
        p = self._tsv("run_accession\tcountry\tcollection_date\nSRR1\tHaiti\t\n")
        m = backfill_metadata.load_mapping(p)
        self.assertIsNone(m["SRR1"][0])


if __name__ == "__main__":
    unittest.main()
