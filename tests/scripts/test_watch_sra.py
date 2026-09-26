import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))
import watch_sra  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


class TestExtractRunAccessions(unittest.TestCase):
    """Regression coverage for the real bug this session found: NCBI
    nests run accessions inside a "Runs" Item as escaped XML text, not
    an Item literally named "Run". The old code looked for "Run"
    (singular), which never matched -- silently, since .iter() just
    finds nothing rather than raising -- so the watcher returned zero
    new accessions on every scheduled run."""

    def _docsum(self, runs_item_name, runs_text):
        xml = f'''<eSummaryResult>
  <DocSum>
    <Id>12345</Id>
    <Item Name="{runs_item_name}" Type="String">{runs_text}</Item>
  </DocSum>
</eSummaryResult>'''
        return ET.fromstring(xml)

    def test_extracts_from_correctly_named_runs_item(self):
        root = self._docsum("Runs", '&lt;Run acc="SRR12345678" total_spots="100"/&gt;')
        self.assertEqual(watch_sra.extract_run_accessions(root), ["SRR12345678"])

    def test_wrong_item_name_yields_nothing_not_an_exception(self):
        # this is exactly the pre-fix bug, kept as a regression guard:
        # a wrong Item name must not raise, but it must also not find
        # anything -- proving the current code depends on getting the
        # name right, rather than "any Item happens to work"
        root = self._docsum("Run", '&lt;Run acc="SRR12345678" total_spots="100"/&gt;')
        self.assertEqual(watch_sra.extract_run_accessions(root), [])

    def test_multiple_runs_in_one_experiment(self):
        root = self._docsum(
            "Runs",
            '&lt;Run acc="SRR001" total_spots="1"/&gt;'
            '&lt;Run acc="SRR002" total_spots="2"/&gt;',
        )
        self.assertEqual(watch_sra.extract_run_accessions(root), ["SRR001", "SRR002"])

    def test_multiple_docsums(self):
        xml = '''<eSummaryResult>
  <DocSum><Id>1</Id><Item Name="Runs" Type="String">&lt;Run acc="SRR001"/&gt;</Item></DocSum>
  <DocSum><Id>2</Id><Item Name="Runs" Type="String">&lt;Run acc="SRR002"/&gt;</Item></DocSum>
</eSummaryResult>'''
        root = ET.fromstring(xml)
        self.assertEqual(watch_sra.extract_run_accessions(root), ["SRR001", "SRR002"])

    def test_no_docsums_at_all(self):
        root = ET.fromstring("<eSummaryResult></eSummaryResult>")
        self.assertEqual(watch_sra.extract_run_accessions(root), [])

    def test_runs_item_present_but_empty(self):
        root = self._docsum("Runs", "")
        self.assertEqual(watch_sra.extract_run_accessions(root), [])

    def test_accepts_err_and_drr_prefixes(self):
        root = self._docsum("Runs", '&lt;Run acc="ERR11684929"/&gt;&lt;Run acc="DRR000001"/&gt;')
        self.assertEqual(
            watch_sra.extract_run_accessions(root), ["ERR11684929", "DRR000001"])


class TestBuildNewItems(unittest.TestCase):
    def test_haiti_tier_takes_priority_over_global_for_same_accession(self):
        # an accession appearing in both geo and recent results must only
        # be added once, tagged haiti (geo is processed first)
        items = watch_sra.build_new_items(
            geo_accs=["SRR1"], recent_accs=["SRR1", "SRR2"], seen=set(), max_new=10)
        self.assertEqual(items, [
            {"accession": "SRR1", "priority": "haiti"},
            {"accession": "SRR2", "priority": "global"},
        ])

    def test_already_seen_accessions_are_excluded_from_both_tiers(self):
        items = watch_sra.build_new_items(
            geo_accs=["SRR1", "SRR2"], recent_accs=["SRR2", "SRR3"],
            seen={"SRR2"}, max_new=10)
        self.assertEqual(items, [
            {"accession": "SRR1", "priority": "haiti"},
            {"accession": "SRR3", "priority": "global"},
        ])

    def test_max_new_caps_total_across_both_tiers(self):
        items = watch_sra.build_new_items(
            geo_accs=["SRR1", "SRR2", "SRR3"], recent_accs=["SRR4", "SRR5"],
            seen=set(), max_new=2)
        self.assertEqual(len(items), 2)
        self.assertEqual([i["accession"] for i in items], ["SRR1", "SRR2"])

    def test_empty_inputs_yield_empty_list(self):
        self.assertEqual(watch_sra.build_new_items([], [], set(), 10), [])

    def test_duplicate_within_same_tier_only_added_once(self):
        # esummary can return the same accession twice (e.g. a chunking
        # edge case); the priority-tagged list must not carry duplicates
        items = watch_sra.build_new_items(
            geo_accs=["SRR1", "SRR1"], recent_accs=[], seen=set(), max_new=10)
        self.assertEqual(items, [{"accession": "SRR1", "priority": "haiti"}])


class TestLoadSeen(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_file_returns_empty_set_not_exception(self):
        self.assertEqual(watch_sra.load_seen(Path(self.tmp.name) / "nope.txt"), set())

    def test_comments_and_blank_lines_ignored(self):
        p = Path(self.tmp.name) / "seen.txt"
        write(p, "# header comment\nSRR1\n\nSRR2\n# trailing\n")
        self.assertEqual(watch_sra.load_seen(p), {"SRR1", "SRR2"})

    def test_whitespace_stripped(self):
        p = Path(self.tmp.name) / "seen.txt"
        write(p, "  SRR1  \nSRR2\t\n")
        self.assertEqual(watch_sra.load_seen(p), {"SRR1", "SRR2"})


if __name__ == "__main__":
    unittest.main()
