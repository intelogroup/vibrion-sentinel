import json
import sys
import tempfile
import unittest
import urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest import mock

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


class TestHttpRetry(unittest.TestCase):
    """The 2026-09-27 scheduled watcher died on a single NCBI HTTP 429.
    These tests pin the retry contract: 429/5xx are retried with
    exponential backoff, Retry-After is honored, anything else raises
    immediately, and the loop gives up after MAX_ATTEMPTS."""

    def _ok(self, body=b'{"esearchresult": {}}'):
        resp = mock.MagicMock()
        resp.__enter__.return_value = resp
        resp.read.return_value = body
        resp.headers = {}
        return resp

    def _err(self, code, retry_after=None):
        hdrs = {}
        if retry_after is not None:
            hdrs["Retry-After"] = str(retry_after)
        return urllib.error.HTTPError("https://x", code, "err", hdrs, None)

    def _fake_urlopen(self, behaviors):
        it = iter(behaviors)

        def fake(req, timeout=None):
            b = next(it)
            if isinstance(b, Exception):
                raise b
            return b

        return fake

    def _patched(self, behaviors):
        return (
            mock.patch("urllib.request.urlopen",
                       side_effect=self._fake_urlopen(behaviors)),
            mock.patch.object(watch_sra.time, "sleep"),
        )

    def test_429_retried_then_succeeds_with_exponential_backoff(self):
        urlopen_p, sleep_p = self._patched(
            [self._err(429), self._err(429), self._ok()])
        with urlopen_p, sleep_p as slp:
            body, _ = watch_sra._http_get("https://x")
        self.assertEqual(json.loads(body), {"esearchresult": {}})
        self.assertEqual(slp.call_count, 2)
        d1 = slp.call_args_list[0][0][0]
        d2 = slp.call_args_list[1][0][0]
        self.assertTrue(1.6 <= d1 <= 2.4, d1)  # 2s +/-20% jitter
        self.assertTrue(3.2 <= d2 <= 4.8, d2)  # 4s +/-20% jitter

    def test_retry_after_header_honored(self):
        urlopen_p, sleep_p = self._patched(
            [self._err(429, retry_after=5), self._ok()])
        with urlopen_p, sleep_p as slp:
            watch_sra._http_get("https://x")
        d = slp.call_args[0][0]
        self.assertTrue(4.0 <= d <= 6.0, d)  # 5s +/-20% jitter

    def test_gives_up_after_max_attempts(self):
        urlopen_p, sleep_p = self._patched(
            [self._err(429)] * (watch_sra.MAX_ATTEMPTS + 1))
        with urlopen_p, sleep_p as slp:
            with self.assertRaises(urllib.error.HTTPError):
                watch_sra._http_get("https://x")
        self.assertEqual(slp.call_count, watch_sra.MAX_ATTEMPTS - 1)

    def test_400_not_retried(self):
        urlopen_p, sleep_p = self._patched([self._err(400), self._ok()])
        with urlopen_p as uo, sleep_p as slp:
            with self.assertRaises(urllib.error.HTTPError):
                watch_sra._http_get("https://x")
        self.assertEqual(uo.call_count, 1)
        slp.assert_not_called()

    def test_503_retried(self):
        urlopen_p, sleep_p = self._patched([self._err(503), self._ok()])
        with urlopen_p, sleep_p as slp:
            body, _ = watch_sra._http_get("https://x")
        self.assertEqual(json.loads(body), {"esearchresult": {}})
        self.assertEqual(slp.call_count, 1)


class TestApiKeyParam(unittest.TestCase):
    def test_key_added_when_env_set(self):
        with mock.patch.dict(watch_sra.os.environ, {"NCBI_API_KEY": "abc123"}):
            p = watch_sra._params({"db": "sra"}, True)
        self.assertEqual(p["api_key"], "abc123")
        self.assertEqual(p["retmode"], "json")

    def test_no_key_by_default(self):
        with mock.patch.dict(watch_sra.os.environ, {}, clear=True):
            p = watch_sra._params({"db": "sra"}, False)
        self.assertNotIn("api_key", p)
        self.assertNotIn("retmode", p)

    def test_blank_key_ignored(self):
        with mock.patch.dict(watch_sra.os.environ, {"NCBI_API_KEY": "   "}):
            self.assertNotIn("api_key", watch_sra._params({}, True))


if __name__ == "__main__":
    unittest.main()
