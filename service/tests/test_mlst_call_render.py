"""MLST tri-state rendering: the portal must SHOW the typing outcome.

The pipeline's assembly block carries mlst_call (assigned / unassigned /
not_interpretable) because the bare ST value is ambiguous when null.
These tests pin that each state renders as a visible, localized line in
the lineage table -- an unassigned ST is a result a reader can see, in
all three report languages, with the raw allele profile next to it.

Style mirrors test_phase3.py (fixtures + render_html + i18n catalog).
Run: python -m unittest discover -s service/tests -t .
"""

from __future__ import annotations

import copy
import html as html_mod
import json
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from service.app.reports import build_context, render_html  # noqa: E402
from service.app.reports import i18n as i18n_mod  # noqa: E402

FIXTURES = REPO_ROOT / "service" / "tests" / "fixtures"
PROFILE = ["51", "41", "4", "79", "72", "1", "83"]


def load_fixture(name: str) -> dict:
    with open(FIXTURES / name) as f:
        return json.load(f)


class MlstCallRenderTestCase(unittest.TestCase):
    def setUp(self):
        self.rep = load_fixture("report_qc_pass.json")

    def _with_assembly(self, **assembly):
        rep = copy.deepcopy(self.rep)
        rep["assembly"] = assembly
        return rep

    def _html(self, report, lang="en"):
        return html_mod.unescape(
            render_html(report, "JOB-1", "org-a", lang=lang))

    def test_unassigned_renders_localized_with_profile(self):
        rep = self._with_assembly(
            mlst_st=None, mlst_call="unassigned",
            mlst_scheme="vcholerae", mlst_alleles=PROFILE,
            vibecheck_lineage=None, amr_genes=[])
        for lang in ("en", "fr", "ht"):
            rendered = self._html(rep, lang)
            s = i18n_mod.STRINGS[lang]
            self.assertIn(s["mlst_unassigned"], rendered, lang)
            self.assertIn(s["label_mlst_profile"], rendered, lang)
            self.assertIn("51-41-4-79-72-1-83", rendered, lang)

    def test_unassigned_context_values(self):
        rep = self._with_assembly(
            mlst_st=None, mlst_call="unassigned",
            mlst_scheme="vcholerae", mlst_alleles=PROFILE,
            vibecheck_lineage=None, amr_genes=[])
        ctx = build_context(rep, "JOB-1", "org-a", "en")
        self.assertEqual(ctx["mlst_call"], "unassigned")
        self.assertEqual(ctx["mlst_st"],
                         i18n_mod.STRINGS["en"]["mlst_unassigned"])
        self.assertEqual(ctx["mlst_profile"], "51-41-4-79-72-1-83")

    def test_assigned_renders_st_value(self):
        rep = self._with_assembly(
            mlst_st=69, mlst_call="assigned", mlst_scheme="vcholerae",
            mlst_alleles=["7", "5", "3"], vibecheck_lineage="T12",
            amr_genes=[])
        rendered = self._html(rep)
        self.assertIn("ST69", rendered)
        ctx = build_context(rep, "JOB-1", "org-a", "en")
        self.assertEqual(ctx["mlst_call"], "assigned")

    def test_not_interpretable_renders_localized(self):
        rep = self._with_assembly(
            mlst_st=None, mlst_call="not_interpretable",
            mlst_scheme="vcholerae", mlst_alleles=["7", "5", "3"],
            vibecheck_lineage=None, amr_genes=[])
        for lang in ("en", "fr", "ht"):
            rendered = self._html(rep, lang)
            self.assertIn(
                i18n_mod.STRINGS[lang]["mlst_not_interpretable"],
                rendered, lang)

    def test_legacy_report_without_mlst_call_still_renders_st(self):
        # The shipped fixture predates mlst_call: the bare ST renders as
        # before, and no profile row appears.
        rendered = self._html(self.rep)
        self.assertIn("69", rendered)
        self.assertNotIn("51-41-4-79-72-1-83", rendered)
        ctx = build_context(self.rep, "JOB-1", "org-a", "en")
        self.assertEqual(ctx["mlst_profile"], "")


if __name__ == "__main__":
    unittest.main()
