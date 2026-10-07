"""Phase 3 tests: report rendering (3 languages), worker hardening,
human auth + dashboard, batch backends, banned-phrase CI check.

Conventions mirror the earlier suites: unittest, fakes, no network, no
credentials. Run: /tmp/p3-venv/bin/python -m unittest discover -s service/tests -t .
"""

from __future__ import annotations

import hashlib
import html
import json
import os
import re
import sys
import tempfile
import threading
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from fastapi.testclient import TestClient  # noqa: E402

from service import worker  # noqa: E402
from service.app import batch as batch_mod  # noqa: E402
from service.app.auth import hash_key  # noqa: E402
from service.app.batch import LocalBackend, run_to_outcome  # noqa: E402
from service.app.config import Settings  # noqa: E402
from service.app.human_auth import (  # noqa: E402
    SESSION_COOKIE,
    hash_password,
    verify_password,
)
from service.app.main import create_app  # noqa: E402
from service.app.reports import build_context, render_html, render_pdf  # noqa: E402
from service.app.reports import i18n as i18n_mod  # noqa: E402
from service.app.storage import MemoryStorage  # noqa: E402
from service.tests.fakes import FakeDB  # noqa: E402

FIXTURES = REPO_ROOT / "service" / "tests" / "fixtures"


def load_fixture(name: str) -> dict:
    with open(FIXTURES / name) as f:
        return json.load(f)


def make_fastq(n_reads: int = 50) -> bytes:
    recs = []
    for i in range(n_reads):
        seq = "ACGT" * 25
        recs.append(f"@read{i}\n{seq}\n+\n{'I' * 100}\n")
    return "".join(recs).encode()


# ======================================================================
# Report rendering
# ======================================================================
class ReportRenderTestCase(unittest.TestCase):
    def setUp(self):
        self.rep_pass = load_fixture("report_qc_pass.json")
        self.rep_fail = load_fixture("report_qc_fail.json")

    def _html(self, report, lang="en"):
        return render_html(report, "JOB-1", "org-a", lang=lang,
                           uploaded_at="2026-10-07", completed_at="2026-10-07")

    # -- RUO + disclaimer in every language --------------------------------
    def test_ruo_label_and_disclaimer_present_all_languages(self):
        for lang in ("en", "fr", "ht"):
            html = self._html(self.rep_pass, lang)
            # RUO label appears in header AND footer (at least twice).
            ruo = i18n_mod.RUO_LABEL[lang]
            self.assertGreaterEqual(html.count(ruo), 2,
                                    f"RUO label missing/once in {lang}")
            self.assertIn(i18n_mod.DISCLAIMER[lang], html,
                          f"disclaimer missing in {lang}")

    def test_translator_review_markers_on_non_english(self):
        for lang in ("fr", "ht"):
            html = self._html(self.rep_pass, lang)
            self.assertIn("TRANSLATOR-REVIEW", html,
                          f"no review marker in {lang}")
            self.assertIn(i18n_mod.REVIEW_BANNERS[lang], html,
                          f"no visible review banner in {lang}")
        html_en = self._html(self.rep_pass, "en")
        self.assertNotIn("TRANSLATOR-REVIEW", html_en)

    def test_unknown_lang_falls_back_to_english(self):
        html = self._html(self.rep_pass, "de")
        self.assertIn("Genomic surveillance report", html)
        self.assertNotIn("TRANSLATOR-REVIEW", html)

    # -- QC verdict first and loud ------------------------------------------
    def test_qc_pass_renders_results_plainly(self):
        html = self._html(self.rep_pass, "en")
        self.assertIn("PASS", html)
        self.assertNotIn('<div class="qc-fail-zone">', html)
        self.assertNotIn('<div class="not-interpretable">', html)
        # Key results present.
        self.assertIn("T12", html)
        self.assertIn("69", html)  # MLST ST
        self.assertIn("floR", html)
        # Gene-only disclaimer present.
        self.assertIn(i18n_mod.STRINGS["en"]["amr_disclaimer"], html)

    def test_qc_fail_greys_out_results(self):
        html = self._html(self.rep_fail, "en")
        self.assertIn("FAIL", html)
        self.assertIn('<div class="qc-fail-zone">', html)
        self.assertIn('<div class="not-interpretable">', html)
        self.assertIn(i18n_mod.STRINGS["en"]["not_interpretable"], html)

    def test_qc_fail_all_languages(self):
        for lang in ("fr", "ht"):
            html = self._html(self.rep_fail, lang)
            self.assertIn('<div class="qc-fail-zone">', html, lang)
            self.assertIn(i18n_mod.STRINGS[lang]["not_interpretable"], html, lang)

    def test_header_carries_sample_org_versions(self):
        html = self._html(self.rep_pass, "en")
        self.assertIn("JOB-1", html)
        self.assertIn("org-a", html)
        self.assertIn("0.2.0", html)
        self.assertIn("nanopore", html)
        self.assertIn("hac", html)
        self.assertIn("amrfinderplus", html)

    def test_missing_keys_render_dash_not_crash(self):
        html = render_html({}, "JOB-X", "org-a", lang="en")
        self.assertIn("—", html)
        html2 = render_html({"qc": {"status": "pass"}}, "JOB-X", "org-a", lang="ht")
        self.assertIn("TRANSLATOR-REVIEW", html2)

    # -- PDF -----------------------------------------------------------------
    def test_pdf_generates(self):
        pdf = render_pdf(self.rep_pass, "JOB-1", "org-a", lang="en")
        self.assertTrue(pdf.startswith(b"%PDF"), pdf[:20])
        self.assertGreater(len(pdf), 1000)

    def test_pdf_french(self):
        pdf = render_pdf(self.rep_fail, "JOB-2", "org-a", lang="fr")
        self.assertTrue(pdf.startswith(b"%PDF"))


# ======================================================================
# Banned phrases (CI check)
# ======================================================================
class BannedPhrasesTestCase(unittest.TestCase):
    MANDATED = (
        list(i18n_mod.RUO_LABEL.values())
        + list(i18n_mod.DISCLAIMER.values())
        # The AMR table's gene-only restatement carries the mandated concept
        # ("not a susceptibility result") in all three languages.
        + [i18n_mod.STRINGS[lang]["amr_disclaimer"] for lang in ("en", "fr", "ht")]
    )

    def _read(self, path: str) -> str:
        with open(path, encoding="utf-8") as f:
            return f.read()

    def _corpus(self) -> str:
        """Everything customer-facing the templates/i18n can emit.

        Templates (which carry {{ }} placeholders) plus the i18n string
        catalog plus the *rendered* report in all three languages — the
        rendered output is what carries the mandated disclaimer verbatim.
        """
        parts: list[str] = []
        for root, _dirs, files in os.walk(REPO_ROOT / "service" / "app"):
            for fn in files:
                if fn.endswith(".html"):
                    parts.append(self._read(os.path.join(root, fn)))
        for lang, strings in i18n_mod.STRINGS.items():
            parts.extend(strings.values())
        parts.extend(i18n_mod.REVIEW_BANNERS.values())
        report = load_fixture("report_qc_pass.json")
        for lang in ("en", "fr", "ht"):
            parts.append(render_html(report, "JOB-1", "org-a", lang=lang))
        return "\n".join(parts)

    def test_no_banned_phrases(self):
        phrases_path = REPO_ROOT / "service" / "tests" / "banned_phrases.txt"
        phrases = [
            ln.strip()
            for ln in self._read(str(phrases_path)).splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        self.assertTrue(phrases, "banned phrase list is empty")
        corpus = self._corpus()
        # Jinja2 autoescape turns ' into &#39; etc. in rendered output —
        # unescape before stripping so mandated strings match verbatim.
        corpus = html.unescape(corpus)
        # Strip the mandated strings first: the verbatim RUO label and the
        # extended surveillance disclaimers legitimately contain forbidden
        # words (e.g. "diagnose"/"treat"/"patient" inside a prohibition).
        for mandated in self.MANDATED:
            corpus = corpus.replace(mandated, "")
        hits = []
        for phrase in phrases:
            pat = re.compile(r"\b" + re.escape(phrase) + r"\b", re.IGNORECASE)
            m = pat.search(corpus)
            if m:
                hits.append(
                    f"{phrase!r} near: ...{corpus[max(0, m.start()-60):m.end()+60]}..."
                )
        self.assertEqual(hits, [], "banned phrases found in templates/i18n")

    def test_mandated_strings_are_actually_present(self):
        # Guard against the exemption list drifting out of sync: every
        # mandated string must actually occur in the rendered corpus, so the
        # strip in test_no_banned_phrases is never vacuous.
        corpus = self._corpus()
        for s in self.MANDATED:
            self.assertIn(s, corpus, f"mandated string not rendered: {s[:60]}...")


# ======================================================================
# Batch backends
# ======================================================================
class BatchBackendTestCase(unittest.TestCase):
    def test_local_ok(self):
        with tempfile.TemporaryDirectory() as d:
            be = LocalBackend()
            rc, log = be.run(d, ["echo", "hello"], timeout_hours=1)
            self.assertEqual(rc, 0)
            self.assertIn("hello", log)

    def test_local_timeout_kills(self):
        with tempfile.TemporaryDirectory() as d:
            be = LocalBackend()
            try:
                be.run(d, ["sleep", "30"], timeout_hours=0.0003)  # ~1s
                self.fail("expected JobTimeoutError")
            except batch_mod.JobTimeoutError as e:
                self.assertIn("wall clock", str(e))

    def test_local_missing_binary_is_infra(self):
        with tempfile.TemporaryDirectory() as d:
            be = LocalBackend()
            try:
                be.run(d, ["/nonexistent/binary-xyz"], timeout_hours=1)
                self.fail("expected InfraError")
            except batch_mod.InfraError:
                pass

    def test_run_to_outcome_classification(self):
        with tempfile.TemporaryDirectory() as d:
            be = LocalBackend()
            self.assertEqual(run_to_outcome(be, d, ["true"], 1)[0], "ok")
            self.assertEqual(run_to_outcome(be, d, ["false"], 1)[0], "pipeline_error")
            self.assertEqual(
                run_to_outcome(be, d, ["sleep", "30"], 0.0003)[0], "timeout"
            )
            self.assertEqual(
                run_to_outcome(be, d, ["/nonexistent/binary-xyz"], 1)[0],
                "infra_error",
            )


# ======================================================================
# Worker hardening: concurrency, timeout, retry discipline
# ======================================================================
class WorkerHardeningTestCase(unittest.TestCase):
    def setUp(self):
        self.db = FakeDB()
        self.storage = MemoryStorage()
        self.tmp = tempfile.mkdtemp(prefix="sentinel-p3-worker-")
        self.settings = Settings()
        self.settings.worker_workdir = os.path.join(self.tmp, "work")
        self.settings.snakemake_cores = 2
        self.settings.repo_dir = str(REPO_ROOT)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _seed(self, n=1, **over):
        payload = make_fastq()
        ids = []
        for _ in range(n):
            job = self.db.create_job(
                {
                    "org_id": "org-a",
                    "filename": "run.fastq.gz",
                    "size_bytes": len(payload),
                    "platform": "nanopore",
                    "basecaller_model": "hac",
                    "tier": "lite",
                    "status": "queued",
                    "sha256": hashlib.sha256(payload).hexdigest(),
                    **over,
                }
            )
            job = self.db.update_job(
                job["id"], {"r2_key": f"intake/org-a/{job['id']}/run.fastq.gz"}
            )
            self.storage.create(job["r2_key"])
            self.storage.append(job["r2_key"], payload)
            self.storage.finalize(job["r2_key"])
            ids.append(job["id"])
        return ids

    def _ok_backend(self, seen=None, report_extra=None):
        lock = threading.Lock()

        class Stub(batch_mod.BatchBackend):
            def run(self, workdir, cmd, timeout_hours):
                outdir = os.path.join(workdir, "out")
                os.makedirs(outdir, exist_ok=True)
                cfg_path = [c for c in cmd if c.endswith("config.yaml")][0]
                with open(cfg_path) as f:
                    cfg = f.read()
                # find which job by workdir name
                jid = os.path.basename(workdir)
                if seen is not None:
                    with lock:
                        seen.append(jid)
                report = {
                    "version": "0.2.0",
                    "accession": jid,
                    "qc": {"status": "pass", "reasons": []},
                }
                report.update(report_extra or {})
                with open(os.path.join(outdir, "report.json"), "w") as f:
                    json.dump(report, f)
                return 0, "stub ok"

        return Stub()

    def test_two_workers_never_double_claim(self):
        ids = self._seed(3)
        seen: list[str] = []

        def loop():
            while worker.process_one_job(
                self.db, self.storage, self.settings,
                backend=self._ok_backend(seen),
            ):
                pass

        threads = [threading.Thread(target=loop) for _ in range(2)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=60)
            self.assertFalse(t.is_alive(), "worker thread hung")

        self.assertEqual(sorted(seen), sorted(ids), "each job claimed exactly once")
        for jid in ids:
            self.assertEqual(self.db.jobs[jid]["status"], "done")

    def test_timeout_marks_failed(self):
        (jid,) = self._seed(1)

        class SlowBatch(batch_mod.BatchBackend):
            def run(self, workdir, cmd, timeout_hours):
                raise batch_mod.JobTimeoutError("command exceeded 0.001h wall clock")

        self.assertTrue(
            worker.process_one_job(self.db, self.storage, self.settings,
                                   backend=SlowBatch())
        )
        job = self.db.jobs[jid]
        self.assertEqual(job["status"], "failed")
        self.assertIn("timeout", job["status_reason"].lower())

    def test_infra_error_retries_exactly_3_then_fails(self):
        (jid,) = self._seed(1)

        class BoomStorage(MemoryStorage):
            def download_to_file(self, key, path):
                raise RuntimeError("R2 503")

        storage = BoomStorage()
        for i in range(1, 5):
            self.assertTrue(
                worker.process_one_job(self.db, storage, self.settings,
                                       backend=self._ok_backend()),
                f"call {i} should have processed the re-queued job",
            )
            job = self.db.jobs[jid]
            if i < 4:
                self.assertEqual(job["status"], "queued", f"after retry {i}")
                self.assertEqual(job.get("retry_count"), i)
            else:
                self.assertEqual(job["status"], "failed")
                self.assertIn("retries exhausted", job["status_reason"])
                self.assertEqual(job.get("retry_count"), 3)
        # No more work left.
        self.assertFalse(
            worker.process_one_job(self.db, storage, self.settings,
                                   backend=self._ok_backend())
        )

    def test_pipeline_error_never_retries(self):
        (jid,) = self._seed(1)

        class BadBatch(batch_mod.BatchBackend):
            def run(self, workdir, cmd, timeout_hours):
                return 1, "snakemake: FATAL something broke"

        self.assertTrue(
            worker.process_one_job(self.db, self.storage, self.settings,
                                   backend=BadBatch())
        )
        job = self.db.jobs[jid]
        self.assertEqual(job["status"], "failed")
        self.assertIn("FATAL", job["status_reason"])
        self.assertEqual(job.get("retry_count") or 0, 0)
        self.assertFalse(
            worker.process_one_job(self.db, self.storage, self.settings,
                                   backend=BadBatch())
        )

    def test_qc_fail_is_a_result_not_an_error(self):
        (jid,) = self._seed(1)
        backend = self._ok_backend(
            report_extra={"qc": {"status": "fail", "reasons": ["mean_depth low"]}}
        )
        self.assertTrue(
            worker.process_one_job(self.db, self.storage, self.settings,
                                   backend=backend)
        )
        job = self.db.jobs[jid]
        self.assertEqual(job["status"], "done")
        self.assertEqual(job["report"]["qc"]["status"], "fail")

    def test_concurrency_setting_flows(self):
        self.assertEqual(Settings().worker_concurrency, 2)
        os.environ["WORKER_CONCURRENCY"] = "4"
        try:
            self.assertEqual(Settings().worker_concurrency, 4)
        finally:
            del os.environ["WORKER_CONCURRENCY"]


# ======================================================================
# Report endpoints (HTML / PDF / tree), API key + session auth
# ======================================================================
class ReportEndpointTestCase(unittest.TestCase):
    def setUp(self):
        self.db = FakeDB()
        self.storage = MemoryStorage()
        self.settings = Settings()
        self.key_a = "sk_test_rep_a"
        self.key_b = "sk_test_rep_b"
        self.db.insert_api_key("org-a", hash_key(self.key_a), "a", "admin")
        self.db.insert_api_key("org-b", hash_key(self.key_b), "b", "admin")
        # A done job with a report in org-a.
        job = self.db.create_job(
            {
                "org_id": "org-a",
                "filename": "run.fastq.gz",
                "platform": "nanopore",
                "basecaller_model": "hac",
                "tier": "assembly",
                "status": "done",
                "report": load_fixture("report_qc_pass.json"),
                "created_at": "2026-10-07T10:00:00+00:00",
                "updated_at": "2026-10-07T11:00:00+00:00",
            }
        )
        self.job_id = job["id"]
        # A queued job (no report yet).
        q = self.db.create_job(
            {"org_id": "org-a", "filename": "x.fastq.gz",
             "platform": "illumina", "status": "queued"}
        )
        self.queued_id = q["id"]
        self.app = create_app(self.db, self.storage, self.settings)
        self.client = TestClient(self.app)

    def _h(self, key):
        return {"Authorization": f"Bearer {key}"}

    def test_html_report_with_api_key(self):
        r = self.client.get(f"/jobs/{self.job_id}/report",
                            headers=self._h(self.key_a))
        self.assertEqual(r.status_code, 200, r.text[:200])
        self.assertIn("text/html", r.headers["content-type"])
        self.assertIn("For Research Use Only", r.text)
        self.assertIn("T12", r.text)

    def test_html_report_french(self):
        r = self.client.get(f"/jobs/{self.job_id}/report?lang=fr",
                            headers=self._h(self.key_a))
        self.assertEqual(r.status_code, 200)
        self.assertIn("TRANSLATOR-REVIEW", r.text)
        self.assertIn("Rapport de surveillance génomique", r.text)

    def test_pdf_report(self):
        r = self.client.get(f"/jobs/{self.job_id}/report.pdf",
                            headers=self._h(self.key_a))
        self.assertEqual(r.status_code, 200, r.text[:200])
        self.assertIn("application/pdf", r.headers["content-type"])
        self.assertTrue(r.content.startswith(b"%PDF"))

    def test_report_not_ready_404(self):
        r = self.client.get(f"/jobs/{self.queued_id}/report",
                            headers=self._h(self.key_a))
        self.assertEqual(r.status_code, 404)

    def test_cross_org_report_404(self):
        r = self.client.get(f"/jobs/{self.job_id}/report",
                            headers=self._h(self.key_b))
        self.assertEqual(r.status_code, 404)
        r = self.client.get(f"/jobs/{self.job_id}/report.pdf",
                            headers=self._h(self.key_b))
        self.assertEqual(r.status_code, 404)

    def test_no_auth_401(self):
        r = self.client.get(f"/jobs/{self.job_id}/report")
        self.assertEqual(r.status_code, 401)

    def test_tree_pending_without_auspice_url(self):
        r = self.client.get(f"/jobs/{self.job_id}/tree",
                            headers=self._h(self.key_a))
        self.assertEqual(r.status_code, 200)
        body = r.json()
        self.assertIsNone(body["tree"])
        self.assertIn("Auspice", body["followup"])

    def test_tree_redirects_when_configured(self):
        self.settings.auspice_base_url = "https://auspice.example.org"
        app = create_app(self.db, self.storage, self.settings)
        client = TestClient(app)
        r = client.get(f"/jobs/{self.job_id}/tree",
                       headers=self._h(self.key_a), follow_redirects=False)
        self.assertEqual(r.status_code, 302)
        self.assertIn("auspice.example.org", r.headers["location"])
        self.assertIn(self.job_id, r.headers["location"])


# ======================================================================
# Human auth: signup / login / dashboard / invites / revocation
# ======================================================================
class HumanAuthTestCase(unittest.TestCase):
    def setUp(self):
        self.db = FakeDB()
        self.storage = MemoryStorage()
        self.settings = Settings()
        self.app = create_app(self.db, self.storage, self.settings)
        self.client = TestClient(self.app)

    def _signup(self, client, email="ada@example.org", password="correct-horse-99"):
        r = client.post("/auth/signup",
                        json={"email": email, "password": password})
        self.assertEqual(r.status_code, 201, r.text)
        return r.json()

    def test_signup_creates_user_org_key_and_session(self):
        body = self._signup(self.client)
        self.assertIn("org_id", body)
        self.assertTrue(body["api_key"].startswith("sk_"))
        self.assertIn("session_token", body)
        # Cookie set.
        self.assertIn(SESSION_COOKIE, self.client.cookies)
        # /auth/me works.
        r = self.client.get("/auth/me")
        self.assertEqual(r.status_code, 200)
        me = r.json()
        self.assertEqual(me["email"], "ada@example.org")
        self.assertEqual(me["orgs"][0]["role"], "admin")
        # The issued API key works too.
        r2 = TestClient(self.app)
        r = r2.get("/jobs", headers={"Authorization": f"Bearer {body['api_key']}"})
        self.assertEqual(r.status_code, 200)

    def test_signup_rejects_bad_input(self):
        r = self.client.post("/auth/signup",
                             json={"email": "not-an-email", "password": "correct-horse-99"})
        self.assertEqual(r.status_code, 400)
        r = self.client.post("/auth/signup",
                             json={"email": "a@b.co", "password": "short"})
        self.assertEqual(r.status_code, 400)

    def test_signup_duplicate_email_409(self):
        self._signup(self.client)
        r = self.client.post("/auth/signup",
                             json={"email": "ada@example.org",
                                   "password": "another-good-pw"})
        self.assertEqual(r.status_code, 409)

    def test_login_wrong_password_401(self):
        self._signup(self.client)
        c2 = TestClient(self.app)
        r = c2.post("/auth/login",
                    json={"email": "ada@example.org", "password": "wrong-password-1"})
        self.assertEqual(r.status_code, 401)
        r = c2.post("/auth/login",
                    json={"email": "nope@example.org", "password": "wrong-password-1"})
        self.assertEqual(r.status_code, 401)

    def test_passwords_are_argon2_not_plaintext(self):
        self._signup(self.client)
        user = self.db.get_user_by_email("ada@example.org")
        self.assertTrue(user["password_hash"].startswith("$argon2"))
        self.assertNotIn("correct-horse-99", user["password_hash"])
        self.assertTrue(verify_password(user["password_hash"], "correct-horse-99"))
        self.assertFalse(verify_password(user["password_hash"], "nope"))

    def test_logout_revokes_session(self):
        self._signup(self.client)
        r = self.client.post("/auth/logout")
        self.assertEqual(r.status_code, 204)
        r = self.client.get("/auth/me")
        self.assertEqual(r.status_code, 401)

    def test_change_password_revokes_all_sessions(self):
        self._signup(self.client)
        # Second session via login on another client.
        c2 = TestClient(self.app)
        r = c2.post("/auth/login",
                    json={"email": "ada@example.org",
                          "password": "correct-horse-99"})
        self.assertEqual(r.status_code, 200)
        # Change password from the first client.
        r = self.client.post("/auth/change-password", json={
            "current_password": "correct-horse-99",
            "new_password": "even-better-horse-42",
        })
        self.assertEqual(r.status_code, 200, r.text)
        # Both sessions dead.
        self.assertEqual(self.client.get("/auth/me").status_code, 401)
        self.assertEqual(c2.get("/auth/me").status_code, 401)
        # New password works.
        c3 = TestClient(self.app)
        r = c3.post("/auth/login",
                    json={"email": "ada@example.org",
                          "password": "even-better-horse-42"})
        self.assertEqual(r.status_code, 200)

    def test_dashboard_lists_org_jobs_with_links(self):
        body = self._signup(self.client)
        org_id = body["org_id"]
        job = self.db.create_job(
            {"org_id": org_id, "filename": "run.fastq.gz",
             "platform": "nanopore", "status": "done",
             "report": load_fixture("report_qc_pass.json")}
        )
        r = self.client.get("/dashboard")
        self.assertEqual(r.status_code, 200, r.text[:200])
        self.assertIn("text/html", r.headers["content-type"])
        self.assertIn("For Research Use Only", r.text)
        self.assertIn(job["id"][:8], r.text)
        self.assertIn(f"/jobs/{job['id']}/report", r.text)
        self.assertIn(f"/jobs/{job['id']}/report.pdf", r.text)
        self.assertIn(f"/jobs/{job['id']}/tree", r.text)

    def test_dashboard_requires_login(self):
        c = TestClient(self.app)
        r = c.get("/dashboard", follow_redirects=False)
        self.assertEqual(r.status_code, 302)
        self.assertIn("/login", r.headers["location"])

    def test_dashboard_report_link_works_with_session(self):
        body = self._signup(self.client)
        job = self.db.create_job(
            {"org_id": body["org_id"], "filename": "run.fastq.gz",
             "platform": "nanopore", "status": "done",
             "report": load_fixture("report_qc_fail.json")}
        )
        # No API key — the session cookie alone opens the report.
        r = self.client.get(f"/jobs/{job['id']}/report")
        self.assertEqual(r.status_code, 200, r.text[:200])
        self.assertIn("FAIL", r.text)

    def test_session_cannot_see_other_org_report(self):
        a = self._signup(self.client)  # client has A's session
        # B's org + job, via a second client.
        cb = TestClient(self.app)
        b = self._signup(cb, email="bob@example.org")
        job = self.db.create_job(
            {"org_id": b["org_id"], "filename": "run.fastq.gz",
             "platform": "nanopore", "status": "done",
             "report": load_fixture("report_qc_pass.json")}
        )
        r = self.client.get(f"/jobs/{job['id']}/report")
        self.assertEqual(r.status_code, 404)

    def test_invite_flow(self):
        a = self._signup(self.client)
        # A (admin) invites carol@example.org.
        r = self.client.post(
            "/org/invites",
            headers={"X-Org-Id": a["org_id"]},
            json={"email": "carol@example.org", "role": "member"},
        )
        self.assertEqual(r.status_code, 201, r.text)
        token = r.json()["invite_token"]
        # Carol signs up and accepts.
        cc = TestClient(self.app)
        self._signup(cc, email="carol@example.org")
        r = cc.post("/auth/accept-invite", json={"token": token})
        self.assertEqual(r.status_code, 200, r.text)
        self.assertEqual(r.json()["org_id"], a["org_id"])
        r = cc.get("/auth/me")
        orgs = {o["org_id"]: o["role"] for o in r.json()["orgs"]}
        self.assertEqual(orgs[a["org_id"]], "member")
        # Token is single-use.
        r = cc.post("/auth/accept-invite", json={"token": token})
        self.assertEqual(r.status_code, 404)

    def test_invite_requires_admin(self):
        a = self._signup(self.client)
        # Carol joins as viewer via an admin-made invite, then tries inviting.
        r = self.client.post(
            "/org/invites",
            headers={"X-Org-Id": a["org_id"]},
            json={"email": "dave@example.org", "role": "viewer"},
        )
        token = r.json()["invite_token"]
        cd = TestClient(self.app)
        self._signup(cd, email="dave@example.org")
        cd.post("/auth/accept-invite", json={"token": token})
        r = cd.post(
            "/org/invites",
            headers={"X-Org-Id": a["org_id"]},
            json={"email": "eve@example.org", "role": "viewer"},
        )
        self.assertEqual(r.status_code, 403)

    def test_api_key_invite_path(self):
        body = self._signup(self.client)
        key = body["api_key"]
        c = TestClient(self.app)
        r = c.post(
            "/org/invites",
            headers={"Authorization": f"Bearer {key}"},
            json={"email": "frank@example.org", "role": "viewer"},
        )
        self.assertEqual(r.status_code, 201, r.text)
