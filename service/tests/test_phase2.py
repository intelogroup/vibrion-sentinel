"""Phase 2 multi-tenancy tests: roles, key management, envelope encryption,
rotation, deletion workflows, migration.

Conventions mirror service/tests/test_tus.py: unittest, no network, no
credentials. Run: python -m unittest discover -s service/tests -t .
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from fastapi.testclient import TestClient  # noqa: E402

from service.app import tus  # noqa: E402
from service.app.auth import hash_key  # noqa: E402
from service.app.config import Settings  # noqa: E402
from service.app.crypto import (  # noqa: E402
    DecryptionError,
    EncryptedStorage,
    KeyProvider,
    generate_dek,
    load_kek,
    unwrap_dek,
    wrap_dek,
)
from service.app.main import create_app  # noqa: E402
from service.app.storage import MemoryStorage  # noqa: E402
from service.migrate_encrypt import discover_orgs, migrate, plan_org  # noqa: E402
from service.tests.fakes import FakeDB  # noqa: E402


def make_fastq(n_reads: int = 20) -> bytes:
    recs = []
    for i in range(n_reads):
        seq = "ACGT" * 25
        recs.append(f"@read{i}\n{seq}\n+\n{'I' * 100}\n")
    return "".join(recs).encode()


TEST_KEK = bytes(range(32))  # deterministic test KEK (never used in prod)


class Phase2Base(unittest.TestCase):
    def make_app(self, encrypted=True):
        self.db = FakeDB()
        self.inner = MemoryStorage()
        if encrypted:
            self.kp = KeyProvider(self.db, TEST_KEK)
            self.storage = EncryptedStorage(self.inner, self.kp)
        else:
            self.kp = None
            self.storage = self.inner
        self.settings = Settings()
        # keys: admin / member / viewer on org-a, admin on org-b
        self.admin_key = "sk_admin_a"
        self.member_key = "sk_member_a"
        self.viewer_key = "sk_viewer_a"
        self.other_admin_key = "sk_admin_b"
        self.db.insert_api_key("org-a", hash_key(self.admin_key), "admin", "admin")
        self.db.insert_api_key("org-a", hash_key(self.member_key), "member", "member")
        self.db.insert_api_key("org-a", hash_key(self.viewer_key), "viewer", "viewer")
        self.db.insert_api_key("org-b", hash_key(self.other_admin_key), "b-admin", "admin")
        self.app = create_app(self.db, self.storage, self.settings, self.kp)
        self.client = TestClient(self.app)
        return self.client

    def _auth(self, key):
        return {"Authorization": f"Bearer {key}"}

    def _create_upload(self, payload, meta, key=None):
        h = self._auth(key or self.admin_key)
        h["Upload-Length"] = str(len(payload))
        h["Upload-Metadata"] = tus.encode_metadata(meta)
        return self.client.post("/uploads", headers=h)

    def _upload_id(self, resp):
        self.assertEqual(resp.status_code, 201, resp.text)
        return resp.headers["Location"].rsplit("/", 1)[-1]

    def _full_upload(self, payload=None, key=None, meta=None):
        payload = payload if payload is not None else make_fastq()
        meta = meta or {"filename": "run.fastq", "platform": "illumina"}
        uid = self._upload_id(self._create_upload(payload, meta, key))
        r = self.client.patch(
            f"/uploads/{uid}",
            headers={**self._auth(key or self.admin_key), "Upload-Offset": "0"},
            content=payload,
        )
        self.assertEqual(r.status_code, 204, r.text)
        return uid

    def _raw_bytes(self, r2_key):
        """Ciphertext bytes as stored (bypassing decryption)."""
        return bytes(self.inner._bufs[r2_key])


# ------------------------------------------------------------------ roles
class RolesTest(Phase2Base):
    def setUp(self):
        self.make_app()

    def test_viewer_cannot_upload(self):
        r = self._create_upload(make_fastq(), {"filename": "x.fastq"}, key=self.viewer_key)
        self.assertEqual(r.status_code, 403, r.text)

    def test_viewer_can_read(self):
        uid = self._full_upload()
        r = self.client.get("/jobs", headers=self._auth(self.viewer_key))
        self.assertEqual(r.status_code, 200)
        self.assertTrue(any(j["id"] == uid for j in r.json()["jobs"]))
        r = self.client.get(f"/jobs/{uid}", headers=self._auth(self.viewer_key))
        self.assertEqual(r.status_code, 200)

    def test_viewer_cannot_patch_or_cancel(self):
        uid = self._upload_id(
            self._create_upload(make_fastq(), {"filename": "x.fastq"})
        )
        r = self.client.patch(
            f"/uploads/{uid}",
            headers={**self._auth(self.viewer_key), "Upload-Offset": "0"},
            content=b"@x\nAC\n+\nII\n",
        )
        self.assertEqual(r.status_code, 403)
        r = self.client.delete(f"/uploads/{uid}", headers=self._auth(self.viewer_key))
        self.assertEqual(r.status_code, 403)

    def test_member_can_upload_but_not_manage_keys(self):
        uid = self._full_upload(key=self.member_key)
        self.assertTrue(uid)
        r = self.client.get("/org/keys", headers=self._auth(self.member_key))
        self.assertEqual(r.status_code, 403)
        r = self.client.post(
            "/org/keys", headers=self._auth(self.member_key), json={"name": "x"}
        )
        self.assertEqual(r.status_code, 403)

    def test_member_cannot_rotate_dek(self):
        r = self.client.post(
            "/org/keys/rotate-dek", headers=self._auth(self.member_key)
        )
        self.assertEqual(r.status_code, 403)

    def test_viewer_cannot_delete_job(self):
        uid = self._full_upload()
        r = self.client.delete(f"/jobs/{uid}", headers=self._auth(self.viewer_key))
        self.assertEqual(r.status_code, 403)

    def test_cross_org_stays_404_not_403(self):
        uid = self._full_upload()  # org-a job
        # org-b admin (right role, wrong org) sees 404, not 403.
        for method, path in [
            ("get", f"/jobs/{uid}"),
            ("delete", f"/jobs/{uid}"),
        ]:
            r = getattr(self.client, method)(path, headers=self._auth(self.other_admin_key))
            self.assertEqual(r.status_code, 404, f"{method} {path}: {r.text}")


# ---------------------------------------------------------- key management
class KeyManagementTest(Phase2Base):
    def setUp(self):
        self.make_app()

    def test_create_key_returns_plaintext_once(self):
        r = self.client.post(
            "/org/keys",
            headers=self._auth(self.admin_key),
            json={"name": "lab-uploader", "role": "member"},
        )
        self.assertEqual(r.status_code, 201, r.text)
        body = r.json()
        self.assertTrue(body["api_key"].startswith("sk_"))
        self.assertEqual(body["role"], "member")
        # The new key works...
        r2 = self.client.get("/jobs", headers=self._auth(body["api_key"]))
        self.assertEqual(r2.status_code, 200)
        # ...but the plaintext is never listed again.
        r3 = self.client.get("/org/keys", headers=self._auth(self.admin_key))
        listed = {k["id"]: k for k in r3.json()["keys"]}
        self.assertIn(body["id"], listed)
        self.assertNotIn("api_key", listed[body["id"]])
        self.assertEqual(len(listed[body["id"]]["key_hash_prefix"]), 10)

    def test_create_key_bad_role_400(self):
        r = self.client.post(
            "/org/keys",
            headers=self._auth(self.admin_key),
            json={"name": "x", "role": "superuser"},
        )
        self.assertEqual(r.status_code, 400)

    def test_revoke_key_fails_closed_everywhere(self):
        r = self.client.post(
            "/org/keys",
            headers=self._auth(self.admin_key),
            json={"name": "temp", "role": "member"},
        )
        kid, plaintext = r.json()["id"], r.json()["api_key"]
        # Revoke.
        d = self.client.delete(f"/org/keys/{kid}", headers=self._auth(self.admin_key))
        self.assertEqual(d.status_code, 204)
        # Fails closed on every route class.
        for method, path, extra in [
            ("post", "/uploads", {"Upload-Length": "10"}),
            ("get", "/jobs", {}),
            ("get", "/org/keys", {}),
            ("get", "/org/audit", {}),
        ]:
            h = self._auth(plaintext)
            h.update(extra)
            rr = getattr(self.client, method)(path, headers=h)
            self.assertEqual(rr.status_code, 401, f"{method} {path}: {rr.text}")
        # Double revoke is 409, not silent.
        d2 = self.client.delete(f"/org/keys/{kid}", headers=self._auth(self.admin_key))
        self.assertEqual(d2.status_code, 409)
        # Audit trail has both entries.
        a = self.client.get("/org/audit", headers=self._auth(self.admin_key)).json()["audit"]
        actions = {(e["action"], e["target"]) for e in a}
        self.assertIn(("key.created", kid), actions)
        self.assertIn(("key.revoked", kid), actions)

    def test_revoke_other_org_key_is_404(self):
        # org-b admin's key id, addressed by org-a admin: 404 (oracle rule).
        b_keys = self.db.list_api_keys("org-b")
        r = self.client.delete(
            f"/org/keys/{b_keys[0]['id']}", headers=self._auth(self.admin_key)
        )
        self.assertEqual(r.status_code, 404)


# --------------------------------------------------------------- encryption
class EncryptionTest(Phase2Base):
    def setUp(self):
        self.make_app(encrypted=True)

    def test_ciphertext_on_disk_proven(self):
        payload = make_fastq()
        uid = self._full_upload(payload)
        job = self.db.get_job("org-a", uid)
        raw = self._raw_bytes(job["r2_key"])
        # PROVEN ciphertext: not the plaintext, doesn't look like FASTQ.
        self.assertNotEqual(raw, payload)
        self.assertFalse(raw.startswith(b"@"))
        self.assertGreater(len(raw), len(payload))  # GCM tag appended
        # Sidecar marks the scheme explicitly.
        sidecar = json.loads(self._raw_bytes(job["r2_key"] + ".meta.json"))
        self.assertEqual(sidecar["encryption"], "aes256-gcm-v1")
        self.assertEqual(sidecar["dek_version"], 1)
        self.assertEqual(sidecar["sha256_plaintext"],
                         __import__("hashlib").sha256(payload).hexdigest())

    def test_round_trip_byte_identical(self):
        payload = make_fastq(50)
        uid = self._full_upload(payload)
        job = self.db.get_job("org-a", uid)
        with tempfile.TemporaryDirectory() as d:
            dest = os.path.join(d, "out.fastq")
            self.storage.download_to_file(job["r2_key"], dest)
            with open(dest, "rb") as f:
                self.assertEqual(f.read(), payload)

    def test_tampered_ciphertext_fails_closed(self):
        payload = make_fastq()
        uid = self._full_upload(payload)
        job = self.db.get_job("org-a", uid)
        # Flip a byte in the stored ciphertext.
        buf = self.inner._bufs[job["r2_key"]]
        buf[100] ^= 0xFF
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(DecryptionError):
                self.storage.download_to_file(job["r2_key"], os.path.join(d, "x"))

    def test_missing_sidecar_fails_closed_no_plaintext_fallback(self):
        payload = make_fastq()
        uid = self._full_upload(payload)
        job = self.db.get_job("org-a", uid)
        del self.inner._bufs[job["r2_key"] + ".meta.json"]
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(DecryptionError):
                self.storage.download_to_file(job["r2_key"], os.path.join(d, "x"))

    def test_wrap_unwrap_round_trip(self):
        dek = generate_dek()
        wrapped = wrap_dek(dek, TEST_KEK)
        self.assertEqual(unwrap_dek(wrapped, TEST_KEK), dek)
        # Wrong KEK -> authentication failure, not garbage.
        with self.assertRaises(DecryptionError):
            unwrap_dek(wrapped, bytes(range(1, 33)))

    def test_load_kek_fail_closed(self):
        old = os.environ.get("SENTINEL_KEK")
        try:
            os.environ.pop("SENTINEL_KEK", None)
            with self.assertRaises(RuntimeError):
                load_kek()
            os.environ["SENTINEL_KEK"] = "not-base64!!!"
            with self.assertRaises(RuntimeError):
                load_kek()
            import base64
            os.environ["SENTINEL_KEK"] = base64.b64encode(b"short").decode()
            with self.assertRaises(RuntimeError):
                load_kek()
            os.environ["SENTINEL_KEK"] = base64.b64encode(TEST_KEK).decode()
            self.assertEqual(load_kek(), TEST_KEK)
        finally:
            if old is None:
                os.environ.pop("SENTINEL_KEK", None)
            else:
                os.environ["SENTINEL_KEK"] = old


# ----------------------------------------------------------------- rotation
class RotationTest(Phase2Base):
    def setUp(self):
        self.make_app(encrypted=True)

    def test_rotate_new_uploads_v2_old_still_decrypt(self):
        p1 = make_fastq(10)
        uid1 = self._full_upload(p1)
        job1 = self.db.get_job("org-a", uid1)

        r = self.client.post("/org/keys/rotate-dek", headers=self._auth(self.admin_key))
        self.assertEqual(r.status_code, 200, r.text)
        self.assertEqual(r.json()["dek_version"], 2)

        p2 = make_fastq(15)
        uid2 = self._full_upload(p2)
        job2 = self.db.get_job("org-a", uid2)

        m1 = json.loads(self._raw_bytes(job1["r2_key"] + ".meta.json"))
        m2 = json.loads(self._raw_bytes(job2["r2_key"] + ".meta.json"))
        self.assertEqual(m1["dek_version"], 1)
        self.assertEqual(m2["dek_version"], 2)

        # Both decrypt (mixed versions handled).
        with tempfile.TemporaryDirectory() as d:
            d1 = os.path.join(d, "a.fastq")
            d2 = os.path.join(d, "b.fastq")
            self.storage.download_to_file(job1["r2_key"], d1)
            self.storage.download_to_file(job2["r2_key"], d2)
            with open(d1, "rb") as f:
                self.assertEqual(f.read(), p1)
            with open(d2, "rb") as f:
                self.assertEqual(f.read(), p2)

        # Audit logged the rotation.
        a = self.client.get("/org/audit", headers=self._auth(self.admin_key)).json()["audit"]
        self.assertTrue(any(
            e["action"] == "dek.rotated" and e["detail"].get("dek_version") == 2
            for e in a
        ))

    def test_rotate_without_dek_starts_at_v1(self):
        # Fresh org (no DEK yet): rotate creates v1, then v2... first call
        # creates v1 implicitly and rotates to v2? No: rotate() on empty
        # starts the chain at 1, then the new version is 2? Check semantics:
        # rotate() with no current key -> nxt=1 -> stores v1 -> returns 1.
        self.db.insert_api_key("org-c", hash_key("sk_c"), "c", "admin")
        r = self.client.post(
            "/org/keys/rotate-dek", headers=self._auth("sk_c")
        )
        self.assertEqual(r.status_code, 200)
        self.assertEqual(r.json()["dek_version"], 1)


# ----------------------------------------------------------------- deletion
class DeletionTest(Phase2Base):
    def setUp(self):
        self.make_app(encrypted=True)

    def test_delete_job_removes_bytes_row_and_audits(self):
        payload = make_fastq()
        uid = self._full_upload(payload)
        job = self.db.get_job("org-a", uid)
        r2_key = job["r2_key"]
        # Mirror a sentinel_runs row as the worker would.
        self.db.mirror_sentinel_run({"accession": uid, "priority": "upload"})

        r = self.client.delete(f"/jobs/{uid}", headers=self._auth(self.member_key))
        self.assertEqual(r.status_code, 204, r.text)

        self.assertIsNone(self.db.get_job("org-a", uid))
        self.assertEqual(self.storage.list_prefix(f"intake/org-a/"), [])
        # Mirror row kept, tombstoned.
        self.assertEqual(len(self.db.mirrored), 1)
        self.assertIsNotNone(self.db.mirrored[0].get("deleted_at"))
        # Audit entry.
        a = self.client.get("/org/audit", headers=self._auth(self.admin_key)).json()["audit"]
        self.assertTrue(any(
            e["action"] == "job.deleted" and e["target"] == uid for e in a
        ))

    def test_delete_job_other_org_404(self):
        uid = self._full_upload()
        r = self.client.delete(f"/jobs/{uid}", headers=self._auth(self.other_admin_key))
        self.assertEqual(r.status_code, 404)
        self.assertIsNotNone(self.db.get_job("org-a", uid))

    def test_purge_wrong_confirm_400_untouched(self):
        uid = self._full_upload()
        job = self.db.get_job("org-a", uid)
        n_keys_before = len(self.db.list_api_keys("org-a"))
        r = self.client.request(
            "DELETE", "/org/data",
            headers=self._auth(self.admin_key), json={"confirm": "org-b"},
        )
        self.assertEqual(r.status_code, 400, r.text)
        # Nothing touched.
        self.assertIsNotNone(self.db.get_job("org-a", uid))
        self.assertTrue(self.storage.list_prefix("intake/org-a/"))
        self.assertEqual(len(self.db.list_api_keys("org-a")), n_keys_before)

    def test_purge_correct_confirm_zeroes_everything(self):
        uid1 = self._full_upload(make_fastq(5))
        uid2 = self._full_upload(make_fastq(7))
        r = self.client.request(
            "DELETE", "/org/data",
            headers=self._auth(self.admin_key), json={"confirm": "org-a"},
        )
        self.assertEqual(r.status_code, 200, r.text)
        body = r.json()
        self.assertEqual(body["org_id"], "org-a")
        self.assertGreaterEqual(body["objects_deleted"], 4)  # 2 objects + 2 sidecars
        self.assertEqual(body["jobs_deleted"], 2)
        # Zero bytes, zero rows.
        self.assertEqual(self.storage.list_prefix("intake/org-a/"), [])
        self.assertEqual(self.db.list_jobs("org-a"), [])
        # All keys revoked: old admin key now 401s.
        r2 = self.client.get("/jobs", headers=self._auth(self.admin_key))
        self.assertEqual(r2.status_code, 401)
        # Other org untouched.
        self.assertEqual(self.db.list_api_keys("org-b")[0]["revoked_at"], None)
        # Audit: the purge itself is logged (written before keys die).
        a = [e for e in self.db.audit_rows
             if e["org_id"] == "org-a" and e["action"] == "org.purged"]
        self.assertEqual(len(a), 1)
        self.assertEqual(a[0]["detail"]["jobs_deleted"], 2)

    def test_purge_requires_admin(self):
        r = self.client.request(
            "DELETE", "/org/data",
            headers=self._auth(self.member_key), json={"confirm": "org-a"},
        )
        self.assertEqual(r.status_code, 403)

    def test_audit_pagination_and_admin_only(self):
        for i in range(3):
            self.client.post(
                "/org/keys", headers=self._auth(self.admin_key),
                json={"name": f"k{i}", "role": "viewer"},
            )
        r = self.client.get(
            "/org/audit?limit=2&offset=1", headers=self._auth(self.admin_key)
        )
        self.assertEqual(r.status_code, 200)
        body = r.json()
        self.assertEqual(len(body["audit"]), 2)
        self.assertEqual(body["limit"], 2)
        self.assertEqual(body["offset"], 1)
        # Member cannot read the audit log.
        r2 = self.client.get("/org/audit", headers=self._auth(self.member_key))
        self.assertEqual(r2.status_code, 403)


# ------------------------------------------------------- worker integration
class WorkerDecryptTest(unittest.TestCase):
    """The worker's exact download path against EncryptedStorage."""

    def setUp(self):
        self.db = FakeDB()
        self.inner = MemoryStorage()
        self.kp = KeyProvider(self.db, TEST_KEK)
        self.storage = EncryptedStorage(self.inner, self.kp)
        self.tmp = tempfile.mkdtemp(prefix="sentinel-p2-worker-")
        self.settings = Settings()
        self.settings.worker_workdir = os.path.join(self.tmp, "work")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _seed_encrypted_queued_job(self, payload: bytes) -> dict:
        import hashlib

        key = "intake/org-a/JOB1/run.fastq"
        self.storage.create(key)
        self.storage.append(key, payload)
        self.storage.finalize(key)
        job = self.db.create_job({
            "org_id": "org-a", "filename": "run.fastq",
            "size_bytes": len(payload), "platform": "illumina",
            "tier": "lite", "status": "queued", "r2_key": key,
            "sha256": hashlib.sha256(payload).hexdigest(),
        })
        return job

    def test_worker_decrypts_before_pipeline(self):
        from service import worker as worker_mod

        payload = make_fastq(30)
        job = self._seed_encrypted_queued_job(payload)
        seen = {}

        def stub_run_pipeline(workdir, cmd):
            # The pipeline input must be byte-identical plaintext.
            cfg_path = cmd[cmd.index("--configfile") + 1]
            import yaml
            with open(cfg_path) as f:
                cfg = yaml.safe_load(f)
            with open(cfg["reads_r1"], "rb") as f:
                seen["input"] = f.read()
            # Fake a pipeline report.
            outdir = cfg["outdir"]
            os.makedirs(outdir, exist_ok=True)
            with open(os.path.join(outdir, "report.json"), "w") as f:
                json.dump({"version": "0.2.0", "reads": {"raw": 1}}, f)
            return 0, "ok"

        self.assertTrue(
            worker_mod.process_one_job(self.db, self.storage, self.settings,
                                       run_pipeline=stub_run_pipeline)
        )
        self.assertEqual(seen["input"], payload)
        done = self.db.get_job("org-a", job["id"])
        self.assertEqual(done["status"], "done")
        self.assertEqual(done["report"]["version"], "0.2.0")

    def test_worker_tampered_object_fails_job(self):
        from service import worker as worker_mod

        payload = make_fastq(10)
        job = self._seed_encrypted_queued_job(payload)
        buf = self.inner._bufs[job["r2_key"]]
        buf[50] ^= 0xFF  # tamper with ciphertext
        self.assertTrue(
            worker_mod.process_one_job(self.db, self.storage, self.settings,
                                       run_pipeline=lambda w, c: (0, "ok"))
        )
        failed = self.db.get_job("org-a", job["id"])
        self.assertEqual(failed["status"], "failed")
        self.assertIn("decryption failed", failed["status_reason"])


# ----------------------------------------------------------------- migration
class MigrationTest(unittest.TestCase):
    def setUp(self):
        self.db = FakeDB()
        self.raw = MemoryStorage()  # Phase 1 plaintext, written directly
        self.kp = KeyProvider(self.db, TEST_KEK)
        self.enc = EncryptedStorage(self.raw, self.kp)

    def _plant_plaintext(self, org, name, payload, status="done"):
        key = f"intake/{org}/job-{name}/{name}.fastq"
        self.raw.create(key)
        self.raw.append(key, payload)
        self.raw.finalize(key)
        job = self.db.create_job({
            "org_id": org, "filename": f"{name}.fastq",
            "size_bytes": len(payload), "platform": "illumina",
            "tier": "lite", "status": status, "r2_key": key,
        })
        return job, key

    def test_dry_run_changes_nothing(self):
        p = make_fastq()
        _, key = self._plant_plaintext("org-a", "s1", p)
        report = migrate(self.db, self.raw, self.enc, dry_run=True)
        self.assertEqual(report["orgs"]["org-a"]["would_encrypt"], [key])
        self.assertEqual(report["orgs"]["org-a"]["encrypted"], [])
        # Bytes untouched, no sidecar.
        self.assertEqual(bytes(self.raw._bufs[key]), p)
        self.assertNotIn(key + ".meta.json", self.raw._bufs)

    def test_real_run_encrypts_and_skips(self):
        p1, p2 = make_fastq(5), make_fastq(8)
        _, k1 = self._plant_plaintext("org-a", "s1", p1)
        _, k2 = self._plant_plaintext("org-a", "s2", p2, status="uploading")  # partial: skip
        # Pre-encrypted object: skipped as already done.
        self.enc.create("intake/org-a/job-s3/s3.fastq")
        self.enc.append("intake/org-a/job-s3/s3.fastq", make_fastq(3))
        self.enc.finalize("intake/org-a/job-s3/s3.fastq")

        report = migrate(self.db, self.raw, self.enc, org_ids=["org-a"])
        org = report["orgs"]["org-a"]
        self.assertEqual(org["encrypted"], [k1])
        self.assertEqual(org["skipped_partial"], [k2])
        self.assertEqual(org["already_encrypted"], ["intake/org-a/job-s3/s3.fastq"])

        # k1 is now ciphertext with a sidecar, and decrypts back.
        self.assertNotEqual(bytes(self.raw._bufs[k1]), p1)
        meta = json.loads(bytes(self.raw._bufs[k1 + ".meta.json"]))
        self.assertEqual(meta["encryption"], "aes256-gcm-v1")
        with tempfile.TemporaryDirectory() as d:
            dest = os.path.join(d, "out.fastq")
            self.enc.download_to_file(k1, dest)
            with open(dest, "rb") as f:
                self.assertEqual(f.read(), p1)
        # k2 (partial) untouched.
        self.assertEqual(bytes(self.raw._bufs[k2]), p2)
        self.assertNotIn(k2 + ".meta.json", self.raw._bufs)

    def test_discover_orgs(self):
        self._plant_plaintext("org-a", "s1", make_fastq(2))
        self._plant_plaintext("org-b", "s1", make_fastq(2))
        self.assertEqual(discover_orgs(self.raw), ["org-a", "org-b"])


if __name__ == "__main__":
    unittest.main()
