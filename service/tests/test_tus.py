"""tus protocol tests: round-trip with simulated interruption/resume,
fail-closed validation, auth isolation, cancel, checksum.

Conventions mirror tests/sentinel_lite/: unittest, no network, no
credentials. Run: python -m unittest discover -s service/tests -t .
"""

from __future__ import annotations

import hashlib
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "workflow" / "sentinel_lite"))

from fastapi.testclient import TestClient  # noqa: E402

import report_lib  # noqa: E402  (pipeline-side validation, for the parity test)
from service.app import tus  # noqa: E402
from service.app.auth import hash_key  # noqa: E402
from service.app.config import Settings  # noqa: E402
from service.app.main import create_app  # noqa: E402
from service.app.storage import MemoryStorage  # noqa: E402
from service.app.validation import validate_upload_params  # noqa: E402
from service.tests.fakes import FakeDB  # noqa: E402


def make_fastq(n_reads: int = 50) -> bytes:
    recs = []
    for i in range(n_reads):
        seq = "ACGT" * 25
        recs.append(f"@read{i}\n{seq}\n+\n{'I' * 100}\n")
    return "".join(recs).encode()


class TusTestCase(unittest.TestCase):
    def setUp(self):
        self.db = FakeDB()
        self.storage = MemoryStorage()
        self.settings = Settings()
        # org A and org B keys
        self.key_a = "sk_test_org_a"
        self.key_b = "sk_test_org_b"
        self.db.insert_api_key("org-a", hash_key(self.key_a), "a")
        self.db.insert_api_key("org-b", hash_key(self.key_b), "b")
        self.app = create_app(self.db, self.storage, self.settings)
        self.client = TestClient(self.app)

    def _headers_a(self, extra=None):
        h = {"Authorization": f"Bearer {self.key_a}"}
        if extra:
            h.update(extra)
        return h

    def _create(self, payload: bytes, meta: dict, headers=None, key=None):
        h = {"Authorization": f"Bearer {key or self.key_a}"}
        h["Upload-Length"] = str(len(payload))
        h["Upload-Metadata"] = tus.encode_metadata(meta)
        if headers:
            h.update(headers)
        return self.client.post("/uploads", headers=h)

    def _upload_id(self, resp) -> str:
        self.assertEqual(resp.status_code, 201, resp.text)
        return resp.headers["Location"].rsplit("/", 1)[-1]

    # -- happy path ------------------------------------------------------
    def test_create_returns_201_location_and_zero_offset(self):
        r = self._create(b"", {"filename": "x.fastq.gz", "platform": "illumina"})
        uid = self._upload_id(r)
        self.assertEqual(r.headers["Upload-Offset"], "0")
        head = self.client.head(f"/uploads/{uid}", headers=self._headers_a())
        self.assertEqual(head.status_code, 200)
        self.assertEqual(head.headers["Upload-Offset"], "0")

    def test_chunked_round_trip_with_interruption_and_resume(self):
        payload = make_fastq(200)
        meta = {
            "filename": "run.fastq.gz",
            "platform": "nanopore",
            "basecaller_model": "hac",
            "tier": "assembly",
        }
        uid = self._upload_id(self._create(payload, meta))

        # Chunk 1 of 3, then "the connection drops" (client just stops).
        c1, c2, c3 = payload[:1000], payload[1000:3000], payload[3000:]
        r = self.client.patch(
            f"/uploads/{uid}",
            content=c1,
            headers=self._headers_a({"Upload-Offset": "0",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 204, r.text)
        self.assertEqual(r.headers["Upload-Offset"], str(len(c1)))
        r = self.client.patch(
            f"/uploads/{uid}",
            content=c2,
            headers=self._headers_a({"Upload-Offset": str(len(c1)),
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 204, r.text)

        # Resume: HEAD reports the offset; client continues from there.
        head = self.client.head(f"/uploads/{uid}", headers=self._headers_a())
        self.assertEqual(head.headers["Upload-Offset"], str(len(c1) + len(c2)))

        r = self.client.patch(
            f"/uploads/{uid}",
            content=c3,
            headers=self._headers_a({"Upload-Offset": str(len(c1) + len(c2)),
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 204, r.text)
        self.assertEqual(r.headers["Upload-Offset"], str(len(payload)))

        # Byte-identical round-trip, sha256 recorded, job queued.
        job = self.db.jobs[uid]
        self.assertEqual(job["status"], "queued")
        self.assertEqual(job["upload_offset"], len(payload))
        self.assertEqual(job["sha256"], hashlib.sha256(payload).hexdigest())
        self.assertEqual(job["platform"], "nanopore")
        self.assertEqual(job["basecaller_model"], "hac")
        self.assertEqual(job["tier"], "assembly")
        got = self.client.get(f"/jobs/{uid}", headers=self._headers_a()).json()
        self.assertEqual(got["status"], "queued")

    def test_offset_mismatch_is_409(self):
        payload = make_fastq(10)
        uid = self._upload_id(self._create(payload, {"platform": "illumina"}))
        r = self.client.patch(
            f"/uploads/{uid}",
            content=b"@x",
            headers=self._headers_a({"Upload-Offset": "999",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 409)

    def test_chunk_overrun_is_400(self):
        payload = make_fastq(10)
        uid = self._upload_id(self._create(payload, {"platform": "illumina"}))
        r = self.client.patch(
            f"/uploads/{uid}",
            content=payload + b"EXTRA",
            headers=self._headers_a({"Upload-Offset": "0",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 400)

    # -- fail-closed ------------------------------------------------------
    def test_nanopore_without_basecaller_is_400(self):
        r = self._create(b"", {"platform": "nanopore", "filename": "x.fastq"})
        self.assertEqual(r.status_code, 400, r.text)

    def test_bad_basecaller_is_400(self):
        r = self._create(b"", {"platform": "nanopore", "basecaller_model": "ultra"})
        self.assertEqual(r.status_code, 400, r.text)

    def test_assembly_on_illumina_is_400(self):
        r = self._create(b"", {"platform": "illumina", "tier": "assembly"})
        self.assertEqual(r.status_code, 400, r.text)

    def test_unknown_platform_is_400(self):
        r = self._create(b"", {"platform": "ont"})
        self.assertEqual(r.status_code, 400, r.text)

    def test_missing_upload_length_is_400(self):
        h = self._headers_a({"Upload-Metadata": tus.encode_metadata({"platform": "illumina"})})
        r = self.client.post("/uploads", headers=h)
        self.assertEqual(r.status_code, 400)

    def test_size_cap_is_413(self):
        self.settings.max_upload_bytes = 10
        r = self._create(make_fastq(10), {"platform": "illumina"})
        self.assertEqual(r.status_code, 413)

    def test_non_fastq_first_bytes_rejected(self):
        payload = b"definitely not fastq" * 10
        uid = self._upload_id(self._create(payload, {"platform": "illumina"}))
        r = self.client.patch(
            f"/uploads/{uid}",
            content=payload,
            headers=self._headers_a({"Upload-Offset": "0",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 400)
        self.assertEqual(self.db.jobs[uid]["status"], "failed")

    def test_checksum_mismatch_fails_job(self):
        payload = make_fastq(20)
        meta = {"platform": "illumina", "checksum": "0" * 64}  # wrong on purpose
        uid = self._upload_id(self._create(payload, meta))
        r = self.client.patch(
            f"/uploads/{uid}",
            content=payload,
            headers=self._headers_a({"Upload-Offset": "0",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 400, r.text)
        job = self.db.jobs[uid]
        self.assertEqual(job["status"], "failed")
        self.assertIn("checksum", job["status_reason"])

    def test_checksum_match_passes(self):
        payload = make_fastq(20)
        meta = {"platform": "illumina",
                "checksum": hashlib.sha256(payload).hexdigest()}
        uid = self._upload_id(self._create(payload, meta))
        r = self.client.patch(
            f"/uploads/{uid}",
            content=payload,
            headers=self._headers_a({"Upload-Offset": "0",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(r.status_code, 204, r.text)
        self.assertEqual(self.db.jobs[uid]["status"], "queued")

    def test_validation_parity_with_pipeline(self):
        """The API must reject exactly what report_lib.validate_platform_config
        rejects. If the pipeline's rules change, this fails and
        service/app/validation.py must be updated."""
        matrix = [
            ({}, True),
            ({"platform": "illumina"}, True),
            ({"platform": "illumina", "tier": "lite"}, True),
            ({"platform": "nanopore"}, False),
            ({"platform": "nanopore", "basecaller_model": "hac"}, True),
            ({"platform": "nanopore", "basecaller_model": "hac", "tier": "assembly"}, True),
            ({"platform": "nanopore", "basecaller_model": "ultra"}, False),
            ({"platform": "illumina", "tier": "assembly"}, False),
            ({"platform": "ont"}, False),
            ({"tier": "turbo"}, False),
        ]
        for params, ok in matrix:
            with self.subTest(params=params):
                try:
                    validate_upload_params(params)
                    api_ok = True
                except ValueError:
                    api_ok = False
                try:
                    report_lib.validate_platform_config(params)
                    pipe_ok = True
                except ValueError:
                    pipe_ok = False
                self.assertEqual(api_ok, pipe_ok, f"parity break on {params}")
                self.assertEqual(api_ok, ok, f"wrong verdict on {params}")

    # -- auth --------------------------------------------------------------
    def test_no_key_is_401(self):
        r = self.client.post("/uploads", headers={"Upload-Length": "0"})
        self.assertEqual(r.status_code, 401)

    def test_wrong_key_is_401(self):
        r = self.client.post(
            "/uploads",
            headers={"Authorization": "Bearer sk_wrong", "Upload-Length": "0"},
        )
        self.assertEqual(r.status_code, 401)

    def test_revoked_key_is_401(self):
        self.db.revoke_key(hash_key(self.key_a))
        r = self.client.get("/jobs", headers=self._headers_a())
        self.assertEqual(r.status_code, 401)

    def test_org_isolation_no_oracle(self):
        """Org B gets 404 (not 403) on org A's resources; B's list is empty."""
        payload = make_fastq(5)
        uid = self._upload_id(self._create(payload, {"platform": "illumina"}))
        b = {"Authorization": f"Bearer {self.key_b}"}
        self.assertEqual(self.client.get(f"/jobs/{uid}", headers=b).status_code, 404)
        self.assertEqual(self.client.head(f"/uploads/{uid}", headers=b).status_code, 404)
        r = self.client.patch(
            f"/uploads/{uid}", content=b"@x",
            headers={**b, "Upload-Offset": "0",
                     "Content-Type": "application/offset+octet-stream"},
        )
        self.assertEqual(r.status_code, 404)
        self.assertEqual(self.client.delete(f"/uploads/{uid}", headers=b).status_code, 404)
        self.assertEqual(self.client.get("/jobs", headers=b).json(), {"jobs": []})

    # -- cancel / retry -------------------------------------------------------
    def test_delete_cancels_and_discards_bytes(self):
        payload = make_fastq(50)
        uid = self._upload_id(self._create(payload, {"platform": "illumina"}))
        job = self.db.jobs[uid]
        r2_key = job["r2_key"]
        self.client.patch(
            f"/uploads/{uid}", content=payload[:500],
            headers=self._headers_a({"Upload-Offset": "0",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertTrue(self.storage.exists(r2_key))
        r = self.client.delete(f"/uploads/{uid}", headers=self._headers_a())
        self.assertEqual(r.status_code, 204)
        self.assertEqual(self.db.jobs[uid]["status"], "cancelled")
        self.assertFalse(self.storage.exists(r2_key))

    def test_retry_requeues_failed_job(self):
        payload = make_fastq(5)
        meta = {"platform": "illumina", "checksum": "0" * 64}
        uid = self._upload_id(self._create(payload, meta))
        self.client.patch(
            f"/uploads/{uid}", content=payload,
            headers=self._headers_a({"Upload-Offset": "0",
                                     "Content-Type": "application/offset+octet-stream"}),
        )
        self.assertEqual(self.db.jobs[uid]["status"], "failed")
        r = self.client.post(f"/jobs/{uid}/retry", headers=self._headers_a())
        self.assertEqual(r.status_code, 200)
        self.assertEqual(self.db.jobs[uid]["status"], "queued")

    def test_retry_non_failed_is_409(self):
        uid = self._upload_id(self._create(b"", {"platform": "illumina"}))
        r = self.client.post(f"/jobs/{uid}/retry", headers=self._headers_a())
        self.assertEqual(r.status_code, 409)

    def test_jobs_list_newest_first(self):
        u1 = self._upload_id(self._create(b"", {"platform": "illumina"}))
        u2 = self._upload_id(self._create(b"", {"platform": "illumina"}))
        jobs = self.client.get("/jobs", headers=self._headers_a()).json()["jobs"]
        self.assertEqual([j["id"] for j in jobs], [u2, u1])


if __name__ == "__main__":
    unittest.main()
