"""Object storage abstraction for upload intake.

The tus PATCH handler appends chunks through this interface and never holds
the whole file in app memory. Two implementations:

- MemoryStorage: in-memory (tests only).
- R2Storage: S3-compatible (Cloudflare R2) via boto3. Chunks are staged to
  a local disk file per upload and PUT to R2 once, at completion — so a
  multi-GB upload never sits in app memory, and a crashed API process can
  resume from the staged bytes (the DB row's upload_offset is the source of
  truth). This is a deliberate Phase 1 simplification over S3 multipart
  streaming; the interface supports swapping in a multipart implementation
  later without touching the API or worker code.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
from typing import Optional


class Storage:
    """Chunk-append object store keyed by the intake R2 key."""

    def create(self, key: str) -> None:
        """Initialize an empty object for a new upload."""
        raise NotImplementedError

    def append(self, key: str, data: bytes) -> None:
        raise NotImplementedError

    def finalize(self, key: str) -> None:
        """Mark the object complete (R2: single PUT of the staged bytes)."""
        raise NotImplementedError

    def abort(self, key: str) -> None:
        """Discard a partial upload."""
        raise NotImplementedError

    def delete(self, key: str) -> None:
        raise NotImplementedError

    def size(self, key: str) -> int:
        raise NotImplementedError

    def sha256(self, key: str) -> str:
        """Hex digest over the bytes received so far."""
        raise NotImplementedError

    def download_to_file(self, key: str, dest_path: str) -> None:
        raise NotImplementedError

    def list_prefix(self, prefix: str) -> list[str]:
        """All keys starting with prefix (for org purge / migration)."""
        raise NotImplementedError


class MemoryStorage(Storage):
    """In-memory store. Tests only — never holds more than test fixtures."""

    def __init__(self) -> None:
        self._bufs: dict[str, bytearray] = {}
        self._finalized: set[str] = set()

    def create(self, key: str) -> None:
        self._bufs[key] = bytearray()

    def append(self, key: str, data: bytes) -> None:
        self._bufs[key].extend(data)

    def finalize(self, key: str) -> None:
        if key not in self._bufs:
            raise KeyError(key)
        self._finalized.add(key)

    def abort(self, key: str) -> None:
        self._bufs.pop(key, None)
        self._finalized.discard(key)

    def delete(self, key: str) -> None:
        self.abort(key)

    def size(self, key: str) -> int:
        return len(self._bufs[key])

    def sha256(self, key: str) -> str:
        return hashlib.sha256(bytes(self._bufs[key])).hexdigest()

    def download_to_file(self, key: str, dest_path: str) -> None:
        with open(dest_path, "wb") as f:
            f.write(bytes(self._bufs[key]))

    def list_prefix(self, prefix: str) -> list[str]:
        return sorted(k for k in self._bufs if k.startswith(prefix))

    # test introspection
    def exists(self, key: str) -> bool:
        return key in self._bufs


class R2Storage(Storage):
    """R2 via boto3, with local disk staging.

    Layout on R2: the exact `r2_key` the API assigned at creation
    (intake/{org_id}/{job_id}/{filename}). Staging lives under `stage_dir`
    (default: the system temp dir); one file per in-flight upload.
    """

    def __init__(
        self,
        endpoint_url: str,
        access_key_id: str,
        secret_access_key: str,
        bucket: str,
        stage_dir: Optional[str] = None,
    ):
        import boto3  # lazy: tests never import boto3

        self._bucket = bucket
        self._stage_dir = stage_dir or tempfile.gettempdir()
        os.makedirs(self._stage_dir, exist_ok=True)
        self._s3 = boto3.client(
            "s3",
            endpoint_url=endpoint_url,
            aws_access_key_id=access_key_id,
            aws_secret_access_key=secret_access_key,
        )

    def _stage_path(self, key: str) -> str:
        # Flat staging name: sha of the key (keys contain slashes).
        name = "sentinel-upload-" + hashlib.sha256(key.encode()).hexdigest()
        return os.path.join(self._stage_dir, name)

    def create(self, key: str) -> None:
        path = self._stage_path(key)
        open(path, "wb").close()

    def append(self, key: str, data: bytes) -> None:
        with open(self._stage_path(key), "ab") as f:
            f.write(data)

    def finalize(self, key: str) -> None:
        path = self._stage_path(key)
        with open(path, "rb") as f:
            self._s3.upload_fileobj(f, self._bucket, key)
        os.remove(path)

    def abort(self, key: str) -> None:
        path = self._stage_path(key)
        if os.path.exists(path):
            os.remove(path)

    def delete(self, key: str) -> None:
        self._s3.delete_object(Bucket=self._bucket, Key=key)

    def size(self, key: str) -> int:
        path = self._stage_path(key)
        if os.path.exists(path):
            return os.path.getsize(path)
        head = self._s3.head_object(Bucket=self._bucket, Key=key)
        return int(head["ContentLength"])

    def sha256(self, key: str) -> str:
        h = hashlib.sha256()
        path = self._stage_path(key)
        if os.path.exists(path):
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(1 << 20), b""):
                    h.update(chunk)
            return h.hexdigest()
        # Fall back to streaming the finalized object (rare; completion
        # normally hashes the staged file before finalize()).
        obj = self._s3.get_object(Bucket=self._bucket, Key=key)
        try:
            for chunk in obj["Body"].iter_chunks(1 << 20):
                h.update(chunk)
        finally:
            obj["Body"].close()
        return h.hexdigest()

    def download_to_file(self, key: str, dest_path: str) -> None:
        tmp = dest_path + ".part"
        self._s3.download_file(self._bucket, key, tmp)
        shutil.move(tmp, dest_path)

    def list_prefix(self, prefix: str) -> list[str]:
        paginator = self._s3.get_paginator("list_objects_v2")
        keys: list[str] = []
        for page in paginator.paginate(Bucket=self._bucket, Prefix=prefix):
            for obj in page.get("Contents", []):
                keys.append(obj["Key"])
        return sorted(keys)
