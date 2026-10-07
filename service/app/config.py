"""Service configuration, from environment only. No secrets baked anywhere."""

from __future__ import annotations

import os
from dataclasses import dataclass, field


def _int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except ValueError:
        return default


@dataclass
class Settings:
    # Supabase (PostgREST; same project as sentinel_runs). Service-role key:
    # the API and the worker both bypass RLS with it.
    supabase_url: str = field(default_factory=lambda: os.environ.get("SUPABASE_URL", ""))
    supabase_service_role_key: str = field(
        default_factory=lambda: os.environ.get("SUPABASE_SERVICE_ROLE_KEY", "")
    )

    # R2 (S3-compatible).
    r2_endpoint_url: str = field(default_factory=lambda: os.environ.get("R2_ENDPOINT_URL", ""))
    r2_access_key_id: str = field(default_factory=lambda: os.environ.get("R2_ACCESS_KEY_ID", ""))
    r2_secret_access_key: str = field(
        default_factory=lambda: os.environ.get("R2_SECRET_ACCESS_KEY", "")
    )
    r2_bucket: str = field(default_factory=lambda: os.environ.get("R2_BUCKET", ""))

    # Upload policy.
    max_upload_bytes: int = field(default_factory=lambda: _int("MAX_UPLOAD_BYTES", 10_000_000_000))

    # Pipeline invocation (worker).
    repo_dir: str = field(default_factory=lambda: os.environ.get("REPO_DIR", "/app"))
    pipeline_refdir: str = field(
        default_factory=lambda: os.environ.get("PIPELINE_REFDIR", "/app/data/references")
    )
    snakemake_cores: int = field(default_factory=lambda: _int("SNAKEMAKE_CORES", 4))
    worker_workdir: str = field(default_factory=lambda: os.environ.get("WORKER_WORKDIR", "/tmp/sentinel-work"))
    worker_poll_interval: int = field(default_factory=lambda: _int("WORKER_POLL_INTERVAL", 10))

    @property
    def snakefile(self) -> str:
        return os.path.join(self.repo_dir, "workflow", "sentinel_lite", "Snakefile")

    def require_db(self) -> None:
        if not self.supabase_url or not self.supabase_service_role_key:
            raise RuntimeError("SUPABASE_URL and SUPABASE_SERVICE_ROLE_KEY must be set")

    def require_r2(self) -> None:
        missing = [
            n
            for n, v in (
                ("R2_ENDPOINT_URL", self.r2_endpoint_url),
                ("R2_ACCESS_KEY_ID", self.r2_access_key_id),
                ("R2_SECRET_ACCESS_KEY", self.r2_secret_access_key),
                ("R2_BUCKET", self.r2_bucket),
            )
            if not v
        ]
        if missing:
            raise RuntimeError(f"missing R2 configuration: {', '.join(missing)}")
