"""Upload metadata validation.

Mirrors workflow/sentinel_lite/report_lib.validate_platform_config so the
API rejects at creation time exactly what the pipeline would refuse at parse
time (fail-closed both ends). A parity test in service/tests/test_tus.py
asserts both accept/reject the same matrix of inputs; if the pipeline's
rules change, that test fails and this module must be updated.
"""

from __future__ import annotations

VALID_PLATFORMS = ("illumina", "nanopore")
VALID_BASECALLERS = ("fast", "hac", "sup")
VALID_TIERS = ("lite", "assembly")


def validate_upload_params(params: dict[str, str]) -> tuple[str, str | None, str]:
    """(platform, basecaller_model, tier) from tus Upload-Metadata, or raise ValueError."""
    platform = params.get("platform", "illumina")
    if platform not in VALID_PLATFORMS:
        raise ValueError(f"unknown platform {platform!r}; expected one of {VALID_PLATFORMS}")
    tier = params.get("tier", "lite")
    if tier not in VALID_TIERS:
        raise ValueError(f"unknown tier {tier!r}; expected one of {VALID_TIERS}")
    basecaller = params.get("basecaller_model")
    if platform == "nanopore":
        if not basecaller:
            raise ValueError(
                f"platform=nanopore requires basecaller_model (one of {VALID_BASECALLERS}); "
                "refusing upload"
            )
        if basecaller not in VALID_BASECALLERS:
            raise ValueError(
                f"unknown basecaller_model {basecaller!r}; expected one of {VALID_BASECALLERS}"
            )
    if tier == "assembly" and platform != "nanopore":
        raise ValueError("tier=assembly is nanopore-only in Phase 0")
    return platform, basecaller, tier
