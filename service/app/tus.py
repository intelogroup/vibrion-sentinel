"""tus 1.0.0 protocol helpers (creation core only).

Wire format notes (https://tus.io/protocols/resumable-upload):
- Upload-Metadata: comma-separated `key base64value` pairs. Keys are
  ASCII, values are base64 (may be empty).
- Upload-Offset / Upload-Length: decimal integers.
"""

from __future__ import annotations

import base64
import binascii


def parse_metadata(header_value: str | None) -> dict[str, str]:
    """Parse an Upload-Metadata header into {key: decoded_value}."""
    out: dict[str, str] = {}
    if not header_value:
        return out
    for pair in header_value.split(","):
        pair = pair.strip()
        if not pair:
            continue
        key, _, b64 = pair.partition(" ")
        key = key.strip()
        if not key:
            continue
        b64 = b64.strip()
        if not b64:
            out[key] = ""
            continue
        try:
            out[key] = base64.b64decode(b64, validate=True).decode("utf-8", "replace")
        except (binascii.Error, ValueError):
            raise ValueError(f"Upload-Metadata: value for {key!r} is not valid base64")
    return out


def encode_metadata(params: dict[str, str]) -> str:
    """Inverse of parse_metadata (used by tests/clients)."""
    return ",".join(
        f"{k} {base64.b64encode(v.encode()).decode()}" for k, v in params.items()
    )


def parse_int_header(value: str | None, name: str) -> int:
    if value is None:
        raise ValueError(f"missing required header {name}")
    try:
        n = int(value)
    except ValueError:
        raise ValueError(f"header {name} must be a decimal integer, got {value!r}")
    if n < 0:
        raise ValueError(f"header {name} must be >= 0, got {n}")
    return n
