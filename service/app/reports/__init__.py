"""Report rendering: report.json v0.2.0 -> HTML (Jinja2) -> PDF (WeasyPrint).

All strings come from service/app/reports/i18n.py. Non-English strings are
wrapped in <!-- TRANSLATOR-REVIEW --> at render time (visible in source) and
fr/ht reports carry a visible review banner. The template reads NOTHING
outside the report dict — every value is defensive (.get with fallbacks), so
a report missing optional keys still renders.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, select_autoescape

from .i18n import DISCLAIMER, REVIEW_BANNERS, RUO_LABEL, get_strings, mark_review

_TEMPLATES = Path(__file__).parent / "templates"
_env = Environment(
    loader=FileSystemLoader(str(_TEMPLATES)),
    autoescape=select_autoescape(default_for_string=True, default=False),
)

_MISSING = "—"


def _v(report: dict[str, Any], *path: str, default: Any = _MISSING) -> Any:
    cur: Any = report
    for p in path:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(p, default)
        if cur is default:
            return default
    return cur if cur not in (None, "") else default


def _num(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _fmt1(value: float) -> str:
    return f"{value:.1f}"


def _confidence(report: dict[str, Any], sm: dict[str, str]) -> tuple[str, list[str]]:
    """Operational confidence tier for one sample's report.

    Adapted from the Phase 5 FP/FN plan (P1 item 4) to the per-sample
    service report: High needs QC pass + >=30x mean depth + >=95% consensus
    called + <=1% non-V. cholerae reads. Anything below the high bar is
    Provisional; a QC-failed (or >5% contaminated) sample is Low / do not
    act. The tier is explicitly labeled operational, never a standard.
    """
    qc = report.get("qc", {})
    qc = qc if isinstance(qc, dict) else {}
    if str(qc.get("status", "")).lower() != "pass":
        return "low", [sm["reason_qc_fail"]]

    reasons: list[str] = []
    depth = _num((report.get("mapping") or {}).get("mean_depth"))
    if depth is not None and depth < 30:
        reasons.append(sm["reason_depth"].format(depth=_fmt1(depth)))
    called = _num(report.get("consensus_called_pct"))
    if called is not None and called < 95:
        reasons.append(sm["reason_called"].format(called=_fmt1(called)))
    frac = _num((report.get("classification") or {}).get("v_cholerae_fraction"))
    contam = (1.0 - frac) * 100.0 if frac is not None else None
    if contam is not None and contam > 1.0:
        reasons.append(sm["reason_contamination"].format(pct=_fmt1(contam)))
    if not reasons:
        return "high", []
    if contam is not None and contam > 5.0:
        return "low", reasons
    return "provisional", reasons


def build_context(
    report: dict[str, Any],
    job_id: str,
    org_id: str,
    lang: str = "en",
    uploaded_at: Any = None,
    completed_at: Any = None,
) -> dict[str, Any]:
    s, needs_review = get_strings(lang)
    # Every non-English UI string gets the review marker.
    sm = {k: mark_review(v, needs_review) for k, v in s.items()}

    qc = report.get("qc", {}) if isinstance(report.get("qc"), dict) else {}
    qc_status = str(qc.get("status", "")).lower()
    qc_pass = qc_status == "pass"
    reasons = qc.get("reasons") or []

    loci_raw = report.get("surveillance_loci") or {}
    call_map = {
        "present": sm["present"],
        "partial": sm["partial"],
        "absent": sm["absent"],
    }
    loci = [
        (locus, call_map.get(str(info.get("call", "")).lower(), str(info.get("call", _MISSING))))
        for locus, info in sorted(loci_raw.items())
        if isinstance(info, dict)
    ]

    assembly = report.get("assembly") or {}
    amr_genes = assembly.get("amr_genes") or []
    amr_rows = [(g, sm["detected"]) for g in amr_genes]

    tool_versions = report.get("tool_versions") or report.get("versions") or {}
    if isinstance(tool_versions, dict):
        tv = "; ".join(f"{k} {v}" for k, v in sorted(tool_versions.items()))
    else:
        tv = str(tool_versions)
    tv = tv or _MISSING

    frac = _v(report, "classification", "v_cholerae_fraction")
    try:
        frac = f"{float(frac) * 100:.1f}%" if frac != _MISSING else _MISSING
    except (TypeError, ValueError):
        pass

    tier, tier_reasons = _confidence(report, sm)

    return {
        "lang": lang if lang in ("en", "fr", "ht") else "en",
        "s": sm,
        "review_banner": (
            mark_review(REVIEW_BANNERS[lang], True) if lang in REVIEW_BANNERS else ""
        ),
        "disclaimer": mark_review(DISCLAIMER.get(lang, DISCLAIMER["en"]), needs_review),
        "ruo_label": mark_review(RUO_LABEL.get(lang, RUO_LABEL["en"]), needs_review),
        "job_id": job_id,
        "org_id": org_id,
        "uploaded_at": uploaded_at or _MISSING,
        "completed_at": completed_at or _MISSING,
        "pipeline_version": _v(report, "version"),
        "reference": _v(report, "reference"),
        "platform": _v(report, "platform"),
        "basecaller_model": _v(report, "basecaller_model"),
        "qc_pass": qc_pass,
        "qc_reasons": [str(r) for r in reasons],
        "confidence_tier": tier,
        "confidence_label": sm[f"confidence_{tier}"],
        "confidence_reasons": tier_reasons,
        "limitations": [sm[f"limitation_{i}"] for i in range(1, 6)],
        "v_cholerae_reads": _v(report, "classification", "v_cholerae_reads"),
        "v_cholerae_fraction": frac,
        "mean_depth": _v(report, "mapping", "mean_depth"),
        "breadth_pct": _v(report, "mapping", "breadth_pct"),
        "snps_vs_7pet": _v(report, "variants", "snps_vs_7pet"),
        "mlst_st": assembly.get("mlst_st", _MISSING),
        "vibecheck_lineage": assembly.get("vibecheck_lineage", _MISSING),
        "loci": loci,
        "amr_genes": amr_rows,
        "tool_versions": tv,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    }


def render_html(
    report: dict[str, Any],
    job_id: str,
    org_id: str,
    lang: str = "en",
    uploaded_at: Any = None,
    completed_at: Any = None,
) -> str:
    ctx = build_context(report, job_id, org_id, lang, uploaded_at, completed_at)
    return _env.get_template("report.html").render(ctx)


def render_pdf(
    report: dict[str, Any],
    job_id: str,
    org_id: str,
    lang: str = "en",
    uploaded_at: Any = None,
    completed_at: Any = None,
) -> bytes:
    """HTML -> PDF via WeasyPrint. Raises RuntimeError with a clear message
    if the system libraries are missing (documented in the Dockerfile)."""
    html = render_html(report, job_id, org_id, lang, uploaded_at, completed_at)
    try:
        from weasyprint import HTML
    except Exception as e:
        raise RuntimeError(
            "WeasyPrint is not usable in this environment "
            f"({e}); install the system deps from the Dockerfile service stage."
        )
    return HTML(string=html).write_pdf()
