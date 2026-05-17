"""Synchronize manuscript and site prose with pipeline JSON outputs."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple


PRIMARY_DETECTIONS: Tuple[Tuple[str, str], ...] = (
    ("NEAR", "NEAR"),
    ("Galileo_1990", "Galileo 1990"),
    ("Cassini", "Cassini"),
    ("Rosetta_2005", "Rosetta 2005"),
)


def _round(value: float | None, digits: int = 2) -> float | None:
    if value is None:
        return None
    return round(float(value), digits)


def _format_percent(value: float | None, digits: int = 1) -> str | None:
    if value is None:
        return None
    return f"{_round(value, digits):.{digits}f}%"


def _format_signed_mm(value: float, digits: int = 2) -> str:
    rounded = _round(value, digits)
    if rounded > 0:
        return f"+{rounded:.{digits}f}"
    return f"{rounded:.{digits}f}"


def _mantissa_exponent(value: float) -> Tuple[float, int]:
    if value == 0.0:
        return 0.0, 0
    v = float(value)
    if not math.isfinite(v):
        raise ValueError("mantissa/exponent formatting requires a finite, non-zero value")
    exp = int(math.floor(math.log10(abs(v))))
    mant = round(v / (10**exp), 1)
    return mant, exp


def _format_bayes_factor_approx(symbol: str, value: float, log10_value: float | None = None) -> str:
    """
    TeX for HTML/manuscript. Use only with a callable `re.sub` replacer (or `str.replace`):
    a string replacer would interpret ``\\approx`` / ``\\times`` as bell/tab escapes.
    """
    return _format_bayes_factor_html(symbol, value, log10_value)


def _format_bayes_factor_html(symbol: str, value: float, log10_value: float | None = None) -> str:
    """Single TeX backslash for HTML/manuscript (safe with `str.replace` and callable `re.sub`)."""
    v = float(value)
    if math.isinf(v) and v > 0 and log10_value is not None and math.isfinite(log10_value):
        exp_i = int(math.floor(log10_value))
        mant = round(10 ** (log10_value - exp_i), 2)
        return f"${symbol} \\approx {mant} \\times 10^{{{exp_i}}}$"
    if not math.isfinite(v):
        return f"${symbol} = \\mathrm{{NaN}}$"
    mant, exp = _mantissa_exponent(v)
    if exp == 0:
        return f"${symbol} \\approx {mant}$"
    return f"${symbol} \\approx {mant} \\times 10^{{{exp}}}$"


def _plain_scientific_times(value: float) -> str:
    """Plain-text scientific notation for README/Zenodo (ASCII × and ^)."""
    x = float(value)
    if not math.isfinite(x) or x == 0.0:
        return str(x)
    e = int(math.floor(math.log10(abs(x))))
    m = round(x / (10**e), 1)
    return f"{m}×10^{e}"


def _format_delta_bic_approx(value: float) -> str:
    return f"$\\Delta$BIC $\\approx {value}$"


_STEP026_METADATA_KEYS = frozenset(
    {
        "ensemble_catalog",
        "ensemble_composition_robustness",
        "full_catalog_model_comparison",
        "full_catalog_skipped",
        "sign_agreement_model_comparison",
        "strict_sign_gate",
        "likelihood_sensitivity_published_uncertainties_only",
    }
)


def _step026_primary_block(step026: Dict[str, Any]) -> Dict[str, Any]:
    """Primary gated comparison (top-level Step 026 row set; n=4 when sign gate is relaxed)."""
    return {key: value for key, value in step026.items() if key not in _STEP026_METADATA_KEYS}


def sign_agreement_model_block(step026: Dict[str, Any]) -> Dict[str, Any]:
    """Sign-agreement-restricted subset (typically n=3) for Section 4.5 robustness prose."""
    if "sign_agreement_model_comparison" in step026:
        return step026["sign_agreement_model_comparison"]
    primary = _step026_primary_block(step026)
    if int(primary.get("n_data", 0)) == 3:
        return primary
    raise KeyError(
        "step026_stable_model_comparison.json missing sign_agreement_model_comparison"
    )


def split_step026_model_blocks(step026: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Return (headline_full_catalog, primary_gated_ensemble) blocks from Step 026.

    ``full_catalog_model_comparison`` is the headline n=9 geometry-spread comparison.
    The primary gated row set lives at the top level (n=4 with ``strict_sign_gate``
    false; n=3 when true). The sign-agreement diagnostic is exposed via
    ``sign_agreement_model_block``.
    """
    if "full_catalog_model_comparison" not in step026:
        raise KeyError(
            "step026_stable_model_comparison.json missing full_catalog_model_comparison"
        )
    return step026["full_catalog_model_comparison"], _step026_primary_block(step026)


def _model_block_to_publication_dict(model_block: Dict[str, Any]) -> Dict[str, Any]:
    """Map one Step 026 comparison block to publication-sync / audit fields."""
    info = model_block["information_criteria"]
    bayes = model_block["bayes_factors"]
    selection = model_block.get("model_selection", {})
    akaike = selection.get("akaike_weights", {})
    return {
        "n_data": int(model_block.get("n_data", 0)),
        "log_likelihoods": {
            "Null": _round(model_block["log_likelihoods"]["Null"]),
            "Anderson": _round(model_block["log_likelihoods"]["Anderson"]),
            "TEP_restricted": _round(model_block["log_likelihoods"]["TEP_restricted"]),
            "TEP_flexible": _round(model_block["log_likelihoods"]["TEP_flexible"]),
        },
        "aic": {
            "Null": _round(info["Null"]["AIC"], 1),
            "Anderson": _round(info["Anderson"]["AIC"], 1),
            "TEP_restricted": _round(info["TEP_restricted"]["AIC"], 1),
            "TEP_flexible": _round(info["TEP_flexible"]["AIC"], 1),
        },
        "bic": {
            "Null": _round(info["Null"]["BIC"], 1),
            "Anderson": _round(info["Anderson"]["BIC"], 1),
            "TEP_restricted": _round(info["TEP_restricted"]["BIC"], 1),
            "TEP_flexible": _round(info["TEP_flexible"]["BIC"], 1),
        },
        "bayes": {
            "anderson_vs_null": _round(bayes["Anderson_vs_Null"], 1),
            "tep_vs_null": _round(bayes["TEP_restricted_vs_Null"], 1),
            "flex_vs_null": _round(bayes["TEP_flexible_vs_Null"], 1),
            "tep_vs_anderson": _round(bayes["TEP_restricted_vs_Anderson"], 1),
            "delta_bic_anderson_null": _round(bayes["delta_BIC_Anderson_vs_Null"], 1),
            "delta_bic_tep_null": _round(bayes["delta_BIC_TEP_restricted_vs_Null"], 1),
            "delta_bic_flex_null": _round(bayes["delta_BIC_TEP_flexible_vs_Null"], 1),
            "delta_bic_tep_anderson": _round(bayes["delta_BIC_TEP_restricted_vs_Anderson"], 1),
            "log10_tep_vs_null": round(float(bayes["log10_BF_TEP_restricted_vs_Null"]), 2),
            "log10_anderson_vs_null": round(float(bayes["log10_BF_Anderson_vs_Null"]), 2),
            "log10_flex_vs_null": round(float(bayes["log10_BF_TEP_flexible_vs_Null"]), 2),
            "log10_tep_vs_anderson": round(float(bayes["log10_BF_TEP_restricted_vs_Anderson"]), 2),
            "tep_vs_null_plain": _plain_scientific_times(float(bayes["TEP_restricted_vs_Null"])),
            "anderson_vs_null_plain": _plain_scientific_times(float(bayes["Anderson_vs_Null"])),
            "flex_vs_null_plain": _plain_scientific_times(float(bayes["TEP_flexible_vs_Null"])),
        },
        "bayes_approx": {
            "B10": _format_bayes_factor_approx(
                "B_{10}",
                bayes["TEP_restricted_vs_Null"],
                bayes.get("log10_BF_TEP_restricted_vs_Null"),
            ),
            "BA0": _format_bayes_factor_approx(
                "B_{A0}",
                bayes["Anderson_vs_Null"],
                bayes.get("log10_BF_Anderson_vs_Null"),
            ),
            "Bf0": _format_bayes_factor_approx(
                "B_{f0}",
                bayes["TEP_flexible_vs_Null"],
                bayes.get("log10_BF_TEP_flexible_vs_Null"),
            ),
            "B_vs_anderson": _format_bayes_factor_approx(
                "B",
                bayes["TEP_restricted_vs_Anderson"],
                bayes.get("log10_BF_TEP_restricted_vs_Anderson"),
            ),
            "dBIC_tep_null": _format_delta_bic_approx(_round(bayes["delta_BIC_TEP_restricted_vs_Null"], 1)),
            "dBIC_anderson_null": _format_delta_bic_approx(_round(bayes["delta_BIC_Anderson_vs_Null"], 1)),
            "dBIC_flex_null": _format_delta_bic_approx(_round(bayes["delta_BIC_TEP_flexible_vs_Null"], 1)),
            "dBIC_tep_anderson": _format_delta_bic_approx(
                _round(bayes["delta_BIC_TEP_restricted_vs_Anderson"], 1)
            ),
        },
        "akaike_weight_tep": _format_percent(100.0 * float(akaike.get("TEP_restricted", 0.0)), 1),
        "headline_bf_vs_anderson_html": _format_bayes_factor_html(
            "B",
            bayes["TEP_restricted_vs_Anderson"],
            bayes.get("log10_BF_TEP_restricted_vs_Anderson"),
        ),
        "bayes_factors_raw": bayes,
    }


def load_expected_publication_values(project_root: Path) -> Dict[str, Any]:
    step026 = json.loads(
        (project_root / "results/step026_stable_model_comparison.json").read_text(encoding="utf-8")
    )
    headline_block, gated_block = split_step026_model_blocks(step026)
    variance = json.loads((project_root / "results/step009_variance_analysis.json").read_text(encoding="utf-8"))
    predictions = json.loads((project_root / "results/step007_tep_predictions.json").read_text(encoding="utf-8"))

    model = _model_block_to_publication_dict(headline_block)
    model_gated = _model_block_to_publication_dict(gated_block)
    stages = variance["variance_decomposition"]["stages"]

    components: Dict[str, Dict[str, float]] = {}
    for mission_key, _ in PRIMARY_DETECTIONS:
        mission = predictions["predictions"][mission_key]["tep_predictions"]
        components[mission_key] = {
            "grad": mission["dv_grad_mm_s"],
            "disf": mission["dv_disf_mm_s"],
            "total": mission["dv_tep_mm_s"],
        }

    return {
        "model": model,
        "model_gated": model_gated,
        "variance": {
            "structural": _format_percent(stages["stage1_structural"]["explained_percent"]),
            "observational": _format_percent(stages["stage2_observational"]["explained_percent"]),
            "environmental": _format_percent(stages["stage3_environmental"]["explained_percent"]),
            "residual": _format_percent(stages["stage4_residual"]["explained_percent"]),
            "total_log10_variance_dex2": (
                None
                if variance["variance_decomposition"].get("total_variance_log10") is None
                else round(float(variance["variance_decomposition"]["total_variance_log10"]), 3)
            ),
            "f10_7_r_abs": round(
                abs(float(stages["stage3_environmental"].get("f10.7_correlation", 0.0))), 3
            ),
            "f10_7_p": (
                None
                if stages["stage3_environmental"].get("f10.7_p_value") is None
                else round(float(stages["stage3_environmental"]["f10.7_p_value"]), 2)
            ),
        },
        "components": components,
    }


def _publication_sources(project_root: Path) -> List[Path]:
    # The manuscript markdown is generated from `site/components/*.html`.
    # Do not write to the generated markdown directly.
    sources = sorted((project_root / "site/components").glob("*.html"))
    return [path for path in sources if path.exists()]


def load_step039_values(project_root: Path) -> Dict[str, Any]:
    payload = json.loads(
        (project_root / "results/step039_flyby_prediction_table.json").read_text(encoding="utf-8")
    )
    return payload


def load_step019_values(project_root: Path) -> Dict[str, Any]:
    path = project_root / "results/step019_3d_field_integration_results.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        "metadata" not in payload
        and "per_flyby" not in payload
        and isinstance(payload, dict)
        and payload
        and all(
            isinstance(v, dict) and "dv_3d_mm_s" in v for v in payload.values()
        )
    ):
        return {
            "metadata": {"trajectory_source": "legacy_flat_payload"},
            "per_flyby": payload,
        }
    return payload


def load_step040_values(project_root: Path) -> Dict[str, Any]:
    payload = json.loads(
        (project_root / "results/step040_cosmographic_shear.json").read_text(encoding="utf-8")
    )
    return payload


def load_step008_values(project_root: Path) -> Dict[str, Any]:
    payload = json.loads((project_root / "results/step008_fitting_results.json").read_text(encoding="utf-8"))
    return payload


def load_step025_values(project_root: Path) -> Dict[str, Any]:
    payload = json.loads(
        (project_root / "results/step025_corrected_uncertainty.json").read_text(encoding="utf-8")
    )
    return payload


def load_step032_values(project_root: Path) -> Dict[str, Any]:
    payload = json.loads(
        (project_root / "results/step032_tep_suppression_analysis.json").read_text(encoding="utf-8")
    )
    return payload


def load_step036_values(project_root: Path) -> Dict[str, Any]:
    payload = json.loads(
        (project_root / "results/step036_final_report.json").read_text(encoding="utf-8")
    )
    return payload


def _raw_class_label(raw_class: str) -> str:
    labels = {
        "true_positive": "True positive",
        "true_null": "True null",
        "raw_tension": "Raw tension",
        "raw_surplus": "Raw surplus",
        "insufficient_data": "Insufficient data",
    }
    return labels.get(raw_class, raw_class.replace("_", " ").title())


def _format_prediction_mm(row: Dict[str, Any]) -> str:
    prediction = row.get("raw_tep_prediction_mm_s")
    if prediction is None:
        return "N/A"
    uncertainty = row.get("raw_tep_uncertainty_mm_s")
    if uncertainty is None:
        return f"${_format_signed_mm(prediction, 3)}$"
    return f"${_format_signed_mm(prediction, 3)} \\pm {uncertainty:.3f}$"


def _format_residual_mm(row: Dict[str, Any]) -> str:
    residual = row.get("pooled_beta_residual_mm_s")
    if residual is None:
        return "—"
    return f"${_format_signed_mm(residual, 3)}$"


def _classification_label(row: Dict[str, Any]) -> str:
    raw_class = row.get("raw_classification")
    if raw_class:
        return _raw_class_label(raw_class)
    classification = row.get("classification", "")
    labels = {
        "insufficient_data": "Insufficient data",
        "insufficient_geometry_published_null": "Published null (geometry unavailable)",
        "true_null": "Predicted null",
        "f_od_unavailable": "F_OD unavailable",
    }
    return labels.get(classification, classification.replace("_", " ").title())


def _sync_step039_prediction_row(text: str, row: Dict[str, Any]) -> str:
    flyby = row["flyby"]
    pattern = rf"(<tr>\s*<td>{re.escape(flyby)}</td>)(\s*<td>[^<]+</td>)([\s\S]*?)(</tr>)"
    match = re.search(pattern, text)
    if not match:
        return text

    data_class_cell = match.group(2)
    cos_asym = row.get("cos_asymmetry")
    if cos_asym is None:
        cos_cell = "---"
    elif abs(cos_asym) < 0.01:
        cos_cell = "$\\approx 0$"
    else:
        cos_cell = f"${cos_asym:+.3f}$"

    observed = row.get("observed_dv_mm_s")
    sigma = row.get("observed_sigma_mm_s")
    if observed is None:
        observed_cell = "—"
    elif sigma is None:
        observed_cell = f"${_format_signed_mm(observed, 2)}$"
    else:
        observed_cell = f"${_format_signed_mm(observed, 2)} \\pm {sigma:.2f}$"

    altitude = row.get("altitude_km")
    altitude_cell = f"{altitude:.1f}" if altitude is not None else "—"

    replacement = (
        "<tr>\n"
        f"            <td>{flyby}</td>"
        f"{data_class_cell}\n"
        f"            <td>{altitude_cell}</td>\n"
        f"            <td>{cos_cell}</td>\n"
        f"            <td>{observed_cell}</td>\n"
        f"            <td>{_format_prediction_mm(row)}</td>\n"
        f"            <td>{_format_residual_mm(row)}</td>\n"
        f"            <td>—</td>\n"
        f"            <td>—</td>\n"
        f"            <td>{_classification_label(row)}</td>\n"
        "        </tr>"
    )
    return text[: match.start()] + replacement + text[match.end() :]


def _sync_step039_table(text: str, rows: List[Dict[str, Any]]) -> str:
    pattern = r"(<caption>\s*Table 4:[\s\S]*?<tbody>)([\s\S]*?)(</tbody>)"
    match = re.search(pattern, text)
    if not match:
        return text

    tbody = match.group(2)
    for row in rows:
        tbody = _sync_step039_prediction_row(tbody, row)
    return text[: match.start(2)] + tbody + text[match.end(2) :]


def _latex_beta() -> str:
    return r"$\beta$"


def _sync_step019_table3b_integration(text: str, step019: Dict[str, Any]) -> str:
    """
    Synchronize Table 3b integration-vs-perigee with Step 019 payload.

    Important: Step 019 is diagnostic on idealized toy geometry; do not present it
    as Horizons/SPICE validation. The table is kept numerically consistent with
    `step019_3d_field_integration_results.json`.
    """
    per_flyby = step019.get("per_flyby", {})
    if not isinstance(per_flyby, dict) or not per_flyby:
        return text

    pattern = r"(<caption>\s*Table 3b:[\s\S]*?<thead>[\s\S]*?</thead>\s*<tbody>)([\s\S]*?)(</tbody>)"
    match = re.search(pattern, text)
    if not match:
        return text

    # Prefer stable ordering: primary detections first if present, then remaining.
    preferred = ["NEAR", "Galileo_1990", "Rosetta_2005", "Cassini", "Juno", "Galileo_1992", "Rosetta_2007", "Rosetta_2009", "MESSENGER"]
    keys = [k for k in preferred if k in per_flyby] + [k for k in per_flyby.keys() if k not in preferred]

    def _label(key: str) -> str:
        labels = {
            "Galileo_1990": "Galileo 1990",
            "Galileo_1992": "Galileo 1992",
            "Rosetta_2005": "Rosetta 2005",
            "Rosetta_2007": "Rosetta 2007",
            "Rosetta_2009": "Rosetta 2009",
            "MESSENGER": "MESSENGER 2005",
        }
        return labels.get(key, key)

    rows_out: List[str] = []
    for key in keys:
        row = per_flyby.get(key, {})
        try:
            dv_peri = float(row.get("dv_current_mm_s"))
            dv_int = float(row.get("dv_3d_mm_s"))
        except (TypeError, ValueError):
            continue
        dv_simp = row.get("dv_simplified_mm_s")
        dv_simp_f = float(dv_simp) if dv_simp is not None else float("nan")
        ratio_legacy = dv_int / dv_peri if dv_peri != 0 else float("nan")
        ratio_comp = dv_int / dv_simp_f if dv_simp_f != 0 else float("nan")
        r_leg_cell = "—" if not math.isfinite(ratio_legacy) else f"{ratio_legacy:.2f}"
        r_comp_cell = "—" if not math.isfinite(ratio_comp) else f"{ratio_comp:.2f}"
        simp_cell = f"{dv_simp_f:+.3f}" if math.isfinite(dv_simp_f) else "—"
        rows_out.append(
            "<tr>\n"
            f"            <td>{_label(key)}</td>\n"
            f"            <td>{dv_peri:+.3f}</td>\n"
            f"            <td>{simp_cell}</td>\n"
            f"            <td>{dv_int:+.3f}</td>\n"
            f"            <td>{r_comp_cell}</td>\n"
            f"            <td>{r_leg_cell}</td>\n"
            "            <td>diagnostic</td>\n"
            "            <td>diagnostic</td>\n"
            "        </tr>"
        )

    if not rows_out:
        return text

    new_tbody = "\n".join(rows_out)
    return text[: match.start(2)] + new_tbody + text[match.end(2) :]


def _sync_step040_table8(text: str, step040: Dict[str, Any]) -> str:
    results = step040.get("flyby_results", [])
    if not isinstance(results, list) or not results:
        return text

    # Must match cosmographic caption only — discussion §5.2 also contains "Table 8"
    # which must never be overwritten with Step 040 rows.
    pattern = (
        r"(<caption>\s*Table 8:\s*Cosmographic Modulation Parameters and Residual Ratios"
        r"[\s\S]*?<tbody>)([\s\S]*?)(</tbody>)"
    )
    match = re.search(pattern, text)
    if not match:
        return text

    def _fmt_ratio(val: float) -> str:
        if not math.isfinite(val):
            return "—"
        return f"<strong>{val:.2f}</strong>"

    rows_out: List[str] = []
    for row in results:
        if not row.get("has_3d_vector", False):
            continue
        mission = str(row.get("mission", "")).strip()
        if not mission:
            continue
        r_au = row.get("heliocentric_distance_au")
        v_rad = row.get("radial_velocity_kms")
        cos_theta = row.get("sc_cmb_cos_theta")
        v_sc_proj = row.get("sc_velocity_cmb_proj_kms")
        f_enh = row.get("cmb_disformal_enhancement")
        both = bool(row.get("both_aligned_flag", 0.0) > 0.5)
        obs = row.get("observed_mm_s")
        pred = row.get("predicted_mm_s")
        ratio = row.get("ratio_obs_pred")

        def _num(x, fmt):
            if x is None:
                return "—"
            try:
                return fmt(float(x))
            except (TypeError, ValueError):
                return "—"

        mission_label = mission.replace("_", " ")
        both_cell = "<strong>YES</strong>" if both else "no"
        rows_out.append(
            "<tr>\n"
            f"            <td>{mission_label}</td>\n"
            f"            <td>{_num(r_au, lambda v: f'{v:.3f}')}</td>\n"
            f"            <td>{_num(v_rad, lambda v: f'{v:+.3f}')}</td>\n"
            f"            <td>{_num(cos_theta, lambda v: f'{v:+.3f}')}</td>\n"
            f"            <td>{_num(v_sc_proj, lambda v: f'{v:+.2f}')}</td>\n"
            f"            <td>{_num(f_enh, lambda v: f'{v:.0f}')}</td>\n"
            f"            <td>{both_cell}</td>\n"
            f"            <td>{_num(obs, lambda v: f'{v:.2f}')}</td>\n"
            f"            <td>{_num(pred, lambda v: f'{v:.3f}')}</td>\n"
            f"            <td>{_fmt_ratio(float(ratio)) if ratio is not None else '—'}</td>\n"
            "        </tr>"
        )

    if not rows_out:
        return text
    new_tbody = "\n".join(rows_out)
    return text[: match.start(2)] + new_tbody + text[match.end(2) :]


def _sync_step040_table9b(text: str, step040: Dict[str, Any]) -> str:
    results = step040.get("flyby_results", [])
    reg = step040.get("multivariate_regression", {})
    coeffs = reg.get("coefficients", {}) if isinstance(reg, dict) else {}
    if not isinstance(results, list) or not coeffs:
        return text

    b0 = float(coeffs.get("b0_intercept"))
    b1 = float(coeffs.get("b1_sc_cmb_cos"))
    b2 = float(coeffs.get("b2_earth_cmb_proj"))
    b3 = float(coeffs.get("b3_sc_orbital"))

    pattern = r"(<caption>\s*Table 9b:[\s\S]*?<tbody>)([\s\S]*?)(</tbody>)"
    match = re.search(pattern, text)
    if not match:
        return text

    rows_out: List[str] = []
    for row in results:
        if not row.get("has_3d_vector", False):
            continue
        mission = str(row.get("mission", "")).strip()
        if not mission:
            continue
        obs_ratio = row.get("ratio_obs_pred")
        sc_cos = row.get("sc_cmb_cos_theta")
        earth_proj = row.get("earth_orbital_cmb_proj_kms")
        sc_orb = row.get("sc_orbital_alignment")
        if any(v is None for v in (obs_ratio, sc_cos, earth_proj, sc_orb)):
            continue
        obs_ratio_f = float(obs_ratio)
        pred_ratio = b0 + b1 * float(sc_cos) + b2 * (float(earth_proj) / 30.0) + b3 * float(sc_orb)
        resid = obs_ratio_f - pred_ratio
        mission_label = mission.replace("_", " ")
        rows_out.append(
            "<tr>\n"
            f"            <td>{mission_label}</td>\n"
            f"            <td>{obs_ratio_f:.2f}</td>\n"
            f"            <td>{pred_ratio:.2f}</td>\n"
            f"            <td>{resid:+.2f}</td>\n"
            "        </tr>"
        )

    if not rows_out:
        return text
    # Ensure tbody wrapper is preserved.
    new_tbody = "\n".join(rows_out)
    return text[: match.start(2)] + "\n" + new_tbody + "\n" + text[match.end(2) :]


def _sync_step040_regression_prose(text: str, step040: Dict[str, Any]) -> str:
    """
    Synchronize the multivariate-regression prose block in 4.11.6 with Step 040 JSON.
    """
    reg = step040.get("multivariate_regression", {})
    if not isinstance(reg, dict):
        return text
    coeffs = reg.get("coefficients", {})
    if not isinstance(coeffs, dict):
        return text

    required = ("b0_intercept", "b1_sc_cmb_cos", "b2_earth_cmb_proj", "b3_sc_orbital")
    if not all(k in coeffs for k in required):
        return text

    b0 = float(coeffs["b0_intercept"])
    b1 = float(coeffs["b1_sc_cmb_cos"])
    b2 = float(coeffs["b2_earth_cmb_proj"])
    b3 = float(coeffs["b3_sc_orbital"])
    r2 = float(reg.get("r_squared", float("nan")))
    adj = float(reg.get("adjusted_r_squared", float("nan")))
    std0 = float(reg.get("original_std", float("nan")))
    std1 = float(reg.get("residual_std", float("nan")))
    red = float(reg.get("std_reduction_percent", float("nan")))

    pattern = (
        r"The fitted coefficients are <em>b</em><sub>0</sub> = [^,]+,"
        r"[\s\S]*?<em>b</em><sub>3</sub> = [^\.]+\. "
        r"The model achieves[\s\S]*?significance\."
    )
    replacement = (
        "The fitted coefficients are <em>b</em><sub>0</sub> = "
        f"{b0:+.3f}, <em>b</em><sub>1</sub> = {b1:+.3f}, "
        f"<em>b</em><sub>2</sub> = {b2:+.3f}, <em>b</em><sub>3</sub> = {b3:+.3f}. "
        f"The model achieves <em>R</em><sup>2</sup> = {r2:.3f} and reduces the residual standard "
        f"deviation from {std0:.3f} to {std1:.3f}, a {red:.1f}% reduction. The adjusted "
        f"<em>R</em><sup>2</sup> = {adj:+.3f} indicates that with <em>n</em> = 8 "
        "and four parameters (including intercept), the regression does not explain residual variance "
        "at conventional significance."
    )
    return re.sub(pattern, replacement, text, count=1)


def _sync_step039_prose(text: str, step039: Dict[str, Any]) -> str:
    summary = step039["metadata"]["raw_classification_summary"]
    uncertainty_summary = step039["metadata"].get(
        "raw_uncertainty_aware_classification_summary", {}
    )
    true_null = summary["true_null"]
    raw_tension = summary["raw_tension"]
    true_positive = summary["true_positive"]
    uncertainty_raw_tension = uncertainty_summary.get("raw_tension", raw_tension)
    uncertainty_true_null = uncertainty_summary.get("true_null", true_null)
    juno = next(row for row in step039["rows"] if row["flyby"] == "Juno")
    beta_tex = _latex_beta()

    def _replace_connection_paragraph(_match: re.Match[str]) -> str:
        return (
            "<p>\n"
            "    <strong>Connection to observations:</strong> Step 039 classifies\n"
            f"    fixed-amplitude raw predictions (Step 008 pooled {beta_tex}) with the Step 007 geometry envelope\n"
            f"    ({true_positive} deterministic true positives, {true_null} deterministic true nulls, {raw_tension} deterministic fixed-amplitude warning case). "
            f"After propagating Step 008 random-effects amplitude scatter, the uncertainty-aware raw layer has {uncertainty_raw_tension} raw-tension cases and {uncertainty_true_null} null-compatible cases. Post-OD survival\n"
            "    factors are withheld until mission OD configuration yields defensible\n"
            "    $F_{\\rm OD}$ estimates.\n"
            "</p>"
        )

    def _replace_summary_paragraph(_match: re.Match[str]) -> str:
        return (
            "<p>\n"
            "    <strong>Summary:</strong> At the Step 008 inverse-variance pooled $\\beta$, Step 039\n"
            f"    classifies {true_positive} published anomalies as deterministic fixed-amplitude raw true positives (NEAR, Galileo\n"
            "    1990, Rosetta 2005). "
            f"{true_null} published null or bound cases are consistent with fixed-amplitude Step 039 predictions "
            f"(Step 008 pooled {beta_tex}) under the Step 007 geometry envelope; {raw_tension} deterministic fixed-amplitude warning case remains (Juno). "
            f"When random-effects amplitude scatter is propagated into the prediction uncertainty, {uncertainty_raw_tension} uncertainty-aware raw-tension cases remain. "
            "Cassini and Galileo 1992 are classified as raw true nulls at the refit weighted-mean $\\beta$. Rosetta 2009 is a published null/bound case with insufficient explicit geometry for the Step 039 prediction table, while Stardust, OSIRIS-REx,\n"
            "    and BepiColombo lack public anomaly reports and are not\n"
            "    used in quantitative likelihood.\n"
            "</p>"
        )

    def _replace_key_result_li(_match: re.Match[str]) -> str:
        return (
            "\n        <li>Key result: "
            f"{true_null} published null or bound cases are consistent with "
            f"fixed-amplitude Step 039 predictions (Step 008 pooled {beta_tex}) under the Step 007 geometry envelope; "
            f"{raw_tension} deterministic fixed-amplitude warning case remains (Juno), but {uncertainty_raw_tension} raw-tension cases remain after random-effects prediction uncertainty; Rosetta 2009 is a published null/bound case with insufficient "
            "explicit geometry for the Step 039 table; 3 flybys (Stardust, OSIRIS-REx, BepiColombo) have no public anomaly report</li>"
        )

    text = re.sub(
        r"\s*<li>Key result: \d+ published null or bound cases are consistent with (?:fixed-amplitude Step 039 predictions \(Step 008 pooled \$\\beta\$\)|fixed-amplitude predictions) under the Step 007 geometry envelope; \d+ (?:deterministic fixed-amplitude warning|fixed-amplitude raw-tension|raw-tension) case remains \(Juno\)(?:, but \d+ raw-tension cases? remain after random-effects prediction uncertainty|, but \d+ remain after random-effects prediction uncertainty)?; Rosetta 2009 is a published null/bound case with insufficient explicit geometry for the Step 039 table; \d+ flybys \(Stardust, OSIRIS-REx, BepiColombo\) have no public anomaly report</li>",
        _replace_key_result_li,
        text,
    )

    text = re.sub(
        r"<p>\s*<strong>Connection to observations:</strong> Step 039 classifies[\s\S]*?\$F_\{\\rm OD\}\$ estimates\.\s*</p>",
        _replace_connection_paragraph,
        text,
    )

    text = re.sub(
        r"<p>\s*<strong>Summary:</strong> At the Step 008 inverse-variance pooled \$\\beta\$, Step 039[\s\S]*?used in quantitative likelihood\.\s*</p>",
        _replace_summary_paragraph,
        text,
    )

    juno_prediction = _format_signed_mm(juno["raw_tep_prediction_mm_s"], 2)
    juno_prediction_unc = juno.get("raw_tep_uncertainty_mm_s")
    juno_observed = _format_signed_mm(juno["observed_dv_mm_s"], 2)
    juno_sigma = juno["observed_sigma_mm_s"]
    def _replace_falsifiability_paragraph(_match: re.Match[str]) -> str:
        unc_clause = (
            f" with random-effects prediction uncertainty $\\pm {juno_prediction_unc:.2f}$ mm/s"
            if juno_prediction_unc is not None
            else ""
        )
        return (
            "Juno remains the sole deterministic fixed-amplitude warning case at the refit weighted-mean $\\beta$ "
            f"(${juno_prediction}$ mm/s raw prediction vs. "
            f"${juno_observed} " + r"\pm" + f" {juno_sigma:.2f}$ mm/s observed){unc_clause}. "
            f"After propagating Step 008 random-effects scatter, the uncertainty-aware Step 039 layer has {uncertainty_raw_tension} raw-tension cases. "
            "Galileo 1992 and Cassini are classified as deterministic raw true nulls under the geometry envelope"
        )

    text = re.sub(
        r"Juno remains the sole (?:deterministic fixed-amplitude )?(?:raw-tension|warning) case at the refit weighted-mean \$\\beta\$ \(\$\+0\.\d+\$ mm/s raw prediction vs\. \$0\.00 \\pm 0\.02\$ mm/s observed\)(?: with random-effects prediction uncertainty \$\\pm [\d.]+\$ mm/s)?\. (?:After propagating Step 008 random-effects scatter, the uncertainty-aware Step 039 layer has \d+ raw-tension cases\. )?Galileo 1992 and Cassini are classified as (?:deterministic )?raw true nulls under the geometry envelope",
        _replace_falsifiability_paragraph,
        text,
        flags=re.MULTILINE,
    )
    return text


def _sync_plaintext_bic_sentence(text: str, model_dict: Dict[str, Any]) -> str:
    """
    Rewrite the Unicode-only Bayes/BIC clause in HTML abstract (1_abstract.html)
    to match Step 026 JSON (avoids 10^18 vs 10^71 drift when body uses LaTeX).
    """
    bayes_approx = model_dict.get("bayes_approx", {})
    if not bayes_approx:
        return text
    legacy = (
        "On the same n = 3 gated ensemble as Step 026, BIC-based summaries give B₁₀ ≈ 1.8×10¹⁸ "
        "(TEP restricted vs Null, ΔBIC ≈ 84), B_A0 ≈ 9.9×10¹⁶ (Anderson vs Null, ΔBIC ≈ 78), "
        "and TEP restricted vs Anderson B ≈ 18.2 (ΔBIC ≈ 5.8), indicating positive"
    )
    if legacy not in text:
        return text
    b10 = bayes_approx["B10"]
    ba0 = bayes_approx["BA0"]
    bva = bayes_approx["B_vs_anderson"]
    d0 = bayes_approx["dBIC_tep_null"]
    da = bayes_approx["dBIC_anderson_null"]
    dt = bayes_approx["dBIC_tep_anderson"]
    new_sentence = (
        "On the same n = 3 gated ensemble as Step 026, BIC-derived Bayes summaries give "
        f"{b10} (TEP restricted vs Null, {d0}), {ba0} (Anderson vs Null, {da}), "
        f"and TEP restricted vs Anderson {bva} ({dt}), indicating positive"
    )
    return text.replace(legacy, new_sentence)


def _sync_text(text: str, expected: Dict[str, Any], step039: Dict[str, Any]) -> str:
    """Rewrite static manuscript numerics using `expected` / Step 039 JSON.

    TeX-heavy `re.sub` replacements must use a callable replacer (or doubled escapes
    where historically required); plain replacement strings interpret ``\\a`` / ``\\t``.
    Exception: templates that need ``\\g<n>`` back-references stay as strings.
    """
    text = text.replace("\x07pprox", r"\approx")
    text = re.sub(r"\t+imes", r"\\times", text)
    text = re.sub(r"\t+approx", r"\\approx", text)

    model = expected["model"]
    variance = expected["variance"]
    bayes_approx = model.get("bayes_approx", {})

    replacements = [
        (r"log L = -?\d+\.\d+(?=, AIC = [\d.]+, BIC = 35\.7)", f"log L = {model['log_likelihoods']['Null']}"),
        (r"log L = -?\d+\.\d+(?=, AIC = [\d.]+, BIC = 21\.9)", f"log L = {model['log_likelihoods']['Anderson']}"),
        (r"log L = -?\d+\.\d+(?=, AIC = [\d.]+, BIC = 18\.2)", f"log L = {model['log_likelihoods']['TEP_restricted']}"),
        (r"log L = -?\d+\.\d+(?=, AIC = [\d.]+, BIC = 21\.0)", f"log L = {model['log_likelihoods']['TEP_flexible']}"),
        (r"\$B_\{10\}\s*(?:=|\\approx)\s*[\de.+-]+(?:\s*\\times\s*10\^\{[\de+-]+\})?\$", bayes_approx.get("B10", f"$B_{{10}} = {model['bayes']['tep_vs_null']}")),
        (r"\$B_\{A0\}\s*(?:=|\\approx)\s*[\de.+-]+(?:\s*\\times\s*10\^\{[\de+-]+\})?\$", bayes_approx.get("BA0", f"$B_{{A0}} = {model['bayes']['anderson_vs_null']}")),
        (r"\$B_\{f0\}\s*(?:=|\\approx)\s*[\de.+-]+(?:\s*\\times\s*10\^\{[\de+-]+\})?\$", bayes_approx.get("Bf0", f"$B_{{f0}} = {model['bayes']['flex_vs_null']}")),
        (r"against Anderson gives \$B\s*(?:=|\\approx)\s*[\de.+-]+(?:\s*\\times\s*10\^\{[\de+-]+\})?\$", f"against Anderson gives {bayes_approx.get('B_vs_anderson', f'$B = {model['bayes']['tep_vs_anderson']}')}"),
        (r"Akaike weight for TEP restricted is [\d.]+%", f"Akaike weight for TEP restricted is {model['akaike_weight_tep']}"),
        (r"structural proxy bundle accounts for [\d.]+%", f"structural proxy bundle accounts for {variance['structural']}"),
        (r"observational pipeline effects \(OD filter absorption and systematic uncertainties\) account for [\d.]+%", f"observational pipeline effects (OD filter absorption and systematic uncertainties) account for {variance['observational']}"),
        (r"environmental modulation contributes [\d.]+%", f"environmental modulation contributes {variance['environmental']}"),
        (r"residual \(small-sample statistics, intrinsic scatter, model incompleteness\) accounts for [\d.]+%", f"residual (small-sample statistics, intrinsic scatter, model incompleteness) accounts for {variance['residual']}"),
        (r"tracked structural proxy bundle accounts for [\d.]+% of fitted-β variance", f"tracked structural proxy bundle accounts for {variance['structural']} of fitted-β variance"),
        (r"observational effects account for [\d.]+%", f"observational effects account for {variance['observational']}"),
        (r"the residual accounts for [\d.]+%", f"the residual accounts for {variance['residual']}"),
        # Alternative Step 009 sentence structures found in 5_results.html and 7_conclusion.html
        (r"environmental modulation is assigned [\d.]+%", f"environmental modulation is assigned {variance['environmental']}"),
        (r"and the residual is [\d.]+%", f"and the residual is {variance['residual']}"),
        (r"structural and observational proxies are negligible fractions", f"structural and observational proxies are negligible fractions"),
    ]

    for pattern, replacement in replacements:
        text = re.sub(pattern, lambda _m, rep=replacement: rep, text)

    # Match the full multi-value sentence from 5_results.html line ~1170
    # (backreferences require a callable replacer, not a static string)
    text = re.sub(
        r"(The Step 009 variance decomposition on log₁₀ β assigns) [\d.]+% (to structural proxies), [\d.]+% (to observational pipeline effects), [\d.]+% (to environmental modulation \(F10\.7; sample-limited at n = 3\)), and [\d.]+% (to residual terms)",
        lambda m: (
            f"{m.group(1)} {variance['structural']} {m.group(2)}, "
            f"{variance['observational']} {m.group(3)}, "
            f"{variance['environmental']} {m.group(4)}, and "
            f"{variance['residual']} {m.group(5)}"
        ),
        text,
    )
    # Match a shorter variant without "F10.7; sample-limited at n = 3"
    text = re.sub(
        r"(assigns [\d.]+% to structural proxies), [\d.]+% (to observational pipeline effects), [\d.]+% (to environmental modulation), and [\d.]+% (to residual terms)",
        lambda m: (
            f"{m.group(1)}, {variance['observational']} {m.group(2)}, "
            f"{variance['environmental']} {m.group(3)}, and "
            f"{variance['residual']} {m.group(4)}"
        ),
        text,
    )

    tv_dex = variance.get("total_log10_variance_dex2")
    if tv_dex is not None:
        text = re.sub(
            r"log\u2081\u2080\(\u03b2\) is [\d.]+ dex\u00b2",
            lambda _m: f"log\u2081\u2080(\u03b2) is {tv_dex} dex\u00b2",
            text,
        )

    f10_r = variance.get("f10_7_r_abs")
    f10_p = variance.get("f10_7_p")
    if f10_r is not None and f10_p is not None:
        text = re.sub(
            r"The Step 009 stage-3 environmental diagnostic for solar_activity_f10\.7 lists Pearson \|r\| ≈ [\d.]+ "
            r"\(two-sided p ≈ [\d.]+, n = 3 gated fits\)\.",
            lambda _m: (
                "The Step 009 stage-3 environmental diagnostic for solar_activity_f10.7 lists Pearson "
                f"|r| ≈ {f10_r} (two-sided p ≈ {f10_p}, n = 3 gated fits)."
            ),
            text,
        )
        text = re.sub(
            r"\(F10\.7 correlation \|r\| ≈ [\d.]+,\s*p ≈ [\d.]+ on n = 3\)",
            lambda _m: (
                f"(F10.7 correlation |r| ≈ {f10_r}, p ≈ {f10_p} on n = 3)"
            ),
            text,
        )

    gbf_html = model.get("headline_bf_vs_anderson_html")
    if gbf_html and bayes_approx:
        text = text.replace(r"$B \approx 15536.4$", gbf_html)
    d_tep_null = model["bayes"]["delta_bic_tep_null"]
    d_tep_anderson = model["bayes"]["delta_bic_tep_anderson"]
    d_anderson_null = model["bayes"]["delta_bic_anderson_null"]
    for old, new in (
        (r"$\Delta{\rm BIC}\approx84$", rf"$\Delta{{\rm BIC}}\approx{d_tep_null}$"),
        (r"$\Delta{\rm BIC}\approx 84$", rf"$\Delta{{\rm BIC}}\approx {d_tep_null}$"),
        (r"$\Delta{\rm BIC}\approx5.8$", rf"$\Delta{{\rm BIC}}\approx{d_tep_anderson}$"),
        (r"$\Delta{\rm BIC}\approx 5.8$", rf"$\Delta{{\rm BIC}}\approx {d_tep_anderson}$"),
        (r"$\Delta{\rm BIC}\approx78$", rf"$\Delta{{\rm BIC}}\approx{d_anderson_null}$"),
        (r"$\Delta{\rm BIC}\approx 78$", rf"$\Delta{{\rm BIC}}\approx {d_anderson_null}$"),
    ):
        text = text.replace(old, new)
    gbf_html_dup = model.get("headline_bf_vs_anderson_html")
    if gbf_html_dup:
        text = text.replace(
            r"$B \approx 2.94 \times 10^{911}$ \times 10^{911}$$",
            gbf_html_dup,
        )
        text = text.replace(
            r"$B \approx 2.94 \times 10^{911}$ \times 10^{911}$",
            gbf_html_dup,
        )

    for label, key in (
        ("Null", "Null"),
        ("Anderson", "Anderson"),
        ("TEP restricted", "TEP_restricted"),
        ("TEP flexible", "TEP_flexible"),
    ):
        text = re.sub(
            rf"<strong>{re.escape(label)}:</strong> log L = [-+]?\d+\.\d+, AIC = [-+]?\d+\.\d+, BIC = [-+]?\d+\.\d+",
            lambda _m, lbl=label, k=key: (
                f"<strong>{lbl}:</strong> log L = {model['log_likelihoods'][k]}, "
                f"AIC = {model['aic'][k]}, BIC = {model['bic'][k]}"
            ),
            text,
        )

    if not bayes_approx:
        text = re.sub(
            r"\$B_\{10\} = [\de.+-]+",
            lambda _m: f"$B_{{10}} = {model['bayes']['tep_vs_null']}",
            text,
        )
        text = re.sub(
            r"\$B_\{A0\} = [\de.+-]+",
            lambda _m: f"$B_{{A0}} = {model['bayes']['anderson_vs_null']}",
            text,
        )
        text = re.sub(
            r"\$B_\{f0\} = [\d.]+",
            lambda _m: f"$B_{{f0}} = {model['bayes']['flex_vs_null']}",
            text,
        )

    if bayes_approx:
        # Repair legacy corruption where repeated "\times 10^{...}$" fragments were appended
        # outside the closing "$" of the Bayes-factor math segment.
        for symbol_key, replacement in (
            ("B_{10}", bayes_approx["B10"]),
            ("B_{A0}", bayes_approx["BA0"]),
            ("B_{f0}", bayes_approx["Bf0"]),
        ):
            pattern = rf"(\${re.escape(symbol_key)}[\s\S]*?\$)(?:\s*\\times\s*10\^\{{-?\d+\}}\$)+"
            text = re.sub(pattern, lambda _m, rep=replacement: rep, text)
        # Legacy corruption: finite formatter emitted "inf" when JSON BF was +Infinity
        text = re.sub(
            r"\$B\s*\\approx\s*inf\$\s*\\times\s*10\^\{\d+\}\$",
            lambda _m, s=bayes_approx["B_vs_anderson"]: s,
            text,
        )
    if bayes_approx:
        text = re.sub(
            r"B₁₀ = [\d.]+ \(ΔBIC = [\d.]+, strong evidence\)",
            lambda _m: (
                f"B₁₀ ≈ {model['bayes']['tep_vs_null_plain']} "
                f"(ΔBIC ≈ {model['bayes']['delta_bic_tep_null']}, headline BIC surrogate; not literal posterior odds at n=3)"
            ),
            text,
        )
    else:
        text = re.sub(
            r"B₁₀ = [\d.]+ \(ΔBIC = [\d.]+, strong evidence\)",
            lambda _m: (
                f"B₁₀ = {model['bayes']['tep_vs_null']} (ΔBIC = {model['bayes']['delta_bic_tep_null']}, strong evidence)"
            ),
            text,
        )
    text = re.sub(
        r"TEP restricted vs Anderson empirical B = [\d.]+ \(ΔBIC = [\d.]+, positive evidence\)",
        lambda _m: (
            f"TEP restricted vs Anderson empirical B = {model['bayes']['tep_vs_anderson']} "
            f"(ΔBIC = {model['bayes']['delta_bic_tep_anderson']}, positive evidence)"
        ),
        text,
    )

    def _sync_delta_bic(text_in: str, label: str, value: float) -> str:
        # Match ($\Delta$BIC $\approx 123.4$) after the model row label (not legacy $\Delta\$ split).
        pattern = rf"(<strong>{re.escape(label)}:</strong>[\s\S]*?\(\$\\Delta\$BIC\s*\$\\approx\s*)[\d.]+(?=\$\))"
        return re.sub(pattern, rf"\g<1>{value}", text_in, count=1)

    def _sync_delta_bic_paren(text_in: str, label: str, value: float) -> str:
        pattern = rf"(<strong>{re.escape(label)}:</strong>[\s\S]*?\(\$\\Delta\\$BIC\s*\$\\approx\s*)[\d.]+(?=\$\))"
        return re.sub(pattern, rf"\g<1>{value}", text_in, count=1)

    for label, value in (
        ("TEP restricted vs Null", model["bayes"]["delta_bic_tep_null"]),
        ("Anderson vs Null", model["bayes"]["delta_bic_anderson_null"]),
        ("TEP flexible vs Null", model["bayes"]["delta_bic_flex_null"]),
        ("TEP restricted vs Anderson", model["bayes"]["delta_bic_tep_anderson"]),
    ):
        text = _sync_delta_bic(text, label, value)
        text = _sync_delta_bic_paren(text, label, value)
    def _delta_bic_replacement(value: float) -> str:
        return _format_delta_bic_approx(value)

    def _replace_anderson_vs_null(_match: re.Match[str]) -> str:
        return (
            "<strong>Anderson vs Null:</strong> "
            f"$B_{{A0}} = {model['bayes']['anderson_vs_null']} "
            + _delta_bic_replacement(model["bayes"]["delta_bic_anderson_null"])
        )

    def _replace_tep_flexible_vs_null(_match: re.Match[str]) -> str:
        return (
            "<strong>TEP flexible vs Null:</strong> "
            f"$B_{{f0}} = {model['bayes']['flex_vs_null']} "
            + _delta_bic_replacement(model["bayes"]["delta_bic_flex_null"])
        )

    text = re.sub(
        r"\(\$\\Delta\$BIC = [\d.]+\)(?= — <em>strong evidence</em>)",
        lambda _match: _delta_bic_replacement(model["bayes"]["delta_bic_tep_null"]),
        text,
        count=1,
    )
    text = re.sub(
        r"\(\$\\Delta\$BIC = [\d.]+\)(?= — <em>positive evidence</em>)",
        lambda _match: _delta_bic_replacement(model["bayes"]["delta_bic_tep_anderson"]),
        text,
        count=1,
    )
    text = re.sub(
        r"<strong>Anderson vs Null:</strong> \$B_\{A0\} = [\de.+-]+ \(\$\\Delta\$BIC = [\d.]+\)",
        _replace_anderson_vs_null,
        text,
    )
    text = re.sub(
        r"<strong>TEP flexible vs Null:</strong> \$B_\{f0\} = [\d.]+ \(\$\\Delta\$BIC = [\d.]+\)",
        _replace_tep_flexible_vs_null,
        text,
    )

    for pattern, value in (
        (r"\$\\Delta\$BIC = [\d.]+(?= — \*strong evidence\*)", model["bayes"]["delta_bic_tep_null"]),
        (r"\$\\Delta\$BIC = [\d.]+(?= — \*positive evidence\*)", model["bayes"]["delta_bic_tep_anderson"]),
    ):
        text = re.sub(pattern, lambda _match, replacement=value: f"$\\Delta$BIC = {replacement}", text)

    for mission_key, label in PRIMARY_DETECTIONS:
        component = expected["components"][mission_key]
        # Allow optional $...$ wrappers around numeric cells (LaTeX-style tables).
        markdown_row = (
            rf"\| {re.escape(label)} \| (\$?)([+-]?\d+\.\d+)(\$?) \| (\$?)([+-]?\d+\.\d+)(\$?) \| (\$?)([+-]?\d+\.\d+)(\$?) \|"
        )
        markdown_replacement = (
            f"| {label} | {_format_signed_mm(component['grad'])} | "
            f"{_format_signed_mm(component['disf'])} | {_format_signed_mm(component['total'])} |"
        )
        text = re.sub(markdown_row, lambda _m, rep=markdown_replacement: rep, text)

        html_row = (
            rf"(<td>{re.escape(label)}</td>\s*<td>)(?:\$)?([+-]?\d+\.\d+)(?:\$)?(</td>\s*<td>)(?:\$)?"
            rf"([+-]?\d+\.\d+)(?:\$)?(</td>\s*<td>)(?:\$)?([+-]?\d+\.\d+)(?:\$)?(</td>)"
        )
        html_replacement = (
            rf"\g<1>{_format_signed_mm(component['grad'])}\g<3>"
            rf"{_format_signed_mm(component['disf'])}\g<5>"
            rf"{_format_signed_mm(component['total'])}\g<7>"
        )
        text = re.sub(html_row, html_replacement, text, flags=re.IGNORECASE)

    text = _sync_step039_table(text, step039.get("rows", []))
    text = _sync_step039_prose(text, step039)
    text = _sync_plaintext_bic_sentence(text, model)
    return text


def _sync_step025_uncertainty(text: str, step025: Dict[str, Any]) -> str:
    """Synchronize uncertainty budget claims with Step 025 JSON."""
    budget = step025.get("corrected_uncertainty_budget", {})
    total_rel = budget.get("total_relative_uncertainty")
    sys_rel = budget.get("systematic_relative_uncertainty")
    het_rel = budget.get("heterogeneity_relative_uncertainty")
    stat_rel = budget.get("statistical_relative_uncertainty")
    vc = budget.get("variance_contributions", {})
    vc_sys = vc.get("systematic")
    vc_het = vc.get("heterogeneity")
    vc_stat = vc.get("statistical")

    # Sync interpretation prose in 7_conclusion.html, 5_results.html, and 6_discussion.html
    if total_rel is not None:
        text = re.sub(
            r"total relative uncertainty is [\d.]+%",
            lambda _m: f"total relative uncertainty is {total_rel*100:.1f}%",
            text,
        )
        text = re.sub(
            r"total relative uncertainty of [\d.]+%",
            lambda _m: f"total relative uncertainty of {total_rel*100:.1f}%",
            text,
        )
        text = re.sub(
            r"total relative uncertainty [\d.]+% on",
            lambda _m: f"total relative uncertainty {total_rel*100:.1f}% on",
            text,
        )
    if sys_rel is not None:
        text = re.sub(
            r"systematic uncertainty contributing [\d.]+%",
            lambda _m: f"systematic uncertainty contributing {sys_rel*100:.1f}%",
            text,
        )
        text = re.sub(
            r"systematic uncertainty \([\d.]+%\)",
            lambda _m: f"systematic uncertainty ({sys_rel*100:.1f}%)",
            text,
        )
    if het_rel is not None:
        text = re.sub(
            r"heterogeneity contributing [\d.]+%",
            lambda _m: f"heterogeneity contributing {het_rel*100:.1f}%",
            text,
        )
        text = re.sub(
            r"heterogeneity \([\d.]+%\)",
            lambda _m: f"heterogeneity ({het_rel*100:.1f}%)",
            text,
        )
        text = re.sub(
            r"heterogeneity \([\d.]+% relative uncertainty",
            lambda _m: f"heterogeneity ({het_rel*100:.1f}% relative uncertainty",
            text,
        )

    # Sync "X% of variance" wording in 6_discussion.html
    if vc_het is not None:
        text = re.sub(
            r"heterogeneity \([\d.]+% relative uncertainty; [\d.]+% of variance\)",
            lambda _m: f"heterogeneity ({het_rel*100:.1f}% relative uncertainty; {vc_het*100:.1f}% of variance)",
            text,
        )
    if vc_sys is not None:
        text = re.sub(
            r"systematic component \([\d.]+% relative uncertainty; [\d.]+% of variance\)",
            lambda _m: f"systematic component ({sys_rel*100:.1f}% relative uncertainty; {vc_sys*100:.1f}% of variance)",
            text,
        )

    # Sync the "Total Relative Uncertainty" <ul> block in 5_results.html
    if all(v is not None for v in (stat_rel, sys_rel, het_rel, total_rel)):
        text = re.sub(
            r"(<p>\s*<strong>Total Relative Uncertainty:</strong>\s*</p>\s*<ul>)([\s\S]*?)(</ul>)",
            lambda _m: (
                f"{_m.group(1)}\n"
                f"    <li>Statistical: {stat_rel*100:.2f}%</li>\n"
                f"    <li>Systematic: {sys_rel*100:.2f}%</li>\n"
                f"    <li>Heterogeneity: {het_rel*100:.2f}%</li>\n"
                f"    <li>Total: {total_rel*100:.2f}%</li>\n"
                f"{_m.group(3)}"
            ),
            text,
        )

    # Sync the "Variance Contributions" <ul> block in 5_results.html
    if all(v is not None for v in (vc_stat, vc_sys, vc_het)):
        text = re.sub(
            r"(<p>\s*<strong>Variance Contributions:</strong>\s*</p>\s*<ul>)([\s\S]*?)(</ul>)",
            lambda _m: (
                f"{_m.group(1)}\n"
                f"    <li>Statistical: {vc_stat*100:.1f}%</li>\n"
                f"    <li>Systematic: {vc_sys*100:.1f}%</li>\n"
                f"    <li>Heterogeneity: {vc_het*100:.1f}%</li>\n"
                f"{_m.group(3)}"
            ),
            text,
        )
    return text


def _sync_step032_spearman(text: str, step032: Dict[str, Any]) -> str:
    """Synchronize Spearman correlation claims with Step 032 JSON."""
    enhanced = step032.get("enhanced_evidence", {})
    spearman = enhanced.get("altitude_anomaly_spearman", {})
    rho = spearman.get("rho")
    p_val = spearman.get("p_value")
    n = spearman.get("n")

    if rho is not None and n is not None:
        # Match LaTeX-style Spearman claim (e.g., 7_conclusion.html) with explicit p-value
        if p_val is not None:
            text = re.sub(
                r"Spearman\s*\$?\\rho\$?\s*=\s*[+-]?[\d.]+\s*\(not\s+significant,\s*\$?p\s*=?\s*[\d.]+\$?,\s*\$?n\s*=\s*\d+\$?\)",
                lambda _m: (
                    f"Spearman $\\rho$ = {rho:.2f} (not significant, "
                    f"$p$ = {p_val:.2f}, $n$ = {n})"
                ),
                text,
            )
            # Match legacy hardcoded format like 'Spearman $\rho = -0.85$, $p = 0.004$'
            text = re.sub(
                r"Spearman\s*\$?\\rho\$?\s*=\s*[+-]?[\d.]+\$?\s*,\s*\$?p\s*=?\s*[\d.]+\$?",
                lambda _m: (
                    f"Spearman $\\rho$ = {rho:.2f} (not significant, "
                    f"$p$ = {p_val:.2f}, $n$ = {n})"
                ),
                text,
            )
        # Match Unicode-rho claim without p-value (e.g., 9_reproducibility.html)
        text = re.sub(
            r"(Key result:\s*)?Altitude-anomaly\s+correlation\s+\$?\\rho\$?\s*=\s*[+-]?[\d.]+\s*\(not\s+significant,\s*\$?n\s*=\s*\d+\$?\)",
            lambda _m: (
                f"{_m.group(1) or ''}Altitude-anomaly correlation $\\rho$ = {rho:.2f} "
                f"(not significant, $n$ = {n})"
            ),
            text,
        )
        # Also match plain Unicode rho (not preceded by dollar signs)
        text = re.sub(
            r"Altitude-anomaly\s+correlation\s+ρ\s*=\s*[+-]?[\d.]+\s*\(not\s+significant,\s*n\s*=\s*\d+\)",
            lambda _m: (
                f"Altitude-anomaly correlation $\\rho$ = {rho:.2f} "
                f"(not significant, $n$ = {n})"
            ),
            text,
        )
    return text


def _sync_step036_ppn(text: str, step036: Dict[str, Any]) -> str:
    """Synchronize PPN margin claims with Step 036 JSON."""
    summary = step036.get("summary", {})
    ppn = summary.get("ppn_analysis", {})
    solar = ppn.get("solar_screened", {})

    worst_margin = ppn.get("ppn_worst_fit_safety_margin")
    margin_sun_surface = solar.get("margin_sun_surface")
    margin_sun_path = solar.get("margin_sun_path_cassini")
    max_gamma = ppn.get("max_fitted_gamma_deviation")

    if worst_margin is not None:
        text = re.sub(
            r"factor of about [\d.]+(?=\s*\(Section 4\.6\.1a\))",
            lambda _m: f"factor of about {round(float(worst_margin)):.0f}",
            text,
        )
        text = re.sub(
            r"worst-fit safety margin of [\d.]+x",
            lambda _m: f"worst-fit safety margin of {float(worst_margin):.1f}x",
            text,
        )
    if margin_sun_path is not None:
        text = re.sub(
            r"margin of about [\d.]+(?= for the largest gated coupling along the Cassini radio path)",
            lambda _m: f"margin of about {float(margin_sun_path):.1f}",
            text,
        )
        text = re.sub(
            r"margin of about [\d.]+(?= rather than the much larger Earth-screened factors)",
            lambda _m: f"margin of about {float(margin_sun_path):.1f}",
            text,
        )
    if margin_sun_surface is not None:
        text = re.sub(
            r"margin of about [\d.]+(?= at the solar surface)",
            lambda _m: f"margin of about {float(margin_sun_surface):.1f}",
            text,
        )
    if max_gamma is not None:
        # Update worst-fit gamma deviation claim if present
        text = re.sub(
            r"\|\$?\\gamma\s*-\s*1\$?\|\s*\\approx\s*\$?1\.0\s*\\times\s*10\^\{-6\}\$?",
            lambda _m: f"$|\\gamma - 1| \\approx {max_gamma:.1e}$",
            text,
        )
    return text


def _format_tex_scientific(value: float, sigfigs: int = 3) -> str:
    if value == 0.0:
        return "0"
    exponent = int(math.floor(math.log10(abs(float(value)))))
    mantissa = float(value) / (10**exponent)
    digits = max(sigfigs - 1, 0)
    mantissa_rounded = round(mantissa, digits)
    mantissa_fmt = f"{mantissa_rounded:.{digits}f}" if digits > 0 else f"{mantissa_rounded:.0f}"
    return f"{mantissa_fmt} \\times 10^{{{exponent}}}"


def _sync_step008_beta_claims(text: str, step008: Dict[str, Any]) -> str:
    stats = step008.get("overall_analysis", {}).get("beta_statistics", {})
    robust = step008.get("robustness_analysis", {})
    loo = robust.get("leave_one_out", {}).get("leave_one_out_results", {})
    bootstrap = robust.get("bootstrap", {})

    if not stats or not bootstrap or not loo:
        raise RuntimeError("Step 008 values missing required fields for publication sync.")

    # Older `re.sub(..., repl_string, ...)` passes corrupted TeX when `repl` contains `\approx`
    # (`\a` → BEL) or `\times` (`\t` → TAB). Callable replacers avoid replacement-template parsing.
    text = text.replace("\x07pprox", "\\approx")

    beta_min = float(stats["min"])
    beta_max = float(stats["max"])
    beta_w = float(stats["weighted_mean"])
    beta_w_sigma = float(stats["weighted_uncertainty"])

    factor_span = beta_max / beta_min if beta_min != 0.0 else float("inf")

    beta_min_tex = _format_tex_scientific(beta_min, 3)
    beta_max_tex = _format_tex_scientific(beta_max, 3)
    beta_w_tex = _format_tex_scientific(beta_w, 3)
    beta_w_sigma_tex = _format_tex_scientific(beta_w_sigma, 3)

    stability = float(robust["leave_one_out"]["stability_coefficient"])
    stability_fmt = f"{stability:.3f}"

    ci95_lo = float(bootstrap["ci_95_lower"])
    ci95_hi = float(bootstrap["ci_95_upper"])
    boot_median = float(bootstrap["bootstrap_median"])
    boot_median_tex = _format_tex_scientific(boot_median, 3)
    ci95_lo_tex = _format_tex_scientific(ci95_lo, 3)
    ci95_hi_tex = _format_tex_scientific(ci95_hi, 3)

    loo_near = float(loo["NEAR"]["beta_without_this"])
    loo_gal = float(loo["Galileo_1990"]["beta_without_this"])
    loo_ros = float(loo["Rosetta_2005"]["beta_without_this"])

    loo_near_tex = _format_tex_scientific(loo_near, 3)
    loo_gal_tex = _format_tex_scientific(loo_gal, 3)
    loo_ros_tex = _format_tex_scientific(loo_ros, 3)

    # Results section canonical paragraph (site/components/5_results.html)
    text = re.sub(
        r"The three gated fitted \$β\$ values span a factor of about 4\.0 \(\$5\.03 \\times 10\^{-4}\$ to\s*\$2\.02 \\times 10\^{-3}\$\),",
        lambda _m: (
            f"The three gated fitted $β$ values span a factor of about {factor_span:.1f} "
            f"(${beta_min_tex}$ to ${beta_max_tex}$),"
        ),
        text,
    )
    text = text.replace(
        "$β = 5.65 \\times 10^{-4} \\pm 2.79 \\times 10^{-5}$",
        f"$β = {beta_w_tex} \\pm {beta_w_sigma_tex}$",
    )
    text = re.sub(
        r"\(stability coefficient [\d.]+\s*&lt;\s*0\.5\)",
        lambda _m: f"(stability coefficient {stability_fmt} &lt; 0.5)",
        text,
    )

    # Abstract (site/components/1_abstract.html): update the span + weighted mean fragment.
    abstract_pattern = (
        r"span β = 5\.03×10⁻⁴ to 2\.02×10⁻³ \(factor 4\.0\), "
        r"with inverse-variance weighted mean β = 5\.65×10⁻⁴ ± 2\.79×10⁻⁵\."
    )
    abstract_replacement = (
        f"span β = ${beta_min_tex}$ to ${beta_max_tex}$ (factor {factor_span:.1f}), "
        f"with inverse-variance weighted mean β = ${beta_w_tex} \\pm {beta_w_sigma_tex}$."
    )
    text = re.sub(
        r"span β = \$[^\$]+\$\s+to \$[^\$]+\$\s*\(factor [\d.]+\),\s*"
        r"with inverse-variance weighted mean β = \$[^\$]+\$\.",
        lambda _: abstract_replacement,
        text,
    )
    # Repair any prior missing closing $ in the abstract sentence.
    text = text.replace(
        f"with inverse-variance weighted mean β = ${beta_w_tex} \\pm {beta_w_sigma_tex}.",
        f"with inverse-variance weighted mean β = ${beta_w_tex} \\pm {beta_w_sigma_tex}$.",
    )
    text = re.sub(abstract_pattern, lambda _: abstract_replacement, text)
    # If an earlier sync wrote scientific e-notation, normalize it to TeX.
    text = text.replace(
        f"span β = {beta_min:.2e} to {beta_max:.2e} (factor {factor_span:.1f}), with inverse-variance weighted mean β = {beta_w:.2e} ± {beta_w_sigma:.2e}.",
        abstract_replacement,
    )

    # Discussion limitation paragraph (site/components/6_discussion.html): update the variance/robustness summary.
    text = re.sub(
        r"The fitted β values span 5\.03×10⁻⁴ to 2\.02×10⁻³ across the three primary ensemble fits—a factor of 4\.0\.",
        lambda _m: (
            f"The fitted β values span ${beta_min_tex}$ to ${beta_max_tex}$ across the three primary ensemble fits—a factor of {factor_span:.1f}."
        ),
        text,
    )
    text = text.replace(
        f"The fitted β values span {beta_min:.2e} to {beta_max:.2e} across the three primary ensemble fits—a factor of {factor_span:.1f}.",
        f"The fitted β values span ${beta_min_tex}$ to ${beta_max_tex}$ across the three primary ensemble fits—a factor of {factor_span:.1f}.",
    )
    text = re.sub(
        r"Cross-validation indicates model stability \(stability coefficient [\d.]+\s*&lt;\s*0\.5\)\.",
        f"Cross-validation indicates model stability (stability coefficient {stability_fmt} &lt; 0.5).",
        text,
    )
    text = re.sub(
        r"The inverse-variance weighted mean β = 5\.65×10⁻⁴ ± 2\.79×10⁻⁵ is representative across flyby geometries\.",
        lambda _: f"The inverse-variance weighted mean β = ${beta_w_tex} \\pm {beta_w_sigma_tex}$ is representative across flyby geometries.",
        text,
    )
    text = re.sub(
        r"The inverse-variance weighted mean β = [\d.eE+-]+ ± [\d.eE+-]+ is representative across flyby geometries\.",
        lambda _: (
            f"The inverse-variance weighted mean β = ${beta_w_tex} \\pm {beta_w_sigma_tex}$ is representative "
            "across flyby geometries."
        ),
        text,
    )

    # Robustness verification paragraph is easiest to sync as a whole.
    robustness_replacement = (
        "<p>Robustness verification: Step 008 parametric bootstrap ($10^4$ draws) yields median "
        f"$\\beta \\approx {boot_median_tex}$ with 95% interval $[{ci95_lo_tex},\\,{ci95_hi_tex}]$, "
        "and leave-one-out recomputations "
        f"${loo_near_tex}$ (without NEAR), ${loo_gal_tex}$ (without Galileo 1990), and ${loo_ros_tex}$ (without Rosetta 2005). "
        f"The stability coefficient $\\approx {stability_fmt}$ is below the 0.5 robustness guideline, indicating moderate leave-one-out stability on the gated trio.</p>"
    )
    text = re.sub(
        r"<p>Robustness verification:[\s\S]*?</p>",
        lambda _: robustness_replacement,
        text,
    )

    # Conclusion (site/components/7_conclusion.html): update the summary estimate + robustness coefficient.
    text = text.replace(
        "mean $\\beta \\approx 5.65\\times10^{-4} \\pm 2.79\\times10^{-5}$",
        f"mean $\\beta \\approx {beta_w_tex} \\pm {beta_w_sigma_tex}$",
    )
    text = re.sub(
        r"show (?:moderate|high) stability \(coefficient \$\\approx [\d.]+\$\), with NEAR as the dominant lever\.",
        lambda _m: (
            f"show moderate stability (coefficient $\\approx {stability_fmt}$), "
            "with NEAR as the dominant lever."
        ),
        text,
    )
    text = re.sub(
        r"Bootstrap resampling \(\$10\^4\$ draws\) and leave-one-out recomputations in Step 008 show (?:moderate|high) stability \(coefficient \$\\approx [\d.]+\$\)\.",
        lambda _m: (
            "Bootstrap resampling ($10^4$ draws) and leave-one-out recomputations in Step 008 show moderate stability "
            f"(coefficient $\\approx {stability_fmt}$)."
        ),
        text,
    )

    text = re.sub(
        r"The stability coefficient \(relative standard deviation of LOO estimates\s+"
        r"divided by their mean\) is [\d.]+, indicating (?:moderate|high) robustness",
        lambda _m: (
            "The stability coefficient (relative standard deviation of LOO estimates "
            f"divided by their mean) is {stability_fmt}, indicating moderate robustness"
        ),
        text,
    )
    text = re.sub(
        r"The stability coefficient is \$[\d.]+\$, indicating (?:moderate|high) robustness",
        lambda _m: f"The stability coefficient is ${stability_fmt}$, indicating moderate robustness",
        text,
    )
    text = re.sub(
        r"Leave-one-out analysis on the gated trio shows moderate stability \(coefficient \$\\approx [\d.]+\$\)\.",
        lambda _m: (
            f"Leave-one-out analysis on the gated trio shows moderate stability "
            f"(coefficient $\\approx {stability_fmt}$)."
        ),
        text,
    )
    text = re.sub(
        r"Leave-one-out analysis on the gated trio shows high stability \(coefficient \$\\approx [\d.]+\$\)\.",
        lambda _m: (
            f"Leave-one-out analysis on the gated trio shows moderate stability "
            f"(coefficient $\\approx {stability_fmt}$)."
        ),
        text,
    )
    text = re.sub(
        r"Cross-validation reports stability coefficient \$\\approx [\d.]+\$(?:\.[\d.]+\$)*\.",
        lambda _m: f"Cross-validation reports stability coefficient $\\approx {stability_fmt}$.",
        text,
    )

    return text


def _normalize_abstract_mathjax_backslashes(text: str) -> str:
    """
    Inline math in `1_abstract.html` historically used doubled TeX escapes (e.g. `\\\\approx` in the source file);
    normalize to a single command backslash so MathJax/KaTeX parse inline `$...$` segments reliably.
    """
    for cmd in ("approx", "Delta", "times", "pm", "sigma", "alpha", "beta", "gamma", "lambda", "rho", "phi", "theta"):
        text = text.replace("\\\\" + cmd, "\\" + cmd)
    return text


def sync_publication_artifacts(project_root: Path) -> List[Path]:
    expected = load_expected_publication_values(project_root)
    step039 = load_step039_values(project_root)
    step019 = load_step019_values(project_root)
    step040 = load_step040_values(project_root)
    step008 = load_step008_values(project_root)
    step025 = load_step025_values(project_root)
    step032 = load_step032_values(project_root)
    step036 = load_step036_values(project_root)
    updated: List[Path] = []

    for path in _publication_sources(project_root):
        original = path.read_text(encoding="utf-8")
        synced = _sync_text(original, expected, step039)
        synced = _sync_step008_beta_claims(synced, step008)
        synced = _sync_step019_table3b_integration(synced, step019)
        synced = _sync_step040_table8(synced, step040)
        synced = _sync_step040_table9b(synced, step040)
        synced = _sync_step040_regression_prose(synced, step040)
        synced = _sync_step025_uncertainty(synced, step025)
        synced = _sync_step032_spearman(synced, step032)
        synced = _sync_step036_ppn(synced, step036)
        if synced != original:
            path.write_text(synced, encoding="utf-8")
            updated.append(path)

    abstract_path = project_root / "site/components/1_abstract.html"
    if abstract_path.exists():
        abs_text = abstract_path.read_text(encoding="utf-8")
        abs_norm = _normalize_abstract_mathjax_backslashes(abs_text)
        if abs_norm != abs_text:
            abstract_path.write_text(abs_norm, encoding="utf-8")
            if abstract_path not in updated:
                updated.append(abstract_path)

    return updated
