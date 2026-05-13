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


def _round(value: float, digits: int = 2) -> float:
    return round(float(value), digits)


def _format_percent(value: float, digits: int = 1) -> str:
    return f"{_round(value, digits):.{digits}f}%"


def _format_signed_mm(value: float, digits: int = 2) -> str:
    rounded = _round(value, digits)
    if rounded > 0:
        return f"+{rounded:.{digits}f}"
    return f"{rounded:.{digits}f}"


def _mantissa_exponent(value: float) -> Tuple[float, int]:
    if value == 0.0:
        return 0.0, 0
    exp = int(math.floor(math.log10(abs(float(value)))))
    mant = round(float(value) / (10**exp), 1)
    return mant, exp


def _format_bayes_factor_approx(symbol: str, value: float) -> str:
    mant, exp = _mantissa_exponent(value)
    return f"${symbol} \\\\approx {mant} \\\\times 10^{{{exp}}}$"


def _format_delta_bic_approx(value: float) -> str:
    return f"$\\\\Delta$BIC $\\\\approx {value}$"


def load_expected_publication_values(project_root: Path) -> Dict[str, Any]:
    model = json.loads((project_root / "results/step026_stable_model_comparison.json").read_text(encoding="utf-8"))
    variance = json.loads((project_root / "results/step009_variance_analysis.json").read_text(encoding="utf-8"))
    predictions = json.loads((project_root / "results/step007_tep_predictions.json").read_text(encoding="utf-8"))

    info = model["information_criteria"]
    bayes = model["bayes_factors"]
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
        "model": {
            "log_likelihoods": {
                "Null": _round(model["log_likelihoods"]["Null"]),
                "Anderson": _round(model["log_likelihoods"]["Anderson"]),
                "TEP_restricted": _round(model["log_likelihoods"]["TEP_restricted"]),
                "TEP_flexible": _round(model["log_likelihoods"]["TEP_flexible"]),
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
            },
            "bayes_approx": {
                "B10": _format_bayes_factor_approx("B_{10}", bayes["TEP_restricted_vs_Null"]),
                "BA0": _format_bayes_factor_approx("B_{A0}", bayes["Anderson_vs_Null"]),
                "Bf0": _format_bayes_factor_approx("B_{f0}", bayes["TEP_flexible_vs_Null"]),
                "B_vs_anderson": f"$B \\\\approx {_round(bayes['TEP_restricted_vs_Anderson'], 1)}$",
                "dBIC_tep_null": _format_delta_bic_approx(_round(bayes["delta_BIC_TEP_restricted_vs_Null"], 1)),
                "dBIC_anderson_null": _format_delta_bic_approx(_round(bayes["delta_BIC_Anderson_vs_Null"], 1)),
                "dBIC_flex_null": _format_delta_bic_approx(_round(bayes["delta_BIC_TEP_flexible_vs_Null"], 1)),
                "dBIC_tep_anderson": _format_delta_bic_approx(_round(bayes["delta_BIC_TEP_restricted_vs_Anderson"], 1)),
            },
            "akaike_weight_tep": _format_percent(
                100.0 * model["model_selection"]["akaike_weights"]["TEP_restricted"],
                1,
            ),
        },
        "variance": {
            "structural": _format_percent(stages["stage1_structural"]["explained_percent"]),
            "observational": _format_percent(stages["stage2_observational"]["explained_percent"]),
            "environmental": _format_percent(stages["stage3_environmental"]["explained_percent"]),
            "residual": _format_percent(stages["stage4_residual"]["explained_percent"]),
        },
        "components": components,
    }


def _publication_sources(project_root: Path) -> List[Path]:
    sources = [project_root / "15-TEP-EFA-v0.1-Yogyakarta.md"]
    sources.extend(sorted((project_root / "site/components").glob("*.html")))
    return [path for path in sources if path.exists()]


def load_step039_values(project_root: Path) -> Dict[str, Any]:
    payload = json.loads(
        (project_root / "results/step039_flyby_prediction_table.json").read_text(encoding="utf-8")
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
    residual = row.get("universal_beta_residual_mm_s")
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


def _sync_step039_prose(text: str, step039: Dict[str, Any]) -> str:
    summary = step039["metadata"]["raw_classification_summary"]
    true_null = summary["true_null"]
    raw_tension = summary["raw_tension"]
    true_positive = summary["true_positive"]
    juno = next(row for row in step039["rows"] if row["flyby"] == "Juno")
    beta_tex = _latex_beta()

    def _replace_repro_summary(_match: re.Match[str]) -> str:
        return (
            f"{true_null} published null or bound cases are consistent with universal-{beta_tex} predictions "
            f"under the Step 007 geometry envelope; {raw_tension} raw-tension case remains"
        )

    def _replace_step039_block(_match: re.Match[str]) -> str:
        return (
            f"Step 039 classifies universal-{beta_tex} predictions with the Step 007 geometry envelope "
            f"({true_positive} true positives, {true_null} true nulls, {raw_tension} raw-tension case). "
            "Post-OD survival factors are withheld until mission OD configuration yields defensible "
            r"$F_{\rm OD}$ estimates."
        )

    def _replace_table_summary(_match: re.Match[str]) -> str:
        return (
            f"{true_null} published null or bound cases are consistent with universal-{beta_tex} predictions "
            f"under the Step 007 geometry envelope; {raw_tension} raw-tension case remains (Juno). "
            "Cassini and Galileo 1992 are classified as raw true nulls at the refit weighted-mean "
            f"{beta_tex}."
        )

    text = re.sub(
        r"\d+ missions \(MESSENGER, Rosetta 2007\) show published nulls consistent with universal-\$\\beta\$ predictions; \d+ missions \([^)]+\) are raw-tension cases",
        _replace_repro_summary,
        text,
    )
    text = re.sub(
        r"Step 039 identifies raw-tension cases where universal-\$\\beta\$[\s\S]*?likelihood\.",
        _replace_step039_block,
        text,
    )
    text = re.sub(
        r"two published nulls or bounds as raw true nulls\s*\(MESSENGER, Rosetta 2007\), and three cases as raw tensions where the\s*scaled prediction exceeds the \$3\\sigma\$ threshold while the published\s*bound is null or sub-threshold \(Galileo 1992, Juno, and Cassini, where\s*the published \$0\.11\$ mm/s anomaly lies below the \$3\\sigma\$ detection\s*threshold\)",
        _replace_table_summary,
        text,
        flags=re.MULTILINE,
    )
    juno_prediction = _format_signed_mm(juno["raw_tep_prediction_mm_s"], 2)
    juno_observed = _format_signed_mm(juno["observed_dv_mm_s"], 2)
    juno_sigma = juno["observed_sigma_mm_s"]
    def _replace_falsifiability_paragraph(_match: re.Match[str]) -> str:
        return (
            "Juno remains the sole raw-tension case at the refit weighted-mean $\\beta$ "
            f"(${juno_prediction}$ mm/s raw prediction vs. "
            f"${juno_observed} " + r"\pm" + f" {juno_sigma:.2f}$ mm/s observed). "
            "Galileo 1992 and Cassini are classified as raw true nulls under the geometry envelope"
        )

    text = re.sub(
        r"Galileo 1992\s*\(\$\+0\.69\$ mm/s raw prediction vs\. \$0\.00 \\pm 0\.05\$ mm/s observed\) and Juno\s*\(\$\+0\.57\$ mm/s vs\. \$0\.00 \\pm 0\.02\$ mm/s\) remain the strongest null\s*tensions\. Cassini is classified as raw tension because the universal-\$\\beta\$\s*prediction overshoots the published bound at the stated precision",
        _replace_falsifiability_paragraph,
        text,
        flags=re.MULTILINE,
    )
    return text


def _sync_text(text: str, expected: Dict[str, Any], step039: Dict[str, Any]) -> str:
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
        (r"\$B_\{10\}\s*(?:=|\\approx)\s*[\de.+-]+", bayes_approx.get("B10", f"$B_{{10}} = {model['bayes']['tep_vs_null']}")),
        (r"\$B_\{A0\}\s*(?:=|\\approx)\s*[\de.+-]+", bayes_approx.get("BA0", f"$B_{{A0}} = {model['bayes']['anderson_vs_null']}")),
        (r"\$B_\{f0\}\s*(?:=|\\approx)\s*[\de.+-]+", bayes_approx.get("Bf0", f"$B_{{f0}} = {model['bayes']['flex_vs_null']}")),
        (r"against Anderson gives \$B\s*(?:=|\\approx)\s*[\de.+-]+", f"against Anderson gives {bayes_approx.get('B_vs_anderson', f'$B = {model['bayes']['tep_vs_anderson']}')}"),
        (r"Akaike weight for TEP restricted is [\d.]+%", f"Akaike weight for TEP restricted is {model['akaike_weight_tep']}"),
        (r"structural proxy bundle accounts for [\d.]+%", f"structural proxy bundle accounts for {variance['structural']}"),
        (r"observational pipeline effects \(OD filter absorption and systematic uncertainties\) account for [\d.]+%", f"observational pipeline effects (OD filter absorption and systematic uncertainties) account for {variance['observational']}"),
        (r"environmental modulation contributes [\d.]+%", f"environmental modulation contributes {variance['environmental']}"),
        (r"residual \(small-sample statistics, intrinsic scatter, model incompleteness\) accounts for [\d.]+%", f"residual (small-sample statistics, intrinsic scatter, model incompleteness) accounts for {variance['residual']}"),
        (r"tracked structural proxy bundle accounts for [\d.]+% of fitted-β variance", f"tracked structural proxy bundle accounts for {variance['structural']} of fitted-β variance"),
        (r"observational effects account for [\d.]+%", f"observational effects account for {variance['observational']}"),
        (r"the residual accounts for [\d.]+%", f"the residual accounts for {variance['residual']}"),
    ]

    for pattern, replacement in replacements:
        text = re.sub(pattern, replacement, text)

    for label, key in (
        ("Null", "Null"),
        ("Anderson", "Anderson"),
        ("TEP restricted", "TEP_restricted"),
        ("TEP flexible", "TEP_flexible"),
    ):
        text = re.sub(
            rf"<strong>{re.escape(label)}:</strong> log L = [-+]?\d+\.\d+, AIC = [-+]?\d+\.\d+, BIC = [-+]?\d+\.\d+",
            (
                f"<strong>{label}:</strong> log L = {model['log_likelihoods'][key]}, "
                f"AIC = {model['aic'][key]}, BIC = {model['bic'][key]}"
            ),
            text,
        )

    if bayes_approx:
        approx_token = r"(?:=|\\approx|\x07pprox)"
        text = re.sub(rf"\$B_\{{10\}}\s*{approx_token}\s*[\de.+-]+", bayes_approx["B10"], text)
        text = re.sub(rf"\$B_\{{A0\}}\s*{approx_token}\s*[\de.+-]+", bayes_approx["BA0"], text)
        text = re.sub(rf"\$B_\{{f0\}}\s*{approx_token}\s*[\de.+-]+", bayes_approx["Bf0"], text)
    else:
        text = re.sub(r"\$B_\{10\} = [\de.+-]+", f"$B_{{10}} = {model['bayes']['tep_vs_null']}", text)
        text = re.sub(r"\$B_\{A0\} = [\de.+-]+", f"$B_{{A0}} = {model['bayes']['anderson_vs_null']}", text)
        text = re.sub(r"\$B_\{f0\} = [\d.]+", f"$B_{{f0}} = {model['bayes']['flex_vs_null']}", text)

    if bayes_approx:
        # Repair legacy corruption where repeated "\times 10^{...}$" fragments were appended
        # outside the closing "$" of the Bayes-factor math segment.
        for symbol_key, replacement in (
            ("B_{10}", bayes_approx["B10"]),
            ("B_{A0}", bayes_approx["BA0"]),
            ("B_{f0}", bayes_approx["Bf0"]),
        ):
            pattern = rf"(\${re.escape(symbol_key)}[\s\S]*?\$)(?:\s*\\times\s*10\^\{{-?\d+\}}\$)+"
            text = re.sub(pattern, replacement, text)
        text = re.sub(r"\$B\s*\\approx\s*[\d.]+\$+", bayes_approx["B_vs_anderson"], text)
    text = re.sub(
        r"against Anderson gives \$B = [\d.]+",
        f"against Anderson gives {bayes_approx.get('B_vs_anderson', f'$B = {model['bayes']['tep_vs_anderson']}')}",
        text,
    )
    text = re.sub(
        r"Direct comparison of TEP restricted against Anderson gives \$B = [\d.]+",
        f"Direct comparison of TEP restricted against Anderson gives {bayes_approx.get('B_vs_anderson', f'$B = {model['bayes']['tep_vs_anderson']}')}",
        text,
    )
    text = re.sub(
        r"B₁₀ = [\d.]+ \(ΔBIC = [\d.]+, strong evidence\)",
        f"B₁₀ = {model['bayes']['tep_vs_null']} (ΔBIC = {model['bayes']['delta_bic_tep_null']}, strong evidence)",
        text,
    )
    text = re.sub(
        r"TEP restricted vs Anderson empirical B = [\d.]+ \(ΔBIC = [\d.]+, positive evidence\)",
        (
            f"TEP restricted vs Anderson empirical B = {model['bayes']['tep_vs_anderson']} "
            f"(ΔBIC = {model['bayes']['delta_bic_tep_anderson']}, positive evidence)"
        ),
        text,
    )

    def _sync_delta_bic(text_in: str, label: str, value: float) -> str:
        pattern = rf"(<strong>{re.escape(label)}:</strong>[\s\S]*?\$\\Delta\\$BIC\s*\$\\approx\s*)[\d.]+(?=\$)"
        return re.sub(pattern, rf"\g<1>{value}", text_in, count=1)

    def _sync_delta_bic_paren(text_in: str, label: str, value: float) -> str:
        pattern = rf"(<strong>{re.escape(label)}:</strong>[\s\S]*?\(\$\\Delta\\$BIC\s*\$\\approx\s*)[\d.]+(?=\$\))"
        return re.sub(pattern, rf"\g<1>{value}", text_in, count=1)

    for label, value in (
        ("TEP restricted vs Null", model["bayes"]["delta_bic_tep_null"]),
        ("Anderson vs Null", model["bayes"]["delta_bic_anderson_null"]),
        ("TEP flexible vs Null", model["bayes"]["delta_bic_flex_null"]),
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
        text = re.sub(markdown_row, markdown_replacement, text)

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
    return text


def sync_publication_artifacts(project_root: Path) -> List[Path]:
    expected = load_expected_publication_values(project_root)
    step039 = load_step039_values(project_root)
    updated: List[Path] = []

    for path in _publication_sources(project_root):
        original = path.read_text(encoding="utf-8")
        synced = _sync_text(original, expected, step039)
        if synced != original:
            path.write_text(synced, encoding="utf-8")
            updated.append(path)

    return updated
