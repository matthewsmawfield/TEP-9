#!/usr/bin/env python3
"""
Manuscript Numerical Audit Script

Systematically compares numerical claims in the TEP-EFA manuscript HTML components
against the authoritative pipeline JSON outputs. Emits a machine-readable diff
for CI integration.

Usage:
    python scripts/utils/audit_manuscript_numbers.py

Output:
    results/manuscript_audit.json
"""

import json
import re
import os
from pathlib import Path
from collections import defaultdict

# Paths
REPO_ROOT = Path(__file__).resolve().parents[2]
COMPONENTS_DIR = REPO_ROOT / "site" / "components"
RESULTS_DIR = REPO_ROOT / "results"
AUDIT_OUTPUT = RESULTS_DIR / "manuscript_audit.json"

# Tolerance for floating-point comparison
REL_TOL = 0.02    # 2% relative
ABS_TOL_MM_S = 0.001  # 0.001 mm/s absolute for velocity shifts
ABS_TOL_BETA = 1e-7   # for beta values
ABS_TOL_GENERIC = 0.05 # for generic numbers


def load_json(path):
    with open(path, "r") as f:
        return json.load(f)


def parse_scientific(num_str):
    """Parse strings like '5.64 \\times 10^{-4}' or '1.97 \\times 10^{-4}' into float."""
    num_str = num_str.strip()
    # Handle LaTeX scientific notation
    m = re.match(r"([\d.]+)\s*\\times\s*10\^\{-?(\d+)\}", num_str)
    if m:
        base = float(m.group(1))
        exp = -int(m.group(2))
        return base * (10 ** exp)
    # Handle plain scientific like 1.23e-4
    try:
        return float(num_str)
    except ValueError:
        return None


def extract_numbers_from_html(html_content):
    """Extract all numerical claims from HTML, keeping surrounding context."""
    claims = []
    # Remove HTML tags for cleaner parsing but keep line structure
    lines = html_content.split("\n")
    for line_no, line in enumerate(lines, 1):
        text = re.sub(r"<[^>]+>", " ", line)
        # Find floating-point numbers
        for m in re.finditer(r"([+-]?[\d.,]+(?:×10\^-?\d+)?(?:e[+-]?\d+)?)", text):
            raw = m.group(1).replace(",", "")
            # Try to parse as float
            val = parse_number(raw)
            if val is not None and abs(val) < 1e20:
                claims.append({
                    "line": line_no,
                    "raw": raw,
                    "value": val,
                    "context": text.strip()[:200]
                })
    return claims


def parse_number(s):
    s = s.strip()
    # Handle ×10^n notation
    m = re.match(r"([+-]?[\d.]+)\s*×\s*10\^\{?([+-]?\d+)\}?", s)
    if m:
        return float(m.group(1)) * (10 ** float(m.group(2)))
    # Handle plain ×10^n without braces
    m = re.match(r"([+-]?[\d.]+)\s*×\s*10\^([+-]?\d+)", s)
    if m:
        return float(m.group(1)) * (10 ** float(m.group(2)))
    # Handle standard float
    try:
        return float(s)
    except ValueError:
        return None


def matches(val1, val2, rel_tol=REL_TOL, abs_tol=ABS_TOL_GENERIC):
    if val1 is None or val2 is None:
        return False
    if abs(val1 - val2) <= abs_tol:
        return True
    if val1 == 0 or val2 == 0:
        return False
    return abs(val1 - val2) / max(abs(val1), abs(val2)) <= rel_tol


def audit_component(html_path, json_sources):
    """
    Audit a single HTML component against its corresponding JSON sources.
    json_sources is a list of dicts with keys: json_path, extract_fn, section_name
    """
    with open(html_path, "r") as f:
        html = f.read()

    claims = extract_numbers_from_html(html)
    mismatches = []
    matched = []

    for claim in claims:
        val = claim["value"]
        found_match = False
        for src in json_sources:
            try:
                data = load_json(src["json_path"])
            except Exception as e:
                continue
            extracted = src["extract_fn"](data)
            for ext in extracted:
                if matches(val, ext["value"]):
                    matched.append({
                        **claim,
                        "source": src["section_name"],
                        "matched_key": ext.get("key", "")
                    })
                    found_match = True
                    break
            if found_match:
                break

        if not found_match:
            mismatches.append(claim)

    return mismatches, matched


# ---------------------------------------------------------------------------
# JSON extractors
# ---------------------------------------------------------------------------

def extract_step008_overall(data):
    """Extract key overall statistics from step008."""
    out = []
    oa = data.get("overall_analysis", {})
    if "recommended_beta" in oa:
        out.append({"key": "recommended_beta", "value": oa["recommended_beta"]})
    if "recommended_beta_eff" in oa:
        out.append({"key": "recommended_beta_eff", "value": oa["recommended_beta_eff"]})
    if "recommended_uncertainty" in oa:
        out.append({"key": "recommended_uncertainty", "value": oa["recommended_uncertainty"]})
    # Beta statistics
    bs = oa.get("beta_statistics", {})
    for k in ["mean", "std", "weighted_mean", "weighted_uncertainty", "random_effects_mean", "random_effects_uncertainty"]:
        if k in bs:
            out.append({"key": f"beta_statistics.{k}", "value": bs[k]})
    # Per-flyby fits
    for name, fit_data in data.get("individual_fits", {}).items():
        fit = fit_data.get("fit", {})
        for k in ["beta_fitted", "beta_eff", "uncertainty", "snr", "ppn_gamma_deviation"]:
            if k in fit and fit[k] is not None:
                out.append({"key": f"{name}.{k}", "value": fit[k]})
        pred = fit_data.get("tep_predictions", {})
        for k in ["dv_tep_mm_s", "dv_grad_mm_s", "dv_disf_mm_s"]:
            if k in pred:
                out.append({"key": f"{name}.{k}", "value": pred[k]})
        obs = fit_data.get("observed", {})
        for k in ["dv_obs_mm_s", "sigma_mm_s"]:
            if k in obs:
                out.append({"key": f"{name}.{k}", "value": obs[k]})
    return out


def extract_step026(data):
    """Extract key values from step026 model comparison."""
    out = []
    for section in ["log_likelihoods", "information_criteria", "bayes_factors"]:
        sec = data.get(section, {})
        for k, v in sec.items():
            if isinstance(v, (int, float)):
                out.append({"key": f"{section}.{k}", "value": v})

    # Fitted parameters
    fp = data.get("fitted_parameters", {})
    for model, params in fp.items():
        for k, v in params.items():
            if isinstance(v, (int, float)):
                out.append({"key": f"fitted_parameters.{model}.{k}", "value": v})

    # Pre-specified parameters
    pp = data.get("pre_specified_parameters", {})
    for k, v in pp.items():
        if isinstance(v, (int, float)):
            out.append({"key": f"pre_specified_parameters.{k}", "value": v})

    # Full catalog model comparison
    fcmc = data.get("full_catalog_model_comparison", {})
    for section in ["log_likelihoods", "information_criteria", "bayes_factors"]:
        sec = fcmc.get(section, {})
        for k, v in sec.items():
            if isinstance(v, (int, float)):
                out.append({"key": f"full_catalog.{section}.{k}", "value": v})
    fp2 = fcmc.get("fitted_parameters", {})
    for model, params in fp2.items():
        for k, v in params.items():
            if isinstance(v, (int, float)):
                out.append({"key": f"full_catalog.fitted_parameters.{model}.{k}", "value": v})

    # Sign agreement model comparison
    samc = data.get("sign_agreement_model_comparison", {})
    for section in ["log_likelihoods", "information_criteria", "bayes_factors"]:
        sec = samc.get(section, {})
        for k, v in sec.items():
            if isinstance(v, (int, float)):
                out.append({"key": f"sign_agreement.{section}.{k}", "value": v})
    fp3 = samc.get("fitted_parameters", {})
    for model, params in fp3.items():
        for k, v in params.items():
            if isinstance(v, (int, float)):
                out.append({"key": f"sign_agreement.fitted_parameters.{model}.{k}", "value": v})

    return out


def extract_step010(data):
    out = []
    for k, v in data.items():
        if isinstance(v, (int, float)):
            out.append({"key": k, "value": v})
    return out


def extract_step011(data):
    out = []
    phys = data.get("physics", {})
    for k, v in phys.items():
        if isinstance(v, (int, float)):
            out.append({"key": f"physics.{k}", "value": v})
    for name, res in data.get("individual_results", {}).items():
        for k in ["dv_integrated_mm_s", "dv_perigee_mm_s", "difference_mm_s", "ratio_integrated_perigee", "path_length_km", "n_integration_points"]:
            if k in res:
                out.append({"key": f"{name}.{k}", "value": res[k]})
        mod = res.get("modulation_factors", {})
        for k, v in mod.items():
            if isinstance(v, (int, float)):
                out.append({"key": f"{name}.modulation.{k}", "value": v})
    for item in data.get("modulation_summary", []):
        name = item["name"]
        for k, v in item.items():
            if k != "name" and isinstance(v, (int, float)):
                out.append({"key": f"modulation_summary.{name}.{k}", "value": v})
    return out


def extract_step015(data):
    out = []
    ga = data.get("geometry_analysis", {})
    for fb in ga.get("per_flyby", []):
        name = fb["name"]
        for k, v in fb.items():
            if k != "name" and isinstance(v, (int, float)):
                out.append({"key": f"geometry_analysis.{name}.{k}", "value": v})
    for k in ["g_eff_median", "g_eff_mean", "g_eff_std", "g_eff_min", "g_eff_max"]:
        if k in ga:
            out.append({"key": k, "value": ga[k]})
    corr = ga.get("correlations", {})
    for param, vals in corr.items():
        for k, v in vals.items():
            out.append({"key": f"correlations.{param}.{k}", "value": v})
    reg = ga.get("regression_coefficients", {})
    for k, v in reg.items():
        out.append({"key": f"regression.{k}", "value": v})
    return out


def extract_step019(data):
    out = []
    for name, vals in data.get("per_flyby", {}).items():
        for k, v in vals.items():
            out.append({"key": f"{name}.{k}", "value": v})
    return out


def extract_step039(data):
    out = []
    meta = data.get("metadata", {})
    for k, v in meta.items():
        if isinstance(v, (int, float)):
            out.append({"key": f"metadata.{k}", "value": v})
    for item in data.get("full_catalog_raw_likelihood", {}).get("per_flyby", []):
        name = item["flyby"]
        for k, v in item.items():
            if k != "flyby" and isinstance(v, (int, float)):
                out.append({"key": f"per_flyby.{name}.{k}", "value": v})
    for item in data.get("per_flyby", []):
        name = item["flyby"]
        for k, v in item.items():
            if k != "flyby" and isinstance(v, (int, float)):
                out.append({"key": f"per_flyby_table.{name}.{k}", "value": v})
    return out


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def main():
    sources = [
        {
            "json_path": RESULTS_DIR / "step008_fitting_results.json",
            "extract_fn": extract_step008_overall,
            "section_name": "step008_fitting"
        },
        {
            "json_path": RESULTS_DIR / "step026_stable_model_comparison.json",
            "extract_fn": extract_step026,
            "section_name": "step026_model_comparison"
        },
        {
            "json_path": RESULTS_DIR / "step010_ucd_saturation_results.json",
            "extract_fn": extract_step010,
            "section_name": "step010_ucd"
        },
        {
            "json_path": RESULTS_DIR / "step011_trajectory_integration.json",
            "extract_fn": extract_step011,
            "section_name": "step011_trajectory"
        },
        {
            "json_path": RESULTS_DIR / "step015_hierarchical_bayesian_results.json",
            "extract_fn": extract_step015,
            "section_name": "step015_hierarchical"
        },
        {
            "json_path": RESULTS_DIR / "step019_3d_field_integration_results.json",
            "extract_fn": extract_step019,
            "section_name": "step019_3d_integration"
        },
        {
            "json_path": RESULTS_DIR / "step039_flyby_prediction_table.json",
            "extract_fn": extract_step039,
            "section_name": "step039_predictions"
        },
    ]

    all_mismatches = []
    all_matched = []

    components = sorted(COMPONENTS_DIR.glob("*.html"))
    for comp in components:
        mismatches, matched = audit_component(comp, sources)
        for m in mismatches:
            m["component"] = comp.name
            all_mismatches.append(m)
        for m in matched:
            m["component"] = comp.name
            all_matched.append(m)

    # Group mismatches by component for readability
    grouped = defaultdict(list)
    for m in all_mismatches:
        grouped[m["component"]].append(m)

    summary = {
        "total_claims": len(all_mismatches) + len(all_matched),
        "matched": len(all_matched),
        "mismatches": len(all_mismatches),
        "mismatches_by_component": dict(grouped),
        "matched_by_component": defaultdict(list),
    }

    for m in all_matched:
        summary["matched_by_component"][m["component"]].append(m)
    summary["matched_by_component"] = dict(summary["matched_by_component"])

    # Write audit report
    AUDIT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with open(AUDIT_OUTPUT, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"Audit complete: {summary['matched']}/{summary['total_claims']} claims matched")
    print(f"Mismatches found: {summary['mismatches']}")
    print(f"Report written to: {AUDIT_OUTPUT}")


if __name__ == "__main__":
    main()
