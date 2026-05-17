"""
Step: Anderson Empirical Formula Audit

Computes best-fit (A, B) for Anderson et al. (2008) empirical formula:
    Δv = A * cos_asymmetry + B

and generates comparison tables:
1. Anderson predicted vs TEP predicted vs Observed for all flybys
2. Geometry audit: published declinations, computed cos asymmetry, Horizons-derived

Output: results/step_anderson_audit.json
"""
import json
import math
import numpy as np
from pathlib import Path

ROOT = Path(__file__).parent.parent.parent
RESULTS = ROOT / "results"
DATA = ROOT / "data"


def load_catalog():
    with open(RESULTS / "step003_archival_flyby_catalog.json") as f:
        return json.load(f)


def load_tep_predictions():
    with open(RESULTS / "step007_tep_predictions.json") as f:
        return json.load(f)


def fit_anderson(flybys, include_names=None):
    """
    Fit Δv = A * cos_asymmetry + B to flybys with published anomalies
    and non-null cos_asymmetry values.
    If include_names is given, restrict to those mission names.
    """
    x, y, sigma, names = [], [], [], []
    for fb in flybys:
        name = fb["mission_name"]
        if include_names is not None and name not in include_names:
            continue
        dv = fb.get("published_anomaly_mm_s")
        cos_a = fb.get("cos_asymmetry")
        unc = fb.get("published_anomaly_uncertainty_mm_s", 0.05)
        if dv is not None and cos_a is not None and unc is not None:
            x.append(cos_a)
            y.append(dv)
            sigma.append(unc)
            names.append(name)

    x = np.array(x)
    y = np.array(y)
    sigma = np.array(sigma)
    w = 1.0 / sigma**2

    # Weighted least squares: y = A*x + B
    # Design matrix [x, 1]
    X = np.vstack([x, np.ones_like(x)]).T
    W = np.diag(w)

    # (X^T W X)^{-1} X^T W y
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ y
    coeffs = np.linalg.solve(XtWX, XtWy)
    A_fit, B_fit = coeffs

    # Covariance matrix
    cov = np.linalg.inv(XtWX)
    A_err = math.sqrt(cov[0, 0])
    B_err = math.sqrt(cov[1, 1])

    # Predictions and residuals
    y_pred = A_fit * x + B_fit
    residuals = y - y_pred
    chi2 = np.sum(w * residuals**2)
    dof = len(x) - 2
    rms = math.sqrt(np.mean(residuals**2))

    # R^2 (weighted)
    y_mean = np.average(y, weights=w)
    ss_res = np.sum(w * residuals**2)
    ss_tot = np.sum(w * (y - y_mean)**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return {
        "A_fit_mm_s": float(A_fit),
        "B_fit_mm_s": float(B_fit),
        "A_uncertainty_mm_s": float(A_err),
        "B_uncertainty_mm_s": float(B_err),
        "chi2": float(chi2),
        "dof": int(dof),
        "reduced_chi2": float(chi2 / dof) if dof > 0 else None,
        "rms_residual_mm_s": float(rms),
        "r2_weighted": float(r2),
        "n_data_points": len(x),
        "fit_mission_names": names,
    }


def compute_anderson_predictions(flybys, A, B):
    """Compute Anderson formula prediction for each flyby."""
    results = []
    for fb in flybys:
        cos_a = fb.get("cos_asymmetry")
        dv_obs = fb.get("published_anomaly_mm_s")
        if cos_a is not None:
            dv_anderson = A * cos_a + B
        else:
            dv_anderson = None
        results.append({
            "mission_name": fb["mission_name"],
            "cos_asymmetry": cos_a,
            "declination_in_deg": fb.get("declination_in_deg"),
            "declination_out_deg": fb.get("declination_out_deg"),
            "published_anomaly_mm_s": dv_obs,
            "published_anomaly_uncertainty_mm_s": fb.get("published_anomaly_uncertainty_mm_s"),
            "anderson_prediction_mm_s": dv_anderson,
        })
    return results


def compute_declination_cos_asymmetry(flybys):
    """Compute cos(dec_in) - cos(dec_out) from published declinations."""
    results = []
    for fb in flybys:
        dec_in = fb.get("declination_in_deg")
        dec_out = fb.get("declination_out_deg")
        if dec_in is not None and dec_out is not None:
            computed = math.cos(math.radians(dec_in)) - math.cos(math.radians(dec_out))
        else:
            computed = None
        results.append({
            "mission_name": fb["mission_name"],
            "declination_in_deg": dec_in,
            "declination_out_deg": dec_out,
            "cos_asymmetry_published": fb.get("cos_asymmetry"),
            "cos_asymmetry_computed_from_decs": computed,
        })
    return results


def load_horizons_geometry():
    """Load Horizons-derived geometry where available."""
    # Try to load from step007 which has geometry sections
    tep = load_tep_predictions()
    results = {}
    for name, data in tep.get("predictions", {}).items():
        geom = data.get("geometry", {})
        if geom:
            results[name] = {
                "dec_in_deg": geom.get("dec_in_deg"),
                "dec_out_deg": geom.get("dec_out_deg"),
                "cos_dec_asymmetry": geom.get("cos_dec_asymmetry"),
                "declination_source": geom.get("declination_source"),
            }
    return results


def build_geometry_audit_table(flybys, horizons_geom):
    """Build the geometry audit table (Task D)."""
    rows = []
    for fb in flybys:
        name = fb["mission_name"]
        dec_in = fb.get("declination_in_deg")
        dec_out = fb.get("declination_out_deg")
        cos_pub = fb.get("cos_asymmetry")

        # Computed from Anderson declinations
        if dec_in is not None and dec_out is not None:
            cos_computed = math.cos(math.radians(dec_in)) - math.cos(math.radians(dec_out))
        else:
            cos_computed = None

        # Horizons-derived
        hg = horizons_geom.get(name, {})
        dec_in_hor = hg.get("dec_in_deg")
        dec_out_hor = hg.get("dec_out_deg")
        cos_hor = hg.get("cos_dec_asymmetry")

        rows.append({
            "mission_name": name,
            "declination_in_deg_anderson": dec_in,
            "declination_out_deg_anderson": dec_out,
            "cos_asymmetry_published": cos_pub,
            "cos_asymmetry_computed_from_anderson_decs": round(cos_computed, 4) if cos_computed is not None else None,
            "declination_in_deg_horizons": dec_in_hor,
            "declination_out_deg_horizons": dec_out_hor,
            "cos_asymmetry_horizons": cos_hor,
        })
    return rows


def load_step008_fits():
    """Load per-flyby fitted β values from Step 008."""
    path = RESULTS / "step008_fitting_results.json"
    if not path.exists():
        return {}
    with open(path) as f:
        data = json.load(f)
    fits = {}
    for name, entry in data.get("individual_fits", {}).items():
        fit = entry.get("fit", {})
        beta_fitted = fit.get("beta_fitted")
        beta_ref = entry.get("tep_predictions", {}).get("beta_reference", 1e-4)
        dv_ref = entry.get("tep_predictions", {}).get("dv_tep_mm_s")
        fits[name] = {
            "beta_fitted": beta_fitted,
            "beta_ref": beta_ref,
            "dv_ref": dv_ref,
        }
    # Pooled inverse-variance weighted mean for null flybys
    pooled = data.get("pooled_fit", {})
    fits["_pooled_beta"] = pooled.get("weighted_mean", 5.65e-4)
    return fits


def compute_fitted_beta_tep_prediction(name, tep_data, step008_fits):
    """
    Compute TEP prediction at fitted β using the 3/4 power law:
        dv_fitted = dv_ref * (β_fitted / β_ref)^(3/4)
    For flybys without individual fits, use the pooled β.
    """
    tep_pred = tep_data.get("predictions", {}).get(name, {})
    dv_ref = tep_pred.get("tep_predictions", {}).get("dv_tep_mm_s")
    if dv_ref is None:
        return None

    fit_info = step008_fits.get(name, {})
    beta_fitted = fit_info.get("beta_fitted")
    beta_ref = fit_info.get("beta_ref", 1e-4)

    if beta_fitted is not None and beta_ref is not None and beta_ref > 0:
        return dv_ref * (beta_fitted / beta_ref) ** (3 / 4)

    # Fallback: pooled beta for null flybys
    pooled_beta = step008_fits.get("_pooled_beta", 5.65e-4)
    if beta_ref is not None and beta_ref > 0:
        return dv_ref * (pooled_beta / beta_ref) ** (3 / 4)
    return None


def build_anderson_tep_comparison_table(flybys, anderson_preds, tep_data, step008_fits):
    """Build Anderson vs TEP vs Observed table (Tasks B + D)."""
    rows = []
    for ap in anderson_preds:
        name = ap["mission_name"]
        tep_pred = tep_data.get("predictions", {}).get(name, {})
        tep_dv = None
        if tep_pred:
            tep_dv = tep_pred.get("tep_predictions", {}).get("dv_tep_mm_s")

        tep_fitted = compute_fitted_beta_tep_prediction(name, tep_data, step008_fits)

        rows.append({
            "mission_name": name,
            "cos_asymmetry": ap["cos_asymmetry"],
            "observed_dv_mm_s": ap["published_anomaly_mm_s"],
            "observed_uncertainty_mm_s": ap["published_anomaly_uncertainty_mm_s"],
            "anderson_prediction_mm_s": round(ap["anderson_prediction_mm_s"], 4) if ap["anderson_prediction_mm_s"] is not None else None,
            "tep_prediction_beta_ref_mm_s": round(tep_dv, 4) if tep_dv is not None else None,
            "tep_prediction_fitted_beta_mm_s": round(tep_fitted, 4) if tep_fitted is not None else None,
        })
    return rows


def main():
    catalog = load_catalog()
    flybys = catalog["flybys"]
    tep_data = load_tep_predictions()

    # Task B: Anderson fit and predictions (n=6, all with published anomalies + geometry)
    fit_result = fit_anderson(flybys)
    A = fit_result["A_fit_mm_s"]
    B = fit_result["B_fit_mm_s"]

    # Task C: Anderson n=4 sensitivity check (primary detections only)
    primary_names = {"NEAR", "Galileo_1990", "Cassini", "Rosetta_2005"}
    fit_result_n4 = fit_anderson(flybys, include_names=primary_names)
    A_n4 = fit_result_n4["A_fit_mm_s"]
    B_n4 = fit_result_n4["B_fit_mm_s"]

    anderson_preds = compute_anderson_predictions(flybys, A, B)
    anderson_preds_n4 = compute_anderson_predictions(flybys, A_n4, B_n4)

    # Task D: Load Step 008 fitted β values for TEP fitted-β predictions
    step008_fits = load_step008_fits()
    anderson_comparison = build_anderson_tep_comparison_table(flybys, anderson_preds, tep_data, step008_fits)

    # Task D: Geometry audit
    horizons_geom = load_horizons_geometry()
    geometry_audit = build_geometry_audit_table(flybys, horizons_geom)

    # Cassini-specific analysis
    cassini_anderson = next((r for r in anderson_preds if r["mission_name"] == "Cassini"), {})
    cassini_row = next((r for r in anderson_comparison if r["mission_name"] == "Cassini"), {})

    output = {
        "metadata": {
            "step": "anderson_audit",
            "description": "Anderson empirical formula best-fit and geometry audit",
            "timestamp": "2026-05-16T16:00:00Z",
        },
        "anderson_fit": fit_result,
        "anderson_fit_n4_primary_detections_only": fit_result_n4,
        "cassini_sign_analysis": {
            "published_cos_asymmetry": cassini_anderson.get("cos_asymmetry"),
            "observed_anomaly_mm_s": cassini_anderson.get("published_anomaly_mm_s"),
            "anderson_prediction_mm_s": cassini_row.get("anderson_prediction_mm_s"),
            "tep_prediction_mm_s": cassini_row.get("tep_prediction_beta_ref_mm_s"),
            "conclusion": "Both Anderson empirical formula and TEP model predict negative anomaly for Cassini because the published geometry factor cos_asymmetry = -0.088 is negative. The sign mismatch is inherited from the literature geometry, not a TEP-specific failure.",
        },
        "anderson_vs_tep_vs_observed_table": anderson_comparison,
        "geometry_audit_table": geometry_audit,
    }

    out_path = RESULTS / "step_anderson_audit.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Anderson audit written to {out_path}")
    print(f"  n=6 fit: A = {A:.4f} ± {fit_result['A_uncertainty_mm_s']:.4f} mm/s")
    print(f"  n=6 fit: B = {B:.4f} ± {fit_result['B_uncertainty_mm_s']:.4f} mm/s")
    print(f"  n=6 fit: R² (weighted) = {fit_result['r2_weighted']:.4f}")
    print(f"  n=6 fit: RMS residual = {fit_result['rms_residual_mm_s']:.4f} mm/s")
    print(f"  n=6 fit: Reduced χ² = {fit_result['reduced_chi2']:.4f}" if fit_result['reduced_chi2'] else "  n=6 fit: Reduced χ² = N/A")
    print()
    print(f"  n=4 fit (primary only): A = {A_n4:.4f} ± {fit_result_n4['A_uncertainty_mm_s']:.4f} mm/s")
    print(f"  n=4 fit (primary only): B = {B_n4:.4f} ± {fit_result_n4['B_uncertainty_mm_s']:.4f} mm/s")
    print(f"  n=4 fit (primary only): Reduced χ² = {fit_result_n4['reduced_chi2']:.4f}" if fit_result_n4['reduced_chi2'] else "  n=4 fit: Reduced χ² = N/A")
    print()
    print("--- Anderson vs TEP vs Observed ---")
    for row in anderson_comparison:
        name = row["mission_name"]
        obs = row["observed_dv_mm_s"]
        anderson = row["anderson_prediction_mm_s"]
        tep = row["tep_prediction_beta_ref_mm_s"]
        tep_fit = row["tep_prediction_fitted_beta_mm_s"]
        if obs is not None:
            a_str = f"{anderson:+.3f}" if anderson is not None else "N/A"
            t_str = f"{tep:+.3f}" if tep is not None else "N/A"
            tf_str = f"{tep_fit:+.3f}" if tep_fit is not None else "N/A"
            print(f"  {name:20s}: Obs={obs:+.3f}, Anderson={a_str}, TEP(β_ref)={t_str}, TEP(β_fit)={tf_str}")
    print()
    print("--- Geometry Audit ---")
    for row in geometry_audit:
        name = row["mission_name"]
        pub = row["cos_asymmetry_published"]
        comp = row["cos_asymmetry_computed_from_anderson_decs"]
        hor = row["cos_asymmetry_horizons"]
        if pub is not None or comp is not None or hor is not None:
            print(f"  {name:20s}: Pub={pub}, Computed={comp}, Horizons={hor}")


if __name__ == "__main__":
    main()
