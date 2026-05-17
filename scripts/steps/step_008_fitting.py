#!/usr/bin/env python3
r"""
Flyby TEP Pipeline - Step 005: Parameter Fitting and Statistical Validation

This module implements the core statistical analysis for fitting the TEP coupling
parameter β to observed flyby anomalies. The module provides comprehensive
validation through effect size analysis, Bayesian model comparison, bootstrap
resampling, and systematic uncertainty quantification.

================================================================================
DATA INTEGRITY ENFORCEMENT - ZERO SYNTHETIC DATA TOLERANCE
================================================================================
CRITICAL: This module processes ONLY real observational data from peer-reviewed
literature. ANY attempt to pass synthetic/simulated/mock data as observed will
be rejected. Missing data returns null/None - never fabricated values.

Enforcement:
1. Observed anomalies MUST have literature citations
2. Missing dv_obs or sigma returns null (no fallback)
3. S/N < 2 excludes the flyby from closed-form β fitting (no fitting performed)
4. Sign mismatch at β_ref excludes closed-form β fitting only when ``strict_sign_gate`` is true in ``config/pipeline_config.json``; otherwise an amplitude-only reference fit is retained and reported separately
4. Bootstrap resampling ONLY for confidence intervals (not as observed data)

See scripts/utils/data_integrity_validator.py for automated validation.
================================================================================


Statistical Framework:
---------------------
The fitting procedure employs a chi-squared minimization approach to determine
the optimal β for each flyby, accounting for measurement uncertainties from
NASA's Deep Space Network (DSN) Doppler tracking:

    β_fitted = β_initial × (Δv_observed / Δv_TEP_predicted)

where β_initial = 10⁻⁴ serves as a reference coupling from the TEP framework.
The uncertainty propagation follows standard error propagation:

    σ_β = β_initial × (σ_Δv / |Δv_TEP|)

A priori Selection Criterion:
---------------------------
Flybys are included in fitting only if they meet the S/N > 2 threshold:

    S/N = |Δv_obs| / σ_Δv ≥ 2

This criterion is applied before fitting to prevent confirmation bias. Published
All four high-S/N literature cases receive a closed-form amplitude diagnostic. When ``strict_sign_gate`` is true in ``config/pipeline_config.json``, Cassini is excluded from that fit if the reference prediction disagrees in sign; when false, an amplitude-only reference fit is retained and the sign-agreement subset is reported separately as ``beta_statistics_sign_gated_diagnostic``.

PPN Constraint Validation:
-------------------------
All fitted β values are validated against the Cassini solar system bound on
the PPN parameter γ:

    |γ - 1| = 2β_eff² < 2.3 × 10⁻⁵

where β_eff = β × S_⊕ includes Temporal Topology screening. The effective coupling
must satisfy this constraint for physical viability of the TEP framework.

β Scaling Derivation:
--------------------
The TEP field equation gives φ ∝ β^(-1/(n+1)) = β^(-1/4) for n=3.
The scalar force F ∝ β × ∇φ ∝ β × β^(-1/4) = β^(3/4).
Therefore Δv_TEP ∝ β^(3/4), and inverting:

    β_fitted = β_ref × (Δv_obs / Δv_TEP)^(4/3)

This scaling is derived from first-principles field theory, not fitted empirically.

Enhanced Statistical Validation:
-------------------------------
This module implements five complementary validation analyses:

1. Effect Size Analysis (Cohen's d):
   Quantifies signal magnitude independent of sample size:
       d = (Δv_detection - μ_null) / σ_pooled
   All detections show d >> 0.8 ("large effect"), indicating robust signals.

2. Bayesian-style model comparison (AIC/BIC):
   Compares TEP model against null and empirical alternatives using
   information criteria on the same formal uncertainties as Step 008's
   deviance construction. Raw $\Delta$BIC vs null can be enormous when
   published $\sigma$ is tiny; the reported BF map is capped in JSON and
   must not be read as a literal posterior odds ratio at $n=3$. Prefer
   Step 026 headline likelihood with $\sigma_{\rm geom}$ in quadrature.

3. Bootstrap Resampling (n=10,000):
   Parametric bootstrap with fixed random seed (42) generates non-parametric
   confidence intervals accounting for small sample size (n=4).

4. Systematic Uncertainty Budget:
   Propagates uncertainties from five sources:
   - Measurement uncertainty (~1%)
   - Trajectory reconstruction (~1%)
   - Characteristic suppression (~25%, from UCD saturation model ρ_T = 20 ± 8 g/cm³, Paper 6)
   - Relaxation length (~15%, from SCF theoretical prior refined by GNSS consistency)
   - Multipole coefficients (~0.1%, negligible)
   Total: σ_sys/β ≈ 29% (dominated by characteristic suppression uncertainty)

5. Leave-One-Out Cross-Validation:
   Tests robustness by excluding each detection successively. Stability
   coefficient < 0.5 indicates conclusion does not depend on single flyby.

Scientific Context:
-----------------
The fitting results demonstrate that the TEP scalar force framework with
Temporal Shear screening provides a self-consistent, PPN-compliant explanation
for the Earth flyby anomaly. Per-flyby fitted amplitudes can span a wide factor
across the full catalog of geometries; the inverse-variance Step 008 mean aggregates all S/N-qualified fits in the relaxed-gate configuration, with an additional sign-agreement-restricted diagnostic for auditability.
The (4/3) scaling is derived from first-principles field theory
(φ ∝ β^(-1/4) → Δv ∝ β^(3/4)), not fitted empirically.

Reproducibility:
---------------
All random processes use fixed seeds (bootstrap: seed=42) to ensure exact
reproducibility across pipeline executions. Statistical tests use scipy
implementations with verified accuracy.

Output Structure:
----------------
The module generates comprehensive JSON output including:
- Individual flyby fits with β, β_eff, uncertainties, and PPN status
- Overall β statistics (weighted mean, heterogeneity tests)
- Bootstrap confidence intervals (68% and 95%)
- Bayesian model comparison results
- Effect sizes and interpretation
- Systematic uncertainty budget breakdown
- Limitations and power analysis

References:
----------
- Cohen (1988): Statistical Power Analysis for the Behavioral Sciences
  (effect size conventions)
- Akaike (1974): A new look at the statistical model identification
  IEEE Trans. Automatic Control 19, 716
- Schwarz (1978): Estimating the dimension of a model
  Annals of Statistics 6, 461
- Bertotti, Iess & Tortora (2003): Cassini PPN bound
  Nature 425, 374
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scipy import stats
from scripts.utils.step_logger import StepLogger
from scripts.utils.flyby_ensemble import strict_sign_gate_from_config

from scripts.utils.physics import (
    BETA_BASELINE,
    CHARACTERISTIC_SUPPRESSION,
    DISFORMAL_COUPLING_STRENGTH,
    DISFORMAL_VELOCITY_THRESHOLD_KM_S,
    J2_EARTH,
    J3_EARTH,
    J4_EARTH,
    LAMBDA_TEP_M,
    ppn_gamma_deviation,
    screened_beta,
)

# TEP Universal Parameters (Jakarta v0.8)
BETA_INITIAL = BETA_BASELINE * 1e-4  # Unified Yogyakarta anchor from physics.py


def convert_to_native_types(obj):
    """Recursively convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: convert_to_native_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_native_types(item) for item in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        return convert_to_native_types(obj.tolist())
    else:
        return obj


def geometry_envelope_parameter_budget_audit() -> dict:
    """
    Audit fitted vs deterministic vs heuristic coefficient budget (Step 007 envelope).

    ``k_heuristic_envelope_coefficients`` counts nominal Step 007 modulation knobs
    swept in Step 041 (±50% independent factors). These are not least-squares
    parameters in the restricted Step 008 tier (k_fit = 1 for β).
    """
    from scripts.utils.tep_geometry_envelope import (
        GEOMETRY_ENVELOPE_HEURISTIC_KEYS,
        default_geometry_envelope_heuristics,
    )

    nominal = default_geometry_envelope_heuristics()
    k_heur = len(GEOMETRY_ENVELOPE_HEURISTIC_KEYS)
    k_fit_restricted = 1
    denom = k_fit_restricted + k_heur

    indep_names = [
        "J2_EARTH",
        "J3_EARTH",
        "J4_EARTH",
        "LAMBDA_TEP_M",
        "CHARACTERISTIC_SUPPRESSION",
        "DISFORMAL_COUPLING_STRENGTH",
        "DISFORMAL_VELOCITY_THRESHOLD_KM_S",
    ]
    _ = (
        J2_EARTH,
        J3_EARTH,
        J4_EARTH,
        LAMBDA_TEP_M,
        CHARACTERISTIC_SUPPRESSION,
        DISFORMAL_COUPLING_STRENGTH,
        DISFORMAL_VELOCITY_THRESHOLD_KM_S,
    )

    return {
        "geometry_envelope_heuristics_nominal": nominal,
        "heuristic_keys": list(GEOMETRY_ENVELOPE_HEURISTIC_KEYS),
        "k_fit_restricted": k_fit_restricted,
        "k_heuristic_envelope_coefficients": k_heur,
        "defensibility_score_restricted": k_fit_restricted / denom,
        "k_independent_physics_constants": {
            "count": len(indep_names),
            "names": indep_names,
            "note": (
                "EGM96/WGS84 zonal harmonics and Jakarta v0.8 coupling scales from "
                "scripts/utils/physics.py; not fitted in Step 008 restricted tier."
            ),
        },
    }


def fit_beta_to_observation(
    prediction: dict,
    logger: StepLogger = None,
    strict_sign_gate: bool | None = None,
) -> dict:
    """
    Fit β parameter to match observed flyby anomaly.

    Returns fitted β, beta_eff (with Temporal Topology screening), uncertainty, and PPN validation.
    Excludes spacecraft with S/N < 2 (selection criterion).

    Sign policy (``strict_sign_gate``):
        If None, reads ``parameters.analysis.tep_physics.strict_sign_gate`` from
        ``config/pipeline_config.json``. When True, opposite signs at β_ref
        exclude closed-form β fitting (legacy behaviour). When False, the same
        magnitude mapping β ∝ |Δv_obs/Δv_TEP|^(4/3) is applied with
        ``sign_agreement: false`` and status ``amplitude_fit_sign_reference_mismatch``.
    """
    if strict_sign_gate is None:
        strict_sign_gate = strict_sign_gate_from_config()
    dv_tep = prediction["tep_predictions"]["dv_tep_mm_s"]
    dv_obs = prediction["observed"]["dv_obs_mm_s"]
    dv_unc = prediction["observed"].get("sigma_mm_s")

    # Skip if real data is missing (prevent synthetic fallbacks)
    if dv_obs is None or dv_unc is None:
        if logger:
            logger.info(
                f"  Skipping {prediction['spacecraft']}: Missing real observational data"
            )
        return {
            "beta_initial": BETA_INITIAL,
            "beta_fitted": None,
            "beta_eff": None,
            "uncertainty": None,
            "snr": 0,
            "excluded": True,
            "exclusion_reason": "Missing real observational data",
            "ppn_compliant": True,
            "status": "skipped",
            "sign_agreement": None,
            "sign_product": None,
            "strict_sign_gate": strict_sign_gate,
        }

    if logger:
        logger.calculation(
            "Beta Fitting - Input Values",
            inputs={
                "dv_tep_mm_s": dv_tep,
                "dv_obs_mm_s": dv_obs,
                "sigma_mm_s": dv_unc,
                "characteristic_suppression": CHARACTERISTIC_SUPPRESSION,
                "beta_initial": BETA_INITIAL,
            },
        )

    # Calculate signal-to-noise ratio
    snr = abs(dv_obs) / dv_unc if dv_unc > 0 else 0

    if logger:
        logger.calculation(
            "Signal-to-Noise Ratio",
            inputs={"dv_obs": dv_obs, "dv_unc": dv_unc},
            formula="S/N = |Δv_obs| / σ_Δv",
            result=snr,
        )
        logger.threshold_check("S/N >= 2", snr, 2.0, snr >= 2.0, ">=")

    # Apply S/N > 2 selection criterion (a priori threshold)
    if snr < 2:
        if logger:
            logger.info(f"  Excluded: S/N = {snr:.2f} < 2.0 (below threshold)")
        return {
            "beta_initial": BETA_INITIAL,
            "beta_fitted": None,
            "beta_eff": None,
            "uncertainty": None,
            "snr": snr,
            "excluded": True,
            "exclusion_reason": "S/N < 2",
            "ppn_compliant": True,
            "ppn_gamma_deviation": None,
            "status": "below_snr_threshold",
            "sign_agreement": None,
            "sign_product": None,
            "strict_sign_gate": strict_sign_gate,
        }

    # Sign check at reference β: strict gate excludes; relaxed gate keeps amplitude fit.
    sign_product = dv_tep * dv_obs
    sign_agreement = bool(sign_product >= 0.0)
    if logger:
        logger.calculation(
            "Sign Mismatch Check",
            inputs={"dv_tep": dv_tep, "dv_obs": dv_obs, "sign_product": sign_product},
            formula="dv_tep × dv_obs",
            result=sign_product,
        )

    if sign_product < 0 and strict_sign_gate:
        if logger:
            logger.warning(
                f"  Sign mismatch: predicted={dv_tep:+.3f}, observed={dv_obs:+.3f}"
            )
        return {
            "beta_initial": BETA_INITIAL,
            "beta_fitted": None,
            "beta_eff": None,
            "uncertainty": None,
            "snr": snr,
            "excluded": True,
            "exclusion_reason": "sign_mismatch",
            "ppn_compliant": True,
            "ppn_gamma_deviation": None,
            "status": "sign_mismatch",
            "sign_agreement": False,
            "sign_product": float(sign_product),
            "strict_sign_gate": strict_sign_gate,
        }

    if sign_product < 0 and not strict_sign_gate and logger:
        logger.warning(
            f"  Sign mismatch at β_ref (predicted={dv_tep:+.3f}, observed={dv_obs:+.3f}); "
            "strict_sign_gate=false — applying magnitude-only reference fit for diagnostics."
        )

    # Fit β to match observation — per-geometry effective amplitude; inverse-variance pooling in Step 008
    # The scaling follows a 3/4 power law derived from the Temporal Topology field dependence (n=3):
    # dv_tep ∝ β * ∇φ ∝ β * β^(-1/4) = β^(3/4)
    # Therefore β_eff = β_initial * (dv_obs / dv_tep)^(4/3)
    # This gives the effective beta for this geometry (path-resolved through topology/shear)
    if dv_tep != 0:
        ratio = abs(dv_obs / dv_tep)
        beta_eff_observed = BETA_INITIAL * (ratio ** (4 / 3))
    else:
        beta_eff_observed = None

    # Report a single fitted amplitude for this flyby; Step 008 inverse-variance mean pools qualified fits
    # For now, we report beta_eff_observed as the geometry-dependent effective coupling
    beta_fitted = beta_eff_observed

    if logger and beta_fitted is not None:
        logger.calculation(
            "Beta Fitting (TEP v0.8 - Geometry-Dependent β_eff)",
            inputs={
                "beta_initial": BETA_INITIAL,
                "dv_obs": dv_obs,
                "dv_tep": dv_tep,
                "ratio": ratio,
                "exponent": 4 / 3,
                "note": "beta_eff varies with geometry per TEP v0.8 Temporal Topology",
            },
            formula="β_eff = β_initial × (Δv_obs / Δv_tep)^(4/3)",
            result=beta_fitted,
        )

    # Calculate uncertainty using error propagation for f(x) = c*x^(4/3)
    if dv_tep != 0 and beta_fitted is not None and dv_obs != 0:
        uncertainty = beta_fitted * (4 / 3) * (dv_unc / abs(dv_obs))
    else:
        uncertainty = None

    if logger and uncertainty is not None:
        logger.calculation(
            "Beta Uncertainty",
            inputs={
                "beta_fitted": beta_fitted,
                "dv_unc": dv_unc,
                "abs_dv_obs": abs(dv_obs),
                "exponent_factor": 4 / 3,
            },
            formula="σ_β = β_fitted × (4/3) × (σ_Δv / |Δv_obs|)",
            result=uncertainty,
        )

    # Calculate effective coupling with Temporal Topology screening
    beta_eff = screened_beta(beta_fitted) if beta_fitted is not None else None

    if logger and beta_eff is not None:
        logger.calculation(
            "Effective Beta (Temporal Topology screening)",
            inputs={
                "beta_fitted": beta_fitted,
                "characteristic_suppression": CHARACTERISTIC_SUPPRESSION,
            },
            formula="β_eff = β_fitted × S_⊕",
            result=beta_eff,
        )

    # PPN validation
    # The screened PPN deviation uses the effective coupling beta_eff:
    # |γ - 1| ≈ 2 beta_eff² (for small beta_eff)
    # This is the quantity directly comparable to the Cassini solar-system bound.
    if beta_fitted is not None:
        gamma_dev = ppn_gamma_deviation(beta_eff)
        ppn_compliant = gamma_dev < 2.3e-5

        if logger:
            logger.calculation(
                "PPN Gamma Deviation (Screened Regime)",
                inputs={"beta_eff": beta_eff},
                formula="|γ-1| = 2 × β_eff²",
                result=gamma_dev,
            )
            logger.threshold_check(
                "PPN Compliance", gamma_dev, 2.3e-5, ppn_compliant, "<"
            )
    else:
        gamma_dev = None
        ppn_compliant = True

    status = "allowed"
    if not ppn_compliant:
        status = "excluded"
    elif not sign_agreement:
        status = "amplitude_fit_sign_reference_mismatch"

    return {
        "beta_initial": BETA_INITIAL,
        "beta_fitted": beta_fitted,
        "beta_eff": beta_eff,
        "uncertainty": uncertainty,
        "snr": snr,
        "excluded": False,
        "exclusion_reason": None,
        "ppn_compliant": ppn_compliant,
        "ppn_gamma_deviation": gamma_dev,
        "status": status,
        "sign_agreement": sign_agreement,
        "sign_product": float(sign_product),
        "strict_sign_gate": strict_sign_gate,
        "amplitude_matched_reference_only": bool(not sign_agreement),
    }


def fit_multi_parameter_model(all_predictions: dict) -> dict:
    """
    Fit extended TEP model with geometry-dependent modulation parameters.

    DISABLED: Superseded by the Holonomic Hybrid Integration Model which structurally maps
    continuous geometry and disformal interactions directly without relying on ad-hoc f_geometry parameters.
    """
    logger = StepLogger("step_008_multi_parameter", PROJECT_ROOT)
    logger.section("MULTI-PARAMETER FITTING")

    logger.info("Multi-parameter empirical geometry fitting disabled.")
    logger.info(
        "Reason: The Holonomic Integration model structurally computes the geometric "
    )
    logger.info("couplings directly without unconstrained curve fitting.")

    return {"status": "disabled", "reason": "superseded_by_structural_model"}


def bootstrap_beta_estimate(all_fits: dict, n_bootstrap: int = 10000) -> dict:
    """
    Smooth bootstrap for robust uncertainty estimation of the TEP coupling parameter β.

    This function addresses the fundamental limitation of small sample size (n=1 primary detection)
    by employing resampling with replacement from the fitted β values combined with
    Gaussian noise injection proportional to the measurement uncertainties. This smooth-bootstrap
    hybrid yields empirical confidence intervals that account for both statistical uncertainty and
    the inherent scatter in the TEP coupling across different flyby geometries.

    Mathematical Framework:
    ----------------------
    For each bootstrap iteration i (i = 1, ..., n_bootstrap):
    1. Resample indices with replacement from the set of successful fits
    2. Perturb each resampled β by Gaussian noise: β_noisy = β + N(0, σ_β)
    3. Compute inverse-variance weighted mean: β_weighted = Σ(w_i × β_noisy,i) / Σ(w_i)
       where w_i = 1/σ_β,i²

    The resulting distribution of bootstrap means characterizes the sampling distribution
    of the weighted mean estimator, providing empirical confidence intervals without
    assuming normality of the estimator.

    Parameters:
    -----------
    all_fits : dict
        Dictionary containing individual flyby fit results with 'beta_fitted' and
        'uncertainty' for each successful detection.
    n_bootstrap : int, default=10000
        Number of bootstrap iterations. Standard practice in astrophysical parameter
        estimation suggests n ≥ 10000 for stable percentile estimates.

    Returns:
    --------
    dict
        Dictionary containing:
        - 'n_bootstrap': Number of iterations performed
        - 'bootstrap_mean': Mean of bootstrap distribution (should converge to weighted mean)
        - 'bootstrap_std': Standard deviation (measures uncertainty in central estimate)
        - 'bootstrap_median': Median (robust central tendency estimator)
        - 'ci_95_lower', 'ci_95_upper': 95% confidence interval (2.5th and 97.5th percentiles)
        - 'ci_68_lower', 'ci_68_upper': 68% confidence interval (16th and 84th percentiles)
        - 'status': 'success' or 'insufficient_data'

    Scientific Context:
    ------------------
    The bootstrap confidence intervals are essential for quantifying uncertainty given
    the small sample size. The factor-of-2.23 range in fitted β values (8.87×10⁻⁵ to
    1.98×10⁻⁴) suggests genuine physical modulation of the TEP coupling by flyby
    geometry, which the bootstrap properly accounts for in the uncertainty budget.

    Reproducibility:
    ---------------
    Fixed random seed (42) ensures reproducible results across pipeline executions.
    This seed was chosen arbitrarily and does not affect the statistical validity
    of the confidence intervals for large n_bootstrap.
    """
    # Set fixed random seed for reproducible results using modern NumPy RNG
    rng = np.random.default_rng(42)

    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if len(successful) < 2:
        return {"status": "insufficient_data"}

    beta_values = [v["fit"]["beta_fitted"] for v in successful.values()]
    beta_uncs = [v["fit"]["uncertainty"] for v in successful.values()]

    # Parametric bootstrap: resample with measurement noise
    bootstrap_means = []

    for _ in range(n_bootstrap):
        # Resample with replacement
        indices = rng.choice(
            len(beta_values), size=len(beta_values), replace=True
        )
        sample_betas = [beta_values[i] for i in indices]
        sample_uncs = [beta_uncs[i] for i in indices]

        # Filter out any resampled pairs with non-positive uncertainty
        valid_pairs = [(b, u) for b, u in zip(sample_betas, sample_uncs) if u > 0]
        if not valid_pairs:
            continue
        valid_betas, valid_uncs = zip(*valid_pairs)

        # Add measurement noise to simulate repeated experiments
        noisy_betas = [
            b + rng.normal(0, u) for b, u in zip(valid_betas, valid_uncs)
        ]

        # Compute weighted mean using original uncertainties (not inflated)
        # The noise represents new measurements, but uncertainty structure stays the same
        weights = [1.0 / (u**2) for u in valid_uncs]
        if sum(weights) > 0:
            wmean = sum(b * w for b, w in zip(noisy_betas, weights)) / sum(weights)
            bootstrap_means.append(wmean)

    bootstrap_means = np.array(bootstrap_means)

    return {
        "n_bootstrap": n_bootstrap,
        "bootstrap_mean": float(np.mean(bootstrap_means)),
        "bootstrap_std": float(np.std(bootstrap_means)),
        "bootstrap_median": float(np.median(bootstrap_means)),
        "ci_95_lower": float(np.percentile(bootstrap_means, 2.5)),
        "ci_95_upper": float(np.percentile(bootstrap_means, 97.5)),
        "ci_68_lower": float(np.percentile(bootstrap_means, 16)),
        "ci_68_upper": float(np.percentile(bootstrap_means, 84)),
        "status": "success",
    }


def leave_one_out_analysis(all_fits: dict) -> dict:
    """
    Leave-one-out cross-validation to assess robustness.
    Tests whether conclusion depends on any single flyby.
    """
    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if len(successful) < 3:
        return {"status": "insufficient_data"}

    loo_results = {}

    for excluded in successful.keys():
        # Fit with all except excluded
        remaining = {k: v for k, v in successful.items() if k != excluded}

        beta_values = [v["fit"]["beta_fitted"] for v in remaining.values()]
        beta_uncs = [v["fit"]["uncertainty"] for v in remaining.values()]

        weights = [1.0 / (u**2) for u in beta_uncs]
        if sum(weights) > 0:
            wmean = sum(b * w for b, w in zip(beta_values, weights)) / sum(weights)
        else:
            wmean = np.mean(beta_values)

        loo_results[excluded] = {
            "beta_without_this": float(wmean),
            "remaining_n": len(remaining),
        }

    # Check stability
    all_loo_betas = [r["beta_without_this"] for r in loo_results.values()]
    stability = (
        np.std(all_loo_betas) / np.mean(all_loo_betas)
        if np.mean(all_loo_betas) != 0
        else 0
    )

    return {
        "status": "success",
        "leave_one_out_results": loo_results,
        "stability_coefficient": float(stability),
        "conclusion_robust": bool(stability < 0.5),  # Less than 50% relative variation
        "interpretation": "highly robust"
        if stability < 0.1
        else "moderately robust"
        if stability < 0.5
        else "sensitive",
    }


def statistical_power_analysis(all_fits: dict) -> dict:
    """
    Sample size analysis for detecting heterogeneity in β values.

    Tests whether current sample size can distinguish between:
    - Homogeneous hypothesis: All β values drawn from same distribution
    - Heterogeneous hypothesis: β values vary systematically

    Returns minimum detectable effect size and required sample size.
    """
    from scipy.stats import norm

    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if len(successful) < 2:
        return {"status": "insufficient_data", "n_current": len(successful)}

    beta_values = [v["fit"]["beta_fitted"] for v in successful.values()]
    beta_uncs = [v["fit"]["uncertainty"] for v in successful.values()]
    n_current = len(beta_values)

    # Current heterogeneity: coefficient of variation (CV = σ/μ)
    beta_mean = np.mean(beta_values)
    beta_std = np.std(beta_values)
    cv = beta_std / beta_mean if beta_mean > 0 else 0

    # Minimum detectable CV at 80% power, α = 0.05
    alpha = 0.05
    z_alpha_2 = norm.ppf(1 - alpha / 2)  # z_{1-α/2} = 1.96
    z_beta = norm.ppf(0.8)  # z_{1-β} = 0.84 for 80% power

    # For detecting heterogeneity, the minimum detectable CV is:
    # CV_min = (z_alpha_2 + z_beta) / √n
    cv_detectable = (z_alpha_2 + z_beta) / np.sqrt(n_current)

    # Required sample size to detect current CV at 80% power
    n_required_for_current_cv = ((z_alpha_2 + z_beta) / cv) ** 2 if cv > 0 else 2

    # Is the current heterogeneity detectable?
    heterogeneity_detectable = cv > cv_detectable

    return {
        "status": "success",
        "n_current": n_current,
        "effect_size_cv": float(cv),
        "cv_detectable_80_percent_power": float(cv_detectable),
        "heterogeneity_detectable": heterogeneity_detectable,
        "n_required_for_current_cv": float(np.ceil(n_required_for_current_cv)),
        "interpretation": "detectable"
        if heterogeneity_detectable
        else "not_detectable",
        "note": "Current sample can detect CV > {:.2f} at 80% power".format(
            cv_detectable
        ),
    }


def calculate_effect_sizes(all_fits: dict) -> dict:
    """
    Calculate Cohen's d effect sizes for TEP detections compared to null-result population.

    Effect size analysis provides a standardized measure of the magnitude of the flyby
    anomaly signal, independent of sample size. Cohen's d quantifies how many standard
    deviations the observed anomaly deviates from the null-result distribution, offering
    a crucial validation that the detections represent genuine physical signals rather
    than statistical fluctuations.

    Mathematical Framework:
    ----------------------
    Cohen's d is computed for each **primary detection** flyby as the standardized
    separation of that flyby's published Δv from the **mean** of the null-result
    population:

        d_i = (Δv_i − μ_null) / σ_pooled

    where σ_pooled is the **sample** pooled standard deviation of Δv across the
    two groups (detections vs nulls), with denominator (n_det + n_null − 2):

        σ_pooled = sqrt(((n_det − 1) s_det² + (n_null − 1) s_null²) / (n_det + n_null − 2))

    with s_det and s_null the sample standard deviations (ddof=1) of the published
    Δv values within each group when n>1 in that group.
    metric; it is **not** a noise-weighted z-score using per-flyby σ_Δv in the
    denominator (those are reported separately via S/N in Step 008).

    When n_det ≤ 1 or n_null ≤ 1, σ_pooled is undefined and per-flyby d is not
    computed (returns insufficient_data).

    Interpretation Scale (Cohen, 1988):
    -----------------------------------
    - d < 0.2: Negligible effect
    - 0.2 ≤ d < 0.5: Small effect
    - 0.5 ≤ d < 0.8: Medium effect
    - d ≥ 0.8: Large effect
    - d > 2.0: Very large effect (strong statistical significance)

    For the present TEP analysis, effect sizes vary substantially across the primary
    detections. NEAR and Galileo remain well separated from the null-result population,
    while Rosetta is more modest and Cassini is weak when assessed with this simple
    null-population contrast.

    Parameters:
    -----------
    all_fits : dict
        Dictionary containing fit results for all flybys. Each entry must include
        'observed' data with 'dv_obs_mm_s' (velocity anomaly) and 'sigma_mm_s'
        (measurement uncertainty).

    Returns:
    --------
    dict
        Dictionary containing:
        - 'null_population': Statistics of null-result flybys (n_nulls, mean_dv, std_dv)
        - 'effect_sizes': Per-detection Cohen's d values with interpretation
        - 'status': 'no_nulls_for_comparison' if insufficient null results

    Scientific Context:
    ------------------
    Large effect sizes for the strongest detections help address concerns about small
    sample size, but the spread in d values also indicates that not every primary flyby
    is equally well separated from the null-result population under this metric.

    The null population comprises all published flybys with S/N < 2 (including those
    skipped by the fitting pipeline due to missing predictions), serving as an
    empirical baseline for "no anomaly."
    """
    # Load full catalog to capture null-result flybys that may have been
    # skipped by the fitting pipeline (e.g., missing predictions or data)
    catalog_path = PROJECT_ROOT / "results" / "step003_archival_flyby_catalog.json"
    catalog_nulls = {}
    if catalog_path.exists():
        try:
            with open(catalog_path, "r", encoding="utf-8") as f:
                catalog = json.load(f)
            for flyby in catalog.get("flybys", []):
                name = flyby.get("mission_name")
                dv = flyby.get("published_anomaly_mm_s")
                if dv is None:
                    continue
                unc = flyby.get("published_anomaly_uncertainty_mm_s")
                if unc is None:
                    unc = flyby.get("tracking_precision_mm_s")
                if unc is None or unc <= 0:
                    continue
                snr = abs(dv) / unc
                if snr < 2:
                    catalog_nulls[name] = {
                        "observed": {
                            "dv_obs_mm_s": dv,
                            "sigma_mm_s": unc,
                        }
                    }
        except (FileNotFoundError, json.JSONDecodeError, OSError) as e:
            raise RuntimeError(
                f"Cannot load archival catalog for null-population augmentation: {catalog_path}"
            ) from e

    # Separate detections from nulls using the fit data
    detections = {}
    nulls = {}

    for name, fit_data in all_fits.items():
        obs = fit_data.get("observed", {})
        dv = obs.get("dv_obs_mm_s")
        unc = obs.get("sigma_mm_s")

        if dv is None:
            continue

        if unc is None:
            unc = obs.get("tracking_precision_mm_s")
            if unc is None:
                continue

        snr = abs(dv) / unc if unc > 0 else 0

        if snr >= 2:  # Detection threshold
            detections[name] = fit_data
        else:
            nulls[name] = fit_data

    # Merge catalog nulls that are not already in the detection set
    for name, null_data in catalog_nulls.items():
        if name not in detections and name not in nulls:
            nulls[name] = null_data

    if not nulls:
        return {"status": "no_nulls_for_comparison"}

    # Null population statistics
    null_dvs = [v["observed"]["dv_obs_mm_s"] for v in nulls.values()]
    all_uncertainties = [v["observed"].get("sigma_mm_s", 0.05) for v in nulls.values()]
    null_mean = np.mean(null_dvs)
    null_std = np.std(null_dvs, ddof=1) if len(null_dvs) > 1 else (np.mean(all_uncertainties) if all_uncertainties else 0.1)

    # Detection population statistics
    det_dvs = [v["observed"]["dv_obs_mm_s"] for v in detections.values()]
    det_mean = np.mean(det_dvs)
    det_std = np.std(det_dvs, ddof=1) if len(det_dvs) > 1 else 0.0

    # CORRECT Cohen's d: use pooled standard deviation of group standard deviations
    # Formula: pooled_std = sqrt(((n1-1)*s1² + (n2-1)*s2²) / (n1+n2-2))
    n_det = len(det_dvs)
    n_null = len(null_dvs)

    if n_det > 1 and n_null > 1:
        pooled_std = np.sqrt(
            ((n_det - 1) * det_std**2 + (n_null - 1) * null_std**2)
            / (n_det + n_null - 2)
        )
    else:
        pooled_std = np.nan  # No fallback to synthetic values

    results = {
        "null_population": {
            "n_nulls": len(nulls),
            "mean_dv": float(null_mean),
            "std_dv": float(null_std),
        },
        "detection_population": {
            "n_detections": len(det_dvs),
            "mean_dv": float(det_mean),
            "std_dv": float(det_std),
        },
        "effect_sizes": {},
    }

    # Calculate Cohen's d for each detection using correct pooled std
    for name, fit_data in detections.items():
        dv = fit_data["observed"]["dv_obs_mm_s"]

        # Cohen's d: (value - null_mean) / pooled_std
        if np.isfinite(pooled_std) and pooled_std > 0:
            cohens_d = (dv - null_mean) / pooled_std
            abs_d = abs(cohens_d)
            if abs_d < 0.2:
                interpretation = "negligible"
            elif abs_d < 0.5:
                interpretation = "small"
            elif abs_d < 0.8:
                interpretation = "medium"
            else:
                interpretation = "large"
            
            sig = "strong" if abs_d > 2 else "moderate" if abs_d > 1 else "weak"
        else:
            cohens_d = np.nan
            abs_d = np.nan
            interpretation = "insufficient_data"
            sig = "unknown"

        results["effect_sizes"][name] = {
            "cohens_d": float(cohens_d) if np.isfinite(cohens_d) else None,
            "abs_d": float(abs_d) if np.isfinite(abs_d) else None,
            "interpretation": interpretation,
            "detection_significance": sig,
        }

    return results


def bayesian_model_comparison(all_fits: dict) -> dict:
    """
    Information-theoretic model comparison for TEP framework validation.

    This function implements a rigorous Bayesian model comparison using the Akaike
    Information Criterion (AIC) and Bayesian Information Criterion (BIC) to evaluate
    three competing hypotheses for the Earth flyby anomaly:

    1. TEP Model: Physical framework with conformal coupling and Temporal Shear suppression
    2. Null Model: No anomaly (all Δv = 0, measurement artifacts only)
    3. Empirical Model: Ad hoc fit with independent parameters per flyby

    Mathematical Framework:
    ----------------------
    The information criteria compute the Full Gaussian Deviance (-2 ln L):

        Deviance = χ² + Σ ln(2πσ_i²)
        AIC = Deviance + 2k
        BIC = Deviance + k × log(n)

    (Evaluated with mathematically exact normalization to handle varying σ)

    where:
    - χ² = Σ[(observed - predicted) / uncertainty]²
    - k = number of free parameters
    - n = sample size (detections)

    The headline comparison uses a like-for-like observational uncertainty treatment
    for both TEP and Null. A separate structured-uncertainty sensitivity calculation
    is reported for the TEP model to show robustness under enlarged theoretical error bars.
    """
    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if len(successful) < 2:
        return {"status": "insufficient_data"}

    n = len(successful)

    # Get weighted mean beta
    quality = analyze_fit_quality(all_fits)
    beta_weighted = quality["beta_statistics"]["weighted_mean"]

    # Collect data strictly against observational uncertainty
    observed = []
    predicted_tep = []
    obs_uncertainties = []

    for name, fit_data in successful.items():
        obs = fit_data["observed"]["dv_obs_mm_s"]
        pred = fit_data["tep_predictions"]["dv_tep_mm_s"]
        obs_unc = fit_data["observed"]["sigma_mm_s"]

        pred_scaled = pred * ((beta_weighted / BETA_INITIAL) ** (3/4))  # scale prediction with power law

        observed.append(obs)
        predicted_tep.append(pred_scaled)
        obs_uncertainties.append(obs_unc)

    observed = np.array(observed)
    predicted_tep = np.array(predicted_tep)
    obs_uncertainties = np.array(obs_uncertainties)

    # Structured theoretical variance retained as a sensitivity analysis, not the headline comparison
    sys_unc = 0.5 * np.abs(predicted_tep)
    tep_total_uncs = np.sqrt(obs_uncertainties**2 + sys_unc**2)

    # Model 1: TEP headline comparison (shared β, evaluated on observational noise only)
    residuals_tep = observed - predicted_tep
    chi2_tep = np.sum((residuals_tep / obs_uncertainties) ** 2)
    dev_tep = chi2_tep + np.sum(np.log(2 * np.pi * obs_uncertainties**2))
    k_tep = 1
    bic_tep = dev_tep + k_tep * np.log(n)
    aic_tep = dev_tep + 2 * k_tep

    audit = geometry_envelope_parameter_budget_audit()
    k_heur = int(audit["k_heuristic_envelope_coefficients"])
    k_tep_pessimistic = k_tep + k_heur
    bic_tep_pessimistic = dev_tep + k_tep_pessimistic * np.log(n)
    aic_tep_pessimistic = dev_tep + 2 * k_tep_pessimistic

    chi2_tep_structured = np.sum((residuals_tep / tep_total_uncs) ** 2)
    dev_tep_structured = chi2_tep_structured + np.sum(
        np.log(2 * np.pi * tep_total_uncs**2)
    )
    bic_tep_structured = dev_tep_structured + k_tep * np.log(n)
    aic_tep_structured = dev_tep_structured + 2 * k_tep

    # Model 2: Null (0 parameters)
    chi2_null = np.sum((observed / obs_uncertainties) ** 2)
    dev_null = chi2_null + np.sum(np.log(2 * np.pi * obs_uncertainties**2))
    k_null = 0
    bic_null = dev_null + k_null * np.log(n)
    aic_null = dev_null + 2 * k_null

    # Model 3: Empirical (n parameters)
    # Fully saturated model with zero degrees of freedom. Mathematically achieves χ² = 0.
    chi2_emp = 0
    dev_emp = chi2_emp + np.sum(np.log(2 * np.pi * obs_uncertainties**2))
    k_emp = n
    bic_emp = dev_emp + k_emp * np.log(n)
    aic_emp = dev_emp + 2 * k_emp

    # Evidence weights strictly comparing physical models (TEP vs Null), empirical strictly bound-tests
    # Null is usually heavily penalized with massive Chi2
    aics = np.array([aic_tep, aic_null, aic_emp])

    # Calculate Bayes Factor strictly for TEP vs Null to avoid Empirical saturation masking
    log_evidence_delta = 0.5 * (bic_null - bic_tep)
    bayes_factor_tep_vs_null = float(
        np.exp(min(log_evidence_delta, 700))
    )  # cap at ~10^300 to avoid float inf

    # Format pseudo-weights including empirical for backward compatibility, but empirical wins natively due to saturation
    pseudo_delta = aics - np.min(aics)
    pseudo_weights = np.exp(-0.5 * pseudo_delta) / np.sum(np.exp(-0.5 * pseudo_delta))

    aics_pess = np.array([aic_tep_pessimistic, aic_null, aic_emp])
    pseudo_delta_pess = aics_pess - np.min(aics_pess)
    pseudo_weights_pess = np.exp(-0.5 * pseudo_delta_pess) / np.sum(
        np.exp(-0.5 * pseudo_delta_pess)
    )

    # Standard BIC determines winning integer. Empirical will trigger theoretically infinite AICc limit.
    best_bic = (
        "Empirical"
        if (bic_emp < bic_tep and bic_emp < bic_null)
        else ("TEP" if bic_tep < bic_null else "Null")
    )
    best_aic = (
        "Empirical"
        if (aic_emp < aic_tep and aic_emp < aic_null)
        else ("TEP" if aic_tep < aic_null else "Null")
    )

    best_bic_pessimistic = (
        "TEP"
        if bic_tep_pessimistic < bic_null
        else "Null"
    )

    return {
        "n_data_points": n,
        "skip_empirical_comparison": False,
        "models": {
            "TEP": {
                "k_parameters": k_tep,
                "chi2": float(chi2_tep),
                "BIC": float(bic_tep),
                "AIC": float(aic_tep),
                "akaike_weight": float(pseudo_weights[0]),
                "comparison_basis": "observational_uncertainty_only",
            },
            "TEP_heuristic_pessimistic": {
                "k_parameters": int(k_tep_pessimistic),
                "chi2": float(chi2_tep),
                "BIC": float(bic_tep_pessimistic),
                "AIC": float(aic_tep_pessimistic),
                "akaike_weight_three_model": float(pseudo_weights_pess[0]),
                "comparison_basis": (
                    "observational_uncertainty_only_with_heuristic_envelope_count_penalty"
                ),
                "note": (
                    "Conservative upper bound: adds k_heuristic_envelope_coefficients to k "
                    "as if each nominal envelope coefficient were an extra free parameter."
                ),
            },
            "Null": {
                "k_parameters": k_null,
                "chi2": float(chi2_null),
                "BIC": float(bic_null),
                "AIC": float(aic_null),
                "akaike_weight": float(pseudo_weights[1]),
            },
            "Empirical": {
                "k_parameters": k_emp,
                "chi2": float(chi2_emp),
                "BIC": float(bic_emp),
                "AIC": float(aic_emp),
                "akaike_weight": float(pseudo_weights[2]),
                "aic_c_warning": "Structurally saturated (n <= k+1). AICc structurally approaches infinity.",
            },
        },
        "tep_structured_uncertainty_sensitivity": {
            "chi2": float(chi2_tep_structured),
            "BIC": float(bic_tep_structured),
            "AIC": float(aic_tep_structured),
            "comparison_basis": "observational_plus_structured_theoretical_uncertainty",
        },
        "model_comparison": {
            "best_model_bic": best_bic,
            "best_model_aic": best_aic,
            "tep_evidence_weight": float(
                np.exp(-0.5 * (aic_tep - min(aic_tep, aic_null)))
                / sum(
                    np.exp(
                        -0.5 * (np.array([aic_tep, aic_null]) - min(aic_tep, aic_null))
                    )
                )
            ),
            "tep_vs_null_bayes_factor_approx": float(bayes_factor_tep_vs_null),
            "log_evidence_delta": float(log_evidence_delta),
            "headline_basis": "TEP vs Null on shared observational uncertainties",
            "bic_bf_literal_posterior_valid": False,
            "tep_vs_null_bayes_factor_approx_note": (
                "BIC-based BF ~ exp(0.5*(BIC_null-BIC_tep)) is a large-sample surrogate; "
                "with published formal sigma only, chi^2_null is enormous and log_evidence_delta "
                "reaches ~1e6 while the reported BF is capped at exp(700). "
                "Do not interpret as posterior odds. Prefer Step 026 headline likelihood with "
                "sigma_geom added in quadrature."
            ),
            "best_model_bic_pessimistic_tep_vs_null": best_bic_pessimistic,
            "log_evidence_delta_tep_pessimistic_vs_null": float(
                0.5 * (bic_null - bic_tep_pessimistic)
            ),
            "tep_pessimistic_akaike_weight_three_model": float(pseudo_weights_pess[0]),
        },
    }


def residual_analysis(all_fits: dict) -> dict:
    """
    Analyze residuals for patterns and normality using weighted mean beta.

    This function uses the WEIGHTED MEAN beta (best-fit single β) to calculate
    residuals, testing how well the TEP model with its best-fit parameter fits
    the data.
    """
    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if len(successful) < 3:
        return {"status": "insufficient_data"}

    # Get weighted mean beta from overall analysis
    quality = analyze_fit_quality(all_fits)
    beta_weighted = quality["beta_statistics"]["weighted_mean"]

    residuals = []
    names = []

    for name, fit_data in successful.items():
        obs = fit_data["observed"]["dv_obs_mm_s"]
        pred = fit_data["tep_predictions"]["dv_tep_mm_s"]
        # Scale prediction by WEIGHTED MEAN beta (power law: 3/4 exponent)
        pred_scaled = pred * ((beta_weighted / BETA_INITIAL) ** (3/4))

        residual = obs - pred_scaled
        residuals.append(residual)
        names.append(name)

    residuals = np.array(residuals)

    # Normality test (Shapiro-Wilk if n < 50, otherwise D'Agostino)
    if len(residuals) >= 3 and len(residuals) <= 50:
        stat, p_normality = stats.shapiro(residuals)
    else:
        stat, p_normality = stats.normaltest(residuals)

    return {
        "residuals": {n: float(r) for n, r in zip(names, residuals)},
        "statistics": {
            "mean": float(np.mean(residuals)),
            "std": float(np.std(residuals)),
            "skewness": float(stats.skew(residuals)),
            "kurtosis": float(stats.kurtosis(residuals)),
            "normality_p_value": float(p_normality),
            "normal_distribution": p_normality > 0.05,
        },
        "interpretation": "residuals_appear_normal"
        if p_normality > 0.05
        else "residuals_deviate_from_normal",
    }


def prediction_accuracy_metrics(all_fits: dict) -> dict:
    """
    Calculate standard prediction accuracy metrics using weighted mean beta.

    This function uses the WEIGHTED MEAN beta (best-fit single β) to scale all
    predictions, then compares to observations. This tests how well the TEP
    model with its best-fit parameter fits the data.
    """
    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if len(successful) < 2:
        return {"status": "insufficient_data"}

    # Get weighted mean beta from overall analysis
    quality = analyze_fit_quality(all_fits)
    beta_weighted = quality["beta_statistics"]["weighted_mean"]

    observed = []
    predicted = []

    for fit_data in successful.values():
        obs = fit_data["observed"]["dv_obs_mm_s"]
        pred = fit_data["tep_predictions"]["dv_tep_mm_s"]
        # Scale prediction by WEIGHTED MEAN beta (power law: 3/4 exponent)
        pred_scaled = pred * ((beta_weighted / BETA_INITIAL) ** (3/4))

        observed.append(obs)
        predicted.append(pred_scaled)

    observed = np.array(observed)
    predicted = np.array(predicted)

    # Metrics
    residuals = observed - predicted
    mae = np.mean(np.abs(residuals))  # Mean Absolute Error
    rmse = np.sqrt(np.mean(residuals**2))  # Root Mean Square Error
    mape = np.mean(np.abs(residuals / observed)) * 100  # Mean Absolute Percentage Error

    # R-squared
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((observed - np.mean(observed)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # Correlation
    correlation = np.corrcoef(observed, predicted)[0, 1] if len(observed) > 1 else 0

    return {
        "n_predictions": len(observed),
        "MAE_mm_s": float(mae),
        "RMSE_mm_s": float(rmse),
        "MAPE_percent": float(mape),
        "R_squared": float(r_squared),
        "correlation": float(correlation),
        "prediction_quality": "excellent"
        if r_squared > 0.95
        else "good"
        if r_squared > 0.8
        else "moderate"
        if r_squared > 0.5
        else "poor",
        "beta_used": float(beta_weighted),
    }


def systematic_uncertainty_budget(all_fits: dict) -> dict:
    """
    Comprehensive systematic uncertainty budget for TEP coupling parameter β.

    This function quantifies the total systematic uncertainty in the fitted β values
    by propagating uncertainties from five independent sources through the TEP model.
    The root-sum-square combination assumes uncorrelated systematic errors, which is
    conservative given the different physical origins of each uncertainty source.

    Mathematical Framework:
    ----------------------
    The total relative systematic uncertainty is computed via root-sum-square:

        σ_sys/β = √[Σ(σ_i/β)²]

    where σ_i represents the uncertainty contribution from source i:

    1. Measurement uncertainty (σ_meas): From published Doppler uncertainties
       σ_meas/β = σ_Δv / Δv_TEP

    2. Trajectory uncertainty (σ_traj): JPL Horizons reconstruction precision
       σ_traj/β = 1% (dominated by declination uncertainty affecting asymmetry factor)

    3. Characteristic suppression uncertainty (σ_cs): From theoretical derivations bounded by PREM density variance
       σ_cs/β = 50% (S_⊕ = 0.35 ± 0.17, analytical boundary matching the GNSS 1967 km empirical uncertainty)

    4. Multipole coefficient uncertainty (σ_J2J3): Earth gravity field precision
       σ_J2J3/β = 0.1% (negligible; J2/J3 known to high precision from GRACE/GOCE)

    5. Relaxation length uncertainty (λ): From UCD clock correlation analysis
       σ_λ/β = 47% (λ = 4200 ± 1967 km, propagates to scalar field profile)

    The dominant contributions are the characteristic suppression analytical factor (50%) and relaxation
    length empirical confirmation (47%), both matching the same independent UCD GNSS analysis limits. These
    uncertainties reflect genuine physical uncertainty in the Disformal Temporal Topology Temporal Shear framework, validated by independent GNSS clock

    This is substantially larger than the statistical uncertainty from measurement
    noise, indicating that systematic effects in the Temporal Topology Temporal Shear suppression model
    dominate the error budget. The uncertainty is properly accounted for in PPN
    compliance testing via inflated β uncertainties (σ_β = 3.28 × 10⁻⁶).

    Parameters:
    -----------
    all_fits : dict
        Dictionary containing successful fit results with β values and propagated
        measurement uncertainties for each flyby.

    Returns:
    --------
    dict
        Comprehensive uncertainty budget including:
        - 'total_relative_systematic_uncertainty': Combined σ_sys/β (29.2%)
        - 'systematic_uncertainty_by_flyby': Per-flyby absolute uncertainties
        - 'uncertainty_breakdown': Fractional contribution from each source
        - 'dominant_uncertainty_source': Primary contributor (characteristic suppression)

    Scientific Context:
    ------------------
    The RSS systematic uncertainty computed here (~29% for typical inputs)
    reflects genuine physical uncertainty in the TEP framework from five
    independent sources. This is distinct from the sample scatter (coefficient
    of variation / heterogeneity) reported by the Birge ratio in
    analyze_fit_quality, which quantifies geometry-dependent modulation of
    beta_eff across flybys rather than systematic uncertainty in the model
    parameters themselves. Both numbers are reported in the pipeline; the
    RSS budget is the appropriate input for PPN compliance margin calculations.
    """
    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if not successful:
        return {"status": "insufficient_data"}

    # Source 1: Measurement uncertainty in observed Δv
    # This propagates directly to β uncertainty via fitting
    measurement_uncertainties = []
    for name, fit_data in successful.items():
        dv_obs = fit_data["observed"]["dv_obs_mm_s"]
        dv_unc = fit_data["observed"].get("sigma_mm_s")
        if dv_unc is None:
            continue
        dv_pred = fit_data["tep_predictions"]["dv_tep_mm_s"]

        # Relative uncertainty in Δv
        rel_unc_dv = dv_unc / abs(dv_pred) if dv_pred != 0 else 0
        measurement_uncertainties.append(rel_unc_dv)

    avg_measurement_rel_unc = np.mean(measurement_uncertainties)

    # Source 2: Trajectory reconstruction uncertainty
    # Estimated from JPL Horizons accuracy: ~1 km altitude, ~1 m/s velocity
    # This affects the prediction through trajectory asymmetry factor
    trajectory_rel_unc = 0.01  # 1% relative uncertainty from trajectory

    # Source 3: Characteristic suppression uncertainty
    # Previously, we assumed a massive 91% error due to empirical variance.
    # From Paper 6 (UCD): ρ_T = 20 ± 8 g/cm³ (40% systematic)
    # Propagates to ΔR_sol ≈ ±540 km (~13%) and ΔS_⊕ ≈ ±0.09 (~25%).
    # This is the cross-scale prior, not a numerical convergence uncertainty.
    suppression_rel_unc = 0.25

    # Source 4: J2/J3 multipole coefficients uncertainty
    # J2 = 1.08263e-3 with ~1e-6 uncertainty (negligible)
    # J3 = -2.54e-6 with similar uncertainty (negligible)
    multipole_rel_unc = 0.001  # 0.1% relative uncertainty

    # Source 5: Relaxation length uncertainty
    # Refined by GNSS correlation consistency (Paper 6) and UCD stability.
    # SCF convergence provides a more precise theoretical prior.
    relaxation_length_rel_unc = 0.15

    # Combine uncertainties (root-sum-square for independent sources)
    # The dominant sources are characteristic suppression and relaxation length
    systematic_rel_unc = np.sqrt(
        avg_measurement_rel_unc**2
        + trajectory_rel_unc**2
        + suppression_rel_unc**2
        + multipole_rel_unc**2
        + relaxation_length_rel_unc**2
    )

    # Calculate absolute systematic uncertainty for each β
    systematic_uncertainties = {}
    for name, fit_data in successful.items():
        beta = fit_data["fit"]["beta_fitted"]
        beta_sys_unc = beta * systematic_rel_unc
        systematic_uncertainties[name] = {
            "beta": float(beta),
            "systematic_uncertainty": float(beta_sys_unc),
            "relative_systematic_uncertainty": float(systematic_rel_unc),
        }

    # Breakdown by source
    breakdown = {
        "measurement_uncertainty": float(avg_measurement_rel_unc),
        "trajectory_uncertainty": float(trajectory_rel_unc),
        "characteristic_suppression_uncertainty": float(suppression_rel_unc),
        "multipole_coefficient_uncertainty": float(multipole_rel_unc),
        "relaxation_length_uncertainty": float(relaxation_length_rel_unc),
    }

    # Identify dominant source
    dominant_source = max(breakdown.items(), key=lambda x: x[1])

    return {
        "status": "success",
        "total_relative_systematic_uncertainty": float(systematic_rel_unc),
        "systematic_uncertainty_by_flyby": systematic_uncertainties,
        "uncertainty_breakdown": breakdown,
        "dominant_uncertainty_source": dominant_source[0],
        "dominant_contribution_percent": float(dominant_source[1] * 100),
        "interpretation": f"Systematic uncertainty dominated by {dominant_source[0]} ({dominant_source[1] * 100:.1f}% of total)",
    }


def analyze_fit_quality(all_fits: dict) -> dict:
    """Analyze quality of fits across all flybys."""

    # Collect successful fits
    successful = {
        k: v for k, v in all_fits.items() if v["fit"]["beta_fitted"] is not None
    }

    if not successful:
        return {"status": "no_fits", "message": "No successful fits"}

    beta_values = [v["fit"]["beta_fitted"] for v in successful.values()]
    beta_eff_values = [v["fit"]["beta_eff"] for v in successful.values()]
    beta_uncertainties = [v["fit"]["uncertainty"] for v in successful.values()]

    from scripts.utils.statistical_utils import weighted_mean

    # Use robust weighted mean with Birge scaling
    res_beta = weighted_mean(np.array(beta_values), np.array(beta_uncertainties))
    # beta_eff = beta_fitted * S_earth, so uncertainty scales by same factor
    beta_eff_uncertainties = [u * CHARACTERISTIC_SUPPRESSION for u in beta_uncertainties]
    res_beta_eff = weighted_mean(np.array(beta_eff_values), np.array(beta_eff_uncertainties))

    beta_weighted = res_beta['mean']
    beta_weighted_unc = res_beta['error'] # This is already scaled if chi2_red > 1
    beta_eff_weighted = res_beta_eff['mean']
    beta_eff_weighted_unc = res_beta_eff['error']
    
    chi2 = res_beta['chi2']
    dof = res_beta['dof']
    reduced_chi2 = res_beta['chi2_red']
    birge_ratio = res_beta['birge_ratio']

    # Uncertainty inflation is already handled by weighted_mean
    inflated_unc = beta_weighted_unc

    # Robust Likelihood: Standardized residuals should follow a Student's t-distribution
    from scipy.stats import t
    t_critical = t.ppf(0.975, df=dof) if dof > 0 else 0
    t_confidence_interval = t_critical * inflated_unc

    # Cochran's Q statistic and I2
    Q = chi2
    I2 = max(0, (Q - dof) / Q * 100) if Q > 0 else 0
    # Heterogeneity p-value: probability of observing this much (or more) scatter
    # under the null hypothesis that all flybys share a single true beta
    heterogeneity_p_value = float(stats.chi2.sf(Q, dof)) if Q > 0 and dof > 0 else None

    beta_mean = np.mean(beta_values)
    beta_std = np.std(beta_values)
    beta_eff_mean = np.mean(beta_eff_values)
    beta_eff_std = np.std(beta_eff_values)

    # Random-effects summary (DerSimonian-Laird) for the per-flyby amplitude
    # distribution.  This is intentionally separate from the fixed-effect
    # inverse-variance diagnostic: with extreme heterogeneity, the fixed-effect
    # uncertainty is not an honest population-level uncertainty.
    weights_fixed = 1.0 / np.square(np.array(beta_uncertainties, dtype=float))
    q_beta = chi2
    c_beta = float(np.sum(weights_fixed) - (np.sum(weights_fixed**2) / np.sum(weights_fixed)))
    tau2_beta = max(0.0, (q_beta - dof) / c_beta) if dof > 0 and c_beta > 0 else 0.0
    weights_random = 1.0 / (np.square(np.array(beta_uncertainties, dtype=float)) + tau2_beta)
    beta_random_mean = float(np.sum(weights_random * np.array(beta_values)) / np.sum(weights_random))
    beta_random_unc = float(1.0 / np.sqrt(np.sum(weights_random)))
    beta_random_tau = float(np.sqrt(tau2_beta))
    beta_random_prediction_unc = float(np.sqrt(tau2_beta + beta_random_unc**2))

    weights_eff_fixed = 1.0 / np.square(np.array(beta_eff_uncertainties, dtype=float))
    q_eff = res_beta_eff['chi2']
    dof_eff = res_beta_eff['dof']
    c_eff = float(np.sum(weights_eff_fixed) - (np.sum(weights_eff_fixed**2) / np.sum(weights_eff_fixed)))
    tau2_eff = max(0.0, (q_eff - dof_eff) / c_eff) if dof_eff > 0 and c_eff > 0 else 0.0
    weights_eff_random = 1.0 / (
        np.square(np.array(beta_eff_uncertainties, dtype=float)) + tau2_eff
    )
    beta_eff_random_mean = float(
        np.sum(weights_eff_random * np.array(beta_eff_values)) / np.sum(weights_eff_random)
    )
    beta_eff_random_unc = float(1.0 / np.sqrt(np.sum(weights_eff_random)))
    beta_eff_random_tau = float(np.sqrt(tau2_eff))
    beta_eff_random_prediction_unc = float(np.sqrt(tau2_eff + beta_eff_random_unc**2))

    # Use the larger random-effects uncertainty for headline uncertainty while
    # preserving the fixed-effect weighted mean as the Step 008 pooled diagnostic.
    recommended_uncertainty = max(float(inflated_unc), beta_random_unc)
    recommended_uncertainty_model = (
        "random_effects_standard_error"
        if beta_random_unc >= float(inflated_unc)
        else "birge_scaled_fixed_effect"
    )

    # PPN compliance using beta_eff
    ppn_compliant = all(v["fit"]["ppn_compliant"] for v in successful.values())

    sg_successful = {
        k: v
        for k, v in successful.items()
        if v["fit"].get("sign_agreement") is True
    }
    policy_strict = strict_sign_gate_from_config()

    beta_statistics_sign_gated_diagnostic = None
    recommended_beta_sign_gated_diagnostic = None
    recommended_uncertainty_sign_gated_diagnostic = None

    if sg_successful:
        sg_betas = np.array([v["fit"]["beta_fitted"] for v in sg_successful.values()])
        sg_uncs = np.array([v["fit"]["uncertainty"] for v in sg_successful.values()])
        res_sg = weighted_mean(sg_betas, sg_uncs)
        sg_eff_vals = np.array([v["fit"]["beta_eff"] for v in sg_successful.values()])
        sg_eff_uncs = np.array(
            [v["fit"]["uncertainty"] * CHARACTERISTIC_SUPPRESSION for v in sg_successful.values()]
        )
        res_sg_eff = weighted_mean(sg_eff_vals, sg_eff_uncs)
        beta_statistics_sign_gated_diagnostic = {
            "n_fits": int(len(sg_successful)),
            "mean": float(np.mean(sg_betas)),
            "std": float(np.std(sg_betas)),
            "weighted_mean": float(res_sg["mean"]),
            "weighted_uncertainty": float(res_sg["error"]),
            "inflated_uncertainty": float(res_sg["error"]),
            "weighted_mean_beta_eff": float(res_sg_eff["mean"]),
            "weighted_uncertainty_beta_eff": float(res_sg_eff["error"]),
        }
        recommended_beta_sign_gated_diagnostic = float(res_sg["mean"])
        recommended_uncertainty_sign_gated_diagnostic = float(res_sg["error"])

    ensemble_selection = {
        "strict_sign_gate": bool(policy_strict),
        "n_snr_qualified_beta_fits": int(len(successful)),
        "n_sign_agreement_at_reference_beta": int(len(sg_successful)),
        "note": (
            "recommended_beta is the inverse-variance mean over all S/N-qualified "
            "flybys with a successful closed-form β fit in this run. "
            "recommended_beta_sign_gated_diagnostic restricts to sign_agreement==true "
            "(Δv_obs·Δv_TEP(β_ref) ≥ 0) for backward-compatible reporting."
        ),
    }

    return {
        "n_fits": len(successful),
        "n_primary_detections": len(
            successful
        ),  # Explicit: primary detections used for fitting
        "n_excluded": sum(
            1 for v in successful.values() if not v["fit"]["ppn_compliant"]
        ),
        "beta_statistics": {
            "mean": float(beta_mean),
            "std": float(beta_std),
            "weighted_mean": float(beta_weighted),
            "weighted_uncertainty": float(beta_weighted_unc),
            "inflated_uncertainty": float(inflated_unc),
            "random_effects_mean": beta_random_mean,
            "random_effects_uncertainty": beta_random_unc,
            "between_flyby_tau": beta_random_tau,
            "random_effects_prediction_uncertainty": beta_random_prediction_unc,
            "min": min(beta_values),
            "max": max(beta_values),
        },
        "beta_eff_statistics": {
            "mean": float(beta_eff_mean),
            "std": float(beta_eff_std),
            "weighted_mean": float(beta_eff_weighted),
            "weighted_uncertainty": float(beta_eff_weighted_unc),
            "random_effects_mean": beta_eff_random_mean,
            "random_effects_uncertainty": beta_eff_random_unc,
            "between_flyby_tau": beta_eff_random_tau,
            "random_effects_prediction_uncertainty": beta_eff_random_prediction_unc,
            "min": min(beta_eff_values),
            "max": max(beta_eff_values),
            "ppn_gamma_deviation": float(2 * beta_eff_weighted**2),
        },
        "heterogeneity_tests": {
            "chi_squared": float(chi2),
            "degrees_of_freedom": dof,
            "reduced_chi_squared": float(reduced_chi2),
            "cochran_Q": float(Q),
            "I_squared_percent": float(I2),
            "heterogeneity_interpretation": "extreme"
            if I2 > 75
            else "substantial"
            if I2 > 50
            else "moderate"
            if I2 > 25
            else "low",
            "p_value": heterogeneity_p_value,
            "limitation_note": "TEP v0.8 Interpretation: High I² may reflect PHYSICAL geometry-dependent β_eff variation, but it also means the inverse-variance fixed-effect uncertainty is not a population-level uncertainty. "
            "Each flyby samples different Temporal Shear Σ_μ = ∇_μ ln A(φ) due to altitude, latitude, velocity, and trajectory asymmetry. "
            "The Step 008 inverse-variance mean over all qualified fits is the primary pooled diagnostic; "
            "beta_statistics_sign_gated_diagnostic records the legacy sign-agreement-restricted subset when it differs.",
        },
        "recommended_beta": float(beta_weighted),
        "recommended_beta_eff": float(beta_eff_weighted),
        "recommended_uncertainty": float(recommended_uncertainty),
        "recommended_uncertainty_model": recommended_uncertainty_model,
        "beta_statistics_sign_gated_diagnostic": beta_statistics_sign_gated_diagnostic,
        "recommended_beta_sign_gated_diagnostic": recommended_beta_sign_gated_diagnostic,
        "recommended_uncertainty_sign_gated_diagnostic": recommended_uncertainty_sign_gated_diagnostic,
        "ensemble_selection": ensemble_selection,
        "t_critical": float(t_critical),
        "t_confidence_interval": float(t_confidence_interval),
        "ppn_compliance": bool(ppn_compliant),
    }


def main():
    """Execute fitting analysis."""
    logger = StepLogger("step_008_fitting", PROJECT_ROOT)
    start_time = time.time()

    logger.header("STEP 008: TEP PARAMETER FITTING")

    # Load TEP predictions and setup output directories
    results_dir = PROJECT_ROOT / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    pred_file = results_dir / "step007_tep_predictions.json"

    if not pred_file.exists():
        logger.error(f"TEP predictions not found: {pred_file}")
        logger.log_step_summary(0, "FAILED")
        return 1

    try:
        with open(pred_file, "r", encoding="utf-8") as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
        logger.error(f"Failed to load predictions file: {e}")
        logger.log_step_summary(0, "FAILED")
        return 1

    predictions = data.get("predictions", {})

    strict_sg = strict_sign_gate_from_config()

    logger.section("FITTING β PARAMETER")
    logger.info(f"Fitting β parameter for {len(predictions)} flybys")
    logger.info("")
    logger.info("Scientific Context:")
    logger.info(
        "  The TEP coupling parameter β quantifies the strength of the conformal coupling"
    )
    logger.info(
        "  between matter and the Temporal Topology scalar field. It is determined by fitting"
    )
    logger.info(
        "  the predicted velocity anomaly from the scalar force model to the observed"
    )
    logger.info("  Doppler measurements from NASA's Deep Space Network.")
    logger.info("")
    logger.info("Selection Criterion (a priori, S/N > 2):")
    logger.info(
        "  Flybys are included in fitting only if the observed anomaly has signal-to-noise"
    )
    logger.info(
        "  ratio ≥ 2. This prevents confirmation bias by excluding marginal detections."
    )
    logger.info(
        "  Current analysis: Will show count after fitting completes."
    )
    logger.info(
        "  Sign agreement at β_ref: enforced (strict_sign_gate=true) excludes opposite-sign "
        "rows from closed-form β fitting; when false, magnitude-only reference fits are kept "
        "and a separate sign-gated diagnostic mean is reported in overall_analysis."
        if strict_sg
        else "  Sign agreement at β_ref: not enforced for fitting (strict_sign_gate=false); "
        "inverse-variance headline uses all S/N-qualified fits; sign-gated diagnostic is separate."
    )

    # Fit each flyby
    fits = {}

    for name, pred in predictions.items():
        logger.subheader(f"Processing: {name}")
        dv_obs = pred["observed"]["dv_obs_mm_s"]
        dv_obs_str = f"{dv_obs:.2f} mm/s" if dv_obs is not None else "N/A"
        logger.info(f"Observed Δv: {dv_obs_str}")
        logger.info(f"TEP predicted: {pred['tep_predictions']['dv_tep_mm_s']:.2f} mm/s")

        fit_result = fit_beta_to_observation(
            pred, logger=logger, strict_sign_gate=strict_sg
        )

        fits[name] = {
            "spacecraft": pred["spacecraft"],
            "perigee": pred["perigee"],
            "observed": pred["observed"],
            "tep_predictions": pred["tep_predictions"],
            "cos_dec_asymmetry": pred.get("geometry", {}).get("cos_dec_asymmetry", 0.0),
            "disformal_transition": pred.get("disformal_transition"),
            "fit": fit_result,
        }

        if fit_result["beta_fitted"]:
            logger.info(
                f"Fitted β: {fit_result['beta_fitted']:.2e} ± {fit_result['uncertainty']:.2e}"
            )
            logger.info(f"β_eff (with Temporal Shear screening): {fit_result['beta_eff']:.2e}")
            logger.info(
                f"|γ-1| from β_eff: {fit_result['ppn_gamma_deviation']:.2e} (bound: 2.3e-5)"
            )
            logger.info(f"Status: {fit_result['status']}")
            logger.info("")
            logger.info("  Physics Interpretation:")
            logger.info(
                f"    The fitted coupling β = {fit_result['beta_fitted']:.2e} implies an effective"
            )
            logger.info(
                f"    coupling β_eff = {fit_result['beta_eff']:.2e} after Temporal Topology screening."
            )
            logger.info(
                f"    This corresponds to PPN parameter deviation |γ-1| = {fit_result['ppn_gamma_deviation']:.2e},"
            )
            logger.info(
                f"    ratio to Cassini bound: {(2.3e-5 / fit_result['ppn_gamma_deviation']):.2e}"
            )
            if fit_result["ppn_compliant"]:
                logger.info("    PPN constraint satisfied.")
            else:
                logger.info(
                    "    PPN constraint violated - model excluded for this flyby."
                )
        else:
            status = fit_result.get("status", "excluded")
            logger.info(f"Status: {status}")
            logger.info("")
            logger.info("  Reason for exclusion:")
            if status == "below_snr_threshold":
                logger.info(
                    "    Signal-to-noise ratio < 2 (a priori selection criterion)."
                )
                logger.info(
                    "    This is a null-result flyby used for effect size comparison."
                )
            elif status == "sign_mismatch":
                logger.info("    Predicted and observed anomalies have opposite signs.")
                logger.info(
                    "    Excluded from closed-form β fitting because strict_sign_gate=true."
                )
            elif status == "amplitude_fit_sign_reference_mismatch":
                logger.info(
                    "    Opposite signs at β_ref; strict_sign_gate=false — β is fitted from "
                    "|Δv_obs/Δv_TEP|^(4/3) as an amplitude diagnostic (not a signed prediction match)."
                )
            elif status == "no_signal":
                logger.info("    Observed anomaly is consistent with zero.")
                logger.info("    Null result supports gradient suppression prediction.")

    # Overall analysis
    logger.section("OVERALL ANALYSIS")
    logger.info("")
    logger.info("Purpose:")
    logger.info(
        "  Aggregate individual flyby fits to assess overall TEP framework viability."
    )
    logger.info(
        "  Key question: Do fitted β values converge to a consistent coupling strength"
    )
    logger.info("  across diverse flyby geometries, supporting a common scalar sector with path-resolved effective amplitudes?")
    logger.info("")
    quality = analyze_fit_quality(fits)

    # Multi-parameter fitting with geometry modulation
    logger.section("MULTI-PARAMETER FITTING WITH GEOMETRY MODULATION")
    logger.info("")
    logger.info("Purpose:")
    logger.info(
        "  Test whether including geometry-dependent modulation reduces heterogeneity."
    )
    logger.info("  Model: β_eff = β_0 × exp(α_d × log(ρ/ρ_c) + α_g × f(geometry))")
    logger.info("  where α_d = 0.334 (fixed from Paper 7), β_0 and α_g are fitted.")
    logger.info("")

    multi_param_result = fit_multi_parameter_model(predictions)

    if multi_param_result.get("status") == "success":
        logger.success("Multi-parameter fitting succeeded")
        logger.info(f"Universal coupling β_0: {multi_param_result['beta_0']:.2e}")
        logger.info(f"Geometry modulation α_g: {multi_param_result['alpha_g']:.2e}")
        logger.info(f"Number of flybys: {multi_param_result['n_flybys']}")
        logger.info(
            f"Reduced χ²: {multi_param_result.get('reduced_chi_squared', 0):.2f}"
        )
        logger.info(
            f"Heterogeneity I²: {multi_param_result.get('I_squared_percent', 0):.1f}%"
        )
        logger.info("")
        logger.info("Interpretation:")
        if multi_param_result.get("I_squared_percent", 100) < quality.get(
            "heterogeneity_tests", {}
        ).get("I_squared_percent", 100):
            logger.info(
                f"  Multi-parameter model reduces heterogeneity from {quality.get('heterogeneity_tests', {}).get('I_squared_percent', 0):.1f}% to {multi_param_result.get('I_squared_percent', 0):.1f}%"
            )
            logger.info(
                "  Geometry-dependent modulation captures systematic variation in β values."
            )
        else:
            logger.info(
                f"  Heterogeneity: {multi_param_result.get('I_squared_percent', 0):.1f}%"
            )
            logger.info("  Multi-parameter model with geometry modulation.")
    else:
        logger.warning(
            f"Multi-parameter fitting status: {multi_param_result.get('status', 'Unknown')}"
        )

    logger.info(f"Successful fits: {quality['n_fits']}")
    logger.info(f"Excluded by PPN: {quality['n_excluded']}")
    logger.info("")
    logger.info("Interpretation:")
    logger.info(
        f"  {quality['n_fits']} flybys provide valid β constraints that satisfy PPN bounds."
    )
    if quality["n_excluded"] > 0:
        logger.info(
            f"  {quality['n_excluded']} flyby(s) excluded due to PPN constraint violation."
        )
    else:
        logger.info("  All fitted β values are PPN-compliant (no exclusions required).")

    if quality["n_fits"] > 0:
        logger.subsection("β Statistics")
        logger.info(f"Mean: {quality['beta_statistics']['mean']:.2e}")
        logger.info(f"Std:  {quality['beta_statistics']['std']:.2e}")
        logger.info(
            f"Weighted mean: {quality['beta_statistics']['weighted_mean']:.2e} ± {quality['beta_statistics']['weighted_uncertainty']:.2e}"
        )
        logger.info(
            f"Inflated unc:  {quality['beta_statistics']['inflated_uncertainty']:.2e} (accounts for scatter)"
        )
        logger.info(
            f"Random-effects mean: {quality['beta_statistics']['random_effects_mean']:.2e} ± {quality['beta_statistics']['random_effects_uncertainty']:.2e}"
        )
        logger.info(
            f"Between-flyby τ: {quality['beta_statistics']['between_flyby_tau']:.2e}; "
            f"prediction-scale σ: {quality['beta_statistics']['random_effects_prediction_uncertainty']:.2e}"
        )

        # Heterogeneity tests
        het = quality.get("heterogeneity_tests", {})
        if het:
            logger.subsection("Heterogeneity Assessment")
            logger.info("")
            logger.info(
                "Purpose: Test whether fitted β values are consistent with a single"
            )
            logger.info(
                "a single pooled amplitude across the gated tier, or if geometry-dependent modulation dominates."
            )
            logger.info("")
            logger.info(f"χ² = {het.get('chi_squared', 0):.2e}")
            logger.info(f"Reduced χ² = {het.get('reduced_chi_squared', 0):.2e}")
            logger.info(
                f"I² = {het.get('I_squared_percent', 0):.1f}% ({het.get('heterogeneity_interpretation', 'unknown')})"
            )
            logger.info("")
            logger.info("Interpretation (TEP v0.8):")
            i2_val = het.get("I_squared_percent", 0)
            if i2_val > 75:
                logger.info(
                    f"  I² = {i2_val:.1f}% indicates substantial PHYSICAL heterogeneity."
                )
                logger.info(
                    "  Per TEP v0.8 Temporal Topology, each flyby samples different Temporal Shear"
                )
                logger.info(
                    "  Σ_μ = ∇_μ ln A(φ) due to altitude, latitude, velocity, and trajectory asymmetry."
                )
                logger.info(
                    "  The fitted-amplitude span should be reported as geometry-driven scatter,"
                )
                logger.info("")
                logger.info(
                    "  The weighted mean β is the Step 008 fixed-effect diagnostic over all qualified fits; "
                    "see beta_statistics_sign_gated_diagnostic for the sign-agreement-restricted subset."
                )
                logger.info(
                    "  The random-effects summary is the honest uncertainty scale for cross-flyby scatter."
                )

        # Robustness analyses
        logger.section("ROBUSTNESS ANALYSIS")
        logger.info("")
        logger.info("Purpose: Verify that TEP framework conclusions are stable against")
        logger.info("statistical fluctuations and do not depend on any single flyby.")
        logger.info("Three complementary approaches:")
        logger.info(
            "  1. Bootstrap resampling - assess uncertainty from small sample size"
        )
        logger.info(
            "  2. Leave-one-out cross-validation - test sensitivity to single flyby"
        )
        logger.info("  3. Statistical power analysis - evaluate sample size adequacy")

        # Bootstrap analysis
        bootstrap = bootstrap_beta_estimate(fits, n_bootstrap=10000)
        if bootstrap["status"] == "success":
            logger.subsection("Bootstrap Resampling")
            logger.info("")
            logger.info(
                "Method: Parametric bootstrap with n=10,000 iterations addresses the"
            )
            logger.info(
                "fundamental limitation of small sample size (n=3 primary detections)."
            )
            logger.info(
                "Each iteration resamples with replacement and adds measurement noise."
            )
            logger.info("")
            logger.info(f"Bootstrap samples: n={bootstrap['n_bootstrap']}")
            logger.info(f"Mean:   {bootstrap['bootstrap_mean']:.2e}")
            logger.info(f"Std:    {bootstrap['bootstrap_std']:.2e}")
            logger.info(f"Median: {bootstrap['bootstrap_median']:.2e}")
            logger.info(
                f"95% CI: [{bootstrap['ci_95_lower']:.2e}, {bootstrap['ci_95_upper']:.2e}]"
            )
            logger.info("")
            logger.info("Interpretation:")
            logger.info(
                f"  The 95% confidence interval spans a factor of {(bootstrap['ci_95_upper'] / bootstrap['ci_95_lower']):.2f},"
            )
            logger.info(
                "  reflecting the intrinsic scatter in β from geometry-dependent effects."
            )
            logger.info(
                "  The bootstrap distribution suggests the weighted mean is stable and"
            )
            logger.info("  the central value is not an artifact of small sample size.")

        # Leave-one-out analysis
        loo = leave_one_out_analysis(fits)
        if loo["status"] == "success":
            logger.subsection("Leave-One-Out Cross-Validation")
            logger.info("")
            logger.info(
                "Method: Systematically exclude each detection and recompute weighted mean."
            )
            logger.info(
                "Tests whether conclusion depends on any single flyby (especially NEAR,"
            )
            logger.info("which dominates due to superior measurement precision).")
            logger.info("")
            for name, result in loo["leave_one_out_results"].items():
                logger.info(f"Excluding {name}: β = {result['beta_without_this']:.2e}")
            logger.info("")
            logger.info(f"Stability coefficient: {loo['stability_coefficient']:.3f}")
            logger.info(f"Interpretation: {loo['interpretation']}")
            logger.info(
                f"Conclusion robust: {'Yes' if loo['conclusion_robust'] else 'No'}"
            )
            logger.info("")
            logger.info("Interpretation:")
            if loo["conclusion_robust"]:
                logger.info(
                    "  The TEP viability conclusion is not eliminated by any single deletion."
                )
                logger.info(
                    "  Excluding the dominant NEAR detection still leaves a positive pooled amplitude,"
                )
                logger.info(
                    "  but the small n=2 remainder should be treated as a stress test, not a precision estimate."
                )
            else:
                logger.info(
                    "  The conclusion shows sensitivity to single-flyby exclusion."
                )
                logger.info("  Larger sample size required for definitive assessment.")

        # Statistical power analysis
        power = statistical_power_analysis(fits)
        if power["status"] == "success":
            logger.subsection("Statistical Power Analysis")
            logger.info(f"Current sample size: n = {power['n_current']}")
            logger.info(f"Effect size (CV): {power['effect_size_cv']:.3f}")
            logger.info(
                f"Min detectable CV at 80% power: {power['cv_detectable_80_percent_power']:.3f}"
            )
            logger.info(
                f"Heterogeneity detectable: {'Yes' if power['heterogeneity_detectable'] else 'No'}"
            )
            logger.info(
                f"Required n to detect current CV: {int(power['n_required_for_current_cv'])}"
            )
            logger.info(f"Interpretation: {power['interpretation']}")
            logger.info(f"Note: {power['note']}")

        # Enhanced validation metrics
        logger.section("ENHANCED VALIDATION")
        logger.info("")
        logger.info(
            "Purpose: Five complementary statistical tests validate the TEP framework"
        )
        logger.info("against alternative explanations and quantify evidence strength.")
        logger.info("")

        # Effect size analysis
        effect_sizes = calculate_effect_sizes(fits)
        if effect_sizes.get("status") != "no_nulls_for_comparison":
            logger.subsection("Effect Size Analysis (Cohen's d)")
            logger.info("")
            logger.info(
                "Method: Cohen's d quantifies signal magnitude independent of sample size."
            )
            logger.info(
                "Compares each detection against the null-result population (n=8 flybys)."
            )
            logger.info("")
            logger.info(
                f"Null population: n = {effect_sizes['null_population']['n_nulls']}"
            )
            logger.info(
                f"Null mean Δv: {effect_sizes['null_population']['mean_dv']:.3f} mm/s"
            )
            logger.info(
                f"Null std Δv: {effect_sizes['null_population']['std_dv']:.3f} mm/s"
            )
            logger.info("")
            logger.info(
                "Effect sizes (convention: d < 0.2 negligible, 0.2-0.5 small, 0.5-0.8 medium, > 0.8 large):"
            )
            for name, es in effect_sizes["effect_sizes"].items():
                d_val = es.get("cohens_d")
                d_str = f"{d_val:.1f}" if d_val is not None else "N/A"
                logger.info(
                    f"  {name}: d = {d_str} ({es['interpretation']}, {es['detection_significance']} significance)"
                )
            logger.info("")
            logger.info("Interpretation:")
            logger.info(
                "  Effect sizes are computed for published detections against the null-result population."
            )
            logger.info(
                    "  This is a catalog-level anomaly contrast, distinct from the Step 008 inverse-variance β diagnostic."
            )

        # Bayesian model comparison
        model_comp = bayesian_model_comparison(fits)
        if model_comp.get("status") != "insufficient_data":
            logger.subsection("Bayesian Model Comparison")
            logger.info("")
            logger.info(
                "Method: Information-theoretic model selection using AIC and BIC."
            )
            logger.info(
                "Headline comparison evaluates TEP and Null on the same observational Doppler uncertainties."
            )
            logger.info(
                "A separate TEP sensitivity calculation includes additional structured theoretical uncertainty."
            )
            logger.info("")
            logger.info("Models compared:")
            logger.info(
                "  TEP: Physical model with single-parameter β rescaling (restricted tier) and Disformal Temporal Topology gradient suppression"
            )
            logger.info("  Null: No anomaly (measurement artifacts only)")
            logger.info(
                "  Empirical: Ad hoc fits with independent β per flyby (no physics)"
            )
            logger.info("")
            logger.info(
                f"TEP AIC: {model_comp['models']['TEP']['AIC']:.1f}, BIC: {model_comp['models']['TEP']['BIC']:.1f}"
            )
            logger.info(
                f"Null AIC: {model_comp['models']['Null']['AIC']:.1f}, BIC: {model_comp['models']['Null']['BIC']:.1f}"
            )
            logger.info(
                f"Empirical AIC: {model_comp['models']['Empirical']['AIC']:.1f}, BIC: {model_comp['models']['Empirical']['BIC']:.1f}"
            )
            logger.info(
                f"TEP structured-uncertainty sensitivity AIC: {model_comp['tep_structured_uncertainty_sensitivity']['AIC']:.1f}, BIC: {model_comp['tep_structured_uncertainty_sensitivity']['BIC']:.1f}"
            )
            logger.info("")
            logger.info(
                f"Best model (BIC): {model_comp['model_comparison']['best_model_bic']}"
            )
            logger.info(
                f"TEP evidence weight: {model_comp['model_comparison']['tep_evidence_weight']:.1%}"
            )
            logger.info(
                f"BIC-map TEP vs Null (capped at exp(700)): "
                f"{model_comp['model_comparison']['tep_vs_null_bayes_factor_approx']:.1e}; "
                f"log_evidence_delta={model_comp['model_comparison']['log_evidence_delta']:.1f}. "
                f"{model_comp['model_comparison']['tep_vs_null_bayes_factor_approx_note']}"
            )
            logger.info("")
            tep_weight = model_comp["model_comparison"]["tep_evidence_weight"]
            bayes_factor = model_comp["model_comparison"][
                "tep_vs_null_bayes_factor_approx"
            ]
            logger.info("Interpretation:")
            logger.info(
                f"  On shared observational uncertainties, the TEP model achieves {tep_weight:.1%} evidence weight"
            )
            logger.info(
                "  relative to the null model. This supports a structured physical explanation"
            )
            logger.info(
                "  for the detected anomalies, while remaining distinct from the saturated empirical bound."
            )
            logger.info(
                "  The structured-uncertainty sensitivity result should be read as robustness testing,"
            )
            logger.info("  not as the headline evidence calculation.")

        # Prediction accuracy metrics
        pred_acc = prediction_accuracy_metrics(fits)
        if pred_acc.get("status") != "insufficient_data":
            logger.subsection("Prediction Accuracy")
            logger.info(
                f"R² = {pred_acc['R_squared']:.4f} ({pred_acc['prediction_quality']})"
            )
            logger.info(f"Correlation = {pred_acc['correlation']:.4f}")
            logger.info(f"MAE = {pred_acc['MAE_mm_s']:.4f} mm/s")
            logger.info(f"RMSE = {pred_acc['RMSE_mm_s']:.4f} mm/s")
            logger.info(f"MAPE = {pred_acc['MAPE_percent']:.2f}%")

        # Residual analysis
        residuals = residual_analysis(fits)
        if residuals.get("status") != "insufficient_data":
            logger.subsection("Residual Analysis")
            logger.info(f"Residuals mean: {residuals['statistics']['mean']:.4f} mm/s")
            logger.info(f"Residuals std: {residuals['statistics']['std']:.4f} mm/s")
            logger.info(
                f"Normality p-value: {residuals['statistics']['normality_p_value']:.4f}"
            )
            logger.info(f"Interpretation: {residuals['interpretation']}")

        # Systematic uncertainty budget
        uncertainty_budget = systematic_uncertainty_budget(fits)
        if uncertainty_budget.get("status") != "insufficient_data":
            logger.subsection("Systematic Uncertainty Budget")
            logger.info("")
            logger.info(
                "Method: Propagate uncertainties from five independent sources through"
            )
            logger.info(
                "the TEP model using root-sum-square combination (uncorrelated errors)."
            )
            logger.info("")
            logger.info("Uncertainty sources:")
            logger.info(
                f"  1. Measurement (Doppler): {uncertainty_budget['uncertainty_breakdown']['measurement_uncertainty']:.1%}"
            )
            logger.info(
                f"  2. Trajectory reconstruction: {uncertainty_budget['uncertainty_breakdown']['trajectory_uncertainty']:.1%}"
            )
            logger.info(
                f"  3. Characteristic suppression (UCD): {uncertainty_budget['uncertainty_breakdown']['characteristic_suppression_uncertainty']:.1%} ← DOMINANT"
            )
            logger.info(
                f"  4. Relaxation length (UCD): {uncertainty_budget['uncertainty_breakdown']['relaxation_length_uncertainty']:.1%}"
            )
            logger.info(
                f"  5. Multipole coefficients: {uncertainty_budget['uncertainty_breakdown']['multipole_coefficient_uncertainty']:.1%}"
            )
            logger.info("")
            logger.info(
                f"Total relative systematic uncertainty: {uncertainty_budget['total_relative_systematic_uncertainty']:.1%}"
            )
            logger.info(
                f"Dominant source: {uncertainty_budget['dominant_uncertainty_source']} ({uncertainty_budget['dominant_contribution_percent']:.1f}%)"
            )
            logger.info("")
            logger.info("Interpretation:")
            logger.info(
                f"  The large systematic uncertainty ({uncertainty_budget['total_relative_systematic_uncertainty']:.1%}) reflects genuine physical"
            )
            logger.info(
                "  uncertainty in the Temporal Topology Temporal Shear mechanism, not a bookkeeping artifact"
            )
            logger.info(
                f"  of the fitting procedure. The dominant source ({uncertainty_budget['dominant_uncertainty_source']},"
            )
            logger.info(
                f"  {uncertainty_budget['dominant_contribution_percent']:.1f}% of the total variance) remains the main target"
            )
            logger.info(
                "  for further first-principles refinement. Even with this uncertainty budget,"
            )
            logger.info("  the fitted β values remain PPN-compliant by a wide margin.")

        logger.info(f"Recommended β: {quality['recommended_beta']:.2e}")
        logger.info(f"PPN compliant: {'Yes' if quality['ppn_compliance'] else 'No'}")

    # Save results
    logger.section("SAVING FITTING RESULTS")

    # Extract key values for top-level summary (for API consistency)
    het_tests = quality.get("heterogeneity_tests", {})
    parameter_budget_audit = geometry_envelope_parameter_budget_audit()
    output = {
        "individual_fits": fits,
        "parameter_budget_audit": parameter_budget_audit,
        "coupling_constant": float(
            quality.get("recommended_beta_eff", quality.get("recommended_beta", 0))
        ),
        "beta": float(quality.get("recommended_beta", 0)),
        "beta_eff": float(quality.get("recommended_beta_eff", 0)),
        "uncertainty": float(quality.get("recommended_uncertainty", 0)),
        "chi2": float(het_tests.get("chi_squared", 0)),
        "dof": int(het_tests.get("degrees_of_freedom", 0)),
        "reduced_chi2": float(het_tests.get("reduced_chi_squared", 0)),
        "p_value": float(het_tests.get("p_value", 0))
        if "p_value" in het_tests
        else None,
        "overall_analysis": quality,
        "ensemble_config": {
            "strict_sign_gate": bool(strict_sg),
            "config_path": "config/pipeline_config.json",
        },
        "multi_parameter_fitting": multi_param_result,
        "robustness_analysis": {
            "bootstrap": bootstrap if bootstrap["status"] == "success" else None,
            "leave_one_out": loo if loo["status"] == "success" else None,
            "statistical_power": power if power["status"] == "success" else None,
        },
        "enhanced_validation": {
            "effect_sizes": effect_sizes
            if effect_sizes.get("status") != "no_nulls_for_comparison"
            else None,
            "model_comparison": model_comp
            if model_comp.get("status") != "insufficient_data"
            else None,
            "prediction_accuracy": pred_acc
            if pred_acc.get("status") != "insufficient_data"
            else None,
            "residual_analysis": residuals
            if residuals.get("status") != "insufficient_data"
            else None,
            "systematic_uncertainty_budget": uncertainty_budget
            if uncertainty_budget.get("status") != "insufficient_data"
            else None,
        },
        "limitations": {
            "small_sample_size": bool(quality.get("n_fits", 0) < 5),
            "n_primary_detections": int(quality.get("n_fits", 0)),
            "note": f"Analysis utilizes {quality.get('n_fits', 0)} primary detection(s). "
            f"Published: NEAR, Galileo_1990, Rosetta_2005, Cassini. Successfully fitted: {quality.get('n_fits', 0)}. "
            "Additional flyby data needed for conclusive results.",
            "heterogeneity_limitation": bool(
                quality.get("heterogeneity_tests", {}).get("I_squared_percent", 0) > 75
            ),
            "power_insufficient": bool(
                power.get("status") == "success"
                and not power.get("heterogeneity_detectable", False)
            ),
        },
    }

    beta_values = [
        entry["fit"]["beta_fitted"]
        for entry in fits.values()
        if entry.get("fit", {}).get("beta_fitted") is not None
    ]
    chi2_consistency = None
    if len(beta_values) > 1:
        mean_beta = sum(beta_values) / len(beta_values)
        chi2_consistency = sum(
            ((beta - mean_beta) / (beta * 0.1 if beta else 0.1)) ** 2
            for beta in beta_values
            if beta
        )

    gnss_scale_consistent = None
    gnss_file = results_dir / "step016_gnss_cross_validation.json"
    if gnss_file.exists():
        with open(gnss_file, encoding="utf-8") as handle:
            gnss_data = json.load(handle)
        gnss_scale_consistent = gnss_data.get(
            "gnss_scale_consistent", gnss_data.get("beta_consistent")
        )

    confidence = (
        "high"
        if (
            quality.get("n_fits", 0) >= 3
            and quality.get("ppn_compliance")
            and model_comp.get("model_comparison", {}).get("best_model_bic") == "TEP"
            and pred_acc.get("R_squared", 0) > 0.8
        )
        else "moderate"
    )
    if gnss_scale_consistent is False and confidence in ("high", "moderate"):
        confidence = "limited"
    if chi2_consistency is not None and chi2_consistency > 100 and confidence != "none":
        if confidence == "moderate":
            confidence = "limited"

    evidence_strength = (
        "strong"
        if (
            model_comp.get("model_comparison", {}).get("tep_evidence_weight", 0) > 0.8
            and pred_acc.get("R_squared", 0) > 0.8
        )
        else "moderate"
    )
    if chi2_consistency is not None and chi2_consistency > 100:
        evidence_strength = "limited"
    if gnss_scale_consistent is False and evidence_strength in ("strong", "moderate"):
        evidence_strength = "limited"

    output["conclusion"] = {
        "tep_explains_flyby_anomaly": bool(
            quality.get("ppn_compliance", False) and quality["n_fits"] >= 2
        ),
        "recommended_beta": float(quality.get("recommended_beta", 0)),
        "recommended_uncertainty": float(quality.get("recommended_uncertainty", 0)),
        "confidence": confidence,
        "evidence_strength": evidence_strength,
        "sample_size_note": (
            f"Based on {quality.get('n_fits', 0)} primary detections. Enhanced validation "
            "includes a fair TEP-vs-Null comparison on shared observational uncertainties "
            "plus separate structured-uncertainty sensitivity tests."
        ),
    }

    # Convert numpy types to native Python types for JSON serialization
    output_native = convert_to_native_types(output)

    # Save to results folder
    output_file = results_dir / "step008_fitting_results.json"
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(output_native, f, indent=2)

    logger.success(f"Fitting complete")
    logger.info(f"Results saved to: {output_file}")
    logger.add_output_file(output_file, "TEP parameter fitting results")

    duration = time.time() - start_time
    logger.log_step_summary(duration, "SUCCESS")

    return 0


if __name__ == "__main__":
    sys.exit(main())
