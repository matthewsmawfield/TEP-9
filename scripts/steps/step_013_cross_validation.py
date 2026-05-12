#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 019: Comprehensive Cross-Validation

This module implements multiple cross-validation strategies to test TEP model
robustness and predictive power.

Validation Strategies:
----------------------
1. Leave-One-Out (LOO): Exclude each of the 4 primary detections, predict
   using the remaining 3. Uses step008 fitted β values directly so that the
   disformal sign correction and S/N filtering are already applied.
2. Bootstrap: Resample the 4 fitted betas (with Gaussian noise on their
   uncertainties) to produce empirical confidence intervals.
3. Altitude-Stratified: Group by perigee altitude regime.

Key fix over previous version
------------------------------
The previous implementation:
  - Loaded ALL flybys from step007_tep_predictions.json including
    Rosetta_2007 (S/N=0.4, below the a-priori S/N>2 threshold).
  - Used a naive ratio fit beta = 1e-4 * (dv_obs / dv_pred) applied
    to training sets that mixed primary detections with sub-threshold
    flybys, corrupting the median β.
  - Did not enforce β > 0 physically.

The corrected implementation:
  - Loads primary detections (excluded=False, β>0) from step008.
  - LOO uses arithmetic mean of the 3 remaining fitted betas to predict
    the left-out flyby (consistent with manuscript Section 4.5).
  - Bootstrap resamples the 4 fitted betas with measurement noise.
  - Enforces β > 0 at all stages.
"""

import numpy as np
import json
from pathlib import Path
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


class CrossValidator:
    """
    Comprehensive cross-validation for TEP model.

    Loads validated primary detections from step008_fitting_results.json
    (which already applies the Holonomic Hybrid disformal coupling and
    S/N > 2 selection). The raw step007 predictions (at beta_ref=1e-4)
    are used solely to scale predictions for the left-out flyby.
    """

    def __init__(self):
        self.logger = StepLogger("step_013_cross_validation", PROJECT_ROOT)

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def load_primary_detections(self):
        """
        Load the 4 primary detections from step008.

        Returns a list of dicts with keys:
          name, beta_fitted, beta_unc, dv_obs, dv_unc, dv_pred_ref,
          beta_ref, altitude_km
        where dv_pred_ref is the step007 raw prediction at beta_ref=1e-4.
        """
        fit_file = PROJECT_ROOT / 'results' / 'step008_fitting_results.json'
        pred_file = PROJECT_ROOT / 'results' / 'step007_tep_predictions.json'

        try:
            with open(fit_file) as f:
                fit_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            logger.error(f"Failed to load fitting data: {e}")
            return None
        
        try:
            with open(pred_file) as f:
                pred_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            logger.error(f"Failed to load predictions data: {e}")
            return None

        detections = []
        for name, entry in fit_data['individual_fits'].items():
            fit = entry['fit']
            # Only include primary detections that are fitted with β > 0
            if fit.get('excluded', True):
                continue
            beta = fit.get('beta_fitted')
            if beta is None or beta <= 0:
                continue

            # Step007 raw prediction at beta_ref = 1e-4
            step4 = pred_data['predictions'].get(name, {})
            dv_pred_ref = step4.get('tep_predictions', {}).get('dv_tep_mm_s', None)
            beta_ref = step4.get('tep_predictions', {}).get('beta_initial', 1e-4)  # 1e-4

            detections.append({
                'name': name,
                'beta_fitted': beta,
                'beta_unc': fit.get('uncertainty', beta * 0.1),
                'dv_obs': entry['observed']['dv_obs_mm_s'],
                'dv_unc': entry['observed'].get('sigma_mm_s', entry['observed'].get('dv_unc_mm_s', 0.1)),
                'dv_pred_ref': dv_pred_ref,  # at beta_ref
                'beta_ref': beta_ref,
                'altitude_km': entry['perigee']['altitude_km'],
            })

        return detections

    # ------------------------------------------------------------------
    # Core metric calculation
    # ------------------------------------------------------------------

    def _predict_with_beta(self, beta_loo, flyby):
        """
        Scale the step007 raw prediction by the LOO β.

        dv_pred(β) = dv_pred_ref * (β / β_ref)^(3/4)
        The TEP velocity shift follows a 3/4 power law in β because the
        scalar force ∝ β * ∇φ ∝ β * β^(-1/4) = β^(3/4).
        """
        if flyby['dv_pred_ref'] is None or flyby['beta_ref'] == 0:
            return None
        return flyby['dv_pred_ref'] * ((beta_loo / flyby['beta_ref']) ** 0.75) if flyby['beta_ref'] != 0 else 0.0

    # ------------------------------------------------------------------
    # Leave-One-Out CV
    # ------------------------------------------------------------------

    def leave_one_out_cv(self, detections):
        """
        Leave-one-out cross-validation on the 4 primary detections.

        For each held-out flyby:
          - Compute the arithmetic mean β from the remaining 3 fitted betas.
            (Arithmetic mean used for consistency with manuscript §4.5.)
          - Enforce β_loo > 0 (it always should be since all individual
            betas are positive).
          - Scale the step007 raw prediction to get dv_pred at β_loo.
          - Record residual vs observed anomaly.
        """
        self.logger.section("LEAVE-ONE-OUT CROSS-VALIDATION")

        results = []

        for i, test in enumerate(detections):
            train = [d for j, d in enumerate(detections) if j != i]

            # Arithmetic mean β from remaining 3 (manuscript §4.5)
            beta_loo = float(np.mean([d['beta_fitted'] for d in train]))
            beta_loo = max(beta_loo, 0.0)  # enforce positivity (redundant but safe)

            # Predict the held-out flyby
            dv_pred = self._predict_with_beta(beta_loo, test)

            if dv_pred is None:
                self.logger.warning(f"Cannot predict {test['name']}: no step007 raw prediction")
                continue

            residual = test['dv_obs'] - dv_pred
            abs_err = abs(residual)
            rel_err = abs_err / abs(test['dv_obs']) * 100 if test['dv_obs'] != 0 else float('inf')
            correct_sign = (dv_pred * test['dv_obs']) > 0

            results.append({
                'left_out': test['name'],
                'beta_loo': beta_loo,
                'dv_obs': test['dv_obs'],
                'dv_pred': dv_pred,
                'residual_mm_s': residual,
                'abs_error_mm_s': abs_err,
                'rel_error_pct': rel_err,
                'correct_sign': correct_sign,
            })

            self.logger.info(f"Left out {test['name']}:")
            self.logger.info(f"  β_loo            = {beta_loo:.3e}")
            self.logger.info(f"  Predicted Δv     = {dv_pred:.4f} mm/s")
            self.logger.info(f"  Observed  Δv     = {test['dv_obs']:.4f} mm/s")
            self.logger.info(f"  |Residual|       = {abs_err:.4f} mm/s ({rel_err:.1f}%)")
            self.logger.info(f"  Sign correct     = {correct_sign}")

        # --- Stability analysis ---
        betas_loo = [r['beta_loo'] for r in results]
        beta_mean_loo = float(np.mean(betas_loo))
        beta_std_loo = float(np.std(betas_loo))
        stability = beta_std_loo / beta_mean_loo if beta_mean_loo > 0 else float('inf')

        sign_accuracy = sum(1 for r in results if r['correct_sign']) / len(results) if results else 0.0

        summary = {
            'n_folds': len(results),
            'beta_loo_mean': beta_mean_loo,
            'beta_loo_std': beta_std_loo,
            'stability_coefficient': float(stability),
            'stability_assessment': (
                'stable' if stability < 0.5 else
                'moderate' if stability < 1.0 else 'unstable'
            ),
            'sign_accuracy': sign_accuracy,
            'fold_results': results,
        }

        self.logger.info(f"\nLOO-CV Summary:")
        self.logger.info(f"  β_loo = {beta_mean_loo:.3e} ± {beta_std_loo:.3e}")
        self.logger.info(f"  Stability coefficient: {stability:.3f} ({summary['stability_assessment']})")
        self.logger.info(f"  Sign accuracy: {sign_accuracy*100:.0f}%")

        return summary

    # ------------------------------------------------------------------
    # Bootstrap
    # ------------------------------------------------------------------

    def bootstrap_validation(self, detections, n_bootstrap=10000):
        """
        Bootstrap resampling of the 4 primary detections.

        Each bootstrap iteration:
          - Resamples with replacement from the 4 fitted betas.
          - Adds Gaussian noise scaled to each beta's uncertainty.
          - Enforces β > 0 before recording (negative draws from noise
            are physically impossible and discarded).
        """
        self.logger.section("BOOTSTRAP VALIDATION")

        if len(detections) < 3:
            return {'status': 'insufficient_data'}

        np.random.seed(42)
        beta_boots = []

        for _ in range(n_bootstrap):
            # Resample with replacement
            indices = np.random.choice(len(detections), size=len(detections), replace=True)
            sample = [detections[j] for j in indices]

            # Add measurement noise and enforce β > 0
            valid_betas = []
            for d in sample:
                noise = np.random.normal(0, d['beta_unc'])
                b = d['beta_fitted'] + noise
                if b > 0:
                    valid_betas.append(b)

            if valid_betas:
                beta_boots.append(float(np.mean(valid_betas)))

        beta_boots = np.array(beta_boots)

        results = {
            'n_bootstrap': n_bootstrap,
            'n_valid_samples': len(beta_boots),
            'beta_mean': float(np.mean(beta_boots)),
            'beta_median': float(np.median(beta_boots)),
            'beta_std': float(np.std(beta_boots)),
            'ci_68': [
                float(np.percentile(beta_boots, 16)),
                float(np.percentile(beta_boots, 84)),
            ],
            'ci_95': [
                float(np.percentile(beta_boots, 2.5)),
                float(np.percentile(beta_boots, 97.5)),
            ],
        }

        self.logger.info(f"Bootstrap Results ({n_bootstrap} samples):")
        self.logger.info(f"  β = {results['beta_mean']:.3e} ± {results['beta_std']:.3e}")
        self.logger.info(f"  68% CI: [{results['ci_68'][0]:.3e}, {results['ci_68'][1]:.3e}]")
        self.logger.info(f"  95% CI: [{results['ci_95'][0]:.3e}, {results['ci_95'][1]:.3e}]")

        return results

    # ------------------------------------------------------------------
    # Altitude-stratified
    # ------------------------------------------------------------------

    def altitude_stratified_cv(self, detections):
        """
        Report fitted β grouped by perigee altitude regime.
        Uses the step008 fitted betas directly (no re-fitting).
        """
        self.logger.section("ALTITUDE-STRATIFIED ANALYSIS")

        thresholds = {'low_altitude': (0, 1000), 'mid_altitude': (1000, 5000)}
        results = {}

        for label, (lo, hi) in thresholds.items():
            subset = [d for d in detections if lo <= d['altitude_km'] < hi]
            if not subset:
                continue
            betas = [d['beta_fitted'] for d in subset]
            results[label] = {
                'n_flybys': len(subset),
                'beta_mean': float(np.mean(betas)),
                'beta_median': float(np.median(betas)),
                'beta_range': [float(min(betas)), float(max(betas))],
                'members': [d['name'] for d in subset],
            }
            self.logger.info(
                f"{label}: β = {results[label]['beta_mean']:.3e} "
                f"(n={len(subset)}, members={[d['name'] for d in subset]})"
            )

        return results

    # ------------------------------------------------------------------
    # Master runner
    # ------------------------------------------------------------------

    def run_validation(self):
        """Execute all cross-validation strategies."""
        self.logger.header("STEP 013: CROSS-VALIDATION AND ROBUSTNESS ANALYSIS")
        self.logger.info(
            "Loading primary detections from step008 (S/N>2, excluded=False, β>0). "
            "Rosetta_2007 (S/N=0.4) and other sub-threshold flybys are correctly "
            "excluded from training data."
        )

        detections = self.load_primary_detections()
        self.logger.info(f"Primary detections loaded: {len(detections)}")
        for d in detections:
            self.logger.info(
                f"  {d['name']:20s}: β={d['beta_fitted']:.3e}, "
                f"dv_obs={d['dv_obs']:.2f} mm/s, alt={d['altitude_km']:.0f} km"
            )

        # Run validations
        loo_results = self.leave_one_out_cv(detections)
        boot_results = self.bootstrap_validation(detections)
        alt_results = self.altitude_stratified_cv(detections)

        # ---------------------------------------------------------------
        # Overall assessment
        # ---------------------------------------------------------------
        self.logger.section("OVERALL VALIDATION SUMMARY")

        beta_fitted_values = [d['beta_fitted'] for d in detections]
        beta_overall_mean = float(np.mean(beta_fitted_values))
        beta_overall_std = float(np.std(beta_fitted_values))
        beta_scatter_ratio = max(beta_fitted_values) / min(beta_fitted_values)

        all_results = {
            'data_source': 'step008_fitting_results.json (primary detections only, S/N>2)',
            'n_primary_detections': len(detections),
            'detection_names': [d['name'] for d in detections],
            'fitted_betas': {d['name']: d['beta_fitted'] for d in detections},
            'fitted_beta_uncertainties': {d['name']: d['beta_unc'] for d in detections},
            'beta_summary': {
                'mean': beta_overall_mean,
                'std': beta_overall_std,
                'min': float(min(beta_fitted_values)),
                'max': float(max(beta_fitted_values)),
                'scatter_ratio': float(beta_scatter_ratio),
                'all_positive': all(b > 0 for b in beta_fitted_values),
            },
            'leave_one_out': loo_results,
            'bootstrap': boot_results,
            'altitude_stratified': alt_results,
            'summary': {
                'n_total_flybys': len(detections),
                'n_with_observations': len(detections),
                'validation_status': 'complete',
            },
        }

        if 'beta_mean' in boot_results:
            self.logger.info(f"Bootstrap β: {boot_results['beta_mean']:.3e} ± {boot_results['beta_std']:.3e}")
        if 'stability_coefficient' in loo_results:
            sc = loo_results['stability_coefficient']
            self.logger.info(f"LOO stability: {sc:.3f} ({loo_results['stability_assessment']})")
        self.logger.info(f"β scatter ratio: {beta_scatter_ratio:.2f}× (all positive: {all_results['beta_summary']['all_positive']})")

        # Save
        output_file = PROJECT_ROOT / 'results' / 'step013_cross_validation.json'
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2)

        self.logger.success(f"Cross-validation complete. Saved to {output_file}")
        return all_results


def main():
    logger = StepLogger("step_013_cross_validation", PROJECT_ROOT)
    start_time = time.time()

    validator = CrossValidator()
    results = validator.run_validation()

    duration = time.time() - start_time
    if results:
        validator.logger.log_step_summary(duration, "SUCCESS")
        return 0
    else:
        validator.logger.log_step_summary(duration, "FAILED")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
