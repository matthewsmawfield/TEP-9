#!/usr/bin/env python3
"""
Step 032: Component-Level Hierarchical Model
==========================================

This module implements a component-level hierarchical model that separates
the gradient and disformal contributions to the TEP prediction, rather than
scaling only the already-cancelled total.

This directly addresses the core issue identified by the investigation:
the hierarchical model should not scale only the already-cancelled dv_tep_ref,
but should treat the physical components separately.

Model:
    Δv_i = a_grad * Δv_grad,i + a_disf * Δv_disf,i + u_r(i) + ε_i

where:
    a_grad = β_0 / 10^-4 * G_traj * S_earth * F_plasma
    a_disf = b_disf / 0.05 * G_disf
    u_r(i) ~ N(0, σ_r) - regime-level random effect
    ε_i ~ Student-t(ν, 0, σ) - residual error

This allows:
1. Cassini's cancellation regime to be properly modeled
2. Different scaling for gradient vs disformal components
3. Regime-specific random effects (not flyby-specific fudge factors)
"""

import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Any

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import DISFORMAL_COUPLING_STRENGTH, BETA_BASELINE

# BETA_INITIAL from physics.py (Unified Yogyakarta anchor)
BETA_INITIAL = BETA_BASELINE * 1e-4  # Convert baseline to actual coupling


class RegimeClassifier:
    """Classify flybys into physical regimes for random effects."""
    
    REGIMES = {
        "low_altitude_gradient": {
            "definition": "h < 800 km, positive asymmetry",
            "examples": ["NEAR"]
        },
        "mid_altitude_enhancement": {
            "definition": "800 km < h < 1500 km, positive asymmetry",
            "examples": ["Galileo_1990"]
        },
        "high_velocity_anti_aligned_cancellation": {
            "definition": "v > 16 km/s, negative asymmetry",
            "examples": ["Cassini"]
        },
        "high_altitude_suppressed": {
            "definition": "h > 2000 km or low asymmetry",
            "examples": ["Rosetta_2005", "Rosetta_2007", "MESSENGER", "Stardust"]
        },
        "modern_od_suppressed": {
            "definition": "modern tracking + empirical accelerations",
            "examples": ["Juno"]
        }
    }
    
    @staticmethod
    def classify(flyby_data: Dict) -> str:
        """Classify a flyby into its physical regime."""
        altitude = flyby_data.get("altitude_km")
        velocity = flyby_data.get("velocity_km_s")
        cos_asymmetry = flyby_data.get("cos_asymmetry")
        
        # Validate required data
        if altitude is None or velocity is None or cos_asymmetry is None:
            raise ValueError(f"Missing required data for classification: altitude={altitude}, velocity={velocity}, cos_asymmetry={cos_asymmetry}")
        
        # Check high-velocity anti-aligned cancellation regime first
        if velocity > 16.0 and cos_asymmetry < 0:
            return "high_velocity_anti_aligned_cancellation"
        
        # Check low-altitude gradient-dominated
        if altitude < 800 and cos_asymmetry > 0:
            return "low_altitude_gradient"
        
        # Check mid-altitude enhancement
        if 800 <= altitude <= 1500 and cos_asymmetry > 0:
            return "mid_altitude_enhancement"
        
        # Default to high-altitude suppressed
        return "high_altitude_suppressed"


class ComponentHierarchicalModel:
    """Component-level hierarchical Bayesian model."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.logger = StepLogger("step_032_component_hierarchical", project_root)
        
    def load_component_data(self) -> Dict:
        """Load component-level data from step_004."""
        step_004_file = self.project_root / "results/step004_tep_predictions.json"
        
        try:
            with open(step_004_file, 'r') as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load flyby data: {e}")
            return None
        
        # Extract component data
        components = {}
        for name, flyby in data["predictions"].items():
            if flyby["observed"]["dv_obs_mm_s"] == 0:
                continue
            
            # Skip if sigma is None
            if flyby["observed"].get("sigma_mm_s") is None:
                self.logger.warning(f"Skipping {name}: sigma_mm_s is None")
                continue
            
            components[name] = {
                "dv_grad_mm_s": flyby["tep_predictions"]["dv_grad_mm_s"],
                "dv_disf_mm_s": flyby["tep_predictions"]["dv_disf_mm_s"],
                "dv_total_mm_s": flyby["tep_predictions"]["dv_tep_mm_s"],
                "dv_obs_mm_s": flyby["observed"]["dv_obs_mm_s"],
                "sigma_mm_s": flyby["observed"]["sigma_mm_s"],
                "altitude_km": flyby["perigee"]["altitude_km"],
                "velocity_km_s": flyby["perigee"]["velocity_km_s"],
                "cos_asymmetry": flyby["geometry"]["cos_dec_asymmetry"]
            }
        
        return components
    
    def fit_component_model(self, components: Dict) -> Dict:
        """
        Fit component-level hierarchical model.
        
        Uses simplified approach (full Bayesian with PyMC would be better but requires more setup).
        This provides initial estimates for the component scaling parameters.
        """
        self.logger.section("COMPONENT-LEVEL HIERARCHICAL MODEL")
        
        # Classify regimes
        regimes = {}
        for name, data in components.items():
            regimes[name] = RegimeClassifier.classify(data)
        
        # Simple least-squares fit for component scaling
        # Δv_obs = a_grad * Δv_grad + a_disf * Δv_disf + ε
        
        grad_values = []
        disf_values = []
        obs_values = []
        sigma_values = []
        
        for name, data in components.items():
            grad_values.append(data["dv_grad_mm_s"])
            disf_values.append(data["dv_disf_mm_s"])
            obs_values.append(data["dv_obs_mm_s"])
            sigma_values.append(data["sigma_mm_s"])
        
        grad_values = np.array(grad_values)
        disf_values = np.array(disf_values)
        obs_values = np.array(obs_values)
        sigma_values = np.array(sigma_values)
        
        # Weighted least squares
        X = np.column_stack([grad_values, disf_values])
        W = np.diag(1.0 / sigma_values**2)
        
        # Solve: (X^T W X) β = X^T W y
        XTWX = X.T @ W @ X
        XTWy = X.T @ W @ obs_values
        
        try:
            beta = np.linalg.solve(XTWX, XTWy)
        except np.linalg.LinAlgError:
            self.logger.warning("Singular matrix, using pseudo-inverse")
            beta = np.linalg.pinv(XTWX) @ XTWy
        
        a_grad_fit, a_disf_fit = beta
        
        # Compute predictions and residuals
        pred_values = a_grad_fit * grad_values + a_disf_fit * disf_values
        residuals = obs_values - pred_values
        
        # Compute regime-level residuals
        regime_residuals = {}
        for regime_name in RegimeClassifier.REGIMES.keys():
            regime_residuals[regime_name] = []
        
        for name, regime in regimes.items():
            if regime in regime_residuals:
                regime_residuals[regime].append(residuals[list(components.keys()).index(name)])
        
        # Compute regime standard deviations
        regime_std = {}
        for regime, res_list in regime_residuals.items():
            if len(res_list) > 1:
                regime_std[regime] = float(np.std(res_list))
            elif len(res_list) == 1:
                regime_std[regime] = 0.0
            else:
                regime_std[regime] = None
        
        # Normalize to baseline (low_altitude_gradient as reference)
        baseline_std = regime_std.get("low_altitude_gradient", 0.1)
        
        results = {
            "component_scaling": {
                "a_grad": float(a_grad_fit),
                "a_disf": float(a_disf_fit),
                "a_grad_normalized": float(a_grad_fit / BETA_INITIAL),
                "a_disf_normalized": float(a_disf_fit / DISFORMAL_COUPLING_STRENGTH)
            },
            "regime_classification": regimes,
            "regime_residuals": {
                regime: {
                    "std": std,
                    "relative_to_baseline": std / baseline_std if baseline_std > 0 and std is not None else None
                }
                for regime, std in regime_std.items()
            },
            "predictions": {},
            "residuals": {}
        }
        
        for name, data in components.items():
            idx = list(components.keys()).index(name)
            results["predictions"][name] = {
                "dv_pred_mm_s": float(pred_values[idx]),
                "dv_grad_contribution": float(a_grad_fit * data["dv_grad_mm_s"]),
                "dv_disf_contribution": float(a_disf_fit * data["dv_disf_mm_s"]),
                "regime": regimes[name]
            }
            results["residuals"][name] = {
                "residual_mm_s": float(residuals[idx]),
                "regime": regimes[name]
            }
        
        # Special Cassini diagnostics
        if "Cassini" in components:
            cassini = components["Cassini"]
            results["cassini_cancellation_diagnostics"] = {
                "dv_grad_mm_s": cassini["dv_grad_mm_s"],
                "dv_disf_mm_s": cassini["dv_disf_mm_s"],
                "dv_total_mm_s": cassini["dv_total_mm_s"],
                "dv_obs_mm_s": cassini["dv_obs_mm_s"],
                "cancellation_ratio": abs(cassini["dv_grad_mm_s"] / cassini["dv_disf_mm_s"]),
                "regime": regimes["Cassini"],
                "interpretation": "Cassini operates in cancellation regime where negative gradient and positive disformal terms partially cancel"
            }
        
        return results
    
    def save_results(self, results: Dict):
        """Save component hierarchical model results."""
        output_dir = self.project_root / "results"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save component hierarchical results
        output_file = output_dir / "step032_component_hierarchical.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Save Cassini-specific diagnostics
        if "cassini_cancellation_diagnostics" in results:
            cassini_file = output_dir / "step032_cassini_cancellation_diagnostics.json"
            with open(cassini_file, 'w') as f:
                json.dump(results["cassini_cancellation_diagnostics"], f, indent=2)
        
        # Save regime random effects summary
        regime_file = output_dir / "step032_regime_random_effects.json"
        regime_data = {
            "regime_residuals": results["regime_residuals"],
            "regime_classification": results["regime_classification"]
        }
        with open(regime_file, 'w') as f:
            json.dump(regime_data, f, indent=2)
        
        self.logger.info(f"Component hierarchical results saved to {output_file}")
        self.logger.info(f"Cassini diagnostics saved to {cassini_file}")
        self.logger.info(f"Regime random effects saved to {regime_file}")
    
    def run(self):
        """Run component-level hierarchical model."""
        self.logger.header("STEP 032: COMPONENT-LEVEL HIERARCHICAL MODEL")
        self.logger.info("Implementing component-level inference to address Cassini cancellation regime")
        
        # Load component data
        components = self.load_component_data()
        self.logger.info(f"Loaded {len(components)} flybys with component data")
        
        # Fit component model
        results = self.fit_component_model(components)
        
        # Log results
        self.logger.subsection("COMPONENT SCALING PARAMETERS")
        self.logger.info(f"a_grad = {results['component_scaling']['a_grad']:.4f}")
        self.logger.info(f"a_disf = {results['component_scaling']['a_disf']:.4f}")
        self.logger.info(f"a_grad_normalized = {results['component_scaling']['a_grad_normalized']:.4f}")
        self.logger.info(f"a_disf_normalized = {results['component_scaling']['a_disf_normalized']:.4f}")
        
        self.logger.subsection("REGIME RESIDUALS")
        for regime, data in results["regime_residuals"].items():
            if data["std"] is not None:
                self.logger.info(f"{regime}: σ = {data['std']:.4f} mm/s")
                if data["relative_to_baseline"] is not None:
                    self.logger.info(f"  Relative to baseline: {data['relative_to_baseline']:.2f}")
        
        if "cassini_cancellation_diagnostics" in results:
            self.logger.subsection("CASSINI CANCELLATION DIAGNOSTICS")
            diag = results["cassini_cancellation_diagnostics"]
            self.logger.info(f"Gradient term: {diag['dv_grad_mm_s']:.4f} mm/s (negative)")
            self.logger.info(f"Disformal term: {diag['dv_disf_mm_s']:.4f} mm/s (positive)")
            self.logger.info(f"Total prediction: {diag['dv_total_mm_s']:.4f} mm/s")
            self.logger.info(f"Observed: {diag['dv_obs_mm_s']:.4f} mm/s")
            self.logger.info(f"Cancellation ratio: {diag['cancellation_ratio']:.4f}")
            self.logger.info(f"Interpretation: {diag['interpretation']}")
        
        # Save results
        self.save_results(results)
        
        self.logger.log_step_summary(0, "SUCCESS")
        return 0


def main():
    """Main entry point."""
    model = ComponentHierarchicalModel(PROJECT_ROOT)
    return model.run()


if __name__ == "__main__":
    sys.exit(main())
