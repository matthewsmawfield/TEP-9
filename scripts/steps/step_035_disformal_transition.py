#!/usr/bin/env python3
"""
Step 035: Disformal Transition Criterion

This module implements the disformal transition criterion Ξ for classifying flybys
into conformal-dominated, mixed, or disformal-dominated regimes.

REVISED CRITERION (Yogyakarta v0.1):
    Ξ_i = (v_i / v_trans)^p × |cos(δ_in) - cos(δ_out)| × (|∇φ_i| / |∇φ_⊕|)^q × sgn(cos(δ_in) - cos(δ_out))

where:
- v_trans: Transition velocity (16.0 km/s)
- v_i: Flyby perigee velocity
- p: Velocity exponent (default 2)
- q: Gradient exponent (default 1)
- ∇φ_⊕: Temporal Shear at Earth's surface
- ∇φ_i: Temporal Shear at flyby altitude
- sgn: Sign of asymmetry (positive for aligned, negative for anti-aligned)

Classification (by |Ξ|):
- |Ξ| < 0.05: Conformal-dominated
- 0.05 ≤ |Ξ| ≤ 0.10: Mixed
- |Ξ| > 0.10: Disformal-dominated

The sign of Ξ indicates aligned (positive) vs anti-aligned (negative) disformal response.
Cassini falls into mixed regime with negative sign (anti-aligned).
"""

import numpy as np
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import RHO_T, LAMBDA_TEP_M, M_PL_GEV, DISFORMAL_VELOCITY_THRESHOLD_KM_S


class DisformalTransitionCriterion:
    """Disformal transition criterion analysis."""
    
    def __init__(self):
        self.logger = StepLogger("step_035_disformal_transition", PROJECT_ROOT)
        
        # Transition velocity (km/s) - from physics.py
        self.v_trans = DISFORMAL_VELOCITY_THRESHOLD_KM_S
        self.v_trans_uncertainty = 1.0  # km/s
        self.p = 2  # Velocity exponent
        self.p_uncertainty = 0.1  # Velocity exponent uncertainty
        self.q = 1  # Gradient exponent
        self.q_uncertainty = 0.1  # Gradient exponent uncertainty
        
        self.logger.info(f"Transition velocity: v_trans = {self.v_trans} ± {self.v_trans_uncertainty} km/s (±{self.v_trans_uncertainty/self.v_trans*100:.0f}% uncertainty)")
        self.logger.info(f"Velocity exponent: p = {self.p} ± {self.p_uncertainty} (±{self.p_uncertainty/self.p*100:.0f}% uncertainty)")
        self.logger.info(f"Gradient exponent: q = {self.q} ± {self.q_uncertainty} (±{self.q_uncertainty/self.q*100:.0f}% uncertainty)")
        self.logger.warning("CRITICAL: v_trans, p, q are HEURISTIC parameters with significant uncertainty")     
        # Temporal Shear at Earth's surface (GeV/m)
        # From first-principles field equation
        self.grad_phi_earth = 2.4e4 / (6371.0 * 1000)  # GeV/m
        
    def load_flyby_data(self):
        """Load flyby data from step005 fitting results."""
        results_file = PROJECT_ROOT / "results" / "step005_fitting_results.json"
        
        if not results_file.exists():
            self.logger.error("Fitting results not found. Run step005 first.")
            return None
        
        try:
            with open(results_file, encoding="utf-8") as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load fitting results: {e}")
            return None
        
        return data
    
    def compute_temporal_shear_at_altitude(self, altitude_km):
        """
        Compute Temporal Shear ∇φ at a given altitude.
        
        From the field equation φ(ρ) = Λ [n Λ^(n+4) M_Pl / (β ρ)]^(1/(n+1)),
        the gradient scales with density as ∇φ ∝ ρ^(-4/3) for n=3.
        
        Using exponential atmosphere: ρ(h) = ρ_0 exp(-h/H) with H ≈ 8.5 km
        """
        # Scale height of Earth's atmosphere (km)
        H = 8.5
        
        # Density at altitude
        rho_altitude = np.exp(-altitude_km / H)
        
        # Temporal Shear scales as ρ^(-4/3)
        grad_phi_altitude = self.grad_phi_earth * rho_altitude ** (-4/3)
        
        return grad_phi_altitude
    
    def compute_asymmetry_factor(self, cos_dec_asymmetry):
        """
        Compute trajectory asymmetry factor G_asym.
        
        From cos(δ) where δ is declination of perigee asymptote.
        """
        return abs(cos_dec_asymmetry) if cos_dec_asymmetry is not None else 0.5
    
    def compute_xi(self, flyby_data):
        """
        Compute disformal transition criterion Ξ for a flyby.
        
        REVISED FORMULA (velocity-activated):
        Ξ = (v / v_trans)^p × |cos(δ_in) - cos(δ_out)| × (|∇φ| / |∇φ_⊕|)^q × sgn(cos(δ_in) - cos(δ_out))
        
        where:
        - p = 2 (velocity exponent)
        - q = 1 (gradient exponent)
        - sgn indicates aligned (positive) vs anti-aligned (negative) disformal response
        
        This makes high velocity INCREASE Ξ (opposite to previous inverted logic).
        """
        # Extract flyby parameters
        altitude_km = flyby_data["perigee"]["altitude_km"]
        velocity_km_s = flyby_data["perigee"]["velocity_km_s"]
        
        # Get asymmetry
        cos_dec_asymmetry = flyby_data.get("cos_dec_asymmetry")
        if cos_dec_asymmetry is None:
            raise ValueError(f"Missing cos_dec_asymmetry for flyby data: {flyby_data.get('name', 'unknown')}")
        asymmetry = cos_dec_asymmetry
        
        # Velocity-activated term (high velocity increases Ξ)
        velocity_term = (velocity_km_s / self.v_trans) ** self.p
        
        # Asymmetry magnitude
        asymmetry_magnitude = abs(asymmetry)
        
        # Gradient ratio (weaker at high altitude)
        grad_phi_altitude = self.compute_temporal_shear_at_altitude(altitude_km)
        gradient_ratio = abs(grad_phi_altitude / self.grad_phi_earth) ** self.q
        
        # Sign of asymmetry (positive for aligned, negative for anti-aligned)
        asymmetry_sign = np.sign(asymmetry) if asymmetry != 0 else 1.0
        
        # Compute Ξ with sign
        xi = velocity_term * asymmetry_magnitude * gradient_ratio * asymmetry_sign
        
        # Compute uncertainty on Ξ (propagate from v_trans, p, q)
        v_trans_rel_uncertainty = self.v_trans_uncertainty / self.v_trans
        p_rel_uncertainty = self.p_uncertainty / self.p
        q_rel_uncertainty = self.q_uncertainty / self.q
        
        # Conservative uncertainty estimate (quadrature sum)
        xi_rel_uncertainty = np.sqrt(
            (self.p * v_trans_rel_uncertainty)**2 +
            (p_rel_uncertainty * np.log(velocity_km_s / self.v_trans))**2 +
            (q_rel_uncertainty * np.log(gradient_ratio))**2
        )
        xi_uncertainty = abs(xi) * xi_rel_uncertainty
        
        return xi, xi_uncertainty
    
    def classify_regime(self, xi):
        """
        Classify flyby regime based on |Ξ|.
        
        Classification thresholds (by |Ξ|):
        - |Ξ| < 0.05: Conformal-dominated
        - 0.05 ≤ |Ξ| ≤ 0.10: Mixed
        - |Ξ| > 0.10: Disformal-dominated
        
        The sign of Ξ indicates:
        - Positive: Aligned disformal response
        - Negative: Anti-aligned disformal response (sign reversal regime)
        """
        xi_abs = abs(xi)
        
        if xi_abs < 0.05:
            return "conformal-dominated"
        elif xi_abs <= 0.10:
            return "mixed"
        else:
            return "disformal-dominated"
    
    def run(self):
        """Run disformal transition criterion analysis."""
        self.logger.section("STEP 035: DISFORMAL TRANSITION CRITERION")
        
        # Load flyby data
        flyby_data = self.load_flyby_data()
        if flyby_data is None:
            return None
        
        # Compute Ξ for each flyby
        xi_results = {}
        
        for name, entry in flyby_data["individual_fits"].items():
            if entry["observed"]["dv_obs_mm_s"] == 0:
                continue
            
            xi, xi_uncertainty = self.compute_xi(entry)
            regime = self.classify_regime(xi)
            
            xi_results[name] = {
                'xi': float(xi),
                'xi_uncertainty': float(xi_uncertainty),
                'regime': regime,
                'altitude_km': float(entry["perigee"]["altitude_km"]),
                'velocity_km_s': float(entry["perigee"]["velocity_km_s"]),
                'cos_asymmetry': float(entry.get("cos_dec_asymmetry", 0.5))
            }
            
            self.logger.info(f"{name}: Ξ = {xi:.3f} ({regime})")
        
        # Check Cassini classification
        if 'Cassini' in xi_results:
            cassini_regime = xi_results['Cassini']['regime']
            self.logger.subsection("CASSINI CLASSIFICATION")
            self.logger.info(f"Cassini regime: {cassini_regime}")
            
            if cassini_regime in ['mixed', 'disformal-dominated']:
                self.logger.info("✓ Cassini correctly classified as mixed/disformal")
            else:
                self.logger.warning("⚠ Cassini not in expected mixed/disformal regime")
        
        # Save results
        output = {
            'model_type': 'disformal_transition_criterion',
            'xi_results': xi_results,
            'v_trans_km_s': self.v_trans,
            'v_trans_uncertainty_km_s': self.v_trans_uncertainty,
            'p': self.p,
            'p_uncertainty': self.p_uncertainty,
            'q': self.q,
            'q_uncertainty': self.q_uncertainty,
            'uncertainty_metadata': {
                'v_trans_uncertainty_fraction': self.v_trans_uncertainty / self.v_trans,
                'p_uncertainty_fraction': self.p_uncertainty / self.p,
                'q_uncertainty_fraction': self.q_uncertainty / self.q,
                'data_source': 'heuristic_disformal_transition_model',
                'calibration_status': 'uncalibrated',
                'recommended_action': 'Calibrate v_trans against theoretical disformal coupling calculations when available'
            }
        }
        
        output_path = PROJECT_ROOT / "results" / "step035_disformal_transition.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(output, f, indent=2)
        
        self.logger.info(f"Disformal transition results saved to {output_path}")
        
        return output


if __name__ == "__main__":
    criterion = DisformalTransitionCriterion()
    criterion.run()
