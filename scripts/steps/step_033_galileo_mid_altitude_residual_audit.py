#!/usr/bin/env python3
"""
Step 033: Galileo 1990 Mid-Altitude Residual Audit
====================================================

This module investigates the systematic underprediction of Galileo 1990,
testing four hypotheses for the residual budget.

Galileo 1990 shows 48% underprediction (observed 3.92 mm/s, predicted 2.02 mm/s).
This step tests hypotheses without adding arbitrary corrections.

Hypotheses:
1. Plasma under-modelled - compare F10.7/Kp/ionospheric density against NEAR and Rosetta
2. Altitude transition shape wrong - compare exponential, logistic, broken-exponential
3. Missing J3/J4/inclination term - recompute with full latitude-dependent harmonics
4. OD survival factor wrong - run Galileo-specific OD absorption simulation

Output: galileo_1990_residual_budget.json
"""

import json
import numpy as np
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger
from scripts.utils.physics import LAMBDA_TEP_M, R_EARTH


class GalileoResidualAuditor:
    """Audits Galileo 1990 residual budget."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.logger = StepLogger("step_033_galileo_residual_audit", project_root)
        
    def load_flyby_data(self):
        """Load flyby data from step004 predictions."""
        step_004_file = self.project_root / "results/step004_tep_predictions.json"
        
        try:
            with open(step_004_file, 'r') as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load flyby data: {e}")
            return None
        
        return data.get("predictions", {})
    
    def hypothesis_1_plasma_undermodelled(self, galileo_data):
        """
        Hypothesis 1: Plasma under-modelled.
        
        Compare F10.7/Kp/ionospheric density against NEAR and Rosetta.
        Galileo 1990 has plasma density 2020 cm³, which should produce non-trivial screening.
        """
        self.logger.subsection("HYPOTHESIS 1: PLASMA UNDER-MODELLED")
        
        # Load plasma data from step015
        step_015_file = self.project_root / "results/step015_plasma_modulation.json"
        
        if not step_015_file.exists():
            self.logger.warning("Plasma modulation data not found. Run step015 first.")
            return {"status": "no_data", "hypothesis": "plasma_undermodelled"}
        
        try:
            with open(step_015_file, 'r') as f:
                plasma_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load plasma modulation data: {e}")
            return {"status": "no_data", "hypothesis": "plasma_undermodelled"}
        
        # Compare plasma densities
        galileo_plasma = plasma_data.get("Galileo_1990", {})
        if not galileo_plasma:
            self.logger.warning("Galileo plasma data not found in step034 output.")
        galileo_density = galileo_plasma.get("plasma_density_cm3")
        galileo_screening = galileo_plasma.get("plasma_screening_factor")
        galileo_sign = galileo_plasma.get("plasma_sign_factor")
        
        near_plasma = plasma_data.get("NEAR", {})
        if not near_plasma:
            self.logger.warning("NEAR plasma data not found in step034 output.")
        near_density = near_plasma.get("plasma_density_cm3")
        near_screening = near_plasma.get("plasma_screening_factor")
        
        rosetta_plasma = plasma_data.get("Rosetta_2005", {})
        if not rosetta_plasma:
            self.logger.warning("Rosetta plasma data not found in step034 output.")
        rosetta_density = rosetta_plasma.get("plasma_density_cm3")
        rosetta_screening = rosetta_plasma.get("plasma_screening_factor")
        
        # Use None if data not available, not defaults
        galileo_screening = galileo_screening if galileo_screening is not None else None
        galileo_sign = galileo_sign if galileo_sign is not None else None
        near_screening = near_screening if near_screening is not None else None
        rosetta_screening = rosetta_screening if rosetta_screening is not None else None
        
        self.logger.info(f"Galileo 1990 plasma density: {galileo_density:.1f} cm³")
        self.logger.info(f"Galileo screening factor: {galileo_screening:.4f}")
        self.logger.info(f"Galileo sign factor: {galileo_sign:.4f}")
        self.logger.info(f"NEAR plasma density: {near_density:.1f} cm³")
        self.logger.info(f"NEAR screening factor: {near_screening:.4f}")
        self.logger.info(f"Rosetta plasma density: {rosetta_density:.1f} cm³")
        self.logger.info(f"Rosetta screening factor: {rosetta_screening:.4f}")
        
        # Hypothesis assessment
        if galileo_density > 1000 and galileo_screening < 0.98:
            # Significant plasma density with non-trivial screening
            assessment = "plasma_contributes_significantly"
            explanation = f"Galileo has high plasma density ({galileo_density:.0f} cm³) with screening factor {galileo_screening:.4f}. Current F_plasma=1.0 ignores this effect."
        else:
            assessment = "plasma_not_significant"
            explanation = f"Galileo plasma density ({galileo_density:.0f} cm³) or screening ({galileo_screening:.4f}) is insufficient to explain residual."
        
        return {
            "status": "complete",
            "hypothesis": "plasma_undermodelled",
            "assessment": assessment,
            "explanation": explanation,
            "galileo_plasma": galileo_plasma,
            "near_plasma": near_plasma,
            "rosetta_plasma": rosetta_plasma
        }
    
    def hypothesis_2_altitude_transition_shape(self, galileo_data):
        """
        Hypothesis 2: Altitude transition shape wrong.
        
        Test three transition shapes:
        - Exponential: g(h) = exp(-h/λ_T)
        - Logistic: g(h) = 1 / (1 + exp[(h-h_c)/w])
        - Broken-exponential: g(h) = exp(-h/λ₁) for h<h_b, A·exp(-(h-h_b)/λ₂) for h≥h_b
        
        Galileo 1990 at 972 km may be in a shoulder region.
        """
        self.logger.subsection("HYPOTHESIS 2: ALTITUDE TRANSITION SHAPE")
        
        altitude = galileo_data["perigee"]["altitude_km"]
        lambda_T = LAMBDA_TEP_M / 1000  # Convert to km
        
        # Exponential (current model)
        g_exp = np.exp(-altitude / lambda_T)
        
        # Logistic (h_c = 800 km, w = 200 km)
        h_c = 800.0
        w = 200.0
        g_log = 1.0 / (1.0 + np.exp((altitude - h_c) / w))
        
        # Broken-exponential (h_b = 800 km, λ₁ = 4000 km, λ₂ = 6000 km)
        h_b = 800.0
        lambda_1 = 4000.0
        lambda_2 = 6000.0
        A = np.exp(-h_b / lambda_1) / np.exp(-h_b / lambda_2)  # Continuity factor
        if altitude < h_b:
            g_broken = np.exp(-altitude / lambda_1)
        else:
            g_broken = A * np.exp(-(altitude - h_b) / lambda_2)
        
        # Compare to NEAR (568 km) and Rosetta (1968 km)
        near_alt = 568.0
        rosetta_alt = 1968.0
        
        g_exp_near = np.exp(-near_alt / lambda_T)
        g_log_near = 1.0 / (1.0 + np.exp((near_alt - h_c) / w))
        g_broken_near = np.exp(-near_alt / lambda_1)  # NEAR is below h_b
        
        g_exp_rosetta = np.exp(-rosetta_alt / lambda_T)
        g_log_rosetta = 1.0 / (1.0 + np.exp((rosetta_alt - h_c) / w))
        g_broken_rosetta = A * np.exp(-(rosetta_alt - h_b) / lambda_2)  # Rosetta is above h_b
        
        self.logger.info(f"Galileo altitude: {altitude} km")
        self.logger.info(f"Exponential g(h): {g_exp:.4f}")
        self.logger.info(f"Logistic g(h): {g_log:.4f}")
        self.logger.info(f"Broken-exponential g(h): {g_broken:.4f}")
        self.logger.info("")
        self.logger.info(f"NEAR ({near_alt} km): exp={g_exp_near:.4f}, log={g_log_near:.4f}, broken={g_broken_near:.4f}")
        self.logger.info(f"Rosetta ({rosetta_alt} km): exp={g_exp_rosetta:.4f}, log={g_log_rosetta:.4f}, broken={g_broken_rosetta:.4f}")
        
        # Hypothesis assessment
        # Check if broken-exponential improves Galileo prediction relative to NEAR/Rosetta
        exp_ratio = g_exp / g_exp_near
        log_ratio = g_log / g_log_near
        broken_ratio = g_broken / g_broken_near
        
        rosetta_exp_ratio = g_exp_rosetta / g_exp_near
        rosetta_log_ratio = g_log_rosetta / g_log_near
        rosetta_broken_ratio = g_broken_rosetta / g_broken_near
        
        # Galileo should have intermediate suppression between NEAR and Rosetta
        # Current prediction: Galileo underpredicted, suggesting too much suppression
        if broken_ratio < exp_ratio:
            assessment = "broken_exponential_improves"
            explanation = f"Broken-exponential gives less suppression at {altitude} km (ratio {broken_ratio:.3f}) vs exponential ({exp_ratio:.3f}). This could reduce Galileo underprediction."
        else:
            assessment = "exponential_adequate"
            explanation = f"Exponential transition (ratio {exp_ratio:.3f}) is comparable to or better than alternatives."
        
        return {
            "status": "complete",
            "hypothesis": "altitude_transition_shape",
            "assessment": assessment,
            "explanation": explanation,
            "galileo_altitude_km": altitude,
            "models": {
                "exponential": {"g": float(g_exp), "ratio_to_near": float(exp_ratio)},
                "logistic": {"g": float(g_log), "ratio_to_near": float(log_ratio)},
                "broken_exponential": {"g": float(g_broken), "ratio_to_near": float(broken_ratio)}
            },
            "near_reference": {
                "altitude_km": near_alt,
                "exponential": float(g_exp_near),
                "logistic": float(g_log_near),
                "broken_exponential": float(g_broken_near)
            },
            "rosetta_reference": {
                "altitude_km": rosetta_alt,
                "exponential": float(g_exp_rosetta),
                "logistic": float(g_log_rosetta),
                "broken_exponential": float(g_broken_rosetta)
            }
        }
    
    def hypothesis_3_missing_harmonics(self, galileo_data):
        """
        Hypothesis 3: Missing J3/J4/inclination term.
        
        Recompute with full latitude-dependent harmonics.
        Current model uses only J2 (equatorial bulge).
        """
        self.logger.subsection("HYPOTHESIS 3: MISSING J3/J4/INCLINATION TERMS")
        
        # Extract geometry parameters
        altitude = galileo_data["perigee"]["altitude_km"]
        perigee_lat_deg = galileo_data.get("perigee_latitude_deg", 0.0)
        perigee_lat_rad = np.radians(perigee_lat_deg)
        
        # Standard Earth gravity field coefficients
        J2 = 1.08263e-3  # Equatorial bulge (already in model)
        J3 = -2.54e-6    # Pear-shaped Earth
        J4 = 1.62e-6     # Higher-order term
        
        # Compute latitude-dependent correction factors
        # J3 contribution: proportional to sin(latitude) * (3*sin^2(latitude) - 1)
        sin_lat = np.sin(perigee_lat_rad)
        j3_factor = sin_lat * (3 * sin_lat**2 - 1)
        
        # J4 contribution: proportional to (35*sin^4(latitude) - 30*sin^2(latitude) + 3)
        j4_factor = 35 * sin_lat**4 - 30 * sin_lat**2 + 3
        
        # Relative magnitudes compared to J2
        j3_relative = (J3 / J2) * j3_factor
        j4_relative = (J4 / J2) * j4_factor
        
        # Total correction factor
        total_correction = 1.0 + j3_relative + j4_relative
        
        # Current prediction
        dv_pred = galileo_data["tep_predictions"]["dv_tep_mm_s"]
        dv_obs = galileo_data["observed"]["dv_obs_mm_s"]
        underprediction = dv_obs - dv_pred
        
        # Corrected prediction with J3/J4
        dv_pred_corrected = dv_pred * total_correction
        underprediction_corrected = dv_obs - dv_pred_corrected
        
        self.logger.info(f"Galileo perigee altitude: {altitude} km")
        self.logger.info(f"Perigee latitude: {perigee_lat_deg:.1f}°")
        self.logger.info(f"J2 = {J2:.2e} (equatorial bulge, already in model)")
        self.logger.info(f"J3 = {J3:.2e} (pear-shaped Earth)")
        self.logger.info(f"J4 = {J4:.2e} (higher-order term)")
        self.logger.info(f"")
        self.logger.info(f"J3 correction factor: {j3_relative:.4f} ({j3_relative*100:.2f}%)")
        self.logger.info(f"J4 correction factor: {j4_relative:.4f} ({j4_relative*100:.2f}%)")
        self.logger.info(f"Total correction factor: {total_correction:.4f} ({(total_correction-1)*100:.2f}%)")
        self.logger.info(f"")
        self.logger.info(f"Current prediction: {dv_pred:.2f} mm/s")
        self.logger.info(f"Corrected prediction: {dv_pred_corrected:.2f} mm/s")
        self.logger.info(f"Observed: {dv_obs:.2f} mm/s")
        self.logger.info(f"")
        self.logger.info(f"Current underprediction: {underprediction:.2f} mm/s ({underprediction/dv_obs*100:.1f}%)")
        self.logger.info(f"Corrected underprediction: {underprediction_corrected:.2f} mm/s ({underprediction_corrected/dv_obs*100:.1f}%)")
        
        # Assessment
        if abs(total_correction - 1.0) < 0.01:
            assessment = "negligible_effect"
            explanation = f"J3/J4 corrections are small ({(total_correction-1)*100:.2f}%). Missing harmonics unlikely to explain underprediction."
        elif abs(total_correction - 1.0) < 0.05:
            assessment = "moderate_effect"
            explanation = f"J3/J4 corrections are moderate ({(total_correction-1)*100:.2f}%). May partially explain underprediction but not sufficient."
        else:
            assessment = "significant_effect"
            explanation = f"J3/J4 corrections are significant ({(total_correction-1)*100:.2f}%). Could contribute to underprediction."
        
        self.logger.info(f"")
        self.logger.info(f"Assessment: {assessment}")
        self.logger.info(f"Explanation: {explanation}")
        
        return {
            "status": "complete",
            "hypothesis": "missing_harmonics",
            "assessment": assessment,
            "explanation": explanation,
            "perigee_latitude_deg": float(perigee_lat_deg),
            "j3_relative": float(j3_relative),
            "j4_relative": float(j4_relative),
            "total_correction": float(total_correction),
            "dv_pred_current": float(dv_pred),
            "dv_pred_corrected": float(dv_pred_corrected),
            "underprediction_current": float(underprediction),
            "underprediction_corrected": float(underprediction_corrected)
        }
    
    def hypothesis_4_od_survival_factor(self, galileo_data):
        """
        Hypothesis 4: OD survival factor wrong.
        
        Run Galileo-specific OD absorption simulation.
        Current F_OD is simplified estimation.
        """
        self.logger.subsection("HYPOTHESIS 4: OD SURVIVAL FACTOR")
        
        # Load current F_OD from step035 (mission-specific OD absorption)
        step_035_file = self.project_root / "results/step035_mission_od_absorption.json"
        
        if not step_035_file.exists():
            self.logger.warning("OD absorption data not found. Run step035 first.")
            return {"status": "no_data", "hypothesis": "od_survival_factor"}
        
        try:
            with open(step_035_file, 'r') as f:
                od_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load OD absorption data: {e}")
            return {"status": "no_data", "hypothesis": "od_survival_factor"}
        
        galileo_od = od_data.get("Galileo_1990", {})
        
        if not galileo_od:
            self.logger.warning("Galileo OD data not found in step035 output.")
            return {"status": "no_data", "hypothesis": "od_survival_factor"}
        
        f_od = galileo_od.get("f_od")
        if f_od is None:
            self.logger.warning("Galileo F_OD value not found in step035 output.")
            return {"status": "no_data", "hypothesis": "od_survival_factor"}
        mission_config = galileo_od.get("mission_config", {})
        
        self.logger.info(f"Galileo F_OD: {f_od:.4f}")
        self.logger.info(f"Mission era: {mission_config.get('era', 'unknown')}")
        self.logger.info(f"Tracking arc: {mission_config.get('tracking_arc_days', 0)} days")
        self.logger.info(f"DSN coverage: {mission_config.get('dsn_coverage', 'unknown')}")
        self.logger.info("")
        
        # Assess impact on Galileo 1990 underprediction
        dv_pred = galileo_data["tep_predictions"]["dv_tep_mm_s"]
        dv_obs = galileo_data["observed"]["dv_obs_mm_s"]
        underprediction = dv_obs - dv_pred
        
        self.logger.info(f"Galileo predicted: {dv_pred:.2f} mm/s")
        self.logger.info(f"Galileo observed: {dv_obs:.2f} mm/s")
        self.logger.info(f"Underprediction: {underprediction:.2f} mm/s ({underprediction/dv_obs*100:.1f}%)")
        self.logger.info("")
        
        # OD survival factor assessment
        if f_od < 0.8:
            assessment = "significant_od_absorption"
            explanation = f"Galileo has low OD survival factor ({f_od:.3f}), indicating strong OD absorption. This could explain part of the underprediction if current F_OD=1.0 overestimates signal recovery."
        elif f_od < 0.95:
            assessment = "moderate_od_absorption"
            explanation = f"Galileo has moderate OD survival factor ({f_od:.3f}), indicating some OD absorption. Current F_OD=1.0 may overestimate signal recovery by {(1-f_od)*100:.1f}%."
        else:
            assessment = "minimal_od_absorption"
            explanation = f"Galileo has high OD survival factor ({f_od:.3f}), indicating minimal OD absorption. OD effects unlikely to explain underprediction."
        
        self.logger.info(f"Assessment: {assessment}")
        self.logger.info(f"Explanation: {explanation}")
        
        return {
            "status": "complete",
            "hypothesis": "od_survival_factor",
            "assessment": assessment,
            "explanation": explanation,
            "f_od": f_od,
            "mission_config": mission_config,
            "underprediction_mm_s": float(underprediction),
            "underprediction_fraction": float(underprediction / dv_obs)
        }
    
    def run(self):
        """Run Galileo 1990 residual audit."""
        self.logger.header("STEP 033: GALILEO 1990 MID-ALTITUDE RESIDUAL AUDIT")
        self.logger.info("Testing hypotheses for Galileo 1990 underprediction")
        self.logger.info("Observed: 3.92 mm/s, Predicted: 2.02 mm/s (48% underprediction)")
        
        # Load flyby data
        flyby_data = self.load_flyby_data()
        
        if "Galileo_1990" not in flyby_data:
            self.logger.error("Galileo 1990 data not found.")
            return 1
        
        galileo_data = flyby_data["Galileo_1990"]
        
        # Test hypotheses
        results = {
            "flyby": "Galileo_1990",
            "observed_dv_mm_s": galileo_data["observed"]["dv_obs_mm_s"],
            "predicted_dv_mm_s": galileo_data["tep_predictions"]["dv_tep_mm_s"],
            "residual_mm_s": galileo_data["observed"]["dv_obs_mm_s"] - galileo_data["tep_predictions"]["dv_tep_mm_s"],
            "underprediction_fraction": 1.0 - galileo_data["tep_predictions"]["dv_tep_mm_s"] / galileo_data["observed"]["dv_obs_mm_s"],
            "hypotheses": {}
        }
        
        # Hypothesis 1: Plasma
        results["hypotheses"]["plasma_undermodelled"] = self.hypothesis_1_plasma_undermodelled(galileo_data)
        
        # Hypothesis 2: Altitude transition
        results["hypotheses"]["altitude_transition_shape"] = self.hypothesis_2_altitude_transition_shape(galileo_data)
        
        # Hypothesis 3: Missing harmonics
        results["hypotheses"]["missing_harmonics"] = self.hypothesis_3_missing_harmonics(galileo_data)
        
        # Hypothesis 4: OD survival
        results["hypotheses"]["od_survival_factor"] = self.hypothesis_4_od_survival_factor(galileo_data)
        
        # Summary
        self.logger.subsection("SUMMARY")
        significant_hypotheses = []
        for name, result in results["hypotheses"].items():
            if result.get("assessment") in ["plasma_contributes_significantly", "broken_exponential_improves"]:
                significant_hypotheses.append(name)
        
        if significant_hypotheses:
            self.logger.info(f"Significant hypotheses: {', '.join(significant_hypotheses)}")
        else:
            self.logger.info("No single hypothesis fully explains the residual.")
            self.logger.info("Combination of effects likely required.")
        
        # Save results
        output_file = self.project_root / "results/galileo_1990_residual_budget.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        self.logger.info(f"Residual budget saved to {output_file}")
        
        self.logger.log_step_summary(0, "SUCCESS")
        return 0


def main():
    """Main entry point."""
    auditor = GalileoResidualAuditor(PROJECT_ROOT)
    return auditor.run()


if __name__ == "__main__":
    sys.exit(main())
