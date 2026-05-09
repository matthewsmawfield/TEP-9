#!/usr/bin/env python3
"""
Step 042: Corrected Systematic Uncertainty Calculation

This module fixes the systematic uncertainty contradiction identified in the review.

Issue from Review:
-----------------
The manuscript claims both:
- "0.03% systematic uncertainty contribution"
- "102% total relative uncertainty"

These are contradictory. The issue is a confusion between:
1. Fractional contribution of measurement uncertainty to total variance
2. Total systematic uncertainty as percentage of β

Fix:
----
1. Clearly distinguish between variance contributions and total relative uncertainty
2. Report both sources correctly without contradiction
3. Include heterogeneity (I²) in the total uncertainty budget
4. Ensure all percentages sum to 100% when appropriate

Corrected Calculation:
----------------------
Total uncertainty = sqrt(σ_statistical² + σ_systematic² + σ_heterogeneity²)

Where:
- σ_statistical: measurement uncertainty (weighted average)
- σ_systematic: from theory parameters (suppression, relaxation length, etc.)
- σ_heterogeneity: from I² metric (geometry-dependent variation)

This resolves the contradiction by properly accounting for all sources.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


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


class CorrectedUncertaintyCalculator:
    """
    Corrected systematic uncertainty calculation that resolves contradictions.
    """
    
    def __init__(self, logger):
        self.logger = logger
    
    def load_fitting_results(self):
        """Load fitting results from step_005."""
        results_file = PROJECT_ROOT / 'results' / 'step005_fitting_results.json'
        
        if not results_file.exists():
            self.logger.error(f"Fitting results not found: {results_file}")
            return None
        
        with open(results_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        return data
    
    def calculate_corrected_uncertainty(self, data):
        """
        Calculate corrected uncertainty budget that resolves the contradiction.
        
        Returns:
            dict: Corrected uncertainty budget with clear distinction between
                  variance contributions and total relative uncertainty
        """
        self.logger.section("CORRECTED UNCERTAINTY CALCULATION")
        
        # Extract data
        individual_fits = data['individual_fits']
        overall_analysis = data['overall_analysis']
        beta_stats = overall_analysis['beta_statistics']
        heterogeneity = overall_analysis['heterogeneity_tests']
        
        # Statistical uncertainty (from measurement errors)
        weighted_mean = beta_stats['weighted_mean']
        weighted_unc = beta_stats['weighted_uncertainty']
        statistical_rel_unc = weighted_unc / weighted_mean if weighted_mean > 0 else 0
        
        # Heterogeneity uncertainty (from I²)
        I2 = heterogeneity['I_squared_percent'] / 100  # Convert to fraction
        std = beta_stats['std']
        mean = beta_stats['mean']
        heterogeneity_rel_unc = std / mean if mean > 0 else 0
        
        # Systematic uncertainty sources (from theory parameters)
        # These are the same as before, but we'll report them correctly
        trajectory_rel_unc = 0.01  # 1% trajectory uncertainty
        suppression_rel_unc = 0.05  # 5% characteristic suppression uncertainty
        multipole_rel_unc = 0.001  # 0.1% multipole uncertainty
        relaxation_length_rel_unc = 0.15  # 15% relaxation length uncertainty
        
        # Combine systematic uncertainties (RSS)
        systematic_rel_unc = np.sqrt(
            trajectory_rel_unc**2 + 
            suppression_rel_unc**2 + 
            multipole_rel_unc**2 + 
            relaxation_length_rel_unc**2
        )
        
        # Total uncertainty (RSS of all sources)
        total_rel_unc = np.sqrt(
            statistical_rel_unc**2 + 
            systematic_rel_unc**2 + 
            heterogeneity_rel_unc**2
        )
        
        # Variance contributions (fraction of total variance)
        total_variance = total_rel_unc**2
        stat_variance_frac = statistical_rel_unc**2 / total_variance if total_variance > 0 else 0
        sys_variance_frac = systematic_rel_unc**2 / total_variance if total_variance > 0 else 0
        het_variance_frac = heterogeneity_rel_unc**2 / total_variance if total_variance > 0 else 0
        
        self.logger.info("Corrected Uncertainty Budget:")
        self.logger.info(f"  Statistical uncertainty: {statistical_rel_unc:.2%}")
        self.logger.info(f"  Systematic uncertainty: {systematic_rel_unc:.2%}")
        self.logger.info(f"  Heterogeneity uncertainty: {heterogeneity_rel_unc:.2%}")
        self.logger.info(f"  Total uncertainty: {total_rel_unc:.2%}")
        self.logger.info("")
        self.logger.info("Variance Contributions:")
        self.logger.info(f"  Statistical: {stat_variance_frac:.1%}")
        self.logger.info(f"  Systematic: {sys_variance_frac:.1%}")
        self.logger.info(f"  Heterogeneity: {het_variance_frac:.1%}")
        self.logger.info(f"  Total: {stat_variance_frac + sys_variance_frac + het_variance_frac:.1%}")
        
        results = {
            'corrected_uncertainty_budget': {
                'statistical_relative_uncertainty': float(statistical_rel_unc),
                'systematic_relative_uncertainty': float(systematic_rel_unc),
                'heterogeneity_relative_uncertainty': float(heterogeneity_rel_unc),
                'total_relative_uncertainty': float(total_rel_unc),
                'variance_contributions': {
                    'statistical': float(stat_variance_frac),
                    'systematic': float(sys_variance_frac),
                    'heterogeneity': float(het_variance_frac),
                    'total': float(stat_variance_frac + sys_variance_frac + het_variance_frac)
                }
            },
            'systematic_breakdown': {
                'trajectory_uncertainty': float(trajectory_rel_unc),
                'characteristic_suppression_uncertainty': float(suppression_rel_unc),
                'multipole_coefficient_uncertainty': float(multipole_rel_unc),
                'relaxation_length_uncertainty': float(relaxation_length_rel_unc),
                'total_systematic': float(systematic_rel_unc)
            },
            'heterogeneity_metrics': {
                'I_squared_percent': float(heterogeneity['I_squared_percent']),
                'heterogeneity_relative_uncertainty': float(heterogeneity_rel_unc)
            },
            'resolution': {
                'issue': 'Contradiction between 0.03% and 102% claims',
                'root_cause': 'Confusion between variance contributions and total relative uncertainty',
                'fix': 'Clearly distinguish between fractional variance contributions and total relative uncertainty',
                'interpretation': f'Total uncertainty is {total_rel_unc:.1%}, with heterogeneity ({het_variance_frac:.1%}) and systematic ({sys_variance_frac:.1%}) as dominant sources'
            }
        }
        
        return results


def main():
    """Execute corrected uncertainty calculation."""
    logger = StepLogger("step_042_corrected_uncertainty", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("STEP 042: CORRECTED SYSTEMATIC UNCERTAINTY CALCULATION")
    
    logger.section("ISSUE IDENTIFICATION")
    logger.info("Original manuscript claims:")
    logger.info("  - '0.03% systematic uncertainty contribution'")
    logger.info("  - '102% total relative uncertainty'")
    logger.info("")
    logger.info("These are contradictory.")
    logger.info("Root cause: Confusion between variance contributions and total relative uncertainty")
    
    calculator = CorrectedUncertaintyCalculator(logger)
    
    # Load fitting results
    data = calculator.load_fitting_results()
    if data is None:
        logger.error("Failed to load fitting results")
        return 1
    
    # Calculate corrected uncertainty
    results = calculator.calculate_corrected_uncertainty(data)
    
    # Save results
    output_file = PROJECT_ROOT / 'results' / 'step042_corrected_uncertainty.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(convert_to_native_types(results), f, indent=2)
    
    logger.success(f"Corrected uncertainty calculation complete. Saved to {output_file}")
    logger.add_output_file(output_file, "Corrected uncertainty budget")
    
    duration = time.time() - start_time
    logger.log_step_summary(duration, "SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
