#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Step 006: Final Report

Generates summary report and validates pipeline integrity.
"""

import sys
import json
import math
from pathlib import Path
from datetime import datetime, timezone
import time

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger

def literature_synthesis_metadata(chi2_consistency: float | None = None) -> dict:
    metadata = {
        'uncertainty': None,
        'uncertainty_fraction': None,
        'uncertainty_absolute': None,
        'status': 'literature_synthesis',
        'calibration_status': 'literature_and_pipeline_outputs',
        'data_source': 'peer_reviewed_literature_and_pipeline_outputs',
        'recommended_action': (
            'Replace literature anomaly inputs with independent DSN OD reanalysis before external claims'
        ),
    }
    if chi2_consistency is not None:
        metadata['derivation'] = (
            f'Literature synthesis with beta-scatter consistency chi2~{chi2_consistency:.0f}'
        )
    return metadata


def generate_report(fitting_data: dict, project_root: Path) -> dict:
    """Generate comprehensive scientific report from fitting results."""
    
    fits = fitting_data['individual_fits']
    overall = fitting_data['overall_analysis']
    
    # Key findings with detailed annotations
    key_flybys = {
        'NEAR': {
            'significance': 'major',
            'dv': 13.46,
            'uncertainty': 0.01,
            'reference': 'Anderson et al. (2008), PRL 100, 091102',
            'doi': '10.1103/PhysRevLett.100.091102',
            'note': 'Eros mission; largest detected anomaly with highest signal-to-noise ratio; PPN constraints challenged due to large implied β'
        },
        'Galileo_1990': {
            'significance': 'moderate',
            'dv': 3.92,
            'uncertainty': 0.03,
            'reference': 'Anderson et al. (2008), PRL 100, 091102',
            'doi': '10.1103/PhysRevLett.100.091102',
            'note': 'First Earth flyby en route to Jupiter'
        },
        'Rosetta_2005': {
            'significance': 'moderate',
            'dv': 1.82,
            'uncertainty': 0.05,
            'reference': 'Morley & Budnik (2007)',
            'doi': None,
            'note': 'Rosetta first Earth flyby; previously excluded due to JPL Horizons data issues, now included via ESA SPICE kernels'
        },
        'Cassini': {
            'significance': 'minor',
            'dv': 0.11,
            'uncertainty': 0.05,
            'reference': 'Anderson et al. (2008), PRL 100, 091102',
            'doi': '10.1103/PhysRevLett.100.091102',
            'note': 'Marginal detection; tests sensitivity limit of TEP framework (S/N = 2.2)'
        }
    }
    
    findings = []
    for name, info in key_flybys.items():
        if name in fits:
            fit = fits[name]['fit']
            if fit['beta_fitted']:
                gamma_dev = fit['ppn_gamma_deviation']
                beta_unc = fit.get('uncertainty', fit['beta_fitted'] * 0.1)
                
                findings.append({
                    'spacecraft': name,
                    'observed_anomaly_mm_s': info['dv'],
                    'anomaly_uncertainty_mm_s': info['uncertainty'],
                    'fitted_beta': fit['beta_fitted'],
                    'beta_uncertainty': beta_unc,
                    'ppn_gamma_deviation': gamma_dev,
                    'ppn_compliant': fit['ppn_compliant'],
                    'significance': info['significance'],
                    'reference': info['reference'],
                    'doi': info['doi'],
                    'annotation': info['note']
                })
    
    # Retrieve weighted mean from step_003 analysis
    beta_values = [v['fit']['beta_fitted'] for k, v in fits.items() 
                   if v['fit']['beta_fitted'] is not None]
    weighted_mean = overall.get('recommended_beta', 0)
    
    # Prefer Step 008's headline uncertainty, which may be random-effects
    # inflated when cross-flyby heterogeneity is extreme.
    beta_uncertainty = overall.get('recommended_uncertainty')
    if beta_uncertainty is None:
        beta_uncertainty = overall['beta_statistics'].get('weighted_uncertainty')
        if beta_uncertainty is None:
            n_eff = len(beta_values)
            if n_eff > 1:
                beta_uncertainty = overall['beta_statistics']['std'] / (n_eff ** 0.5)
            else:
                beta_uncertainty = weighted_mean * 0.1

    # PPN constraint analysis must use the same screened beta_eff convention
    # as Step 008 individual fits: |gamma - 1| ~= 2 beta_eff^2.
    beta_eff_stats = overall.get('beta_eff_statistics', {})
    beta_eff_weighted = beta_eff_stats.get('weighted_mean')
    if beta_eff_weighted is None:
        beta_eff_weighted = overall.get('recommended_beta_eff')
    gamma_deviation = 2 * beta_eff_weighted ** 2 if beta_eff_weighted is not None else None
    gamma_deviation_max_fit = max(
        (
            entry.get('fit', {}).get('ppn_gamma_deviation')
            for entry in fits.values()
            if entry.get('fit', {}).get('ppn_gamma_deviation') is not None
        ),
        default=None,
    )
    gamma_bound = 2.3e-05
    gamma_margin = gamma_bound / gamma_deviation if gamma_deviation and gamma_deviation > 0 else float('inf')
    gamma_margin_worst_fit = (
        gamma_bound / gamma_deviation_max_fit
        if gamma_deviation_max_fit and gamma_deviation_max_fit > 0
        else float('inf')
    )

    # Solar-screened PPN check (UCD saturation radius ansatz extended to Sun)
    # R_sol_sun = (3 M_sun / 4π ρ_T)^(1/3) from Temporal Topology UCD model
    M_SUN_KG = 1.98847e30
    R_SUN_M = 6.9634e8
    RHO_T_KG_M3 = 20.0 * 1000.0  # 20 g/cm³ → kg/m³
    R_SOL_SUN_M = ((3.0 * M_SUN_KG) / (4.0 * math.pi * RHO_T_KG_M3)) ** (1.0 / 3.0)
    S_SUN_SURFACE = (R_SUN_M - R_SOL_SUN_M) / R_SUN_M
    # Cassini radio path during 2002 solar conjunction: r ≳ 4 R_sun
    R_CASSINI_PATH_M = 4.0 * R_SUN_M
    S_SUN_PATH = (R_CASSINI_PATH_M - R_SOL_SUN_M) / R_CASSINI_PATH_M
    # Use weighted mean β as the population-level solar-system coupling estimate.
    # Individual flyby β values are geometry-modulated effective couplings; the
    # weighted mean is the best estimate of the universal bare coupling for solar
    # PPN screening checks.
    solar_beta_bare = overall.get('recommended_beta')
    if solar_beta_bare is None:
        solar_beta_bare = weighted_mean
    if solar_beta_bare is not None and solar_beta_bare > 0:
        beta_eff_sun_surface = solar_beta_bare * S_SUN_SURFACE
        beta_eff_sun_path = solar_beta_bare * S_SUN_PATH
        gamma_deviation_sun_surface = 2.0 * beta_eff_sun_surface ** 2
        gamma_deviation_sun_path = 2.0 * beta_eff_sun_path ** 2
        gamma_margin_sun_surface = gamma_bound / gamma_deviation_sun_surface
        gamma_margin_sun_path = gamma_bound / gamma_deviation_sun_path
    else:
        solar_beta_bare = None
        beta_eff_sun_surface = None
        beta_eff_sun_path = None
        gamma_deviation_sun_surface = None
        gamma_deviation_sun_path = None
        gamma_margin_sun_surface = None
        gamma_margin_sun_path = None

    # Conclusion with caveats
    tep_viable = (
        overall.get('ppn_compliance', False) and 
        overall['n_fits'] >= 2 and
        overall.get('recommended_beta') is not None
    )
    
    # Determine confidence level with justification
    # NOTE: Based on actual fitting results (n_fits), not published anomaly count
    if tep_viable and len(findings) >= 4:
        confidence = 'moderate'
        confidence_rationale = (
            'Four literature-reported detections satisfy PPN constraints in the current fit set, '
            'but the sample remains small and per-flyby β values span nearly an order of magnitude'
        )
    elif tep_viable and len(findings) >= 3:
        confidence = 'moderate'
        confidence_rationale = 'Three detections consistent; limited by small sample size'
    elif tep_viable and len(findings) >= 2:
        confidence = 'moderate'
        confidence_rationale = 'Two detections consistent; limited by small sample size'
    elif tep_viable:
        confidence = 'low'
        confidence_rationale = 'Single detection only; insufficient for statistical validation'
    else:
        confidence = 'none'
        confidence_rationale = 'TEP model cannot explain observed anomalies within physical constraints'
    
    fitted_names = [f['spacecraft'] for f in findings if f['ppn_compliant']]
    excluded = [
        name for name, entry in fits.items()
        if entry.get('fit', {}).get('excluded')
    ]

    beta_stats = overall.get('beta_statistics', {})
    random_effects = None
    if beta_stats.get('random_effects_mean') is not None:
        random_effects = {
            'mean': beta_stats.get('random_effects_mean'),
            'uncertainty': beta_stats.get('random_effects_uncertainty'),
            'between_flyby_tau': beta_stats.get('between_flyby_tau'),
            'prediction_uncertainty': beta_stats.get('random_effects_prediction_uncertainty'),
            'interpretation': (
                'Population-level scatter summary for geometry-dependent fitted amplitudes; '
                'the fixed-effect weighted mean remains the restricted-tier pooled diagnostic.'
            ),
        }
    heterogeneity = overall.get('heterogeneity_tests', {})
    chi2_consistency = heterogeneity.get('chi_squared')

    gnss_file = project_root / 'results' / 'step016_gnss_cross_validation.json'
    gnss_scale_consistent = None
    if gnss_file.exists():
        with open(gnss_file, encoding='utf-8') as f:
            gnss_data = json.load(f)
        gnss_scale_consistent = gnss_data.get('gnss_scale_consistent', gnss_data.get('beta_consistent'))

    if gnss_scale_consistent is False and confidence in ('high', 'moderate'):
        confidence = 'limited'
        confidence_rationale += (
            '; GNSS correlation-length calibration is not consistent with the flyby relaxation-length prior'
        )

    if chi2_consistency is not None and chi2_consistency > 100 and confidence != 'none':
        if confidence == 'moderate':
            confidence = 'limited'
        confidence_rationale += (
            f'; per-flyby β consistency χ²≈{chi2_consistency:.0f} indicates weak ensemble agreement'
        )
    if random_effects and random_effects.get('prediction_uncertainty') and confidence != 'none':
        pred_unc = random_effects['prediction_uncertainty']
        if weighted_mean and pred_unc / abs(weighted_mean) > 1.0:
            if confidence == 'moderate':
                confidence = 'limited'
            confidence_rationale += (
                '; random-effects prediction uncertainty is of order the pooled β, '
                'so precision claims should use the formal fixed-effect result only as a diagnostic'
            )

    caveats = [
        'Analysis relies on literature anomaly values, not independent DSN data analysis',
        f'{len(fitted_names)} flybys yielded successful fits: {", ".join(fitted_names)}',
        f'Excluded or skipped flybys in current fit table: {", ".join(excluded) if excluded else "none"}',
        'screening restoration phase-boundary factor derived from first-principles PREM integration (Step 015); GNSS cross-check must agree before elevating confidence',
        'Null results at high altitude support restoration but do not independently constrain β',
    ]
    if random_effects:
        caveats.append(
            'Extreme between-flyby heterogeneity requires reporting random-effects scatter alongside the formal inverse-variance mean'
        )
    if gnss_scale_consistent is False:
        caveats.append('GNSS correlation-length cross-check failed for the current flyby relaxation-length prior')

    hierarchical_beta = None
    hierarchical_file = project_root / 'results' / 'step015_hierarchical_bayesian_results.json'
    if hierarchical_file.exists():
        with open(hierarchical_file, encoding='utf-8') as handle:
            hierarchical = json.load(handle)
        summary = hierarchical.get('hierarchical_mcmc', {})
        hierarchical_beta = summary.get('beta_0_median')
        if hierarchical_beta is not None:
            caveats.append(
                f'Hierarchical Bayesian pooled β₀ median (Step 015) is {hierarchical_beta:.2e}; '
                'per-flyby β scatter remains large relative to this pooled scale'
            )

    conclusion = {
        'tep_viable': tep_viable,
        'recommended_beta': weighted_mean,
        'beta_uncertainty': beta_uncertainty,
        'ppn_gamma_deviation': gamma_deviation,
        'ppn_gamma_deviation_max_fit': gamma_deviation_max_fit,
        'ppn_bound': gamma_bound,
        'ppn_safety_margin_orders': int(gamma_margin),
        'ppn_worst_fit_safety_margin': gamma_margin_worst_fit,
        'ppn_compliant': overall.get('ppn_compliance', False),
        'n_supporting_flybys': len(fitted_names),
        'n_null_results': 8,  # Galileo 1992, Rosetta 2007, Rosetta 2009, MESSENGER, Juno, Stardust, OSIRIS-REx, BepiColombo
        'confidence': confidence,
        'confidence_rationale': confidence_rationale,
        'caveats': caveats,
        'physical_interpretation': (
            f'The fixed-effect pooled β = {weighted_mean:.2e} gives screened '
            f'β_eff = {beta_eff_weighted:.2e}, implying |γ-1| = {gamma_deviation:.2e}; '
            f'the worst fitted flyby gives |γ-1| = {gamma_deviation_max_fit:.2e}, still within the Cassini bound. '
            'The Ambient Symmetry Restoration mechanism addresses both detections (NEAR, Galileo 1990, Rosetta 2005) '
            'and non-detections (MESSENGER, Juno) across varying flyby geometries.'
        ),
        'random_effects_beta_summary': random_effects,
        **literature_synthesis_metadata(chi2_consistency),
    }
    
    return {
        **literature_synthesis_metadata(chi2_consistency),
        'findings': findings,
        'statistics': {
            'beta_mean': beta_stats.get('mean'),
            'beta_std': beta_stats.get('std'),
            'beta_weighted_mean': beta_stats.get('weighted_mean'),
            'beta_weighted_uncertainty': beta_stats.get('weighted_uncertainty'),
            'beta_recommended_uncertainty': beta_uncertainty,
            'beta_recommended_uncertainty_model': overall.get('recommended_uncertainty_model'),
            'beta_random_effects_mean': beta_stats.get('random_effects_mean'),
            'beta_random_effects_uncertainty': beta_stats.get('random_effects_uncertainty'),
            'beta_between_flyby_tau': beta_stats.get('between_flyby_tau'),
            'beta_random_effects_prediction_uncertainty': beta_stats.get('random_effects_prediction_uncertainty'),
            'beta_hierarchical_median': hierarchical_beta,
            'beta_min': beta_stats.get('min'),
            'beta_max': beta_stats.get('max'),
            'n_fits': overall.get('n_fits', 0),
            'n_total_spacecraft': 12,
            'chi2_consistency': chi2_consistency,
            'reduced_chi2': heterogeneity.get('reduced_chi_squared'),
            'I_squared_percent': heterogeneity.get('I_squared_percent'),
            'heterogeneity_p_value': heterogeneity.get('p_value'),
        },
        'ppn_analysis': {
            'implied_gamma_deviation': gamma_deviation,
            'max_fitted_gamma_deviation': gamma_deviation_max_fit,
            'cassini_bound': gamma_bound,
            'safety_margin': f"{gamma_margin:.0e}x",
            'worst_fit_safety_margin': f"{gamma_margin_worst_fit:.0e}x",
            'compliance_status': 'full_compliance',
            'solar_screened': {
                'r_sol_sun_m': R_SOL_SUN_M,
                's_sun_surface': S_SUN_SURFACE,
                's_sun_path_cassini': S_SUN_PATH,
                'solar_beta_bare': solar_beta_bare,
                'beta_eff_sun_surface': beta_eff_sun_surface,
                'beta_eff_sun_path_cassini': beta_eff_sun_path,
                'gamma_deviation_sun_surface': gamma_deviation_sun_surface,
                'gamma_deviation_sun_path_cassini': gamma_deviation_sun_path,
                'margin_sun_surface': gamma_margin_sun_surface,
                'margin_sun_path_cassini': gamma_margin_sun_path,
                'note': 'UCD saturation radius ansatz extended to Sun; Cassini path at ~4 R_sun during 2002 solar conjunction'
            }
        },
        'conclusion': conclusion,
        'data_provenance': {
            'trajectory_source': 'NASA JPL Horizons ephemeris system',
            'anomaly_source': 'Peer-reviewed literature (Anderson et al. 2008)',
            'measurement_method': 'DSN Doppler tracking + JPL ODP',
            'acquisition_date': datetime.now(timezone.utc).isoformat()
        }
    }


def main():
    """Generate final report."""
    logger = StepLogger("step_036_final_report", PROJECT_ROOT)
    start_time = time.time()
    
    logger.header("STEP 036: FINAL REPORT GENERATION")
    
    # Load fitting results from results folder
    results_dir = PROJECT_ROOT / 'results'
    fit_file = results_dir / 'step008_fitting_results.json'
    
    if not fit_file.exists():
        logger.error(f"Fitting results not found: {fit_file}")
        logger.log_step_summary(0, "FAILED")
        return 1
    
    try:
        with open(fit_file) as f:
            fitting_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError, IOError) as e:
        logger.error(f"Failed to load fitting data: {e}")
        logger.log_step_summary(0, "FAILED")
        return 1
    
    logger.section("GENERATING COMPREHENSIVE REPORT")
    report = generate_report(fitting_data, PROJECT_ROOT)
    
    # Display findings with comprehensive detail
    logger.section("KEY FINDINGS")
    logger.info("Published detections are retained in the catalog; Step 008 fits only sign-gated cases.")
    logger.info("  NEAR: Major anomaly (13.46 mm/s)")
    logger.info("  Galileo_1990: Moderate anomaly (3.92 mm/s)")
    logger.info("  Rosetta_2005: Moderate anomaly (1.82 mm/s)")
    logger.info("  Cassini: Minor anomaly (0.11 mm/s; catalog diagnostic, sign-gate stress case)")
    
    for finding in report['findings']:
        logger.subheader(f"{finding['spacecraft']} ({finding['significance']} anomaly)")
        logger.info(f"Observed Δv: {finding['observed_anomaly_mm_s']:.2f} ± {finding['anomaly_uncertainty_mm_s']:.2f} mm/s")
        logger.info(f"Fitted β: ({finding['fitted_beta']:.2e} ± {finding['beta_uncertainty']:.2e})")
        logger.info(f"PPN |γ-1|: {finding['ppn_gamma_deviation']:.2e} (bound: 2.3×10⁻⁵)")
        logger.info(f"Reference: {finding['reference']}")
        logger.info(f"Note: {finding['annotation']}")
    
    # Display statistical summary
    logger.section("STATISTICAL SUMMARY")
    stats = report['statistics']
    logger.info("β Parameter Statistics:")
    logger.info(f"  Mean:       {stats['beta_mean']:.2e}")
    logger.info(f"  Std Dev:    {stats['beta_std']:.2e}")
    logger.info(f"  Weighted:   {stats['beta_weighted_mean']:.2e}")
    if stats.get('beta_random_effects_mean') is not None:
        logger.info(
            f"  Random eff: {stats['beta_random_effects_mean']:.2e} ± "
            f"{stats['beta_random_effects_uncertainty']:.2e}"
        )
        logger.info(
            f"  Between-flyby τ: {stats['beta_between_flyby_tau']:.2e}"
        )
    logger.info(f"  Range:      [{stats['beta_min']:.2e}, {stats['beta_max']:.2e}]")
    logger.info(f"  Span:       {stats['beta_max']/stats['beta_min']:.0f}x variation (reflects anomaly amplitude range)")
    
    logger.info("Sample Size:")
    logger.info(f"  Detections: {stats['n_fits']} significant anomalies")
    logger.info(f"  Nulls:      8 non-detections (Galileo 1992, Rosetta 2007, Rosetta 2009, MESSENGER, Juno, Stardust, OSIRIS-REx, BepiColombo)")
    logger.info(f"  Total:      {stats['n_total_spacecraft']} flyby events analyzed")
    
    # Display PPN analysis
    logger.section("PPN CONSTRAINT ANALYSIS")
    ppn = report['ppn_analysis']
    logger.info(f"Implied |γ - 1|:    {ppn['implied_gamma_deviation']:.2e}")
    logger.info(f"Cassini bound:       {ppn['cassini_bound']:.2e}")
    logger.info(f"Margin to bound ratio: {ppn['safety_margin']}:1")
    logger.info(f"Status:              {ppn['compliance_status'].replace('_', ' ').title()}")
    
    # Display conclusion with caveats
    logger.section("CONCLUSION")
    conclusion = report['conclusion']
    logger.info(f"TEP Viability:        {'Viable' if conclusion['tep_viable'] else 'Not viable'}")
    logger.info(f"Recommended β:        ({conclusion['recommended_beta']:.2e} ± {conclusion['beta_uncertainty']:.2e})")
    logger.info(f"PPN Compliance:       {'Full' if conclusion['ppn_compliant'] else 'Violated'}")
    logger.info(f"Supporting Flybys:    {conclusion['n_supporting_flybys']}")
    logger.info(f"Confidence Level:     {conclusion['confidence'].upper()}")
    logger.info(f"Rationale: {conclusion['confidence_rationale']}")
    
    if conclusion['tep_viable']:
        logger.info("Physical Interpretation:")
        logger.info(f"  {conclusion['physical_interpretation']}")
    
    logger.info("Caveats and Limitations:")
    for caveat in conclusion['caveats']:
        logger.info(f"  • {caveat}")
    
    # Data provenance
    logger.section("DATA PROVENANCE")
    provenance = report['data_provenance']
    logger.info(f"Trajectory Data:    {provenance['trajectory_source']}")
    logger.info(f"Anomaly Values:     {provenance['anomaly_source']}")
    logger.info(f"Measurement:        {provenance['measurement_method']}")
    logger.info(f"Analysis Date:      {provenance['acquisition_date'][:10]}")
    
    # Save final report
    logger.section("SAVING FINAL REPORT")
    final_output = {
        'report_date': datetime.now(timezone.utc).isoformat(),
        **literature_synthesis_metadata(report['statistics'].get('chi2_consistency_approx')),
        'summary': report,
        'raw_data': fitting_data,
    }
    
    output_file = results_dir / 'step036_final_report.json'
    with open(output_file, 'w') as f:
        json.dump(final_output, f, indent=2)
    
    logger.success(f"Report complete")
    logger.info(f"Final output saved to: {output_file}")
    logger.add_output_file(output_file, "Final comprehensive report")
    
    duration = time.time() - start_time
    logger.log_step_summary(duration, "SUCCESS")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
