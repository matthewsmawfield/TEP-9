# Flyby TEP Pipeline v5.0 (Workflow-Organized)

## Overview

This pipeline analyzes Earth flyby anomalies within the Temporal Equivalence Principle (TEP) framework, providing empirical validation of the universal conformal coupling between matter and the dynamical time field φ. The fitted coupling constant β ≈ 1.67×10⁻³ is consistent across diverse flyby geometries and satisfies solar system PPN constraints.

Key Result: TEP is the preferred explanation for Earth flyby anomalies. The 7.8-fold variance in fitted β values is comprehensively explained through a four-stage decomposition (Step 006).

## Reorganization Summary

Version 5.0 reorganizes the pipeline from chronological addition order to workflow-based analytical phases:

- Fixed duplicate step numbers (007, 010, 017, 031)
- Archived development artifacts (029-031 series)
- Unified variance analysis (consolidates old 5b, 5c, 5d, 22, 28 into Step 006)
- Clear phase-based grouping for easier navigation

## Pipeline Structure

### Workflow Phases

| Phase | Steps | Description |
|-------|-------|-------------|
| 1 | 001-004 | Data Acquisition & Preparation |
| 2 | 004-007 | Core Physics & Variance Analysis |
| 3 | 008-010 | Trajectory & Observational Pipeline |
| 4 | 011-014 | Validation & Robustness |
| 5 | 015-017 | Extended Physics |
| 6 | 018-019 | Model Comparison |
| 7 | 020-023 | DSN Data Framework |
| 8 | 024-026 | Mission-Specific Analysis |
| 9 | 027-028 | Advanced Topics |
| 10 | 029 | Reporting |

### Step Sequence

| Step | Script | Purpose | Phase |
|------|--------|---------|-------|
| 001a | `step_001a_download_spice.py` | Download SPICE kernels (NAIF) | 1 |
| 001b | `step_001b_spice_to_json.py` | Convert SPICE to JSON | 1 |
| 002 | `step_002_archival_data_mining.py` | Expand dataset | 1 |
| 002b | `step_002b_jpl_horizons_fetch.py` | JPL Horizons Data Fetch | 1 |
| 003 | `step_003_dsn_data_ingestion.py` | DSN Data Ingestion | 1 |
| 003b | `step_003_dsn_framework.py` | DSN Framework | 1 |
| 004 | `step_004_tep_model.py` | Enhanced TEP Model | 2 |
| 005 | `step_005_fitting.py` | Parameter Fitting | 2 |
| 006 | `step_007_variance_analysis.py` | Unified Variance Analysis (NEW) | 2 |
| 007 | `step_008_Temporal Shear Suppression_first_principles.py` | Temporal Shear Suppression First-Principles | 2 |
| 008 | `step_009_trajectory_integration.py` | Trajectory Integration | 3 |
| 009 | `step_010_od_filter_simulation.py` | OD Filter Simulation | 3 |
| 010 | `step_011_cross_validation.py` | Cross-Validation | 4 |
| 011 | `step_012_sensitivity_analysis.py` | Sensitivity Analysis | 4 |
| 012 | `step_013_hierarchical_bayesian.py` | Hierarchical Bayesian | 4 |
| 013 | `step_014_gnss_validation.py` | GNSS Validation | 4 |
| 014 | `step_015_plasma_modulation.py` | Plasma Modulation | 5 |
| 015 | `step_016_space_weather.py` | Space Weather Correlation | 5 |
| 016 | `step_016b_3d_field_integration.py` | 3D Field Integration | 5 |
| 017 | `step_017_enhanced_bayesian.py` | Enhanced Bayesian | 5 |
| 018 | `step_018_bayesian_model_comparison.py` | Bayesian Model Comparison | 6 |
| 019 | `step_019_saturation_model.py` | Saturation Model | 6 |
| 020 | `step_021_dsn_processing.py` | DSN Processing Framework | 7 |
| 021 | `step_022_read_trk234.py` | Read TRK-2-34 Format | 7 |
| 022 | `step_023_raw_dsn_reanalysis.py` | Raw DSN Reanalysis | 7 |
| 023 | `step_024_juno_reanalysis.py` | Juno 2013 Reanalysis | 8 |
| 024 | `step_025_pds_search.py` | PDS Search | 8 |
| 025 | `step_026_tep_suppression.py` | TEP Suppression | 8 |
| 026 | `step_027_covariant_holonomy.py` | Covariant Holonomy | 9 |
| 027 | `step_028_cross_corpus_export.py` | Cross-Corpus Export | 9 |
| 028 | `step_029_final_report.py` | Final Report | 10 |
| 029 | `step_020_visualizations.py` | Figure Generation | 10 |

### Supporting Modules

- `enhanced_physics.py` - Non-spherical Earth model (WGS84), PREM density, J2/J3/J4 harmonics, 3D integration
- `run_all.py` - Pipeline orchestration with comprehensive logging and error handling
- `step_logger.py` - Standardized logging across all steps with timestamps and progress tracking

## Usage

### Run Complete Pipeline

```bash
cd /Users/matthewsmawfield/www/TEP-EFA
python scripts/run_all.py
```

This executes all steps in sequence with comprehensive logging. The pipeline will:
1. Download and process trajectory data
2. Compute TEP predictions for all 12 flybys
3. Fit β parameter and validate against PPN constraints
4. Run diagnostic analyses (bootstrap, LOO-CV, sensitivity)
5. Perform model comparison (TEP vs null vs empirical)
6. Generate figures and final report

### Run Individual Steps

```bash
python scripts/steps/step_001_data_ingestion.py
python scripts/steps/step_004_tep_model.py
python scripts/steps/step_005_fitting.py
python scripts/steps/step_005b_diagnostics.py
python scripts/steps/step_005c_enhanced_validation.py
# etc.
```

Each step can be run independently, but steps 005+ require outputs from earlier steps.

## Logging System

The pipeline produces comprehensive logs:

### Console Output
- Real-time progress indicators
- Step-by-step status ([1/7], [2/7], etc.)
- Duration tracking for each step
- Clear success/failure indicators

### Log Files
Location: `logs/pipeline_YYYYMMDD_HHMMSS.log`

Each log entry includes:
- Timestamp (UTC)
- Log level (INFO, SUCCESS, ERROR, WARNING)
- Message content

### Execution Summary
At completion, the pipeline outputs:
- Execution timeline (start/end/duration)
- Step-by-step results table
- Success/failure statistics
- Output file locations

## Output Files

| Output | Location | Description |
|--------|----------|-------------|
| Trajectory data | `data/raw/jpl_horizons/` | JPL Horizons ephemeris and raw responses |
| TEP predictions | `results/step004_tep_predictions.json` | Predicted Δv for modeled flybys |
| Fitting results | `results/step005_fitting_results.json` | β values, uncertainties, PPN checks |
| Variance analysis | `results/step007_variance_analysis.json` | Unified variance decomposition |
| Enhanced validation | `results/step018_enhanced_bayesian.json` | Bayesian robustness analysis |
| Archival catalog | `results/step002_archival_flyby_catalog.json` | Literature flyby dataset with metadata |
| DSN framework | `data/processed/dsn_ingestion_framework.json` | Raw data request templates |
| Suppression analysis | `results/step026_tep_suppression_analysis.json` | OD filtering diagnostics |
| Final report | `results/step006_final_report.json` | Comprehensive summary |
| Figures | `results/step007_figure*.png` | Publication-ready plots |
| Pipeline logs | `logs/pipeline.log` | Detailed execution logs |
| Step logs | `logs/step_*.log` | Individual step logs |

## Weaknesses Addressed

### 1. Small Sample Size (n=3)
Solution: Archival data mining expands dataset to 12+ flybys
- Original: 3 significant detections
- Expanded: 12 cataloged flybys
- Effective weighted sample: n≈8.5

### 2. Model Incompleteness
Solution: Enhanced physics in step_004_tep_model.py
- WGS84 ellipsoid (not spherical)
- J2, J3, J4 gravity harmonics
- PREM dynamic density mapping
- 3D trajectory integration framework
- Proper residual calculation (not scaled)

### 3. Data Provenance/Circularity
Solution: DSN framework in step_003_dsn_data_ingestion.py and step_009_dsn_acquisition.py
- Decoupled trajectory model
- Raw DSN Level-1 data access protocol
- Independent OD pathway (GMAT/Orekit)
- Baseline Newtonian + TEP perturbation
- Formal data request templates for NASA/ESA archives

### 4. Statistical Validation
Solution: Two-step validation process (005b + 005c)
- Diagnostics (005b): Cassini sign analysis, parameter sensitivity, systematic uncertainty budget
- Enhanced Validation (005c): Bayesian model comparison (AIC/BIC), effect sizes (Cohen's d), residual analysis, prediction accuracy metrics
- Results: TEP preferred with 88% evidence weight, 5/5 validation tests passed

### 5. PPN Compliance
Solution: Temporal Shear suppression with independently determined factor
- Screening radius R_sol ≈ 4146 km from UCD soliton physics (Paper 7)
- Characteristic suppression S_⊕ ≈ 0.35 (not tuned to fit flybys)
- All fitted β values satisfy Cassini bound |γ-1| < 2.3×10⁻⁵
- Margin to bound ratio: >10⁹:1

## Configuration

Pipeline parameters are configured in `config/pipeline_config.json`:
- TEP physics parameters (β, Λ, n)
- Earth model settings (geoid, harmonics)
- Integration tolerances
- Statistical method options

## Error Handling

- Pipeline stops on first failure
- All prior step outputs preserved
- Detailed error logging
- Exit codes: 0 = success, 1 = failure
