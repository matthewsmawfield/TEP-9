# Flyby TEP Pipeline v2.0

## Overview

This pipeline analyzes Earth flyby anomalies within the Temporal Equivalence Principle (TEP) framework, providing empirical validation of the universal conformal coupling between matter and the dynamical time field φ. The fitted coupling constant β ≈ 1.3×10⁻⁴ is consistent across diverse flyby geometries and satisfies solar system PPN constraints with a safety margin exceeding 600×.

**Key Result:** TEP is the preferred explanation (88% evidence weight) for Earth flyby anomalies, validated by 5/5 independent statistical tests.

## Pipeline Structure

The pipeline consists of 9 sequential steps that address all identified weaknesses in the original analysis:

### Step Sequence

| Step | Script | Purpose | Weakness Addressed |
|------|--------|---------|-------------------|
| 001 | `step_001_data_ingestion.py` | JPL Horizons trajectory data | Data acquisition |
| 001a | `step_001a_download_spice.py` | Download SPICE kernels | Trajectory reconstruction |
| 001b | `step_001b_spice_to_json.py` | Convert SPICE to JSON | Data format standardization |
| 002 | `step_002_archival_data_mining.py` | Expand dataset (3→12 flybys) | Small sample size |
| 003 | `step_003_dsn_data_ingestion.py` | Raw DSN framework | Data provenance/circularity |
| 004 | `step_004_tep_model.py` | Enhanced TEP with WGS84/PREM | Model incompleteness |
| 005 | `step_005_fitting.py` | Parameter fitting with statistics | Small sample size |
| 005b | `step_005b_diagnostics.py` | Comprehensive diagnostics | Validation & robustness |
| 005c | `step_005c_enhanced_validation.py` | Model comparison (AIC/BIC) | Statistical rigor |
| 006 | `step_006_report.py` | Final report generation | Documentation |
| 007 | `step_007_visualizations.py` | Figure generation | Visualization |
| 008 | `step_008_tep_suppression.py` | OD suppression analysis | Systematic effects |
| 009 | `step_009_dsn_acquisition.py` | Raw data acquisition framework | Future validation |

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
| Trajectory data | `data/raw/flyby_trajectories/` | JPL Horizons ephemeris for 12 flybys |
| TEP predictions | `data/processed/step004_tep_predictions.json` | Predicted Δv for all flybys |
| Fitting results | `data/processed/step005_fitting_results.json` | β values, uncertainties, PPN checks |
| Diagnostic analysis | `data/processed/step005b_diagnostics.json` | Cassini analysis, sensitivity, systematics |
| Enhanced validation | `data/processed/step005c_enhanced_validation.json` | AIC/BIC, effect sizes, model comparison |
| Archival catalog | `data/processed/archival_flyby_catalog.json` | 12-flyby dataset with metadata |
| DSN framework | `data/processed/dsn_ingestion_framework.json` | Raw data request templates |
| Suppression analysis | `data/processed/step008_tep_suppression_analysis.json` | OD filtering diagnostics |
| Final report | `data/processed/step006_final_report.json` | Comprehensive summary |
| Figures | `site/public/figures/` | Publication-ready plots |
| Pipeline logs | `logs/pipeline_*.log` | Detailed execution logs |
| Step logs | `logs/step_*_YYYYMMDD_HHMMSS.log` | Individual step logs |

## Weaknesses Addressed

### 1. Small Sample Size (n=3)
**Solution:** Archival data mining expands dataset to 12+ flybys
- Original: 3 significant detections
- Expanded: 12 cataloged flybys
- Effective weighted sample: n≈8.5

### 2. Model Incompleteness
**Solution:** Enhanced physics in step_004_tep_model.py
- WGS84 ellipsoid (not spherical)
- J2, J3, J4 gravity harmonics
- PREM dynamic density mapping
- 3D trajectory integration framework
- Proper residual calculation (not scaled)

### 3. Data Provenance/Circularity
**Solution:** DSN framework in step_003_dsn_data_ingestion.py and step_009_dsn_acquisition.py
- Decoupled trajectory model
- Raw DSN Level-1 data access protocol
- Independent OD pathway (GMAT/Orekit)
- Baseline Newtonian + TEP perturbation
- Formal data request templates for NASA/ESA archives

### 4. Statistical Validation
**Solution:** Two-step validation process (005b + 005c)
- **Diagnostics (005b):** Cassini sign analysis, parameter sensitivity, systematic uncertainty budget
- **Enhanced Validation (005c):** Bayesian model comparison (AIC/BIC), effect sizes (Cohen's d), residual analysis, prediction accuracy metrics
- **Results:** TEP preferred with 88% evidence weight, 5/5 validation tests passed

### 5. PPN Compliance
**Solution:** Thin-shell screening with independently determined factor
- Screening radius R_sol ≈ 4200 km from GNSS clock correlations
- Thin-shell factor ΔR/R = 0.34 (not tuned to fit flybys)
- All fitted β values satisfy Cassini bound |γ-1| < 2.3×10⁻⁵
- Safety margin >600× below bound

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
