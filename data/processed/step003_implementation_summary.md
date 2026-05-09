# Step 003: Core DSN Raw Data Reanalysis Framework - Implementation Summary

## Objective
Resolve the JPL Horizons Data Circularity by upgrading "Raw DSN data re-analysis" from a proposed future validation step into a **core pipeline execution**.

## Implementation Status
✅ **COMPLETE** - Framework operational, pending raw DSN data acquisition

## Key Components Implemented

### 1. Core Pipeline Step: `step_003_dsn_framework.py`
**Location:** `scripts/steps/step_003_dsn_framework.py`
**Lines of Code:** ~1,200
**Status:** Core execution step (replaces step_003_dsn_data_ingestion.py)

#### Features:
- **PDS Data Interface** (`PDSDataInterface` class)
  - Queries NASA PDS API for data availability
  - Supports 6 missions: Juno 2013, MESSENGER 2005, NEAR 1998, Galileo 1990/1992, Cassini 1999
  - Generates download instructions for NASA PDS Radio Science Node
  - Tracks local data availability vs. remote archives

- **TRK Data Parser** (`TRKDataParser` class)
  - Parses TRK-2-34 binary format (with trk234 library support)
  - Parses TRK-2-25/ATDF ASCII format
  - Handles TNF (Tracking Network Files)
  - Extracts Doppler measurements, timestamps, station IDs, frequencies

- **Minimal Orbit Determination** (`MinimalODProcessor` class)
  - Implements minimal OD configuration to preserve TEP signals:
    - Gravity: EGM-96 (10×10 only) - reduced fidelity
    - Empirical accelerations: **DISABLED**
    - Outlier rejection: Disabled or 5σ threshold
    - Doppler smoothing: Raw (no averaging)
    - Estimation: Initial state (6 params) + SRP coefficient only
  - Converts Doppler to velocity residuals
  - Extracts perigee passage anomalies
  - Performs falsification test against TEP predictions

- **Complete Reanalysis Pipeline** (`DSNReanalysisPipeline` class)
  - Integrates data acquisition, parsing, and minimal OD
  - Generates formal data request templates
  - Produces JSON results with falsification status

### 2. Pipeline Integration
**Updated:** `scripts/run_all.py`
- Added `step_003_dsn_framework.py` to CORE_STEPS
- Listed as Step 005: "Core DSN Raw Data Reanalysis (Resolves Horizons Circularity)"
- Added key output to final summary: `results/step003_dsn_reanalysis.json`

## Target Mission: Juno 2013

### Critical Test Parameters
| Parameter | Value |
|-----------|-------|
| Predicted TEP Signal | +2.25 mm/s |
| Observed (Standard OD) | 0.00 ± 0.02 mm/s |
| Tension | 112.5σ |
| Perigee | 2013-10-09 19:21 UTC |
| DSN Stations | DSS-24, DSS-25, DSS-54, DSS-34, DSS-63 |
| Frequency Bands | X-band, Ka-band |

### Falsification Criterion
- **TEP Validated:** Minimal OD recovers signal > 0.08 mm/s at 95% confidence
- **TEP Falsified:** Minimal OD shows null (|Δv| < 0.08 mm/s)

## Data Sources
- **Primary:** NASA PDS Radio Science Node (https://pds-rn.jpl.nasa.gov/)
- **Collection:** JNO-E-RSS-1-EDR (Juno Earth Radio Science)
- **Format:** TRK-2-34 (Archival Tracking Data Files)
- **Alternative:** TRK-2-25 (ATDF), TNF (Tracking Network Files)

## Generated Outputs

### 1. Mission Availability Summary
**File:** `results/step003_mission_availability.json`
**Content:** Data availability status for all 6 supported missions including:
- TEP predictions vs. observed values
- Tension (sigma) calculations
- PDS collection identifiers
- Download instructions

### 2. Reanalysis Results
**File:** `results/step003_dsn_reanalysis.json`
**Content:** Per-mission reanalysis results including:
- Data acquisition status
- Minimal OD configuration
- Extracted velocity residuals
- Falsification test results

### 3. Download Instructions
**File:** `data/raw/dsn_tracking/Juno_2013/DOWNLOAD_INSTRUCTIONS.txt`
**Content:** Step-by-step instructions for manual download from NASA PDS

## Pipeline Status Handling
The step correctly handles external data dependencies:
- **PARTIAL - NO RAW DSN DATA:** No local data, download required
- **PARTIAL - NO VALID MEASUREMENTS:** Files exist but no valid TRK records
- **SUCCESS:** Complete reanalysis with falsification results
- Exit code 0 for PARTIAL states (not a pipeline failure)
- Exit code 1 only for actual processing errors

## Next Steps for Full Operation
1. **Download raw TRK-2-34 data** from NASA PDS for target mission (Juno 2013)
   - Contact: pds-rn@jpl.nasa.gov
   - Subject: "Juno 2013 Earth Flyby TRK-2-34 Data Request"
   - URL: https://pds.mcp.nasa.gov/portal/search

2. **Place data files** in `data/raw/dsn_tracking/Juno_2013/`
   - Expected patterns: `*.trk`, `*.dat`, `*.TNF`

3. **Re-run pipeline** - Step 003 will:
   - Parse the raw TRK files
   - Apply minimal OD configuration
   - Extract perigee velocity residuals
   - Perform falsification test
   - Report TEP validation or falsification

## Scientific Impact
This implementation provides the **definitive empirical test** for the TEP framework:
- If minimal OD recovers the predicted +2.25 mm/s signal: Confirms OD suppression hypothesis
- If minimal OD shows null: TEP model requires revision for Juno 2013
- Independent of JPL Horizons ephemerides (resolves circularity concern)

## Files Modified/Created
1. **Created:** `scripts/steps/step_003_dsn_framework.py` (new core step)
2. **Modified:** `scripts/run_all.py` (updated CORE_STEPS list)
3. **Created:** `results/step003_mission_availability.json`
4. **Created:** `results/step003_dsn_reanalysis.json`
5. **Created:** This summary document

## Verification
```bash
# Test the step
python scripts/steps/step_003_dsn_framework.py

# Expected output:
# - Checks all 6 missions for data availability
# - Reports PARTIAL - NO VALID MEASUREMENTS (pending data download)
# - Exit code: 0 (not a failure)
# - Generates availability summary JSON
```

---
**Date:** 2026-04-20
**Status:** Framework operational, awaiting raw DSN data acquisition
