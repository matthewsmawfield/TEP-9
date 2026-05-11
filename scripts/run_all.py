#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Run All Steps (Workflow-Organized v5.0)
============================================================

Executes the reorganized 29-step pipeline organized by analytical workflow phases.

Phase 1: Data Acquisition & Preparation (001-004)
  1a. Download SPICE Kernels (NAIF)
  1b. Convert SPICE to JSON
  2.  Archival Data Mining
  2b. JPL Horizons Data Fetch
  3.  DSN Data Ingestion
  3b. DSN Framework (Resolves Horizons Circularity)

Phase 2: Core Physics (004-007)
  4.  Enhanced TEP Model (WGS84, PREM, 3D integration)
  5.  Parameter Fitting
  6.  Unified Variance Analysis (Four-Stage Decomposition)
  7.  Temporal Shear Suppression First-Principles

Phase 3: Trajectory & Observational Pipeline (008-010)
  8.  Trajectory Integration
  9.  OD Filter Simulation (Observational pipeline validation)
  10. Cross-Validation Analysis

Phase 4: Validation & Robustness (011-014)
  11. Sensitivity Analysis
  12. Hierarchical Bayesian Model
  13. GNSS Validation

Phase 5: Extended Physics (015-017)
  14. Plasma Modulation (Cassini sign mismatch)
  15. Space Weather Correlation (Solar activity)
  16. 3D Field Integration
  17. Enhanced Bayesian Analysis

Phase 6: Model Comparison (018-019)
  18. Bayesian Model Comparison
  19. Saturation Model Analysis

Phase 7: DSN Data Framework (020-023)
  20. DSN Processing Framework
  21. Read TRK-2-34 Data Format
  22. Raw DSN Reanalysis

Phase 8: Mission-Specific Analysis (024-026)
  23. Juno 2013 Reanalysis
  24. PDS Search
  25. TEP Suppression Analysis

Phase 9: Advanced Topics (026-028)
  26. Covariant Temporal Shear Impulse
  27. Cross-Corpus Parameter Export

Phase 10: Reporting (028-029)
  28. Final Report Generation
  29. Figure Generation

Key Reorganization:
- Workflow-based grouping instead of chronological addition order
- Fixed duplicate step numbers (007, 010, 017, 031)
- Archived development artifacts (029-031 series)
- Unified variance analysis in Step 006 (consolidates 5b, 5c, 5d, 22, 28)

All steps produce detailed logs with timestamps and execution metrics.
"""

import sys
import subprocess
import time
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Tuple, Dict

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STEPS_DIR = PROJECT_ROOT / 'scripts' / 'steps'
LOGS_DIR = PROJECT_ROOT / 'logs'

# Ensure logs directory exists
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# All pipeline steps (run_all runs everything by default)
# Organized by workflow phases for clarity
CORE_STEPS: List[Tuple[str, str]] = [
    # Phase 1: Data Acquisition & Preparation (001-006)
    ('step_001_download_spice.py', 'Step 001: Download SPICE Kernels (NAIF)'),
    ('step_002_spice_to_json.py', 'Step 002: Convert SPICE to JSON'),
    ('step_003_archival_data_mining.py', 'Step 003: Archival Data Mining'),
    ('step_004_jpl_horizons_fetch.py', 'Step 004: JPL Horizons Data Fetch'),
    ('step_005_dsn_data_ingestion.py', 'Step 005: DSN Data Ingestion'),
    ('step_006_dsn_framework.py', 'Step 006: DSN Framework (Resolves Horizons Circularity)'),
    
    # Phase 2: Core Physics (007-010)
    ('step_007_tep_model.py', 'Step 007: Enhanced TEP Model (WGS84, PREM, 3D)'),
    ('step_008_fitting.py', 'Step 008: Parameter Fitting'),
    ('step_009_variance_analysis.py', 'Step 009: Unified Variance Analysis (Four-Stage)'),
    ('step_010_tep_first_principles.py', 'Step 010: Temporal Shear Suppression First-Principles'),
    
    # Phase 3: Trajectory & Observational Pipeline (011-013)
    ('step_011_trajectory_integration.py', 'Step 011: Trajectory Integration'),
    ('step_012_od_filter_simulation.py', 'Step 012: OD Filter Simulation'),
    ('step_013_cross_validation.py', 'Step 013: Cross-Validation Analysis'),
    
    # Phase 4: Validation & Robustness (014-016)
    ('step_014_sensitivity_analysis.py', 'Step 014: Sensitivity Analysis'),
    ('step_015_hierarchical_bayesian.py', 'Step 015: Hierarchical Bayesian Model'),
    ('step_016_gnss_validation.py', 'Step 016: GNSS Validation'),
    
    # Phase 5: Extended Physics (017-019)
    ('step_017_plasma_modulation.py', 'Step 017: Plasma Modulation (Cassini Sign)'),
    ('step_018_space_weather.py', 'Step 018: Space Weather Correlation'),
    ('step_019_3d_field_integration.py', 'Step 019: 3D Field Integration'),
    
    # Phase 6: Plasma Environment (020-021)
    ('step_020_plasma_environment_reconstruction.py', 'Step 020: Plasma Environment Reconstruction'),
    ('step_021_mission_specific_od_absorption.py', 'Step 021: Mission-Specific OD Absorption'),
    
    # Phase 7: Atmospheric & Thermal Effects (022-023)
    ('step_022_atmospheric_drag_simulation.py', 'Step 022: Atmospheric Drag Simulation'),
    ('step_023_thermal_recoil_modeling.py', 'Step 023: Thermal Recoil Simulation'),
    
    # Phase 8: Statistical Validation (024-026)
    ('step_024_systematic_error_monte_carlo.py', 'Step 024: Systematic Error Monte Carlo'),
    ('step_025_corrected_uncertainty.py', 'Step 025: Corrected Uncertainty Analysis'),
    ('step_026_stable_model_comparison.py', 'Step 026: Stable Model Comparison'),
    
    # Phase 9: Additional Validation (027)
    ('step_027_claim_consistency_audit.py', 'Step 027: Claim Consistency Audit'),
    
    # Phase 10: DSN Data Framework (028-029)
    ('step_028_dsn_processing.py', 'Step 028: DSN Processing Framework'),
    ('step_029_read_trk234.py', 'Step 029: Read TRK-2-34 Data Format'),
    
    # Phase 11: Mission-Specific Analysis (030-032)
    ('step_030_juno_reanalysis.py', 'Step 030: Juno 2013 Reanalysis'),
    ('step_031_pds_search.py', 'Step 031: PDS Search'),
    ('step_032_tep_suppression.py', 'Step 032: TEP Suppression Analysis'),
    
    # Phase 12: Advanced Topics (033-035)
    ('step_033_iri_trajectory_profile.py', 'Step 033: Continuous IRI Trajectory Profiles'),
    ('step_034_covariant_holonomy.py', 'Step 034: Covariant Temporal Shear Impulse'),
    ('step_035_cross_corpus_export.py', 'Step 035: Cross-Corpus Parameter Export'),
    
    # Phase 13: Reporting (036-037)
    ('step_036_final_report.py', 'Step 036: Final Report Generation'),
    ('step_037_visualizations.py', 'Step 037: Figure Generation')
]

DATA_INTEGRITY_REQUIRED_OUTPUTS = [
    'step003_archival_flyby_catalog.json',
    'step007_tep_predictions.json',
    'step008_fitting_results.json',
    'step036_final_report.json',
]


class PipelineLogger:
    """Handles consistent logging for pipeline execution."""
    
    def __init__(self):
        self.start_time = datetime.now(timezone.utc)
        self.log_file = LOGS_DIR / "pipeline.log"
        self.step_results: List[Dict] = []
        
    def _write(self, message: str, level: str = "INFO"):
        """Write message to both console and log file."""
        timestamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        log_line = f"[{timestamp}] [{level:8}] {message}"
        
        # Print to console
        print(message)
        
        # Write to log file (overwrite on first write, then append)
        mode = 'w' if not hasattr(self, '_log_initialized') else 'a'
        with open(self.log_file, mode, encoding='utf-8') as f:
            f.write(log_line + '\n')
        self._log_initialized = True
    
    def header(self, title: str, width: int = 70):
        """Print formatted header."""
        self._write("=" * width)
        self._write(title)
        self._write("=" * width)
    
    def subheader(self, title: str, width: int = 70):
        """Print formatted subheader."""
        self._write("-" * width)
        self._write(title)
        self._write("-" * width)
    
    def info(self, message: str):
        """Log info message."""
        self._write(message, "INFO")
    
    def success(self, message: str):
        """Log success message."""
        self._write(f"✓ {message}", "SUCCESS")
    
    def error(self, message: str):
        """Log error message."""
        self._write(f"✗ {message}", "ERROR")
    
    def warning(self, message: str):
        """Log warning message."""
        self._write(f"⚠ {message}", "WARNING")
    
    def step_start(self, step_num: int, total_steps: int, description: str):
        """Log step start with progress indicator."""
        progress = f"[{step_num}/{total_steps}]"
        self._write("")
        self.header(f"{progress} {description}")
        return time.time()
    
    def step_complete(self, description: str, start_time: float, status: str = "SUCCESS"):
        """Log step completion with duration."""
        duration = time.time() - start_time
        duration_str = f"{duration:.2f}s" if duration < 60 else f"{duration/60:.2f}m"
        
        if status == "SUCCESS":
            self.success(f"{description} completed in {duration_str}")
        else:
            self.error(f"{description} failed after {duration_str}")
        
        return duration
    
    def add_step_result(self, step_name: str, description: str, 
                       status: str, duration: float, exit_code: int = 0):
        """Record step result for summary."""
        self.step_results.append({
            'step': step_name,
            'description': description,
            'status': status,
            'duration_seconds': duration,
            'exit_code': exit_code,
            'timestamp': datetime.now(timezone.utc).isoformat()
        })
    
    def generate_research_audit(self):
        """Generate a professional-grade research audit for reproducibility."""
        import platform
        import hashlib
        import json
        import os

        def get_file_hash(path):
            if not path.exists(): return "MISSING"
            sha256_hash = hashlib.sha256()
            with open(path, "rb") as f:
                for byte_block in iter(lambda: f.read(4096), b""):
                    sha256_hash.update(byte_block)
            return sha256_hash.hexdigest()

        # Gather system telemetry
        telemetry = {
            "os": f"{platform.system()} {platform.release()}",
            "python_version": platform.python_version(),
            "processor": platform.processor(),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "project_root": str(PROJECT_ROOT)
        }

        # Gather script integrity hashes
        integrity = {}
        for filename, _ in CORE_STEPS:
            path = STEPS_DIR / filename
            integrity[filename] = get_file_hash(path)

        # Final audit object
        audit_data = {
            "audit_type": "TEP-EFA Research-Grade Audit",
            "theory_version": "Jakarta v0.8.0",
            "telemetry": telemetry,
            "pipeline_results": self.step_results,
            "script_integrity": integrity
        }

        # Save audit
        audit_dir = PROJECT_ROOT / "results" / "audits"
        audit_dir.mkdir(parents=True, exist_ok=True)
        audit_filename = f"RESEARCH_AUDIT_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        audit_path = audit_dir / audit_filename
        
        with open(audit_path, "w", encoding="utf-8") as f:
            json.dump(audit_data, f, indent=4)
        
        self.info(f"Research audit generated: {audit_path.name}")
        return audit_path

    def final_summary(self, total_steps: int) -> bool:
        """Print detailed final summary and generate research audit."""
        total_duration = (datetime.now(timezone.utc) - self.start_time).total_seconds()
        
        success_count = sum(1 for r in self.step_results if r['status'] == 'SUCCESS')
        fail_count = sum(1 for r in self.step_results if r['status'] == 'FAILED')
        
        self._write("")
        self.header("PIPELINE EXECUTION SUMMARY")
        
        self.subheader("Execution Timeline")
        self.info(f"Pipeline started: {self.start_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
        self.info(f"Pipeline ended:   {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}")
        self.info(f"Total duration:   {total_duration/60:.2f} minutes")
        
        self.subheader("Step-by-Step Results")
        self.info(f"{'Step':<6} {'Description':<50} {'Status':<10} {'Duration':<12}")
        self.info("-" * 80)
        
        for result in self.step_results:
            status_symbol = "✓" if result['status'] == 'SUCCESS' else "✗"
            duration_str = f"{result['duration_seconds']:.2f}s"
            desc = result['description'][:48]
            self.info(f"{result['step']:<6} {desc:<50} {status_symbol} {result['status']:<8} {duration_str:<12}")
        
        self.subheader("Statistics")
        self.info(f"Total steps:     {total_steps}")
        self.info(f"Successful:      {success_count}")
        self.info(f"Failed:          {fail_count}")
        self.info(f"Success rate:    {100*success_count/total_steps:.1f}%")

        # Generate Research Audit
        self.subheader("Research Accountability & Reproducibility")
        self.generate_research_audit()
        
        # Validate Data Integrity - Zero Synthetic Data Tolerance
        self.subheader("Data Integrity Validation")
        self.info("Checking all results for synthetic data contamination...")
        try:
            sys.path.insert(0, str(PROJECT_ROOT / 'scripts' / 'utils'))
            from data_integrity_validator import DataIntegrityValidator
            validator = DataIntegrityValidator(PROJECT_ROOT)
            
            all_clean = True
            for filename in DATA_INTEGRITY_REQUIRED_OUTPUTS:
                json_file = PROJECT_ROOT / 'results' / filename
                if not json_file.exists():
                    self.warning(f"  ⚠ {filename}: Required integrity target not found")
                    all_clean = False
                    continue
                try:
                    with open(json_file, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                    if not validator.validate_results(data, context=json_file.name):
                        all_clean = False
                        self.error(f"  ✗ {json_file.name}: Data contamination detected")
                        for v in validator.get_violations():
                            self.error(f"      - {v}")
                except Exception as e:
                    all_clean = False
                    self.error(f"  ✗ {json_file.name}: Could not validate ({e})")
            
            if all_clean:
                self.success("  ✓ All results files passed data integrity check")
            else:
                self.error("  ✗ DATA INTEGRITY VIOLATIONS DETECTED")
                self.error("  Pipeline results contain synthetic data contamination!")
                return False
        except Exception as e:
            self.warning(f"  ⚠ Could not run data integrity validator: {e}")
        
        if success_count == total_steps:
            self.subheader("Output Locations")
            self.info(f"Log file:         {self.log_file}")
            self.info(f"Results:          {PROJECT_ROOT / 'results'}")
            self.info(f"Key output:       results/step003_dsn_reanalysis.json")
            self.info(f"Key output:       results/step008_fitting_results.json")
            self.info(f"Key output:       results/step036_final_report.json")
            self.info(f"Key figure:       results/step007_figure1_altitude_anomaly.png")
            self._write("")
            self.success("PIPELINE COMPLETED SUCCESSFULLY")
            return True
        else:
            self._write("")
            self.error(f"PIPELINE FAILED - {fail_count} step(s) did not complete")
            return False


def run_step(filename: str, description: str, step_num: int, total_steps: int, 
             logger: PipelineLogger) -> bool:
    """Run a single pipeline step with detailed logging."""
    step_path = STEPS_DIR / filename
    
    if not step_path.exists():
        logger.error(f"File not found: {step_path}")
        logger.add_step_result(f"Step {step_num:03d}", description, "FAILED", 0.0, -1)
        return False
    
    # Record start time
    start_time = logger.step_start(step_num, total_steps, description)
    
    try:
        # Run the step
        result = subprocess.run(
            [sys.executable, str(step_path)],
            cwd=str(PROJECT_ROOT),
            capture_output=False,
            text=True
        )
        
        # Calculate duration and log result
        duration = logger.step_complete(description, start_time, 
                                       "SUCCESS" if result.returncode == 0 else "FAILED")
        
        # Record result
        logger.add_step_result(
            f"Step {step_num:03d}",
            description,
            "SUCCESS" if result.returncode == 0 else "FAILED",
            duration,
            result.returncode
        )
        
        return result.returncode == 0
        
    except Exception as e:
        duration = time.time() - start_time
        logger.error(f"Exception running {description}: {str(e)}")
        logger.add_step_result(f"Step {step_num:03d}", description, "FAILED", duration, -1)
        return False


def main():
    """Execute full pipeline with detailed logging."""
    logger = PipelineLogger()
    
    # Header
    logger.header("FLYBY TEP PIPELINE v3.0 - FULL EXECUTION", 80)
    logger.info(f"Project root: {PROJECT_ROOT}")
    logger.info(f"Steps dir:     {STEPS_DIR}")
    logger.info(f"Log file:      {logger.log_file}")
    
    logger.subheader("Enhancements in v3.0", 80)
    logger.info("  • WGS84 non-spherical Earth model with J2, J3, J4 harmonics")
    logger.info("  • PREM dynamic density mapping")
    logger.info("  • 3D trajectory integration framework")
    logger.info("  • Expanded dataset (12+ flybys with archival mining)")
    logger.info("  • Core raw DSN reanalysis with minimal OD (resolves Horizons circularity)")
    logger.info("  • Hierarchical Bayesian model (addresses extreme heterogeneity)")
    logger.info("  • Plasma modulation (addresses Cassini sign mismatch)")
    logger.info("  • First-principles Temporal Shear Suppression calculation (reduces systematic uncertainty)")
    logger.info("  • Bayesian model comparison with Bayes factors (statistical rigor)")
    
    total_steps = len(CORE_STEPS)
    
    # Run all steps
    logger.subheader("EXECUTING PIPELINE STEPS", 80)
    
    for i, (filename, description) in enumerate(CORE_STEPS, 1):
        if not run_step(filename, description, i, total_steps, logger):
            logger.error("Stopping pipeline due to failure")
            break
    
    # Final summary
    success = logger.final_summary(total_steps)
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
