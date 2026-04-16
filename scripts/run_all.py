#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Run All Steps (Enhanced v3.0)

Executes the complete flyby analysis pipeline with detailed logging:
  1. Download SPICE kernels (NAIF)
  2. Convert SPICE to JSON format
  3. Fetch JPL Horizons data
  4. Data ingestion from JPL Horizons
  5. Archival data mining (sample size expansion)
  6. DSN data framework (addressing circularity)
  7. Enhanced TEP chameleon model (WGS84, PREM, 3D integration)
  8. Parameter fitting with small-sample statistics
  9. Comprehensive diagnostics and validation
  10. Enhanced validation and model comparison
  11. Hierarchical Bayesian model (addresses extreme heterogeneity)
  12. Plasma modulation (addresses Cassini sign mismatch)
  13. First-principles chameleon calculation (reduces systematic uncertainty)
  14. Bayesian model comparison with Bayes factors (statistical rigor)
  15. Enhanced final report
  16. Figure generation
  17. TEP suppression analysis

All steps produce detailed logs with timestamps and execution metrics.
"""

import sys
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Tuple, Dict

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STEPS_DIR = PROJECT_ROOT / 'scripts' / 'steps'
LOGS_DIR = PROJECT_ROOT / 'logs'

# Ensure logs directory exists
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# All pipeline steps (run_all runs everything by default)
CORE_STEPS: List[Tuple[str, str]] = [
    ('step_001a_download_spice.py', 'Step 001a: Download SPICE Kernels'),
    ('step_001b_spice_to_json.py', 'Step 001b: Convert SPICE to JSON'),
    ('step_010_jpl_horizons_fetch.py', 'Step 010: JPL Horizons Data Fetch'),
    ('step_001_data_ingestion.py', 'Step 001: Data Ingestion (JPL Horizons)'),
    ('step_002_archival_data_mining.py', 'Step 002: Archival Data Mining (Sample Expansion)'),
    ('step_003_dsn_data_ingestion.py', 'Step 003: DSN Data Framework (Addressing Circularity)'),
    ('step_004_tep_model.py', 'Step 004: Enhanced TEP Chameleon Model'),
    ('step_005_fitting.py', 'Step 005: Parameter Fitting'),
    ('step_005b_diagnostics.py', 'Step 005b: Comprehensive Diagnostics and Validation'),
    ('step_005c_enhanced_validation.py', 'Step 005c: Enhanced Validation and Model Comparison'),
    ('step_005d_rigorous_statistics.py', 'Step 005d: Rigorous Statistical Tests'),
    ('step_009b_read_trk234_data.py', 'Step 009b: Read TRK-2-34 Data Format'),
    ('step_010_hierarchical_bayesian.py', 'Step 010: Hierarchical Bayesian Model (Addresses Heterogeneity)'),
    ('step_011_gnss_cross_validation.py', 'Step 011: GNSS Cross-Validation'),
    ('step_012_saturation_model.py', 'Step 012: Saturation Model Analysis'),
    ('step_013_dsn_processing.py', 'Step 013: DSN Data Processing Framework'),
    ('step_014_plasma_modulation.py', 'Step 014: Plasma Modulation (Addresses Cassini Sign Mismatch)'),
    ('step_015_chameleon_first_principles.py', 'Step 015: First-Principles Chameleon (Reduces Systematic Uncertainty)'),
    ('step_016_bayesian_model_comparison.py', 'Step 016: Bayesian Model Comparison (Statistical Rigor)'),
    ('step_017_3d_field_integration.py', 'Step 017: 3D Field Integration'),
    ('step_017_trajectory_integration.py', 'Step 017b: Trajectory Integration'),
    ('step_018_enhanced_bayesian.py', 'Step 018: Enhanced Bayesian Analysis'),
    ('step_019_cross_validation.py', 'Step 019: Cross-Validation Analysis'),
    ('step_020_sensitivity_analysis.py', 'Step 020: Sensitivity Analysis'),
    ('step_006_report.py', 'Step 021: Enhanced Final Report'),
    ('step_007_visualizations.py', 'Step 022: Figure Generation'),
    ('step_008_tep_suppression.py', 'Step 023: TEP Suppression Analysis'),
]


class PipelineLogger:
    """Handles consistent logging for pipeline execution."""
    
    def __init__(self):
        self.start_time = datetime.now(timezone.utc)
        self.log_file = LOGS_DIR / f"pipeline_{self.start_time.strftime('%Y%m%d_%H%M%S')}.log"
        self.step_results: List[Dict] = []
        
    def _write(self, message: str, level: str = "INFO"):
        """Write message to both console and log file."""
        timestamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        log_line = f"[{timestamp}] [{level:8}] {message}"
        
        # Print to console
        print(message)
        
        # Write to log file
        with open(self.log_file, 'a') as f:
            f.write(log_line + '\n')
    
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
    
    def final_summary(self, total_steps: int) -> bool:
        """Print detailed final summary."""
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
        
        if success_count == total_steps:
            self.subheader("Output Locations")
            self.info(f"Log file:         {self.log_file}")
            self.info(f"Results:          {PROJECT_ROOT / 'data' / 'processed'}")
            self.info(f"Figures:          {PROJECT_ROOT / 'site' / 'public' / 'figures'}")
            self.info(f"Archival catalog: data/processed/archival_flyby_catalog.json")
            self.info(f"DSN framework:    data/processed/dsn_ingestion_framework.json")
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
    logger.info("  • Raw DSN data access framework (addressing circularity)")
    logger.info("  • Hierarchical Bayesian model (addresses extreme heterogeneity)")
    logger.info("  • Plasma modulation (addresses Cassini sign mismatch)")
    logger.info("  • First-principles chameleon calculation (reduces systematic uncertainty)")
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
