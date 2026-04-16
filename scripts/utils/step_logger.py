#!/usr/bin/env python3
"""
Step Logger - Highly verbose logging for pipeline steps

Provides detailed logging with timestamps, progress tracking,
and output file management for pipeline steps.

Each log entry is prefixed with the step number for easy filtering and tracing.
"""

import logging
import sys
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


def extract_step_number(step_name: str) -> str:
    """Extract step number from step name (e.g., 'step_001a_data_ingestion' -> '001A')."""
    match = re.search(r'step[_-]?(\d+[a-z]?)', step_name.lower())
    return match.group(1).upper() if match else "???"


class StepPrefixFormatter(logging.Formatter):
    """Custom formatter that adds step number prefix to all log messages."""
    
    def __init__(self, step_number: str, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)
        self.step_number = step_number
    
    def format(self, record):
        # Add step number prefix to the message
        record.msg = f"[STEP {self.step_number}] {record.msg}"
        return super().format(record)


class StepLogger:
    """Highly verbose logger for pipeline steps with step number prefixes."""
    
    def __init__(self, step_name: str, project_root: Optional[Path] = None):
        """
        Initialize step logger.
        
        Args:
            step_name: Name of the pipeline step (e.g., "step_001_data_ingestion")
            project_root: Project root directory (defaults to 3 levels up from __file__)
        """
        self.step_name = step_name
        self.step_number = extract_step_number(step_name)
        
        if project_root is None:
            # Default to 3 levels up from utils folder
            project_root = Path(__file__).resolve().parent.parent.parent
        
        self.project_root = project_root
        self.logs_dir = project_root / 'logs'
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        
        # Create log file with timestamp
        timestamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
        self.log_file = self.logs_dir / f"{step_name}_{timestamp}.log"
        
        # Track output files
        self.output_files = []
        self.calculation_counter = 0
        
        # Configure logging
        self._setup_logging()
        
        # Log initialization
        self.logger.info(f"LOGGER INITIALIZED: {step_name} (Step {self.step_number})")
        self.logger.debug(f"Log file: {self.log_file}")
        self.logger.debug(f"Project root: {self.project_root}")
        
    def _setup_logging(self):
        """Configure logging with both file and console handlers with step number prefixes."""
        # Create logger
        self.logger = logging.getLogger(self.step_name)
        self.logger.setLevel(logging.DEBUG)
        
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # File handler - DEBUG level for verbose output with step prefix
        file_handler = logging.FileHandler(self.log_file, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_formatter = StepPrefixFormatter(
            self.step_number,
            fmt='[%(asctime)s UTC] [%(levelname)-8s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)
        
        # Console handler - INFO level with step prefix
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_formatter = StepPrefixFormatter(
            self.step_number,
            fmt='[%(levelname)s] %(message)s'
        )
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
    
    def header(self, title: str, width: int = 80):
        """Log a formatted header."""
        self.logger.info("=" * width)
        self.logger.info(title)
        self.logger.info("=" * width)
    
    def subheader(self, title: str, width: int = 80):
        """Log a formatted subheader."""
        self.logger.info("-" * width)
        self.logger.info(title)
        self.logger.info("-" * width)
    
    def debug(self, message: str):
        """Log debug message (file only)."""
        self.logger.debug(message)
    
    def info(self, message: str):
        """Log info message."""
        self.logger.info(message)
    
    def warning(self, message: str):
        """Log warning message."""
        self.logger.warning(f"⚠ {message}")
    
    def error(self, message: str):
        """Log error message."""
        self.logger.error(f"✗ {message}")
    
    def success(self, message: str):
        """Log success message."""
        self.logger.info(f"✓ {message}")
    
    def section(self, title: str):
        """Log a section header."""
        self.logger.debug("")
        self.logger.debug("=" * 80)
        self.logger.debug(f"SECTION: {title}")
        self.logger.debug("=" * 80)
    
    def subsection(self, title: str):
        """Log a subsection header."""
        self.logger.debug("")
        self.logger.debug(f"--- {title} ---")
    
    def progress(self, message: str):
        """Log progress message."""
        self.logger.info(f"  → {message}")
    
    def data(self, message: str):
        """Log data processing message (file only)."""
        self.logger.debug(f"  DATA: {message}")
    
    def parameter(self, param_name: str, value):
        """Log a parameter value (file only)."""
        self.logger.debug(f"  PARAM: {param_name} = {value}")
    
    def calculation(self, calc_name: str, inputs: dict = None, formula: str = None, 
                    result=None, units: str = None):
        """
        Log a detailed calculation with inputs, formula, and result.
        
        Args:
            calc_name: Name/description of the calculation
            inputs: Dictionary of input parameters {name: value}
            formula: Mathematical formula used
            result: Calculation result
            units: Units of the result
        """
        self.calculation_counter += 1
        calc_id = f"CALC-{self.step_number}-{self.calculation_counter:04d}"
        
        self.logger.debug(f"[{calc_id}] CALCULATION: {calc_name}")
        
        if inputs:
            for name, value in inputs.items():
                self.logger.debug(f"[{calc_id}]   INPUT: {name} = {value}")
        
        if formula:
            self.logger.debug(f"[{calc_id}]   FORMULA: {formula}")
        
        if result is not None:
            units_str = f" {units}" if units else ""
            self.logger.debug(f"[{calc_id}]   RESULT: {result}{units_str}")
        
        return calc_id
    
    def formula(self, name: str, expression: str, variables: dict = None):
        """Log a formula application with variable values."""
        self.logger.debug(f"FORMULA [{name}]:")
        self.logger.debug(f"  Expression: {expression}")
        if variables:
            for var, val in variables.items():
                self.logger.debug(f"    {var} = {val}")
    
    def iteration(self, iteration_num: int, total: int, message: str = ""):
        """Log iteration progress in loops."""
        pct = 100.0 * iteration_num / total if total > 0 else 0
        msg = f" | {message}" if message else ""
        self.logger.debug(f"ITERATION [{iteration_num}/{total}, {pct:.1f}%]{msg}")
    
    def intermediate(self, name: str, value, units: str = None):
        """Log an intermediate calculation value."""
        units_str = f" {units}" if units else ""
        self.logger.debug(f"  INTERMEDIATE: {name} = {value}{units_str}")
    
    def array_data(self, name: str, arr, max_display: int = 5):
        """Log array data (truncated for large arrays)."""
        arr = list(arr) if hasattr(arr, '__iter__') and not isinstance(arr, str) else [arr]
        n = len(arr)
        if n <= max_display:
            self.logger.debug(f"  ARRAY [{name}]: {arr} (n={n})")
        else:
            display = arr[:max_display] + ['...'] + arr[-max_display:]
            self.logger.debug(f"  ARRAY [{name}]: {display} (n={n}, showing first/last {max_display})")
    
    def comparison(self, name: str, expected, actual, tolerance: float = None):
        """Log a comparison between expected and actual values."""
        diff = abs(actual - expected) if isinstance(expected, (int, float)) and isinstance(actual, (int, float)) else None
        match = "MATCH" if (diff is not None and tolerance is not None and diff <= tolerance) else "N/A"
        self.logger.debug(f"  COMPARISON [{name}]: expected={expected}, actual={actual}, diff={diff}, {match}")
    
    def add_output_file(self, file_path: Path, description: str = ""):
        """Register an output file created by this step."""
        self.output_files.append({
            'path': str(file_path),
            'description': description,
            'size_bytes': file_path.stat().st_size if file_path.exists() else 0
        })
        self.logger.info(f"  Output: {file_path.name} ({description})")
    
    def log_output_summary(self):
        """Log summary of all output files."""
        if not self.output_files:
            self.logger.info("No output files generated")
            return
        
        self.logger.info("")
        self.logger.info(f"Generated {len(self.output_files)} output file(s):")
        for output in self.output_files:
            path = Path(output['path'])
            size_kb = output['size_bytes'] / 1024
            desc = f" - {output['description']}" if output['description'] else ""
            self.logger.info(f"  • {path.name}: {size_kb:.2f} KB{desc}")
    
    def log_step_summary(self, duration_seconds: float, status: str = "SUCCESS"):
        """Log final step summary with detailed statistics."""
        duration_str = f"{duration_seconds:.2f}s" if duration_seconds < 60 else f"{duration_seconds/60:.2f}m"
        
        self.logger.info("")
        self.header(f"STEP SUMMARY: {self.step_name}")
        self.logger.info(f"Step Number: {self.step_number}")
        self.logger.info(f"Status: {status}")
        self.logger.info(f"Duration: {duration_str}")
        self.logger.info(f"Calculations Logged: {self.calculation_counter}")
        self.logger.info(f"Output Files: {len(self.output_files)}")
        self.logger.info(f"Log file: {self.log_file}")
        self.log_output_summary()
        self.logger.info("=" * 80)
        self.logger.info(f"STEP {self.step_number} COMPLETE: {status}")
    
    def data_load(self, file_path: Path, description: str = ""):
        """Log data loading operation."""
        self.logger.info(f"Loading data: {file_path.name}")
        self.logger.debug(f"  Full path: {file_path}")
        if description:
            self.logger.debug(f"  Description: {description}")
    
    def data_save(self, file_path: Path, description: str = "", record_count: int = None):
        """Log data saving operation."""
        self.logger.info(f"Saving data: {file_path.name}")
        self.logger.debug(f"  Full path: {file_path}")
        if description:
            self.logger.debug(f"  Description: {description}")
        if record_count is not None:
            self.logger.debug(f"  Records: {record_count}")
    
    def metric(self, name: str, value, units: str = None):
        """Log a metric value."""
        units_str = f" {units}" if units else ""
        self.logger.info(f"METRIC: {name} = {value}{units_str}")
    
    def threshold_check(self, name: str, value, threshold, passed: bool, 
                       operator: str = ">"):
        """Log a threshold check result."""
        status = "PASS" if passed else "FAIL"
        self.logger.info(f"THRESHOLD [{name}]: {value} {operator} {threshold} → {status}")
    
    def cross_reference(self, source_step: str, data_type: str, 
                       source_file: str = None):
        """Log cross-reference to another step's output."""
        msg = f"Cross-reference from {source_step}: {data_type}"
        if source_file:
            msg += f" (from {source_file})"
        self.logger.info(msg)
    
    def get_log_file(self) -> Path:
        """Return the log file path."""
        return self.log_file


def get_step_logger(step_name: str, project_root: Optional[Path] = None) -> StepLogger:
    """
    Convenience function to get a step logger.
    
    Args:
        step_name: Name of the pipeline step
        project_root: Project root directory
    
    Returns:
        StepLogger instance
    """
    return StepLogger(step_name, project_root)
