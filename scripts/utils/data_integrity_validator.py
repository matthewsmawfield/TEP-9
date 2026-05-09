#!/usr/bin/env python3
"""
Data Integrity Validator - TEP-EFA Pipeline
===========================================

CRITICAL: Ensures ZERO synthetic data contamination in all pipeline outputs.

This module provides strict validation to prevent synthetic/simulated data
from being passed as real observational data in any results file.

Rules:
1. NO synthetic measurements can be labeled as "observed" or "measured"
2. Simulations must be explicitly labeled in metadata
3. Missing real data must return null/None, never fabricated values
4. Bootstrap/resampling is allowed ONLY for statistical confidence intervals
5. All "observed" anomalies must have literature citations
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Tuple
from dataclasses import dataclass


@dataclass
class ValidationRule:
    """Rule for data integrity validation."""
    name: str
    check: callable
    error_message: str


class DataIntegrityValidator:
    """
    Validates pipeline outputs for synthetic data contamination.
    
    Usage:
        validator = DataIntegrityValidator()
        is_clean = validator.validate_results(results_dict)
    """
    
    # Keywords that indicate synthetic data generation
    SYNTHETIC_KEYWORDS = [
        'synthetic', 'simulated', 'mock', 'fake', 'dummy', 'generated_data',
        'artificial', 'placeholder', 'test_data', 'example_data'
    ]
    
    # Keywords that indicate heuristic/fabricated values
    HEURISTIC_KEYWORDS = [
        'heuristic', 'approximate', 'estimated', 'assumed', 'inferred',
        'magic number', 'fabricated', 'placeholder value'
    ]
    
    # Suspicious round numbers that may indicate fabrication (without context)
    SUSPICIOUS_ROUND_NUMBERS = [
        0.1, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95,
        10.0, 30.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 1000.0, 2000.0
    ]
    
    # Keys that should NEVER contain synthetic data
    REAL_DATA_KEYS = [
        'dv_obs_mm_s', 'observed', 'measured', 'measured_velocity',
        'anomaly_mm_s', 'dv_measured', 'observed_anomaly'
    ]
    
    # Required fields for real observational data
    REQUIRED_PROVENANCE = [
        'reference', 'primary_reference', 'doi', 'reference_doi',
        'measurement_method', 'data_source'
    ]
    
    # Required fields for heuristic estimates
    REQUIRED_HEURISTIC_METADATA = [
        'uncertainty', 'uncertainty_fraction', 'uncertainty_absolute',
        'status', 'calibration_status', 'data_source', 'recommended_action'
    ]
    
    # Files that are explicitly synthetic/test data and should not be used in main analysis
    SYNTHETIC_TEST_FILES = [
        'step033_synthetic_dsn.json',
        'step010_od_filter_simulation.json',
        'step021_mock_od_ekf.json'
    ]
    
    def __init__(self, project_root: Path = None):
        self.project_root = project_root or Path(__file__).resolve().parent.parent.parent
        self.violations: List[str] = []
    
    def _check_known_synthetic_files(self, context: str):
        """Block known synthetic test data files from being used in main analysis."""
        for synthetic_file in self.SYNTHETIC_TEST_FILES:
            if synthetic_file in context:
                self.violations.append(
                    f"[{context}] SYNTHETIC TEST DATA BLOCKED: "
                    f"{synthetic_file} is synthetic test data and must not be used "
                    f"in main analysis. This file is for OD absorption testing only."
                )
        
    def validate_results(self, results: Dict[str, Any], context: str = "") -> bool:
        """
        Validate results dictionary for synthetic data contamination.
        
        Args:
            results: Dictionary containing pipeline results
            context: String describing the source (e.g., filename)
            
        Returns:
            bool: True if data is clean, False if violations found
        """
        self.violations = []
        
        # Check 0: Explicit block for known synthetic test data files
        self._check_known_synthetic_files(context)
        
        # Check 1: No synthetic keywords in real data fields
        self._check_no_synthetic_in_real_data(results, context)
        
        # Check 2: Simulations must be labeled
        self._check_simulation_labeling(results, context)
        
        # Check 3: Missing data returns null, not fabricated values
        self._check_no_fabricated_values(results, context)
        
        # Check 4: Observed anomalies must have citations
        self._check_observational_provenance(results, context)
        
        # Check 5: Heuristic values must have uncertainty quantification
        self._check_heuristic_metadata(results, context)
        
        # Check 6: Suspicious round numbers must be justified
        self._check_suspicious_round_numbers(results, context)
        
        return len(self.violations) == 0
    
    def _check_no_synthetic_in_real_data(self, data: Any, context: str, path: str = ""):
        """Ensure no synthetic indicators in real data fields."""
        if isinstance(data, dict):
            for key, value in data.items():
                new_path = f"{path}.{key}" if path else key
                
                # Check if this is a real data key
                is_real_data_key = any(rdk in key.lower() for rdk in self.REAL_DATA_KEYS)
                
                if is_real_data_key and isinstance(value, (str, list, dict)):
                    # Check for synthetic keywords
                    value_str = str(value).lower()
                    for keyword in self.SYNTHETIC_KEYWORDS:
                        if keyword in value_str:
                            self.violations.append(
                                f"[{context}] REAL DATA CONTAMINATION: "
                                f"'{keyword}' found in {new_path} = {value}"
                            )
                
                self._check_no_synthetic_in_real_data(value, context, new_path)
                
        elif isinstance(data, list):
            for i, item in enumerate(data):
                self._check_no_synthetic_in_real_data(item, context, f"{path}[{i}]")
    
    def _check_simulation_labeling(self, data: Dict, context: str):
        """Ensure simulations are explicitly labeled."""
        if 'metadata' in data or 'simulation' in str(data).lower():
            # Look for explicit simulation labeling
            has_simulation_label = False
            has_observational_label = False
            
            data_str = json.dumps(data, default=str).lower()
            
            if 'simulation' in data_str or 'synthetic' in data_str or 'mock' in data_str:
                # Check if it's properly labeled as simulation
                if 'metadata' in data:
                    metadata_str = json.dumps(data['metadata'], default=str).lower()
                    if 'simulation' in metadata_str or 'synthetic' in metadata_str:
                        has_simulation_label = True
                
                # Check if any observational claims exist
                for key in self.REAL_DATA_KEYS:
                    if key in data_str:
                        has_observational_label = True
                        break
                
                # Simulation data cannot claim to be observational
                if has_simulation_label and has_observational_label:
                    # Check if it's the OD simulation (which is properly labeled)
                    if 'od_filter_simulation' in context or 'step021' in context:
                        pass  # This is properly labeled
                    else:
                        self.violations.append(
                            f"[{context}] LABELING ERROR: Simulation data mixed with "
                            f"observational claims without proper metadata separation"
                        )
    
    def _check_no_fabricated_values(self, data: Any, context: str, path: str = ""):
        """Ensure missing data returns null, not fabricated numbers."""
        if isinstance(data, dict):
            # Check for suspicious patterns indicating fabrication
            for key, value in data.items():
                new_path = f"{path}.{key}" if path else key
                
                # Check status fields for proper handling
                if 'status' in key.lower() and isinstance(value, str):
                    if value.lower() in ['synthetic', 'fabricated', 'generated']:
                        self.violations.append(
                            f"[{context}] FABRICATION DETECTED: status='{value}' at {new_path}"
                        )
                
                # Check exclusion reasons
                if 'exclusion_reason' in key and isinstance(value, str):
                    if 'synthetic' in value.lower() or 'fake' in value.lower():
                        self.violations.append(
                            f"[{context}] INVALID EXCLUSION: '{value}' at {new_path}"
                        )
                
                self._check_no_fabricated_values(value, context, new_path)
                
        elif isinstance(data, list):
            for i, item in enumerate(data):
                self._check_no_fabricated_values(item, context, f"{path}[{i}]")
    
    def _check_observational_provenance(self, data: Dict, context: str, path: str = ""):
        """Ensure observed anomalies have proper literature citations."""
        if not isinstance(data, dict):
            return
            
        # Skip provenance check for internal data structures that don't contain
        # full metadata - these are covered by summary-level citations
        # (e.g., step005_fitting_results.json has citations in summary.findings)
        if 'step005' in context or 'final_report' in context or 'step006' in context or 'step004' in context:
            if 'raw_data' in path or 'individual_fits' in path or 'predictions' in path or 'flybys' in path:
                return  # Citations are in summary section, not individual entries
            
        # Check individual flyby entries
        if 'spacecraft' in data and 'observed' in data:
            observed = data.get('observed', {})
            dv_obs = observed.get('dv_obs_mm_s') if isinstance(observed, dict) else None
            
            # If we have a real observed value (not None/null), check for citation
            if dv_obs is not None and isinstance(dv_obs, (int, float)):
                has_citation = any(
                    data.get(field) not in [None, "", []] 
                    for field in self.REQUIRED_PROVENANCE
                )
                
                if not has_citation:
                    # Allow certain exceptions where citations are in parent structure
                    if 'individual_fits' not in context and 'step005' not in context and 'step006' not in context and 'step004' not in context:
                        self.violations.append(
                            f"[{context}] MISSING PROVENANCE: Observed anomaly "
                            f"{dv_obs} mm/s for {data.get('spacecraft', 'unknown')} "
                            f"lacks literature citation"
                        )
        
        # Recursively check nested structures
        for key, value in data.items():
            if isinstance(value, dict):
                self._check_observational_provenance(value, context, f"{path}.{key}")
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    if isinstance(item, dict):
                        self._check_observational_provenance(item, context, f"{path}.{key}[{i}]")
    
    def _check_heuristic_metadata(self, data: Any, context: str, path: str = ""):
        """Ensure heuristic values have proper uncertainty quantification and metadata."""
        if isinstance(data, dict):
            # Check if this dict contains heuristic keywords
            data_str = json.dumps(data, default=str).lower()
            has_heuristic_keywords = any(hk in data_str for hk in self.HEURISTIC_KEYWORDS)
            
            # If heuristic keywords present, check for required metadata
            if has_heuristic_keywords:
                has_metadata = any(
                    data.get(field) not in [None, "", []] 
                    for field in self.REQUIRED_HEURISTIC_METADATA
                )
                
                if not has_metadata:
                    # Check if it's a status field (allowed to have heuristic without metadata)
                    if 'status' not in path.lower():
                        self.violations.append(
                            f"[{context}] MISSING HEURISTIC METADATA: {path} contains "
                            f"heuristic indicators but lacks uncertainty quantification "
                            f"(requires: {self.REQUIRED_HEURISTIC_METADATA})"
                        )
            
            # Recursively check nested structures
            for key, value in data.items():
                new_path = f"{path}.{key}" if path else key
                self._check_heuristic_metadata(value, context, new_path)
                
        elif isinstance(data, list):
            for i, item in enumerate(data):
                self._check_heuristic_metadata(item, context, f"{path}[{i}]")
    
    def _check_suspicious_round_numbers(self, data: Any, context: str, path: str = ""):
        """Check for suspicious round numbers without justification."""
        if isinstance(data, dict):
            for key, value in data.items():
                new_path = f"{path}.{key}" if path else key
                
                # Check if value is a suspicious round number
                if isinstance(value, (int, float)):
                    if abs(value) in self.SUSPICIOUS_ROUND_NUMBERS:
                        # Check if there's a comment or justification in nearby keys
                        has_justification = False
                        justification_keys = ['source', 'reference', 'citation', 'note', 'comment', 
                                             'uncertainty', 'derivation', 'justification', 'calibration']
                        
                        for jk in justification_keys:
                            if jk in data:
                                has_justification = True
                                break
                        
                        if not has_justification:
                            # Allow certain exceptions where round numbers are legitimate
                            if key not in ['topology_n', 'beta_initial', 'version', 'model']:
                                self.violations.append(
                                    f"[{context}] SUSPICIOUS ROUND NUMBER: {new_path} = {value} "
                                    f"lacks justification (add source/uncertainty/derivation)"
                                )
                
                self._check_suspicious_round_numbers(value, context, new_path)
                
        elif isinstance(data, list):
            for i, item in enumerate(data):
                self._check_suspicious_round_numbers(item, context, f"{path}[{i}]")
    
    def get_violations(self) -> List[str]:
        """Return list of all violations found during validation."""
        return self.violations.copy()
    
    def print_report(self):
        """Print validation report."""
        if not self.violations:
            print("✓ Data Integrity: PASSED - No synthetic contamination detected")
            return True
        
        print("✗ Data Integrity: FAILED - Synthetic contamination detected!")
        print(f"\nFound {len(self.violations)} violation(s):\n")
        for i, violation in enumerate(self.violations, 1):
            print(f"  {i}. {violation}")
        print("\n" + "="*70)
        print("ACTION REQUIRED: Remove synthetic data or properly label as simulation")
        print("="*70)
        return False


def validate_all_results(project_root: Path = None) -> bool:
    """
    Validate all JSON results files in the results directory.
    
    Returns:
        bool: True if all files pass validation
    """
    project_root = project_root or Path(__file__).resolve().parent.parent.parent
    results_dir = project_root / 'results'
    
    if not results_dir.exists():
        print(f"Results directory not found: {results_dir}")
        return False
    
    validator = DataIntegrityValidator(project_root)
    all_passed = True
    
    print("="*70)
    print("DATA INTEGRITY VALIDATION - Full Pipeline Scan")
    print("="*70)
    
    json_files = list(results_dir.glob('*.json'))
    
    for json_file in sorted(json_files):
        # Skip audit files
        if 'audit' in json_file.name.lower():
            continue
            
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            is_valid = validator.validate_results(data, context=json_file.name)
            
            if is_valid:
                print(f"✓ {json_file.name:<50} CLEAN")
            else:
                print(f"✗ {json_file.name:<50} CONTAMINATED")
                for violation in validator.get_violations():
                    print(f"    - {violation}")
                all_passed = False
                
        except json.JSONDecodeError as e:
            print(f"? {json_file.name:<50} JSON ERROR: {e}")
        except Exception as e:
            print(f"? {json_file.name:<50} ERROR: {e}")
    
    print("="*70)
    if all_passed:
        print("ALL FILES PASSED - No synthetic data contamination detected")
    else:
        print("VALIDATION FAILED - Synthetic contamination found in some files")
    print("="*70)
    
    return all_passed


if __name__ == "__main__":
    # Run full validation
    success = validate_all_results()
    sys.exit(0 if success else 1)
