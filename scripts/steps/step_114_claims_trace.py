#!/usr/bin/env python3
"""Exact result-leaf publication audit, run after the analysis steps.

Matching any nearby numeric value in a cited file is not provenance. The
replacement names each headline source leaf and exercises the real renderer.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.audits.audit_publication import main

logger = StepLogger("step_114_claims_trace")
tee_stdout(logger)
logger.header("Manuscript claims-trace audit")

if __name__ == '__main__':
    rc = main()
    results_dir = Path(__file__).resolve().parents[2] / "results"
    for name in ("step_b78_claims_trace.json", "step_b78_claims_trace.csv"):
        out = results_dir / name
        if out.exists():
            logger.data_save(out, "claims-trace audit report")
    sys.exit(rc)
