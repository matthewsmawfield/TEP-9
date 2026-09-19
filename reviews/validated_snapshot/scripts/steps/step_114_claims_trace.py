#!/usr/bin/env python3
"""Exact result-leaf publication audit, run after the analysis steps.

Matching any nearby numeric value in a cited file is not provenance. The
replacement names each headline source leaf and exercises the real renderer.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.audits.audit_publication import main

if __name__ == '__main__':
    sys.exit(main())
