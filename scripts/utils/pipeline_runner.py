#!/usr/bin/env python3
"""
Shared pipeline runner utility for TEP-3I sub-pipelines.

Provides run_step() and run_pipeline() so that the focused runner scripts
do not each duplicate this logic.

M4 Pro Optimizations:
- Parallel step execution where dependencies allow
- Phase-based execution (ingestion → transform → analysis)
"""

import datetime
import subprocess
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

PROJECT_ROOT = Path(__file__).resolve().parents[2]
STEPS_DIR    = PROJECT_ROOT / "scripts" / "steps"

N_JOBS = multiprocessing.cpu_count()


def _fmt_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    return f"{m}m{s:02d}s"


def run_step(script_name: str, step_idx: int, total: int, label: str = "STEP") -> dict:
    """
    Run a single pipeline step as a subprocess.
    Streams output to both console and log file for verbose real-time monitoring.

    Returns a dict with keys: name, status ("PASS" | "FAIL"), elapsed_s, returncode.
    Never raises — failures are captured in the return dict.
    """
    script_path = STEPS_DIR / script_name
    phase       = f"[{step_idx:>3}/{total}]"

    print(f"\n{'─'*70}")
    print(f" {phase}  {label}: {script_name}")
    print(f"{'─'*70}\n")

    # Create log file for this step
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{script_name.replace('.py', '.log')}"

    print(f"[PIPELINE] Step {step_idx}/{total}: Initializing...")
    print(f"[PIPELINE] Script: {script_path}")
    print(f"[PIPELINE] Log file: {log_file}")
    print(f"[PIPELINE] Starting execution at {datetime.datetime.now().strftime('%H:%M:%S')}")
    print()

    t0 = time.perf_counter()
    
    # Use Popen for real-time streaming to both console and log
    with open(log_file, 'w') as log_fh:
        process = subprocess.Popen(
            [sys.executable, str(script_path)],
            cwd=str(PROJECT_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        
        # Stream output to both console and log file
        for line in process.stdout:
            line = line.rstrip()
            print(f"  {line}")  # Console output with indentation
            log_fh.write(line + "\n")
            log_fh.flush()
    
    process.wait()
    elapsed = time.perf_counter() - t0

    status = "PASS" if process.returncode == 0 else "FAIL"
    if status == "PASS":
        print(f"\n✓  PASS  {script_name}  ({_fmt_elapsed(elapsed)})")
    else:
        print(f"\n✗  FAIL  {script_name}  rc={process.returncode}  ({_fmt_elapsed(elapsed)})")
    print(f"    Log: {log_file}")

    return {
        "name":        script_name,
        "status":      status,
        "elapsed_s":   round(elapsed, 2),
        "returncode":  process.returncode,
    }

    t0 = time.perf_counter()
    
    # Use Popen for real-time streaming to both console and log
    with open(log_file, 'w') as log_fh:
        process = subprocess.Popen(
            [sys.executable, str(script_path)],
            cwd=str(PROJECT_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        
        # Stream output to both console and log file
        for line in process.stdout:
            line = line.rstrip()
            print(f"  {line}")  # Console output with indentation
            log_fh.write(line + "\n")
            log_fh.flush()
    
    process.wait()
    elapsed = time.perf_counter() - t0

    status = "PASS" if process.returncode == 0 else "FAIL"
    if status == "PASS":
        print(f"\n✓  PASS  {script_name}  ({_fmt_elapsed(elapsed)})")
    else:
        print(f"\n✗  FAIL  {script_name}  rc={process.returncode}  ({_fmt_elapsed(elapsed)})")
    print(f"    Log: {log_file}")

    return {
        "name":        script_name,
        "status":      status,
        "elapsed_s":   round(elapsed, 2),
        "returncode":  process.returncode,
    }


def run_pipeline(
    pipeline_name: str,
    steps: list,
    description: str = "",
    stop_on_failure: bool = False,
) -> list:
    """
    Run a list of step scripts and print a summary table.

    Parameters
    ----------
    pipeline_name  : Human-readable name for the header banner.
    steps          : List of script filenames (e.g. "step_001_uncover_load.py").
    description    : Optional one-line description printed in the banner.
    stop_on_failure: If True, stop at first failure (default: continue all).

    Returns
    -------
    List of result dicts (one per step).
    """
    wall_start = time.perf_counter()
    n = len(steps)

    print("╔" + "═" * 68 + "╗")
    print(f"║  TEP-JWST: {pipeline_name:<57}  ║")
    if description:
        print(f"║  {description:<66}  ║")
    print(f"║  Steps: {n:<60}  ║")
    print(f"║  Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'):<58}  ║")
    print("╚" + "═" * 68 + "╝\n")

    results = []
    for i, step in enumerate(steps, start=1):
        r = run_step(step, i, n, label=pipeline_name)
        results.append(r)
        if stop_on_failure and r["status"] != "PASS":
            print(f"\n  ⚠  Pipeline stopped at {step} (stop_on_failure=True)")
            break

    total_elapsed = time.perf_counter() - wall_start
    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_fail = sum(1 for r in results if r["status"] == "FAIL")

    print()
    print("╔" + "═" * 68 + "╗")
    print(f"║  {pipeline_name} — COMPLETE" + " " * (42 - len(pipeline_name)) + "  ║")
    print("╠" + "═" * 68 + "╣")
    for r in results:
        icon = "✓" if r["status"] == "PASS" else "✗"
        t    = _fmt_elapsed(r["elapsed_s"])
        name = r["name"][:52]
        stat = r["status"]
        print(f"║  {icon} {name:<52}  {stat:<8}  {t:>5}  ║")
    print("╠" + "═" * 68 + "╣")
    pct = 100 * n_pass // len(results) if results else 0
    print(f"║  PASS {n_pass}/{len(results)} ({pct}%)   FAIL: {n_fail}   "
          f"Total: {_fmt_elapsed(total_elapsed):<37}  ║")
    print("╚" + "═" * 68 + "╝")

    return results


def run_pipeline_parallel(
    pipeline_name: str,
    steps: list,
    description: str = "",
) -> list:
    """
    Run pipeline with M4 Pro parallel optimizations.
    
    Parallel groups:
    - Phase 1: step_000 (sequential - data fetch)
    - Phase 2: step_001, step_002 (parallel - independent transforms)
    - Phase 3: step_003 (sequential - needs 000 output, parallel network I/O)
    - Phase 4: steps_004-008 (parallel - independent analysis)
    """
    wall_start = time.perf_counter()
    
    print("╔" + "═" * 68 + "╗")
    print(f"║  TEP-3I M4 PRO PIPELINE{'':<40}  ║")
    print(f"║  {pipeline_name:<57}  ║")
    if description:
        print(f"║  {description:<66}  ║")
    print(f"║  Steps: {len(steps):<60}  ║")
    print(f"║  Cores: {N_JOBS:<60}  ║")
    print(f"║  Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'):<58}  ║")
    print("╚" + "═" * 68 + "╝\n")
    
    results = []
    
    # Phase 1: step_000 (must run first)
    print("\n" + "="*70)
    print("PHASE 1: Data Ingestion (Sequential)")
    print("="*70)
    r = run_step("step_000_jpl_horizons_ingestion.py", 1, len(steps), "INGESTION")
    results.append(r)
    if r["status"] != "PASS":
        print("\n✗ Phase 1 failed - aborting pipeline")
        return results
    
    # Phase 2: step_001, step_002 (can run in parallel)
    print("\n" + "="*70)
    print("PHASE 2: Frame Transformations (Parallel)")
    print("="*70)
    phase2_steps = [
        ("step_001_kinematic_transformation.py", 2, "CMB BOOST"),
        ("step_002_anisotropy_analysis.py", 3, "ANISOTROPY"),
    ]
    with ProcessPoolExecutor(max_workers=2) as executor:
        futures = {
            executor.submit(run_step, script, idx, len(steps), label): (script, idx, label)
            for script, idx, label in phase2_steps
        }
        for future in as_completed(futures):
            script, idx, label = futures[future]
            try:
                r = future.result()
                results.append(r)
            except Exception as e:
                print(f"\n✗ Step {label} failed: {e}")
                results.append({"name": script, "status": "FAIL", "elapsed_s": 0, "returncode": -1})
    
    # Phase 3: step_003 (high-res data fetch - uses async I/O)
    print("\n" + "="*70)
    print("PHASE 3: High-Resolution Data (Async I/O)")
    print("="*70)
    r = run_step("step_003_hires_nga_timeseries.py", 4, len(steps), "HIRES NGA")
    results.append(r)
    if r["status"] != "PASS":
        print("\n✗ Phase 3 failed - aborting")
        return results
    
    # Phase 4: steps_004-008 (analysis - run in parallel)
    print("\n" + "="*70)
    print("PHASE 4: TEP Analysis (Parallel)")
    print("="*70)
    
    phase4_steps = [
        ("step_004_tep_frame_discriminator.py", 5, "FRAME DISC"),
        ("step_005_dt_analysis.py", 6, "DT ANALYSIS"),
        ("step_006_comprehensive_tep.py", 7, "COMPREHENSIVE"),
        ("step_007_corrected_tep_physics.py", 8, "TEP PHYSICS"),
        ("step_008_irrefutable_evidence.py", 9, "EVIDENCE"),
    ]
    
    with ProcessPoolExecutor(max_workers=min(5, N_JOBS)) as executor:
        futures = {
            executor.submit(run_step, script, idx, len(steps), label): (script, idx, label)
            for script, idx, label in phase4_steps
        }
        for future in as_completed(futures):
            script, idx, label = futures[future]
            try:
                r = future.result()
                results.append(r)
            except Exception as e:
                print(f"\n✗ Step {label} failed: {e}")
                results.append({"name": script, "status": "FAIL", "elapsed_s": 0, "returncode": -1})
    
    # Sort all results by step name for consistent ordering
    results = sorted(results, key=lambda x: x.get("name", ""))
    
    # Summary
    total_elapsed = time.perf_counter() - wall_start
    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_fail = sum(1 for r in results if r["status"] == "FAIL")
    
    print("\n")
    print("╔" + "═" * 68 + "╗")
    print(f"║  {pipeline_name} — COMPLETE (M4 PRO)" + " " * (26 - len(pipeline_name)) + "  ║")
    print("╠" + "═" * 68 + "╣")
    for r in results:
        icon = "✓" if r["status"] == "PASS" else "✗"
        t = _fmt_elapsed(r["elapsed_s"])
        name = r["name"][:52]
        stat = r["status"]
        print(f"║  {icon} {name:<52}  {stat:<8}  {t:>5}  ║")
    print("╠" + "═" * 68 + "╣")
    pct = 100 * n_pass // len(results) if results else 0
    print(f"║  PASS {n_pass}/{len(results)} ({pct}%)   FAIL: {n_fail}   "
          f"Total: {_fmt_elapsed(total_elapsed):<37}  ║")
    print("╚" + "═" * 68 + "╝")
    
    return results
