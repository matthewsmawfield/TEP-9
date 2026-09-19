"""Shared worker-count policy for the compute-heavy steps.

The heavy steps (raw-astrometry refits, boundary propagation, injection
batteries) parallelize over independent comets with a fork Pool.  The
right default on Apple Silicon is the number of *performance* cores:
the efficiency cores are slower per-thread and reserving one P-core for
the OS/GUI keeps long runs responsive.  On machines without a
perf-level sysctl (Linux, Intel Macs) all cores are uniform, so the
default is cpu_count() - 1.

Every step still honours --workers N on the command line; this helper
only changes the default when the flag is absent.
"""

import os
import subprocess


def _perf_cores():
    """Physical performance-core count on Apple Silicon, else None.

    The perflevel index is not stable across Apple chips: on the M5 Pro
    perflevel0 is the efficiency cluster ("Super") and perflevel1 the
    performance cluster, while on M1/M2 it is the reverse order.  The
    performance cluster always has at least as many cores as the
    efficiency cluster on big.LITTLE Apple Silicon, so take the max.
    """
    best = 0
    for lvl in range(4):
        try:
            out = subprocess.run(
                ["sysctl", "-n", f"hw.perflevel{lvl}.physicalcpu"],
                capture_output=True, text=True, timeout=2)
            if out.returncode == 0:
                best = max(best, int(out.stdout.strip()))
        except Exception:
            pass
    return best or None


def default_workers(reserve=1):
    """Worker count for a CPU-bound fork pool on this machine."""
    n = _perf_cores()
    if n is None:
        n = os.cpu_count() or 1
    return max(1, n - reserve)


def cli_workers(argv, default=None):
    """--workers N override; otherwise the machine default."""
    if "--workers" in argv:
        return max(1, int(argv[argv.index("--workers") + 1]))
    return default if default is not None else default_workers()
