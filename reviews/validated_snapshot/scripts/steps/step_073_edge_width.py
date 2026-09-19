#!/usr/bin/env python3
"""Step 073 -- boundary edge width.

The origin question (Section 6.8) turns on edge morphology: a
topological-defect wall has a sharp threshold while a plasma-sourced
gradient is diffuse.  The transition in the unexplained slip across the
axis cap is fitted with a logistic edge
    dtau(theta) = A * sigmoid((theta0 - theta)/w) + c
against the step-function alternative (w -> 0).  The 10-90% transition
width 4.39w and its bootstrap interval quantify "sharp" in degrees of
axis separation -- and, through the measured slip-vs-radius profile
(step_071), in AU of boundary penetration.

Inputs : results/step_b30_proper_time_slip.csv
Outputs: results/step_b38_edge_width.json
         results/figures/step_b38_edge_width.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import RESULTS, tee_stdout
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_073_edge_width")
tee_stdout(logger)
logger.header("Boundary edge width: logistic fit to the slip transition")

import csv
import json
import numpy as np
from scipy.optimize import curve_fit

SEED  = 20260918
NBOOT = 20000
rng   = np.random.default_rng(SEED)

rows = list(csv.DictReader(open(RESULTS / "step_b30_proper_time_slip.csv")))

def load(sub_rows):
    th, dt, ru = [], [], []
    for r in sub_rows:
        try:
            th.append(float(r["theta"]))
            dt.append(float(r["dtau_unexplained"]))
            ru.append(float(r["dtau_total"]))
        except (ValueError, KeyError):
            continue
    return np.array(th), np.array(dt), np.array(ru)

def logistic(th, A, th0, w, c):
    z = np.clip((th0 - th) / np.maximum(w, 1e-6), -60, 60)
    return A / (1.0 + np.exp(-z)) + c

def fit_edge(th, y):
    """logistic fit + step fit; returns params + RSS comparison."""
    if len(th) < 15:
        return None
    try:
        p0 = [np.percentile(y, 80) - np.percentile(y, 20),
              np.median(th), 8.0, np.percentile(y, 20)]
        bounds = ([-50, 0, 0.2, -50], [50, 180, 60, 50])
        popt, _ = curve_fit(logistic, th, y, p0=p0, bounds=bounds,
                            maxfev=20000)
        rss_l = float(np.sum((y - logistic(th, *popt)) ** 2))
    except Exception:
        return None
    # best step: threshold scanned, two means
    best = (np.inf, None)
    for ths in np.linspace(th.min() + 2, th.max() - 2, 400):
        lo, hi = y[th < ths], y[th >= ths]
        if len(lo) < 3 or len(hi) < 3:
            continue
        rss = np.sum((lo - lo.mean()) ** 2) + np.sum((hi - hi.mean()) ** 2)
        if rss < best[0]:
            best = (float(rss), float(ths))
    rss_s, ths = best
    # bootstrap the logistic width
    A_, th0_, w_, c_ = popt
    idx = rng.integers(0, len(th), (NBOOT, len(th)))
    ws, ths_b, fails = [], [], 0
    for bi in idx:
        try:
            pb, _ = curve_fit(logistic, th[bi], y[bi],
                              p0=[A_, th0_, w_, c_], bounds=bounds,
                              maxfev=10000)
            ws.append(pb[2]); ths_b.append(pb[1])
        except Exception:
            fails += 1
    ws = np.array(ws); ths_b = np.array(ths_b)
    return {"A": float(A_), "theta0": float(th0_), "w": float(w_),
            "c": float(c_),
            "width_10_90": float(4.394449 * w_),
            "w_ci68": [float(np.percentile(ws, 16)),
                       float(np.percentile(ws, 84))],
            "theta0_ci68": [float(np.percentile(ths_b, 16)),
                            float(np.percentile(ths_b, 84))],
            "width10_90_ci68": [float(4.394449 * np.percentile(ws, 16)),
                                float(4.394449 * np.percentile(ws, 84))],
            "boot_ok": int(len(ws)), "boot_fail": int(fails),
            "rss_logistic": rss_l, "rss_step": rss_s,
            "step_threshold": ths,
            "step_within_1se": bool(rss_l <= rss_s * 1.05),
            "n": int(len(th))}

res = {"method": "logistic edge fit to unexplained slip vs axis "
                 "separation; 10-90% width = 4.394w; step-function RSS "
                 "comparison; bootstrap 68% intervals (seed 20260918).",
       "seed": SEED}

for co in ("code", "warsaw", "pooled"):
    sub = [r for r in rows if co == "pooled" or r["cohort"] == co]
    th, dt, _ = load(sub)
    out = fit_edge(th, dt)
    res[f"edge_{co}"] = out
    if out:
        logger.info(f"{co}: theta0={out['theta0']:.1f} deg "
                    f"w={out['w']:.1f} deg  10-90% width="
                    f"{out['width_10_90']:.1f} deg "
                    f"[{out['width10_90_ci68'][0]:.1f},"
                    f"{out['width10_90_ci68'][1]:.1f}] "
                    f"RSS log/step={out['rss_logistic']:.1f}/"
                    f"{out['rss_step']:.1f}")

with open(RESULTS / "step_b38_edge_width.json", "w") as f:
    json.dump(res, f, indent=1)
print(f"wrote {RESULTS / 'step_b38_edge_width.json'}")

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sub = rows
th, dt, _ = load(sub)
fig, ax = plt.subplots(figsize=(6.4, 4.4))
inc = th < 60.0
ax.scatter(th[~inc], dt[~inc], s=20, facecolors="none",
           edgecolors="0.55", label="outside cap")
ax.scatter(th[inc], dt[inc], s=24, c="crimson", label="inside cap")
e = res["edge_pooled"]
if e:
    xs = np.linspace(th.min(), th.max(), 400)
    ax.plot(xs, logistic(xs, e["A"], e["theta0"], e["w"], e["c"]),
            "k-", lw=1.5,
            label=f"logistic: $\\theta_0$={e['theta0']:.0f}°, "
                  f"width={e['width_10_90']:.0f}°")
ax.axvline(60, color="0.4", ls=":", lw=1)
ax.set_xlabel(r"aphelion--axis separation $\theta$ (deg)")
ax.set_ylabel(r"unexplained slip $\delta\tau$ (yr)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("edge morphology: sharp threshold vs graded transition",
             fontsize=10)
fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b38_edge_width.png", dpi=150)
print(f"wrote {FIG / 'step_b38_edge_width.png'}")
