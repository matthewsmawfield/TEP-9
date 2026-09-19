#!/usr/bin/env python3
"""Step 072 -- epoch stability of the implied proper-time slip.

A primordial boundary is inertially fixed: the unexplained slip must be
independent of the perihelion epoch.  A heliospheric/interface boundary
advects with the Local Interstellar Cloud flow (~26 km/s), and a
catalogue quality systematic would relax as orbit determination
improved -- both predict an epoch trend.  The measured per-comet
unexplained slip (step_065 residual after the encounter-budget
regression) is tested against perihelion year, inside and outside the
axis cap, per cohort and pooled, and across the photographic->CCD era
boundary.

Inputs : results/step_b30_proper_time_slip.csv
Outputs: results/step_b37_epoch_stability.json
         results/figures/step_b37_epoch_stability.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import RESULTS, tee_stdout
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_072_epoch_stability")
tee_stdout(logger)
logger.header("Epoch stability of the implied proper-time slip")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu

CAP   = 60.0
SEED  = 20260918
NBOOT = 20000
rng   = np.random.default_rng(SEED)

rows = list(csv.DictReader(open(RESULTS / "step_b30_proper_time_slip.csv")))
recs = []
for r in rows:
    m = re.match(r"[CP]*/(\d{4})", r["desig"])
    if not m:
        continue
    try:
        recs.append(dict(
            desig=r["desig"], cohort=r["cohort"],
            year=int(m.group(1)),
            theta=float(r["theta"]),
            dtau=float(r["dtau_unexplained"]),
            dtau_tot=float(r["dtau_total"]),
            drot=float(r["drot"]),
            q=float(r["q"])))
    except (ValueError, KeyError):
        continue

logger.info(f"comets with parsed perihelion year: n={len(recs)} "
            f"({min(r['year'] for r in recs)}-{max(r['year'] for r in recs)})")

def boot_med(v):
    v = np.asarray(v)
    if len(v) < 4:
        return (float(np.median(v)), np.nan, np.nan)
    idx = rng.integers(0, len(v), (NBOOT, len(v)))
    bs = np.median(v[idx], axis=1)
    lo, hi = np.percentile(bs, [16, 84])
    return (float(np.median(v)), float(lo), float(hi))

def epoch_test(sub):
    th = np.array([r["theta"] for r in sub]); inc = th < CAP
    yr = np.array([r["year"] for r in sub])
    dt = np.array([r["dtau"] for r in sub])
    out = {"n": int(len(sub)), "n_in": int(inc.sum())}
    for tag, m in (("in", inc), ("out", ~inc), ("all", np.ones(len(sub), bool))):
        if m.sum() < 6:
            out[tag] = None; continue
        rho, p = spearmanr(yr[m], dt[m])
        b = np.polyfit(yr[m], dt[m], 1)
        out[tag] = {"rho": float(rho), "p_2sided": float(p),
                    "slope_yr_per_century": float(b[0] * 100.0),
                    "n": int(m.sum())}
    return out

res = {"method": "unexplained slip (step_065 budget-regressed residual) "
                 "vs perihelion year, in/out of the 60-deg cap.  Flat = "
                 "inertially fixed boundary; drift = advecting interface "
                 "or improving-orbit systematic.  Era split at 1950.",
       "cap_deg": CAP, "seed": SEED}

for co in ("code", "warsaw", "pooled"):
    sub = [r for r in recs if co == "pooled" or r["cohort"] == co]
    res[f"epoch_{co}"] = epoch_test(sub)

# era-stratified in-vs-out contrast -- the anomaly is the CONTRAST,
# and any global era-dependent term in the residual cancels in it.
# This is the honest test of whether the signal survives in the
# modern (CCD-era) cohort.
res["era_contrast"] = []
for lo, hi in ((1885, 1950), (1950, 1990), (1990, 2030)):
    ein  = [r["dtau"] for r in recs if r["theta"] < CAP and lo <= r["year"] < hi]
    eout = [r["dtau"] for r in recs if r["theta"] >= CAP and lo <= r["year"] < hi]
    if len(ein) < 4 or len(eout) < 4:
        continue
    U = mannwhitneyu(ein, eout, alternative="greater")
    rec_ = dict(era=f"{lo}-{hi}",
                n_in=len(ein), n_out=len(eout),
                med_in=float(np.median(ein)),
                med_out=float(np.median(eout)),
                contrast=float(np.median(ein) - np.median(eout)),
                p_greater=float(U.pvalue))
    res["era_contrast"].append(rec_)
    logger.info(f"era {lo}-{hi}: in {rec_['med_in']:+.2f} (n={len(ein)}) "
                f"out {rec_['med_out']:+.2f} (n={len(eout)}) "
                f"contrast {rec_['contrast']:+.2f} p={U.pvalue:.4f}")

# absolute-level era split retained for context (both caps drift
# together -> a global era term exists in the residual; the contrast
# above is the controlled statistic)
inc_recs = [r for r in recs if r["theta"] < CAP]
early = [r["dtau"] for r in inc_recs if r["year"] < 1950]
late  = [r["dtau"] for r in inc_recs if r["year"] >= 1950]
if len(early) > 4 and len(late) > 4:
    U = mannwhitneyu(early, late, alternative="two-sided")
    res["era_split_1950_incap"] = {
        "n_early": len(early), "n_late": len(late),
        "med_early": float(np.median(early)),
        "med_late": float(np.median(late)),
        "p_2sided": float(U.pvalue)}
out_recs = [r for r in recs if r["theta"] >= CAP]
early = [r["dtau"] for r in out_recs if r["year"] < 1950]
late  = [r["dtau"] for r in out_recs if r["year"] >= 1950]
if len(early) > 4 and len(late) > 4:
    U = mannwhitneyu(early, late, alternative="two-sided")
    res["era_split_1950_outcap"] = {
        "n_early": len(early), "n_late": len(late),
        "med_early": float(np.median(early)),
        "med_late": float(np.median(late)),
        "p_2sided": float(U.pvalue)}

for co in ("code", "warsaw", "pooled"):
    e = res[f"epoch_{co}"]
    for tag in ("in", "out"):
        if e.get(tag):
            logger.info(f"{co}/{tag}: n={e[tag]['n']} "
                        f"rho={e[tag]['rho']:+.3f} p={e[tag]['p_2sided']:.3f} "
                        f"slope={e[tag]['slope_yr_per_century']:+.2f} yr/cy")

with open(RESULTS / "step_b37_epoch_stability.json", "w") as f:
    json.dump(res, f, indent=1)
print(f"wrote {RESULTS / 'step_b37_epoch_stability.json'}")

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

ax = axes[0]
for co, mk, c_in, c_out in (("code", "o", "crimson", "0.55"),
                            ("warsaw", "s", "darkorange", "0.75")):
    sub = [r for r in recs if r["cohort"] == co]
    yr = np.array([r["year"] for r in sub])
    th = np.array([r["theta"] for r in sub]); inc = th < CAP
    dt = np.array([r["dtau"] for r in sub])
    ax.scatter(yr[~inc], dt[~inc], s=18, marker=mk, facecolors="none",
               edgecolors=c_out, label=f"{co} outside")
    ax.scatter(yr[inc], dt[inc], s=22, marker=mk, c=c_in,
               label=f"{co} inside")
ax.axhline(0, color="k", ls=":", lw=1)
ax.set_xlabel("perihelion year")
ax.set_ylabel(r"unexplained slip $\delta\tau$ (yr)")
ax.legend(frameon=False, fontsize=7, ncol=2)
ax.set_title("slip vs epoch -- fixed wall predicts flat", fontsize=10)

ax = axes[1]
ec = res.get("era_contrast", [])
if ec:
    xs = np.arange(len(ec))
    ax.bar(xs - 0.2, [e["med_in"] for e in ec], width=0.4,
           color="crimson", alpha=0.8, label="inside cap")
    ax.bar(xs + 0.2, [e["med_out"] for e in ec], width=0.4,
           color="0.55", alpha=0.8, label="outside")
    for x, e in zip(xs, ec):
        ax.annotate(f"p={e['p_greater']:.3f}", (x, max(e["med_in"], e["med_out"]) + 0.6),
                    ha="center", fontsize=7)
    ax.axhline(0, color="k", ls=":", lw=1)
    ax.set_xticks(xs)
    ax.set_xticklabels([e["era"] for e in ec])
    ax.set_ylabel(r"median $\delta\tau$ (yr)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("era-stratified in/out contrast", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b37_epoch_stability.png", dpi=150)
print(f"wrote {FIG / 'step_b37_epoch_stability.png'}")
