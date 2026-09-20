#!/usr/bin/env python3
"""Step 078 -- discovery-era stability of the detached-TNO cluster.

If the resident alignment were an artefact of early wide-field
surveys (photographic-era pointing or orbit quality), it should fade
in the modern CCD/LSST-era cohort.  The detached sample
(a > 150, q > 30 AU, condition_code <= 3) is split by first-observation
year; cap membership and the varpi resultant are measured per era,
with the footprint-free in-cap rate and Rayleigh significance.

Inputs : data/raw/sbdb/sbdb_outer_ss.json
Outputs: results/step_b43_tno_era.json
         results/figures/supplementary/step_b43_tno_era.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_078_tno_era")
tee_stdout(logger)
logger.header("Discovery-era stability of the detached-TNO cluster")

import json
import math
import numpy as np
from scipy.stats import binomtest

SEED = 20260918
rng  = np.random.default_rng(SEED)
CAP  = 60.0

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

AX = lv(49.0, -17.0)

def fnum(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")

_d = json.load(open(DATA_RAW / "sbdb" / "sbdb_outer_ss.json"))
sbdb = [dict(zip(_d["fields"], rec)) for rec in _d["data"]]
tnos = []
for r in sbdb:
    a, q, cc = fnum(r, "a"), fnum(r, "q"), fnum(r, "condition_code")
    if not (a > 150 and q > 30 and 0 <= cc <= 3):
        continue
    Om, w, i = fnum(r, "om"), fnum(r, "w"), fnum(r, "i")
    fo = r.get("first_obs") or ""
    try:
        yr = int(str(fo)[:4])
    except ValueError:
        continue
    if not (np.isfinite(Om) and np.isfinite(w) and np.isfinite(i)):
        continue
    tnos.append(dict(yr=yr,
        dir=perih_dir(math.radians(w), math.radians(Om), math.radians(i)),
        varpi=math.radians((Om + w) % 360.0)))

logger.info(f"detached TNOs with first_obs: {len(tnos)}")

ERAS = [(1980, 2005), (2005, 2015), (2015, 2030)]
res = {"method": "detached-TNO sample split by first-observation year; "
                 "per-era cap membership, periapsis resultant toward "
                 "the resident axis, and varpi Rayleigh.  An early-"
                 "survey artefact fades with era; an intrinsic "
                 "alignment persists in the modern cohort.",
       "cap_deg": CAP, "axis": "49,-17", "seed": SEED,
       "eras": [f"{a}-{b}" for a, b in ERAS], "n": len(tnos)}

dirs   = np.array([t["dir"] for t in tnos])
varpi  = np.array([t["varpi"] for t in tnos])
th     = np.degrees(np.arccos(np.clip(dirs @ AX, -1, 1)))
inc    = th < CAP
frac_in = inc.mean()

for lo, hi in ERAS:
    m = np.array([lo <= t["yr"] < hi for t in tnos])
    n = int(m.sum())
    if n < 4:
        res[f"{lo}-{hi}"] = None; continue
    nin = int((inc & m).sum())
    # periapsis-direction projection toward axis within the era
    proj = (dirs[m] @ AX).mean() * math.sqrt(3 * n)
    # varpi Rayleigh within era (full, not cap-restricted)
    R = abs(np.exp(1j * varpi[m]).sum()) / n
    p_ray = math.exp(-n * R * R)
    bt = binomtest(nin, n, frac_in)
    bt_g = binomtest(nin, n, frac_in, alternative="greater")
    res[f"{lo}-{hi}"] = {
        "n": n, "n_in_cap": nin,
        "cap_frac": nin / n, "cap_frac_global": float(frac_in),
        # declared directional prediction (in-cap excess) -> greater
        "p_cap_binom": float(bt_g.pvalue),
        "p_cap_binom_2sided": float(bt.pvalue),
        "axis_proj_z": float(proj),
        "varpi_R": float(R), "varpi_p": float(p_ray)}
    logger.info(f"{lo}-{hi}: n={n} in-cap={nin} ({nin/n:.2f}) "
                f"R={R:.2f} p_ray={p_ray:.4f} binom_p={bt.pvalue:.4f}")

# is first_obs era itself correlated with cap membership?
from scipy.stats import spearmanr
yrs = np.array([t["yr"] for t in tnos])
rho, p = spearmanr(yrs, inc.astype(float))
res["cap_vs_year"] = {"rho": float(rho), "p_2sided": float(p)}
logger.info(f"cap membership vs first_obs year: rho={rho:+.3f} p={p:.3f}")

with open(RESULTS / "step_b43_tno_era.json", "w") as f:
    json.dump(res, f, indent=1)
logger.data_save(RESULTS / 'step_b43_tno_era.json')
# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7.2, 4.4))
labels, vals, ns = [], [], []
for lo, hi in ERAS:
    e = res.get(f"{lo}-{hi}")
    if e:
        labels.append(f"{lo}-{hi}\n(n={e['n']})")
        vals.append(e["cap_frac"]); ns.append(e["n_in_cap"])
ax.bar(range(len(vals)), vals, color="crimson", alpha=0.8, width=0.5)
ax.axhline(frac_in, color="k", ls=":", lw=1,
           label=f"full-sample rate {frac_in:.2f}")
ax.set_xticks(range(len(vals))); ax.set_xticklabels(labels)
ax.set_ylabel("fraction inside 60° cap")
ax.legend(frameon=False, fontsize=8)
ax.set_title("resident alignment by discovery era", fontsize=10)
fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b43_tno_era.png", dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b43_tno_era.png')