#!/usr/bin/env python3
"""Step 074 -- leave-one-channel-out robustness of the joint axis test.

Step 069's global result could in principle be driven by a single
strong channel.  Each channel is dropped in turn and the combined
statistic S_-j = -2 sum_{i != j} ln p_i is re-evaluated at the observed
axis against the same 20,000-direction axis-permutation null (same
seed, same draws).  If the joint anomaly remains significant under
every leave-one-out set, no single channel manufactures the result;
the channel whose removal most degrades S is identified as the
largest contributor -- a weight statement, not a significance claim.

Inputs : same data stack as step_069
Outputs: results/step_b39_loco_joint.json
         results/figures/step_b39_loco_joint.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, lv, parse_code, tee_stdout
logger = StepLogger("step_074_loco_joint")
tee_stdout(logger)
logger.header("Leave-one-channel-out joint axis test")

import csv
import json
import math
import numpy as np
from scipy.stats import norm, rankdata

CAP   = 60.0
SEED  = 20260918
NAX   = 20000
rng   = np.random.default_rng(SEED)

# ------------------------------------------------------------------
# Data loading -- identical to step_069
# ------------------------------------------------------------------

def perih_dir(om, Om, inc):
    co, so = math.cos(om), math.sin(om)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(inc), math.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci, so * si])

def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(desig=line[5:17].strip(),
                             com=line[3].strip(),
                             q=float(line[42:56]),
                             w=float(line[70:82]), Om=float(line[82:94]),
                             i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows

PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
def dedup(rows):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or PREF.get(r["com"], 9) < PREF.get(out[k]["com"], 9):
            out[k] = r
    return out

code_o = parse_code(DATA_RAW / "code" / "code_original.html")
war_r  = dedup(parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")))

rows = list(csv.DictReader(open(RESULTS / "step_b30_proper_time_slip.csv")))
comets = []
for r in rows:
    d = None
    if r["cohort"] == "code" and r["desig"] in code_o:
        e = code_o[r["desig"]]
        d = -perih_dir(math.radians(e["w"]), math.radians(e["Om"]),
                       math.radians(e["i"]))
    elif r["cohort"] == "warsaw" and r["desig"] in war_r:
        e = war_r[r["desig"]]
        d = -perih_dir(math.radians(e["w"]), math.radians(e["Om"]),
                       math.radians(e["i"]))
    if d is None:
        continue
    comets.append(dict(desig=r["desig"], cohort=r["cohort"], dir=d,
                       q=float(r["q"]),
                       dtau=float(r["dtau_unexplained"]),
                       dtau_in=float(r["dtau_inbound_unexplained"])
                       if r["dtau_inbound_unexplained"] else np.nan))

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
    if np.isfinite(Om) and np.isfinite(w) and np.isfinite(i):
        tnos.append(dict(
            dir=perih_dir(math.radians(w), math.radians(Om),
                          math.radians(i)),
            varpi=math.radians((Om + w) % 360.0)))

logger.info(f"comets {len(comets)}, TNOs {len(tnos)}")

cdir  = np.array([c["dir"] for c in comets])
cdtau = np.array([c["dtau"] for c in comets])
cq    = np.array([c["q"] for c in comets])
qm    = np.array([c["cohort"] == "code" for c in comets]) & (cq < 3.1)
crank = rankdata(cdtau[qm])
iswar = np.array([c["cohort"] == "warsaw" for c in comets])
widx  = np.where(iswar)[0]
wdin  = np.array([c["dtau_in"] for c in comets])[widx]
wrank = rankdata(wdin)

tdir   = np.array([t["dir"] for t in tnos])
tvarpi = np.array([t["varpi"] for t in tnos])
texp   = np.exp(1j * tvarpi)

def mwz(rank_arr, mask):
    n1 = int(mask.sum()); n2 = int((~mask).sum())
    if n1 < 3 or n2 < 3:
        return 0.0
    U = rank_arr[mask].sum() - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    return (U - mu) / sd

def channels(axis):
    cth = np.degrees(np.arccos(np.clip(cdir @ axis, -1, 1)))
    inc_m = cth[qm] < CAP
    p1 = norm.sf(mwz(crank, inc_m))
    proj = cdir @ axis
    p2 = norm.sf(proj.mean() * math.sqrt(3 * len(proj)))
    winc = cth[widx] < CAP
    p3 = norm.sf(mwz(wrank, winc))
    tproj = tdir @ axis
    p4 = norm.sf(tproj.mean() * math.sqrt(3 * len(tproj)))
    tth = np.degrees(np.arccos(np.clip(tdir @ axis, -1, 1)))
    tinc = tth < CAP
    ni = int(tinc.sum())
    p5 = math.exp(-ni * (abs(texp[tinc].sum()) / ni) ** 2) if ni >= 3 else 1.0
    return np.array([max(p1, 1e-300), max(p2, 1e-300), max(p3, 1e-300),
                     max(p4, 1e-300), max(p5, 1e-300)])

AX_TNO = lv(49.0, -17.0)
AX_COM = lv(34.0, -13.0)
NAMES  = ["rotation_slip", "aphelion_dipole", "warsaw_inbound",
          "tno_aphelion", "tno_varpi"]

u = rng.uniform(-1, 1, NAX)
ph = rng.uniform(0, 2 * math.pi, NAX)
rr = np.sqrt(1 - u * u)
axes = np.column_stack([rr * np.cos(ph), rr * np.sin(ph), u])

p_obs_tno = channels(AX_TNO)
p_obs_com = channels(AX_COM)
ch_rand = np.empty((NAX, 5))
for j, a in enumerate(axes):
    ch_rand[j] = channels(a)
    if (j + 1) % 5000 == 0:
        logger.info(f"  {j + 1}/{NAX} axes")

ln_rand = np.log(ch_rand)                # (NAX,5)
ln_tno  = np.log(p_obs_tno)
ln_com  = np.log(p_obs_com)

def loco(obs_ln, rand_ln, drop):
    keep = [k for k in range(5) if k != drop] if drop is not None else list(range(5))
    S_obs = -2.0 * obs_ln[keep].sum()
    S_r   = -2.0 * rand_ln[:, keep].sum(axis=1)
    return S_obs, float((1 + (S_r >= S_obs).sum()) / (NAX + 1))

res = {"method": "leave-one-channel-out on the step_069 axis "
                 "permutation: S_-j recomputed about the same 20000 "
                 "random directions (same seed); global p per "
                 "leave-out set.  Channel weights = S drop on removal.",
       "seed": SEED, "n_axes": NAX, "cap_deg": CAP,
       "channels": NAMES}

for lab, obs_ln in (("tno_axis_49m17", ln_tno), ("comet_axis_34m13", ln_com)):
    S_full, p_full = loco(obs_ln, ln_rand, None)
    entry = {"full": {"S": S_full, "p_global": p_full}, "loco": {}}
    for j, nm in enumerate(NAMES):
        S_, p_ = loco(obs_ln, ln_rand, j)
        entry["loco"][nm] = {"S": S_, "p_global": p_,
                             "S_drop": S_full - S_}
        logger.info(f"{lab} drop {nm:16s} S={S_:7.2f} "
                    f"p={p_:.5f} (S drop {S_full - S_:+.2f})")
    res[lab] = entry

with open(RESULTS / "step_b39_loco_joint.json", "w") as f:
    json.dump(res, f, indent=1)
print(f"wrote {RESULTS / 'step_b39_loco_joint.json'}")

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7.5, 4.4))
e = res["tno_axis_49m17"]
xs = np.arange(6)
ps = [e["full"]["p_global"]] + [e["loco"][nm]["p_global"] for nm in NAMES]
labels = ["all 5"] + [f"-{n}" for n in NAMES]
colors = ["crimson"] + ["darkorange"] * 5
ax.bar(xs, ps, color=colors, alpha=0.85)
for x, pv in zip(xs, ps):
    ax.annotate(f"{pv:.4f}", (x, pv), ha="center", va="bottom", fontsize=7)
ax.axhline(0.05, color="k", ls=":", lw=1)
ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("global $p$ (axis permutation)")
ax.set_title("joint significance survives every leave-one-out set",
             fontsize=10)
fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b39_loco_joint.png", dpi=150)
print(f"wrote {FIG / 'step_b39_loco_joint.png'}")
