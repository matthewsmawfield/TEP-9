"""step_069: Joint axis-permutation global significance.

Every channel in this paper is conditioned on a direction on the sky:
the TNO cluster, the comet rotation residual, the implied proper-time
slip, the Warsaw inbound leg, and the aphelion dipole are each
significant about their pre-declared axes.  What has not been computed
is the probability that a single arbitrary direction reproduces the
JOINT anomaly across all populations and catalogues at once.

This step permutes the axis -- the one variable common to every
channel.  For each of N random sky directions the full battery is
recomputed about that direction:

  C1  rotation channel : Mann-Whitney contrast of the CODE matched
      unexplained slip (dtau residual, in-cap vs out-of-cap)
  C2  dipole channel   : mean projection of pooled comet aphelia onto
      the trial direction
  C3  inbound channel  : Mann-Whitney contrast of the Warsaw
      inbound-leg unexplained slip
  C4  resident channel : mean projection of detached-TNO aphelia onto
      the trial direction
  C5  cluster channel  : Rayleigh test of the detached-TNO varpi
      directions inside the trial cap

The five channel p-values are combined by Fisher's statistic
S = -2 sum ln p.  Because every draw uses the same data and changes
only the axis, the null distribution of S is exact for the combined
statistic and requires no independence assumption between channels.
The global p is the fraction of random directions whose combined
statistic reaches the value measured at the observed axis.

Evaluated at both pre-declared axes: the TNO cluster axis
(49,-17) and the Warsaw-discovery comet axis (34,-13).

Inputs
------
results/step_b30_proper_time_slip.csv      (per-comet slip, step_065)
data/raw/code/code_original.html           (CODE aphelion dirs)
data/raw/warsaw/warsaw_tablec.dat          (Warsaw aphelion dirs)
data/raw/sbdb/sbdb_outer_ss.json           (detached TNOs)

Outputs
-------
results/step_b34_joint_axis.json
results/figures/supplementary/step_b34_joint_axis.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, lv, parse_code, tee_stdout
logger = StepLogger("step_069_joint_axis_test")
tee_stdout(logger)
logger.header("Joint axis-permutation global significance")

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
# Per-comet data: aphelion direction + slip channels
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
missed = 0
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
        missed += 1
        continue
    comets.append(dict(desig=r["desig"], cohort=r["cohort"], dir=d,
                       q=float(r["q"]),
                       dtau=float(r["dtau_unexplained"]),
                       dtau_in=float(r["dtau_inbound_unexplained"])
                       if r["dtau_inbound_unexplained"] else np.nan))
logger.info(f"comets with aphelion dirs: {len(comets)} (missed {missed})")

# ------------------------------------------------------------------
# Detached TNOs: aphelion direction + varpi
# ------------------------------------------------------------------

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
                          math.radians(i)),   # periapsis direction
            varpi=math.radians((Om + w) % 360.0)))
logger.info(f"detached TNOs: {len(tnos)}")

# ------------------------------------------------------------------
# Channel machinery -- rank-based for speed
# ------------------------------------------------------------------

cdir  = np.array([c["dir"] for c in comets])
cdtau = np.array([c["dtau"] for c in comets])
cq    = np.array([c["q"] for c in comets])
# CODE matched (deep-plunger) subset -- the pre-declared rotation channel
qm    = np.array([c["cohort"] == "code" for c in comets]) & (cq < 3.1)
crank = rankdata(cdtau[qm])                   # ranks within the channel set
iswar = np.array([c["cohort"] == "warsaw" for c in comets])
widx  = np.where(iswar)[0]
wdin  = np.array([c["dtau_in"] for c in comets])[widx]
wrank = rankdata(wdin)

tdir   = np.array([t["dir"] for t in tnos])
tvarpi = np.array([t["varpi"] for t in tnos])
texp   = np.exp(1j * tvarpi)                  # unit phasors

def mwz(rank_arr, mask):
    """one-sided Mann-Whitney z (in > out) from precomputed ranks."""
    n1 = int(mask.sum()); n2 = int((~mask).sum())
    if n1 < 3 or n2 < 3:
        return 0.0
    U = rank_arr[mask].sum() - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    return (U - mu) / sd

def channels(axis):
    """the five channel p-values about a trial axis."""
    cth = np.degrees(np.arccos(np.clip(cdir @ axis, -1, 1)))
    # C1 CODE matched slip contrast
    inc_m = cth[qm] < CAP
    z1 = mwz(crank, inc_m)
    p1 = norm.sf(z1)
    # C2 dipole mean projection (pooled comets)
    proj = cdir @ axis
    z2 = proj.mean() * math.sqrt(3 * len(proj))
    p2 = norm.sf(z2)
    # C3 Warsaw inbound slip
    winc = cth[widx] < CAP
    z3 = mwz(wrank, winc)
    p3 = norm.sf(z3)
    # C4 TNO aphelion projection
    tproj = tdir @ axis
    z4 = tproj.mean() * math.sqrt(3 * len(tproj))
    p4 = norm.sf(z4)
    # C5 TNO in-cap varpi resultant
    tth = np.degrees(np.arccos(np.clip(tdir @ axis, -1, 1)))
    tinc = tth < CAP
    ni = int(tinc.sum())
    if ni >= 3:
        R = abs(texp[tinc].sum()) / ni
        p5 = math.exp(-ni * R * R)
    else:
        p5 = 1.0
    return np.array([max(p1, 1e-300), max(p2, 1e-300), max(p3, 1e-300),
                     max(p4, 1e-300), max(p5, 1e-300)])

def S(axis):
    return float(-2.0 * np.log(channels(axis)).sum())

AX_TNO = lv(49.0, -17.0)
AX_COM = lv(34.0, -13.0)
S_tno = S(AX_TNO)
S_com = S(AX_COM)
p_tno_channels = channels(AX_TNO)
p_com_channels = channels(AX_COM)
logger.info(f"S(TNO axis 49,-17)   = {S_tno:.2f}  channels={p_tno_channels}")
logger.info(f"S(comet axis 34,-13) = {S_com:.2f}  channels={p_com_channels}")

# ------------------------------------------------------------------
# Random-axis null
# ------------------------------------------------------------------

u = rng.uniform(-1, 1, NAX)
ph = rng.uniform(0, 2 * math.pi, NAX)
rr = np.sqrt(1 - u * u)
axes = np.column_stack([rr * np.cos(ph), rr * np.sin(ph), u])

S_rand = np.empty(NAX)
ch_rand = np.empty((NAX, 5))
for j, a in enumerate(axes):
    p = channels(a)
    ch_rand[j] = p
    S_rand[j] = -2.0 * np.log(p).sum()
    if (j + 1) % 5000 == 0:
        logger.info(f"  {j + 1}/{NAX} axes")

p_global_tno = float((1 + (S_rand >= S_tno).sum()) / (NAX + 1))
p_global_com = float((1 + (S_rand >= S_com).sum()) / (NAX + 1))

# per-channel look-elsewhere at the TNO axis
per_ch = {}
for j, name in enumerate(["rotation_slip", "aphelion_dipole",
                          "warsaw_inbound", "tno_aphelion", "tno_varpi"]):
    obs = p_tno_channels[j]
    per_ch[name] = {
        "p_at_tno_axis": float(obs),
        "frac_axes_better": float(
            (int((ch_rand[:, j] <= obs).sum()) + 1) / (NAX + 1))}

# fraction of random axes whose best-fit cap direction lands near the
# observed axis (sanity: how localized is the joint anomaly)
best = np.argmax(S_rand)
sep_best = math.degrees(math.acos(np.clip(
    axes[best] @ AX_TNO, -1, 1)))

res = {
    "method": "axis permutation: all five channels recomputed about "
              "each of N random sky directions; Fisher S = -2 sum ln p; "
              "global p = fraction of random axes matching the observed "
              "combined statistic.  Null preserves cross-channel "
              "dependence (same data, only the axis varies).",
    "cap_deg": CAP, "seed": SEED, "n_axes": NAX,
    "n_comets": len(comets), "n_warsaw": int(iswar.sum()),
    "n_tno": len(tnos),
    "eval_at_tno_axis_49m17": {
        "S": S_tno, "channel_p": p_tno_channels.tolist(),
        "p_global": p_global_tno},
    "eval_at_comet_axis_34m13": {
        "S": S_com, "channel_p": p_com_channels.tolist(),
        "p_global": p_global_com},
    "per_channel_look_elsewhere": per_ch,
    "best_random_axis": {"S": float(S_rand[best]),
                         "sep_from_tno_deg": sep_best},
    "S_rand_quantiles": {
        "q50": float(np.percentile(S_rand, 50)),
        "q95": float(np.percentile(S_rand, 95)),
        "q99": float(np.percentile(S_rand, 99)),
        "q999": float(np.percentile(S_rand, 99.9))},
}

out = str(RESULTS / "step_b34_joint_axis.json")
json.dump(res, open(out, "w"), indent=1)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7.5, 4.6))
ax.hist(S_rand, bins=60, color="0.55", alpha=0.8, density=True,
        label="random-axis null")
ax.axvline(S_tno, color="crimson", lw=1.5,
           label=f"TNO axis (49°,-17°): S={S_tno:.1f}, p={p_global_tno:.4f}")
ax.axvline(S_com, color="darkorange", lw=1.5, ls="--",
           label=f"comet axis (34°,-13°): S={S_com:.1f}, p={p_global_com:.4f}")
ax.set_xlabel(r"combined statistic $S = -2\sum\ln p_j$ (5 channels)")
ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Joint axis-permutation null vs observed axes", fontsize=10)
fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b34_joint_axis.png", dpi=300)

logger.info(f"p_global at TNO axis   = {p_global_tno:.5f}")
logger.info(f"p_global at comet axis = {p_global_com:.5f}")
for k, v in per_ch.items():
    logger.info(f"  {k}: p={v['p_at_tno_axis']:.4f} "
                f"frac_axes_better={v['frac_axes_better']:.4f}")
logger.data_save(out)
logger.data_save(FIG / "supplementary" / "step_b34_joint_axis.png")