#!/usr/bin/env python3
"""Step 092 -- MPCORB-lineage replication of the detached-TNO
perihelion clustering
====================================================================

The flagship resident result (step 010) is built on JPL SBDB
osculating elements: N = 44 detached objects (a > 150, q > 30,
condition code <= 3) give R_varpi = 0.332, p = 0.0073 under the
conditioned clock-sector null, mean perihelion direction
(lam, beta) ~ (50 deg, -17 deg).  Every downstream claim inherits
that single orbit-fit lineage.  This step repeats the measurement on
the MPCORB catalogue -- the Minor Planet Center's own orbit fits,
computed with different weighting, debiasing and (for many objects)
different solution software -- so a catalogue-lineage artifact in
SBDB would reveal itself as a different axis or a collapsed R.

Sample definitions and statistics are IDENTICAL to step 010:
  detached:    a > 150 AU, q > 30 AU
  extreme:     a > 250 AU, q > 30 AU
  sednoid:     a > 150 AU, q > 50 AU
  detached_q40: a > 150 AU, q > 40 AU
Condition-code quality cut is replaced by the MPCORB-native proxy:
multi-opposition (nopp >= 2) or arc >= 1 yr, with a numbered-only
variant (the most secure subset).  The conditioned null is the same:
fix (a, e, i); draw om, w uniform on the circle; 20000 realizations.

The step additionally performs an object-level cross-match against
the SBDB detached table persisted in step_02_clustering.json:
membership overlap, per-object varpi differences between the two
catalogues, and whether objects unique to either lineage carry the
signal.

Outputs
-------
results/step_b57_mpcorb_clustering.json
results/figures/supplementary/step_b57_mpcorb_clustering.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_092_mpcorb_clustering")
tee_stdout(logger)
logger.header("MPCORB-lineage replication of detached-TNO clustering")

import gzip
import json
import math
import re
import numpy as np

rng = np.random.default_rng(20260918)
N_MC = 20000

# ------------------------------------------------------------------
# MPCORB parser (as step 090)
# ------------------------------------------------------------------

MPCORB = DATA_RAW / "mpc" / "MPCORB.DAT.gz"
MONTH_C = "123456789ABC"
DAY_C = "123456789ABCDEFGHIJKLMNOPQRSTUV"

def load_mpcorb():
    objs = []
    with gzip.open(MPCORB, "rt", errors="replace") as f:
        for i, line in enumerate(f):
            if i < 43 or len(line) < 200:
                continue
            try:
                o = dict(
                    desig=line[0:7].strip(),
                    epoch=line[20:25].strip(),
                    M=float(line[26:35]), w=float(line[37:46]),
                    Om=float(line[48:57]), i=float(line[59:68]),
                    e=float(line[70:79]), a=float(line[92:103]),
                    nopp=int(line[123:126]) if line[123:126].strip() else 0,
                    arc=line[128:136].strip(),
                    name=line[166:194].strip())
            except ValueError:
                continue
            if not o["epoch"].startswith("K"):
                continue
            o["q"] = o["a"] * (1.0 - o["e"])
            o["numbered"] = bool(re.match(r"^\(\d+\)", o["name"]))
            objs.append(o)
    return objs

def prov_desig(name):
    return name.split(")")[-1].strip()

def arc_years(o):
    a = o["arc"]
    if "-" in a:
        try:
            y0, y1 = a.split("-")
            return float(y1) - float(y0)
        except ValueError:
            return 0.0
    try:
        return float(a) / 365.25
    except ValueError:
        return 0.0

# ------------------------------------------------------------------
# Step-010 statistics, verbatim
# ------------------------------------------------------------------

def perih_vec(om, w, i):
    return np.stack([
        np.cos(om) * np.cos(w) - np.sin(om) * np.sin(w) * np.cos(i),
        np.sin(om) * np.cos(w) + np.cos(om) * np.sin(w) * np.cos(i),
        np.sin(w) * np.sin(i)], axis=-1)

def pole_vec(om, i):
    n = np.stack([np.sin(i) * np.sin(om),
                  -np.sin(i) * np.cos(om),
                  np.cos(i)], axis=-1)
    n[n[..., 2] < 0] *= -1
    return n

def stats_of(om, w, i):
    varpi = (om + w) % (2 * np.pi)
    R_varpi = abs(np.exp(1j * varpi).mean())
    R_phat = np.linalg.norm(perih_vec(om, w, i).mean(axis=0))
    R_nhat = np.linalg.norm(pole_vec(om, i).mean(axis=0))
    return R_varpi, R_phat, R_nhat

def conditioned_mc(objs, stat_obs):
    i_rad = np.deg2rad(np.array([o["i"] for o in objs]))
    n = len(objs)
    om_mc = rng.uniform(0, 2 * np.pi, (N_MC, n))
    w_mc = rng.uniform(0, 2 * np.pi, (N_MC, n))
    i_b = np.broadcast_to(i_rad, (N_MC, n))
    varpi = (om_mc + w_mc) % (2 * np.pi)
    Rv = abs(np.exp(1j * varpi).mean(axis=1))
    Rp = np.linalg.norm(perih_vec(om_mc, w_mc, i_b).mean(axis=1), axis=1)
    Rn = np.linalg.norm(pole_vec(om_mc, i_b).mean(axis=1), axis=1)
    return [(float((int((Rv >= stat_obs[0]).sum()) + 1)
                  / (N_MC + 1)), float(Rv.mean())),
            (float((int((Rp >= stat_obs[1]).sum()) + 1)
                  / (N_MC + 1)), float(Rp.mean())),
            (float((int((Rn >= stat_obs[2]).sum()) + 1)
                  / (N_MC + 1)), float(Rn.mean()))]

def boot_ci(ang, n_boot=20000):
    th = np.asarray(ang)
    n = th.size
    idx = rng.integers(0, n, (n_boot, n))
    m = np.angle(np.exp(1j * th[idx]).mean(axis=1))
    c = np.angle(np.exp(1j * th).mean())
    d = np.rad2deg((m - c + np.pi) % (2 * np.pi) - np.pi)
    lo, hi = np.percentile(d, [16, 84])
    return float(lo), float(hi)

def sample_report(objs, label):
    om = np.deg2rad(np.array([o["om"] for o in objs]))
    w = np.deg2rad(np.array([o["w"] for o in objs]))
    i = np.deg2rad(np.array([o["i"] for o in objs]))
    varpi = (om + w) % (2 * np.pi)
    s_obs = stats_of(om, w, i)
    (pv, mv), (pp, mp), (pn, mn) = conditioned_mc(objs, s_obs)
    mean_v = float(np.rad2deg(np.angle(np.exp(1j * varpi).mean())) % 360)
    lo, hi = boot_ci(varpi)
    ph = perih_vec(om, w, i)
    mvec = ph.mean(axis=0); mvec /= np.linalg.norm(mvec)
    lam = float(np.rad2deg(np.arctan2(mvec[1], mvec[0])) % 360)
    bet = float(np.rad2deg(np.arcsin(mvec[2])))
    rep = {
        "label": label, "N": len(objs),
        "varpi_deg": sorted(round(float(np.rad2deg(x)), 2)
                            for x in varpi),
        "varpi_mean_deg": round(mean_v, 1),
        "varpi_mean_ci68": [round(mean_v + lo, 1),
                            round(mean_v + hi, 1)],
        "R_varpi": round(s_obs[0], 4),
        "R_varpi_null_mean": round(mv, 4),
        "p_varpi_conditioned": pv,
        "R_phat": round(s_obs[1], 4),
        "R_phat_null_mean": round(mp, 4),
        "p_phat_conditioned": pp,
        "phat_ecliptic_lam_deg": round(lam, 1),
        "phat_ecliptic_beta_deg": round(bet, 1),
    }
    logger.info(f"{label}: N={len(objs)}  R_varpi={s_obs[0]:.3f} "
                f"(null {mv:.3f}) p={pv:.2e}  "
                f"phat->({lam:.0f},{bet:.0f})")
    return rep

# ------------------------------------------------------------------
# Load MPCORB, build samples
# ------------------------------------------------------------------

logger.info("loading MPCORB ...")
objs_all = load_mpcorb()
logger.info(f"  {len(objs_all)} orbits")

for o in objs_all:
    o["om"] = o["Om"]      # step-010 field names

def secure(o):
    return o["nopp"] >= 2 or arc_years(o) >= 1.0

def sel(a_min, q_min, secure_only=True, numbered=False):
    return [o for o in objs_all
            if o["a"] > a_min and o["q"] > q_min and o["e"] < 1.0
            and (not secure_only or secure(o))
            and (not numbered or o["numbered"])]

samples = {
    "detached_a150_q30": sel(150, 30),
    "extreme_a250_q30": sel(250, 30),
    "sednoid_a150_q50": sel(150, 50),
    "detached_a150_q40": sel(150, 40),
    "detached_a150_q30_numbered": sel(150, 30, numbered=True),
    "detached_a150_q30_allq": sel(150, 30, secure_only=False),
}

out = {"meta": dict(
    catalog="MPCORB.DAT (Minor Planet Center orbit fits; "
            "independent lineage of JPL SBDB)",
    quality="nopp>=2 or arc>=1 yr (condition-code proxy); "
            "numbered-only variant shown",
    null_model="fix (a,e,i); om,w uniform; 20000 realizations "
               "(identical to step 010)",
    reference="step_010_tno_clustering.py / "
              "results/step_02_clustering.json"),
    "n_catalog": len(objs_all), "samples": {}}

for label, s in samples.items():
    if len(s) < 5:
        logger.warning(f"{label}: N={len(s)} -- skipped")
        continue
    out["samples"][label] = sample_report(s, label)

# ------------------------------------------------------------------
# Cross-match against the SBDB detached table
# ------------------------------------------------------------------

logger.subheader("SBDB cross-match")
sbdb_tab = json.load(open(RESULTS / "step_02_clustering.json"))[
    "detached_table"]

def norm_desig(name):
    """'765133 (2013 SL102)' or '(90377) Sedna' -> compact key.
    The provisional designation sits inside parentheses in SBDB
    full_name, after ')' in MPCORB readable names."""
    m = re.search(r"\(([^)]+)\)", name)
    key = m.group(1) if m else name.split(")")[-1]
    return re.sub(r"\s+", "", key.strip())

mpc_by_pd = {}
mpc_by_num = {}
for o in objs_all:
    pd = prov_desig(o["name"])
    if pd:
        mpc_by_pd[re.sub(r"\s+", "", pd)] = o
    mn = re.match(r"^\((\d+)\)", o["name"])
    if mn:
        mpc_by_num[mn.group(1)] = o

def find_mpc(sbdb_name):
    key = norm_desig(sbdb_name)
    if key in mpc_by_pd:
        return mpc_by_pd[key]
    mn = re.match(r"^(\d+)\s", sbdb_name)
    return mpc_by_num.get(mn.group(1)) if mn else None

xrows = []
n_both = 0
dv = []
for t in sbdb_tab:
    key = norm_desig(t["name"])
    mo = find_mpc(t["name"])
    v_mpc = None
    if mo is not None:
        v_mpc = round((mo["om"] + mo["w"]) % 360, 1)
        n_both += 1
        dv.append(abs((v_mpc - t["varpi"] + 180) % 360 - 180))
    xrows.append(dict(sbdb_name=t["name"], sbdb_varpi=t["varpi"],
                      mpcorb_varpi=v_mpc,
                      in_mpcorb_detached=mo in samples["detached_a150_q30"]
                      if mo is not None else False))

sbdb_keys = {norm_desig(t["name"]) for t in sbdb_tab}
sbdb_nums = {m.group(1) for t in sbdb_tab
             if (m := re.match(r"^(\d+)\s", t["name"]))}
def mpc_in_sbdb(o):
    if re.sub(r"\s+", "", prov_desig(o["name"])) in sbdb_keys:
        return True
    mn = re.match(r"^\((\d+)\)", o["name"])
    return bool(mn and mn.group(1) in sbdb_nums)

mpc_only = [prov_desig(o["name"]) for o in samples["detached_a150_q30"]
            if not mpc_in_sbdb(o)]

out["cross_match"] = dict(
    n_sbdb_detached=len(sbdb_tab),
    n_matched_in_mpcorb=n_both,
    median_abs_dvarpi_deg=round(float(np.median(dv)), 2) if dv else None,
    max_abs_dvarpi_deg=round(float(np.max(dv)), 2) if dv else None,
    n_mpcorb_only=len(mpc_only), mpcorb_only=mpc_only,
    rows=xrows)
logger.info(f"SBDB detached N={len(sbdb_tab)}; matched in MPCORB "
            f"{n_both}; median |dvarpi| "
            f"{out['cross_match']['median_abs_dvarpi_deg']} deg; "
            f"MPCORB-only members {len(mpc_only)}: {mpc_only}")

# subsample restricted to cross-matched objects
both = [mo for t in sbdb_tab
        if (mo := find_mpc(t["name"])) is not None
        and mo["a"] > 150 and mo["q"] > 30]
if len(both) >= 5:
    out["samples"]["sbdb_matched_subset"] = sample_report(
        both, "SBDB-matched subset (MPCORB elements)")

# ------------------------------------------------------------------
# Write + figure
# ------------------------------------------------------------------

RES = RESULTS
json.dump(out, open(RES / "step_b57_mpcorb_clustering.json", "w"),
          indent=1, default=float)
logger.data_save(RES / "step_b57_mpcorb_clustering.json")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(13, 4.4))
ax1 = fig.add_subplot(131, projection="polar")
rep = out["samples"]["detached_a150_q30"]
v = np.deg2rad(np.array(rep["varpi_deg"]))
ax1.hist(v, bins=18, range=(0, 2 * np.pi), color="teal", alpha=0.85,
         edgecolor="k", lw=0.4)
mv = np.deg2rad(rep["varpi_mean_deg"])
ax1.plot([mv, mv], [0, ax1.get_ylim()[1]], "r-", lw=2)
ax1.set_theta_zero_location("N")
ax1.set_title(f"MPCORB detached sample (N={rep['N']})\n"
              f"R={rep['R_varpi']}, p={rep['p_varpi_conditioned']:.1e}",
              fontsize=9)

ax2 = fig.add_subplot(132)
sb = np.array([t["varpi"] for t in sbdb_tab])
mp = np.array([x["mpcorb_varpi"] if x["mpcorb_varpi"] is not None
               else np.nan for x in xrows])
ax2.plot([0, 360], [0, 360], "k--", lw=1)
ax2.plot(sb, mp, "o", ms=5, color="darkgreen")
ax2.set_xlabel("$\\varpi_{SBDB}$ (deg)")
ax2.set_ylabel("$\\varpi_{MPCORB}$ (deg)")
ax2.set_title("per-object longitude of perihelion:\n"
              "SBDB fit vs MPC fit", fontsize=9)

ax3 = fig.add_subplot(133)
for lb, c in [("detached_a150_q30", "teal"),
              ("extreme_a250_q30", "navy"),
              ("detached_a150_q30_numbered", "crimson"),
              ("detached_a150_q30_allq", "0.6")]:
    if lb not in out["samples"]:
        continue
    vv = np.array(out["samples"][lb]["varpi_deg"])
    ax3.hist(vv, bins=36, range=(0, 360), histtype="step", color=c,
             lw=1.4, label=f"{lb} (N={out['samples'][lb]['N']})")
ax3.axvline(rep["varpi_mean_deg"], color="r", lw=2)
ax3.set_xlabel("$\\varpi$ (deg)"); ax3.legend(fontsize=6.5)
ax3.set_title("MPCORB $\\varpi$ histograms by sample", fontsize=9)

fig.tight_layout()
FIG = RES / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b57_mpcorb_clustering.png", dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b57_mpcorb_clustering.png')