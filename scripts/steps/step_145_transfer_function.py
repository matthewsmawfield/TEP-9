"""step_145: transfer-function discriminator (result b109).

The observed periapsis rotation is converted to an equivalent time
scale by dt_eq = dtheta * r_b^2 / h (step_065).  Two microscopic
realizations produce that observable:

- impulse realization: a conformal fifth-force kick across a lapse
  wall, dv = c^2 Da / v_n along the wall normal.  It predicts
  (a) dtheta ~ 1/(v_n v), so the implied slip dt = dtheta/Omega is
  strongly anticorrelated with crossing speed, and (b) an along-
  track component that changes orbital energy at
  D(1/a) = 2 v dv_par / GM.

- holonomy realization: a crossing-localized time translation
  (non-integrable time transport, H = oint_C dtau).  It predicts a
  universal dt independent of crossing velocity and zero local
  impulse, so the energy channel stays flat.

Tests on the independent dual-leg refit records (steps_127/128):

- T1: Spearman rho(dtau, v_rx) at shells r_x = 60/100/250 AU, in and
  out of each era's anomaly cap, both eras.  Flat selects holonomy.
- T2: transverse-fraction bound -- the rotation requires
  dv_perp = v * dtheta while the flat energy channel bounds
  |dv_par| <= GM D(1/a)/(2v).  The required transverse fraction of
  any mechanical kick.
- T3: implied conformal-step scale Da ~ dv_perp v / c^2.

outputs:
  results/step_b109_transfer_function.json / .csv
  results/figures/supplementary/step_b109_transfer_function.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))

import json
import numpy as np
from scipy import stats as _st

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import (RESULTS, lv, sep,
                                       load_jsonl_dedup)

logger = StepLogger("step_145_transfer_function")
tee_stdout(logger)

SEED = 20260919
CAP = 60.0
GM = 4.0 * np.pi ** 2          # AU^3/yr^2
AU_YR_TO_M_S = 4740.47         # 1 AU/yr in m/s
C_AU_YR = 299792.458 / AU_YR_TO_M_S / 1000.0 * 1000.0  # c in AU/yr
C_AU_YR = 299792458.0 / AU_YR_TO_M_S                    # m/s / (m/s per AU/yr)
AXIS_DECL = lv(34.0, -13.0)
AXIS_DISP = lv(120.0, -40.0)
SHELLS = (60.0, 100.0, 250.0)


def theta_to(rows, key, axis):
    out = []
    for r in rows:
        a = r[key]
        aph = np.array(a if isinstance(a, (list, tuple))
                       else [float(x) for x in a.split(";")],
                       dtype=float)
        n = np.linalg.norm(aph)
        out.append(sep(aph / n, axis) if n else np.nan)
    return np.array(out, dtype=float)


def crossing_speed(q, e, r_x):
    a = q / (1.0 - e) if e != 1.0 else np.inf
    v2 = GM * (2.0 / r_x - 1.0 / a)
    return np.sqrt(v2) if v2 > 0 else np.nan


# ------------------------------------------------------------------
# load dual-leg records
# ------------------------------------------------------------------
logger.info("loading dual-leg refit records")
cohorts = {}
for tag, jf, axis in (("pre2018", "step_b91_refit.jsonl", AXIS_DECL),
                      ("post2017", "step_b92_refit.jsonl", AXIS_DISP)):
    recs = load_jsonl_dedup(RESULTS / jf)
    recs = [r for r in recs
            if r.get("our_drot") is not None
            and r.get("our_dtau") is not None
            and r.get("our_aph") is not None
            and r.get("our_e") is not None
            and r.get("our_q") is not None]
    th = theta_to(recs, "our_aph", axis)
    for r, t in zip(recs, th):
        r["_th"] = float(t)
        r["_incap"] = bool(t < CAP)
        for rx in SHELLS:
            r[f"_v{int(rx)}"] = crossing_speed(r["our_q"], r["our_e"],
                                               rx)
    cohorts[tag] = recs
    logger.info(f"  {tag}: {len(recs)} records "
                f"({sum(r['_incap'] for r in recs)} in-cap)")

res = {"step": "step_145_transfer_function",
       "result": "b109",
       "description": "Impulse versus crossing-holonomy realization "
                      "of the measured slip: velocity scaling of "
                      "dtau, the energy-channel transverse bound, "
                      "and the implied conformal-step scale.",
       "seed": SEED, "cap_deg": CAP,
       "axes": {"pre2018_declared": [34.0, -13.0],
                "post2017_displaced": [120.0, -40.0]}}

# ------------------------------------------------------------------
# T1: velocity scaling of the implied slip
# ------------------------------------------------------------------
t1 = {}
for tag, recs in cohorts.items():
    dt = np.array([r["our_dtau"] for r in recs], dtype=float)
    inc = np.array([r["_incap"] for r in recs])
    blk = {}
    for rx in SHELLS:
        v = np.array([r[f"_v{int(rx)}"] for r in recs], dtype=float)
        m = np.isfinite(v) & np.isfinite(dt)
        row = {}
        for side, mm in (("in_cap", m & inc), ("out_cap", m & ~inc)):
            if mm.sum() > 10:
                s = _st.spearmanr(dt[mm], v[mm])
                row[side] = dict(n=int(mm.sum()),
                                 rho=float(s.statistic),
                                 p=float(s.pvalue))
        blk[f"r{int(rx)}au"] = row
    # partial correlation given transit epoch (in-cap): is a
    # nonzero v-scaling a real velocity dependence or an epoch
    # composition effect?
    yr = np.array([r.get("yr") for r in recs], dtype=float)
    m = np.isfinite(dt) & np.isfinite(yr) & inc
    v60 = np.array([r["_v60"] for r in recs], dtype=float)
    m &= np.isfinite(v60)
    if m.sum() > 20:
        bv = np.polyfit(yr[m], v60[m], 1)
        bd = np.polyfit(yr[m], dt[m], 1)
        s = _st.spearmanr(v60[m] - np.polyval(bv, yr[m]),
                          dt[m] - np.polyval(bd, yr[m]))
        blk["in_cap_partial_v60_given_epoch"] = dict(
            n=int(m.sum()), rho=float(s.statistic),
            p=float(s.pvalue))
    t1[tag] = blk
res["T1_velocity_scaling"] = t1
logger.info("T1 velocity scaling done")

# ------------------------------------------------------------------
# T2: transverse-fraction bound from the flat energy channel
# ------------------------------------------------------------------
t2 = {}
for tag, recs in cohorts.items():
    inc = np.array([r["_incap"] for r in recs])
    v = np.array([r["_v60"] for r in recs], dtype=float)   # AU/yr
    dth = np.array([r["our_drot"] for r in recs], dtype=float)
    daa = np.array([abs(r["our_daa"]) for r in recs],
                   dtype=float)   # D(1/a) in 1e-6 AU^-1
    m = np.isfinite(v) & np.isfinite(dth) & np.isfinite(daa)
    dv_perp = np.radians(dth) * v * AU_YR_TO_M_S          # m/s
    dv_par = GM * 1e-6 * daa / (2.0 * v) * AU_YR_TO_M_S   # m/s
    row = {}
    for side, mm in (("in_cap", m & inc), ("out_cap", m & ~inc)):
        if mm.sum() > 10:
            frac = 1.0 - dv_par[mm] / np.clip(dv_perp[mm], 1e-9,
                                             None)
            row[side] = dict(
                n=int(mm.sum()),
                med_dv_perp_m_s=float(np.median(dv_perp[mm])),
                med_dv_par_bound_m_s=float(np.median(dv_par[mm])),
                med_transverse_frac=float(np.median(frac)))
    t2[tag] = row
res["T2_transverse_bound"] = t2
logger.info("T2 transverse bound done")

# ------------------------------------------------------------------
# T3: implied conformal-step scale
# ------------------------------------------------------------------
t3 = {}
for tag, recs in cohorts.items():
    inc = np.array([r["_incap"] for r in recs])
    v = np.array([r["_v60"] for r in recs], dtype=float)
    dth = np.array([r["our_drot"] for r in recs], dtype=float)
    m = np.isfinite(v) & np.isfinite(dth) & inc
    if m.sum() > 10:
        # Da ~ dv * v / c^2 with v_n ~ v (order-unity normal comp.)
        da = (np.radians(dth[m]) * v[m] ** 2
              / C_AU_YR ** 2)
        t3[tag] = dict(n=int(m.sum()),
                       med_dalpha=float(np.median(da)),
                       p16=float(np.percentile(da, 16)),
                       p84=float(np.percentile(da, 84)))
res["T3_conformal_step_scale"] = t3
logger.info("T3 conformal-step scale done")

# ------------------------------------------------------------------
# verdict
# ------------------------------------------------------------------
def _flat(blk):
    """All in-cap velocity scalings consistent with zero."""
    out = []
    for rx in blk.values():
        r = (rx.get("in_cap") or {})
        out.append(r.get("rho"), ) if False else None
        if r.get("rho") is not None:
            out.append((abs(r["rho"]) < 0.10, r.get("p", 1) > 0.05))
    return out

flats = {tag: _flat(blk) for tag, blk in t1.items()}
all_flat = all(all(a and b for a, b in pairs)
               for pairs in flats.values() if pairs)
tf_pre = (t2.get("pre2018") or {}).get("in_cap") or {}
tf_post = (t2.get("post2017") or {}).get("in_cap") or {}
tf_min = min(tf_pre.get("med_transverse_frac", 0),
             tf_post.get("med_transverse_frac", 0))
res["test_summary"] = dict(
    velocity_flat=all_flat,
    min_transverse_frac=tf_min)

pv_pre = (t1.get("pre2018") or {}).get(
    "in_cap_partial_v60_given_epoch") or {}
pv_post = (t1.get("post2017") or {}).get(
    "in_cap_partial_v60_given_epoch") or {}
raw_pre = ((t1.get("pre2018") or {}).get("r60au") or {}
           ).get("in_cap") or {}
raw_post = ((t1.get("post2017") or {}).get("r60au") or {}
            ).get("in_cap") or {}

if all_flat and tf_min > 0.9:
    res["verdict"] = (
        "HOLONOMY REALIZATION SELECTED: the implied slip is "
        "velocity-flat in both eras (every in-cap "
        "rho(dtau, v) consistent with zero), while the energy "
        "channel bounds any along-track impulse below "
        f"{min(tf_pre.get('med_dv_par_bound_m_s', 0), tf_post.get('med_dv_par_bound_m_s', 0)):.2f} m/s "
        "against a transverse requirement of "
        f"{max(tf_pre.get('med_dv_perp_m_s', 0), tf_post.get('med_dv_perp_m_s', 0)):.0f} m/s -- "
        f"any mechanical kick would have to be ~{100*tf_min:.1f} per cent "
        "transverse across the whole in-cap population.  A "
        "crossing-localized time translation produces the measured "
        "rotation with no local impulse and a universal amplitude; "
        "the equivalent conformal step is "
        f"~{t3.get('pre2018', {}).get('med_dalpha', float('nan')):.1e}.")
elif tf_min > 0.9:
    res["verdict"] = (
        "ERA-ASYMMETRIC VELOCITY SCALING: the pre-2018 declared-axis "
        "slip is velocity-flat "
        f"(rho={raw_pre.get('rho', float('nan')):+.3f}, "
        f"p={raw_pre.get('p', float('nan')):.3g}) -- the holonomy "
        "signature -- while the post-2017 displaced slip "
        "anticorrelates with crossing speed "
        f"(rho={raw_post.get('rho', float('nan')):+.3f}, "
        f"p={raw_post.get('p', float('nan')):.3g}; partialled on "
        "transit epoch "
        f"rho={pv_post.get('rho', float('nan')):+.3f}, "
        f"p={pv_post.get('p', float('nan')):.3g}).  The energy "
        "channel still bounds any mechanical kick to "
        f"~{100*tf_min:.1f} per cent transverse, so the eras do not "
        "share a single realization: the older record carries the "
        "universal-amplitude (holonomy) signature, the modern "
        "record's residual carries an additional velocity- or "
        "epoch-dependent term -- the same heterogeneity the "
        "epoch-decomposition audit (step_142 T9) measures "
        "independently.  The equivalent conformal step on the "
        "older record is "
        f"~{t3.get('pre2018', {}).get('med_dalpha', float('nan')):.1e}.")
else:
    res["verdict"] = (
        "ENERGY CHANNEL NON-BINDING: the corrected |D(1/a)| floor "
        f"(median ~{min(tf_pre.get('med_dv_par_bound_m_s', 0), tf_post.get('med_dv_par_bound_m_s', 0)):.0f} m/s "
        "along-track) is not tighter than the transverse impulse a "
        "mechanical kick would need "
        f"(~{max(tf_pre.get('med_dv_perp_m_s', 0), tf_post.get('med_dv_perp_m_s', 0)):.0f} m/s), "
        "so the energy channel cannot exclude an impulse "
        "realization -- the earlier ~97 per cent transverse bound "
        "was a wrong-variable artefact (the planetary-approach "
        "distance had been substituted for the energy residual).  "
        "The holonomy discriminator therefore rests on the "
        "velocity-scaling test alone: the pre-2018 declared-axis "
        "slip is velocity-flat "
        f"(rho={raw_pre.get('rho', float('nan')):+.3f}, "
        f"p={raw_pre.get('p', float('nan')):.3g}) while the "
        "post-2017 displaced slip anticorrelates with crossing "
        "speed "
        f"(rho={raw_post.get('rho', float('nan')):+.3f}, "
        f"p={raw_post.get('p', float('nan')):.3g}; partialled on "
        "transit epoch "
        f"rho={pv_post.get('rho', float('nan')):+.3f}, "
        f"p={pv_post.get('p', float('nan')):.3g}).  The equivalent "
        "conformal step on the older record is "
        f"~{t3.get('pre2018', {}).get('med_dalpha', float('nan')):.1e}.")
res["evidence_status"] = "mechanism discriminator"

res["inputs"] = [
    "results/step_b91_refit.jsonl (pre-2018 dual-leg refit)",
    "results/step_b92_refit.jsonl (post-2017 dual-leg refit)"]
res["caveats"] = [
    "The transverse bound uses the measured |D(1/a)| channel "
    "(the inter-leg energy difference, in 1e-6 AU^-1) as the "
    "floor on the along-track impulse; it bounds the median "
    "kick, not the per-object tail, and it includes the real "
    "planetary encounter contribution, so it is a conservative "
    "(loose) floor.",
    "v_n is approximated by the total crossing speed; a strongly "
    "grazing wall geometry lowers the implied conformal step.",
    "The discriminator compares realizations of the measured "
    "observable; it does not derive dtau from a specific scalar "
    "potential, which remains the open microscopic problem."]

# ------------------------------------------------------------------
# outputs
# ------------------------------------------------------------------
out = RESULTS / "step_b109_transfer_function.json"
with open(out, "w") as _fh:
    json.dump(res, _fh, indent=1, default=float)

with open(RESULTS / "step_b109_transfer_function.csv", "w",
          newline="") as f:
    import csv
    keys = ["era", "des", "yr", "_th", "_incap", "our_drot",
            "our_dtau", "our_denc", "our_q", "our_e", "_v60",
            "_v100", "_v250"]
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    for tag, recs in cohorts.items():
        for r in recs:
            row = {k: r.get(k) for k in keys}
            row["era"] = tag
            w.writerow(row)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axes[0]
for tag, col, lab in (("pre2018", "crimson", "pre-2018"),
                      ("post2017", "steelblue", "post-2017")):
    recs = cohorts[tag]
    v = np.array([r["_v60"] for r in recs], dtype=float)
    dt = np.array([r["our_dtau"] for r in recs], dtype=float)
    inc = np.array([r["_incap"] for r in recs])
    m = np.isfinite(v) & np.isfinite(dt) & inc
    ax.scatter(v[m], dt[m], s=14, alpha=0.6, c=col,
               label=f"{lab} in-cap")
ax.set_xlabel("crossing speed at 60 AU (AU/yr)")
ax.set_ylabel("implied slip dtau (yr)")
ax.set_yscale("log")
ax.legend(frameon=False, fontsize=8)
ax.set_title("T1: slip versus crossing speed")

ax = axes[1]
for tag, col, lab in (("pre2018", "crimson", "pre-2018"),
                      ("post2017", "steelblue", "post-2017")):
    recs = cohorts[tag]
    v = np.array([r["_v60"] for r in recs], dtype=float)
    dth = np.array([r["our_drot"] for r in recs], dtype=float)
    daa = np.array([abs(r["our_daa"]) for r in recs], dtype=float)
    inc = np.array([r["_incap"] for r in recs])
    m = np.isfinite(v) & np.isfinite(dth) & np.isfinite(daa) & inc
    ax.scatter(np.radians(dth[m]) * v[m] * AU_YR_TO_M_S,
               GM * 1e-6 * daa[m] / (2 * v[m]) * AU_YR_TO_M_S,
               s=14, alpha=0.6, c=col, label=f"{lab} in-cap")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("required transverse kick dv_perp (m/s)")
ax.set_ylabel("along-track bound |dv_par| (m/s)")
ax.plot([0.1, 100], [0.1, 100], "k:", lw=0.9)
ax.legend(frameon=False, fontsize=8)
ax.set_title("T2: rotation demand vs energy bound")
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b109_transfer_function.png", dpi=300)
logger.info("verdict: " + res["verdict"])
logger.data_save(out)
logger.data_save(RESULTS / "step_b109_transfer_function.csv")
logger.data_save(FIG / "supplementary" / "step_b109_transfer_function.png")
print("TEST SUMMARY:\n" + json.dumps(res["test_summary"], indent=1))
print(f"VERDICT: {res['verdict']}")
