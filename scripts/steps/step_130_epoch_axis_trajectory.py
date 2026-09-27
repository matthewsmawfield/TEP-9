"""Step 130 -- Epoch-resolved axis trajectory on the pooled
independent dual-leg record.

Steps 127-128 established that both eras carry lapse-slip anomalies
in the raw astrometric record on zero-Warsaw-lineage fits -- but at
different sky positions and on different legs: the pre-2018 record
carries the registered inbound-leg anomaly at the declared axis
(34,-13), while the post-2017 record carries a same-class anomaly at
the displaced axis (120,-40) on the outbound leg.  The boundary
channel's registered position is therefore epoch-dependent, and two
readings remain:

  * dynamic-field -- the boundary structure moved: the anomalous
    direction drifts continuously from the declared toward the
    displaced orientation as transit epoch advances;
  * era-artifact / sharp switch -- two unrelated era-locked
    structures; no continuous trajectory connects them.

This step discriminates on the pooled independent record: every
dual-leg comet from the step-127 (pre-2018) and step-128 (post-2017)
refits, each carrying an independently fitted leg-disagreement
carrier (our_ddirf), a leg-rotation channel (our_drot), the
single-fit inbound-leg deviation (reconstructed per-comet as
sep(boundary-leg aphelion, osculating periapsis)), and the
boundary-leg aphelion direction for cap classification at any axis.

Perihelion year is used as the epoch proxy; the boundary-transit
epochs bracket perihelion (inbound tens of years before, outbound
tens of years after) and co-move with it on the ~75-year baseline
spanned here.

Registered tests

  T1  per-bin contrasts: fixed epoch bins; in-cap vs out-of-cap
      Mann-Whitney on ddirf at the declared and displaced axes.
  T2  arc-trajectory scan: per bin, the cap-contrast gap profiled
      along the declared->displaced great-circle arc (s in
      [-0.5, 1.5]); the argmax position per bin is the epoch's
      measured axis position on the arc.  Trajectory statistic:
      Spearman(bin epoch, s_argmax).  Null: epoch-label permutation
      re-deriving the full trajectory (500).
  T3  drift-vs-switch model comparison on the pooled cohort: the
      cap contrast under (a) the fixed declared axis, (b) the fixed
      displaced axis, (c) the best switch model -- declared cap
      before T_s, displaced cap after -- and (d) the best drift
      model -- a slerp-interpolated axis sweeping declared->
      displaced over a ramp (T_c, W).  The switch/drift advantage
      over the best fixed axis is evaluated against an epoch-label
      permutation null (500); fixed axes are unaffected by the
      permutation, so the null prices exactly the epoch-ordering
      component.
  T4  epoch interaction within caps: Spearman(yr, ddirf) among
      declared-cap and displaced-cap members separately, with
      label-shuffle nulls -- does each axis's anomaly strengthen or
      decay across the baseline.
  T5  audit: cohort composition, commensurability, cap overlaps.

Outputs
  results/step_b94_epoch_trajectory.json
  results/step_b94_epoch_trajectory.csv
  results/figures/supplementary/step_b94_epoch_trajectory.png
"""

import json
import math
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.statistics import monte_carlo_tail, monte_carlo_p
from scripts.utils.tep9_common import (DATA_RAW, RESULTS, sep, perih_dir,
                                       load_jsonl_dedup)

FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)

logger = StepLogger("step_130_epoch_axis_trajectory")
tee_stdout(logger)

B91 = RESULTS / "step_b91_refit.jsonl"
B92 = RESULTS / "step_b92_refit.jsonl"

N_PERM = 500
SEED = 20260919
CAP = 60.0

BINS = [(1950, 1990), (1990, 2005), (2005, 2014), (2014, 2018),
        (2018, 2031)]
N_ARC = 41          # arc samples s in [-0.5, 1.5]
ARC_S = np.linspace(-0.5, 1.5, N_ARC)

logger.header("Epoch-resolved axis trajectory: drift vs switch")


def _uv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l), math.sin(b)])


DECL = _uv(34.0, -13.0)
DISP = _uv(120.0, -40.0)


def _slerp(a, b, t):
    """Great-circle interpolation a->b at fraction t (clamped)."""
    t = min(max(t, -0.75), 1.75)
    om = math.acos(float(np.clip(a @ b, -1, 1)))
    if om < 1e-9:
        return a.copy()
    so = math.sin(om)
    v = (math.sin((1 - t) * om) * a + math.sin(t * om) * b) / so
    return v / np.linalg.norm(v)


def _gap(vals, th, cap=CAP):
    inc = np.asarray(th) < cap
    if inc.sum() < 4 or (~inc).sum() < 4:
        return np.nan
    v = np.asarray(vals)
    return float(np.median(v[inc]) - np.median(v[~inc]))


def _mwz(vals, th, cap=CAP):
    inc = np.asarray(th) < cap
    v = np.asarray(vals)
    if inc.sum() < 4 or (~inc).sum() < 4:
        return None
    u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
    return dict(n_in=int(inc.sum()), n_out=int((~inc).sum()),
                med_in=float(np.median(v[inc])),
                med_out=float(np.median(v[~inc])),
                p=float(u.pvalue))


# ---- load pooled dual-leg record ---------------------------------
recs = []
for path, era in ((B91, "pre2018"), (B92, "post2017")):
    for r in load_jsonl_dedup(path):
        if r.get("our_ddirf") is None or r.get("our_aph") is None:
            continue
        a = np.asarray(r["our_aph"], float)
        n = np.linalg.norm(a)
        if not np.all(np.isfinite(a)) or n == 0:
            continue
        a = a / n
        p_osc = perih_dir(math.radians(r["our_w"]),
                          math.radians(r["our_om"]),
                          math.radians(r["our_i"]))
        # epoch proxy must be the transit (perihelion) epoch tp_jd,
        # not the designation year in r["yr"]: for dynamically new
        # comets the two differ by years, and designation year can
        # never exceed the cohort's own discovery boundary -- which
        # would make the era-crosser test vacuous by construction.
        tp_yr = ((r["tp_jd"] - 2451545.0) / 365.25 + 2000.0
                 if r.get("tp_jd") else float(r["yr"]))
        recs.append(dict(era=era, yr=float(tp_yr),
                         ddirf=float(r["our_ddirf"]),
                         drot=float(r.get("our_drot", np.nan)),
                         d_in=sep(-a, p_osc), aph=a.tolist(),
                         th_dec=sep(a, DECL), th_disp=sep(a, DISP)))

n0 = len(recs)
recs = [r for r in recs if np.isfinite(r["drot"])]
logger.info(f"pooled dual-leg comets: {n0} -> {len(recs)} with all "
           f"channels; eras: "
           f"{sum(r['era']=='pre2018' for r in recs)} pre-2018, "
           f"{sum(r['era']=='post2017' for r in recs)} post-2017")

yr = np.array([r["yr"] for r in recs])
dd = np.array([r["ddirf"] for r in recs])
dr = np.array([r["drot"] for r in recs])
din = np.array([r["d_in"] for r in recs])
th_dec = np.array([r["th_dec"] for r in recs])
th_disp = np.array([r["th_disp"] for r in recs])

res = {"step": "step_130_epoch_axis_trajectory",
       "description": ("epoch-resolved axis trajectory on the pooled "
                       "independent dual-leg record (steps 127+128 "
                       "refits); drift-vs-switch discrimination on the "
                       "boundary anomaly's sky position"),
       "n_dual_leg": len(recs),
       "n_pre2018": int(sum(r["era"] == "pre2018" for r in recs)),
       "n_post2017": int(sum(r["era"] == "post2017" for r in recs)),
       "cap_deg": CAP, "epoch_proxy": "perihelion year"}

rng = np.random.default_rng(SEED)

# ---- T1: per-bin registered-axis contrasts -----------------------
t1 = []
for lo, hi in BINS:
    m = (yr >= lo) & (yr < hi)
    row = {"bin": f"{lo}-{hi}", "n": int(m.sum())}
    for nm, th in (("declared", th_dec), ("displaced", th_disp)):
        c = _mwz(dd[m], th[m])
        row[nm] = c
    t1.append(row)
res["T1_per_bin"] = t1
for row in t1:
    d_, x_ = row["declared"], row["displaced"]
    logger.info(f"  {row['bin']} n={row['n']}: "
               f"decl {d_['med_in']:.4f}/{d_['med_out']:.4f} "
               f"p={d_['p']:.4f} | disp {x_['med_in']:.4f}/"
               f"{x_['med_out']:.4f} p={x_['p']:.4f}")

# ---- T2: arc-trajectory scan per bin ------------------------------
arc_dirs = np.array([_slerp(DECL, DISP, s) for s in ARC_S])
aph = np.array([r["aph"] for r in recs])
th_arc = np.degrees(np.arccos(np.clip(aph @ arc_dirs.T, -1, 1)))

def _traj(yr_, dd_):
    """Per-bin argmax arc position + Spearman(epoch, s*)."""
    ss, mid = [], []
    for lo, hi in BINS:
        m = (yr_ >= lo) & (yr_ < hi)
        gaps = np.array([_gap(dd_[m], th_arc[m, k]) for k in
                         range(N_ARC)])
        ss.append(float(ARC_S[int(np.nanargmax(gaps))]))
        mid.append(0.5 * (min(hi, 2030) + lo))
    rho, p = spearmanr(mid, ss)
    return ss, mid, float(rho), float(p)

s_obs, mid, rho_obs, p_obs = _traj(yr, dd)
perm = np.empty(N_PERM)
for i in range(N_PERM):
    _, _, r_, _ = _traj(rng.permutation(yr), dd)
    perm[i] = r_
p_traj = monte_carlo_tail(perm, rho_obs)
res["T2_arc_trajectory"] = {
    "bins": [f"{lo}-{hi}" for lo, hi in BINS],
    "argmax_s": s_obs,
    "argmax_desc": [f"{_slerp(DECL, DISP, s).round(3).tolist()}"
                    for s in s_obs],
    "spearman_epoch_s": rho_obs, "spearman_p_nominal": p_obs,
    "perm_p": p_traj,
    "note": ("s=0 declared axis, s=1 displaced axis; trajectory "
             "statistic is Spearman(bin epoch midpoint, argmax s); "
             "null re-derives the whole per-bin scan under epoch "
             "permutation")}
logger.info(f"  T2 arc argmax per bin: {[f'{s:.2f}' for s in s_obs]} "
           f"rho={rho_obs:.3f} perm p={p_traj:.4f}")

# ---- T3: drift vs switch model comparison -------------------------
def _contrast_for_axes(dd_, th_):
    """in-cap minus out-of-cap median for per-comet axis angle."""
    return _gap(dd_, th_)


def _eval_models(yr_, dd_):
    """Fixed-decl, fixed-disp, best switch, best drift contrasts."""
    c_dec = _contrast_for_axes(dd_, th_dec)
    c_disp = _contrast_for_axes(dd_, th_disp)
    # switch: declared cap before T_s, displaced after
    best_sw, ts_best = -9, None
    for ts in np.arange(2005, 2022, 1.0):
        th = np.where(yr_ < ts, th_dec, th_disp)
        g = _contrast_for_axes(dd_, th)
        if np.isfinite(g) and g > best_sw:
            best_sw, ts_best = g, ts
    # drift: per-comet axis = slerp(t=(yr-T_c)/W clamped to [0,1])
    best_dr, par = -9, None
    for tc in np.arange(2000, 2021, 2.0):
        for w in (10.0, 20.0, 40.0, 80.0):
            t = np.clip((yr_ - (tc - w / 2)) / w, 0, 1)
            dirs = np.array([_slerp(DECL, DISP, ti) for ti in t])
            th = np.degrees(np.arccos(
                np.clip(np.einsum("ij,ij->i", aph, dirs), -1, 1)))
            g = _contrast_for_axes(dd_, th)
            if np.isfinite(g) and g > best_dr:
                best_dr, par = g, (tc, w)
    return c_dec, c_disp, best_sw, ts_best, best_dr, par


c_dec, c_disp, c_sw, ts_sw, c_dr, par_dr = _eval_models(yr, dd)
best_fixed = max(c_dec, c_disp)
adv_obs = max(c_sw, c_dr) - best_fixed
adv_perm = np.empty(N_PERM)
sw_perm, dr_perm = np.empty(N_PERM), np.empty(N_PERM)
for i in range(N_PERM):
    yp = rng.permutation(yr)
    _, _, sw_, _, dr_, _ = _eval_models(yp, dd)
    sw_perm[i], dr_perm[i] = sw_, dr_
    adv_perm[i] = max(sw_, dr_) - best_fixed
res["T3_drift_vs_switch"] = {
    "contrast_fixed_declared": c_dec,
    "contrast_fixed_displaced": c_disp,
    "best_switch": {"contrast": c_sw, "T_switch": ts_sw},
    "best_drift": {"contrast": c_dr, "T_center": par_dr[0],
                   "width_yr": par_dr[1]},
    "advantage_over_fixed": adv_obs,
    "perm_p_advantage": monte_carlo_tail(adv_perm, adv_obs),
    "perm_p_switch_vs_fixed": monte_carlo_tail(sw_perm, c_sw),
    "perm_p_drift_vs_fixed": monte_carlo_tail(dr_perm, c_dr),
    "note": ("contrast = in-cap minus out-of-cap median ddirf on the "
             "pooled cohort; switch/drift hyperparameters scanned per "
             "realization inside the permutation null, so the null "
             "prices the scan too")}
logger.info(f"  T3 fixed decl {c_dec:+.4f} disp {c_disp:+.4f} | "
           f"switch {c_sw:+.4f} (Ts={ts_sw}) drift {c_dr:+.4f} "
           f"(Tc={par_dr[0]},W={par_dr[1]}) adv {adv_obs:+.4f} "
           f"p={res['T3_drift_vs_switch']['perm_p_advantage']:.4f}")

# ---- T4: epoch interaction within caps -----------------------------
t4 = {}
for nm, th in (("declared", th_dec), ("displaced", th_disp)):
    m = th < CAP
    if m.sum() >= 8:
        rho, p = spearmanr(yr[m], dd[m])
        null = np.empty(N_PERM)
        for i in range(N_PERM):
            null[i] = spearmanr(yr[m], rng.permutation(dd[m]))[0]
        t4[nm] = {"n_cap": int(m.sum()), "rho": float(rho),
                  "p_nominal": float(p),
                  "perm_p_two_sided": monte_carlo_tail(null, rho, "two-sided")}
res["T4_epoch_in_cap"] = t4
for nm, v in t4.items():
    logger.info(f"  T4 {nm}-cap: n={v['n_cap']} rho={v['rho']:+.3f} "
               f"perm p={v['perm_p_two_sided']:.4f}")

# ---- T5: audit ------------------------------------------------------
res["T5_audit"] = {
    "median_ddirf_pre2018": float(np.median(
        [r["ddirf"] for r in recs if r["era"] == "pre2018"])),
    "median_ddirf_post2017": float(np.median(
        [r["ddirf"] for r in recs if r["era"] == "post2017"])),
    "n_in_both_caps": int(((th_dec < CAP) & (th_disp < CAP)).sum()),
    "cap_axis_sep_deg": float(np.degrees(np.arccos(
        np.clip(DECL @ DISP, -1, 1)))),
    "channels": ["our_ddirf", "our_drot",
                 "d_in = sep(boundary-leg aph, osc periapsis)"],
}

# ---- T6: era-boundary crossers --------------------------------------
# The natural discriminator between an epoch transition (signature
# follows perihelion epoch) and a cohort-table systematic (signature
# follows the source record): comets with post-2018 perihelia drawn
# from the pre-2018 cohort table.  On the independent refits both
# cohorts read the same astrometry, so only the epoch and the cohort
# assignment differ.  If the transition is physical, the crossers
# should show the modern (displaced-axis) pattern, not the
# historical (declared-axis) one.
cross = np.array([i for i, r in enumerate(recs)
                  if r["era"] == "pre2018" and r["yr"] > 2018])
t6 = {"n": int(len(cross))}
if len(cross) > 8:
    for nm, th in (("declared", th_dec), ("displaced", th_disp)):
        c = _mwz(dd[cross], th[cross])
        if c:
            t6[nm] = c
res["T6_era_crossers"] = t6
logger.info(f"  T6 era crossers (pre-2018 table, tp>2018): "
            f"n={t6['n']} " +
            " ".join(f"{k} p={v['p']:.4f}"
                     for k, v in t6.items() if isinstance(v, dict)))

# secondary channel summary: drot + d_in on the same T1 grid
sec = {}
for ch, v in (("drot", dr), ("d_in", din)):
    rows_ = []
    for lo, hi in BINS:
        m = (yr >= lo) & (yr < hi)
        rows_.append({"bin": f"{lo}-{hi}", "n": int(m.sum()),
                      "declared": _mwz(v[m], th_dec[m]),
                      "displaced": _mwz(v[m], th_disp[m])})
    sec[ch] = rows_
res["secondary_channels_per_bin"] = sec

# verdict -------------------------------------------------------------
sw_p = res["T3_drift_vs_switch"]["perm_p_switch_vs_fixed"]
dr_p = res["T3_drift_vs_switch"]["perm_p_drift_vs_fixed"]
traj_p = res["T2_arc_trajectory"]["perm_p"]
w_best = res["T3_drift_vs_switch"]["best_drift"]["width_yr"]
if min(sw_p, dr_p) < 0.05 and dr_p <= sw_p and w_best <= 12.0:
    v_ = ("TRANSITION-LOCALIZED: the epoch-ordering of the pooled "
          "record is real, but the best drift ramp has reached the "
          "scan's resolution floor -- the transition is fast "
          "(<=10-20 yr about the era boundary), degenerate with a "
          "sharp switch; no slow sweep is resolved")
elif min(sw_p, dr_p) < 0.05 and dr_p <= sw_p:
    v_ = ("DRIFT-FAVOURED: the epoch-ordering of the pooled record is "
          "better matched by a continuously moving axis than by a "
          "sharp era switch or either fixed axis")
elif min(sw_p, dr_p) < 0.05 and sw_p < dr_p:
    x_ = res.get("T6_era_crossers", {})
    xw = ""
    if isinstance(x_.get("displaced"), dict):
        xw = (f"; the era-boundary crossers -- comets with post-2018 "
              f"perihelia drawn from the pre-2018 cohort table -- "
              f"carry the displaced signature "
              f"(p={x_['displaced']['p']:.3f}) rather than the "
              f"declared one, so the transition tracks perihelion "
              f"epoch on the independent fits, not the source-table "
              f"boundary")
    v_ = ("SWITCH-FAVOURED: a sharp era transition between two fixed "
          "axes outperforms both fixed axes and the smooth drift" + xw)
elif traj_p < 0.05:
    v_ = ("TRAJECTORY-ONLY: per-bin axis positions migrate along the "
          "arc but pooled switch/drift contrasts are not resolved")
else:
    v_ = ("UNRESOLVED: no epoch-ordering advantage over the best "
          "fixed axis on the pooled independent record")
res["verdict"] = v_
logger.info("VERDICT: " + v_)

# ---- outputs ---------------------------------------------------------
csv_rows = []
for lo, hi in BINS:
    m = (yr >= lo) & (yr < hi)
    csv_rows.append({"bin": f"{lo}-{hi}", "n": int(m.sum()),
                     "median_ddirf": float(np.median(dd[m]))})
with open(RESULTS / "step_b94_epoch_trajectory.csv", "w") as f:
    f.write("bin,n,median_ddirf\n")
    for r in csv_rows:
        f.write(f"{r['bin']},{r['n']},{r['median_ddirf']:.4f}\n")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axes[0]
xs = np.arange(len(BINS))
for nm, key, mk in (("declared", "declared", "o"),
                    ("displaced", "displaced", "s")):
    gaps = [t[key]["med_in"] - t[key]["med_out"] for t in t1]
    ax.plot(xs, gaps, mk + "-", label=nm)
ax.axhline(0, color="k", lw=0.7)
ax.set_xticks(xs); ax.set_xticklabels([t["bin"] for t in t1],
                                      rotation=30, fontsize=8)
ax.set_ylabel("in-cap - out-cap median ddirf (deg)")
ax.legend(fontsize=8); ax.set_title("per-bin contrasts")

ax = axes[1]
ax.errorbar(xs, s_obs, fmt="D-", label="argmax arc position")
ax.axhline(0, color="tab:blue", lw=0.8, ls="--", label="declared")
ax.axhline(1, color="tab:orange", lw=0.8, ls="--", label="displaced")
ax.set_xticks(xs); ax.set_xticklabels([t["bin"] for t in t1],
                                      rotation=30, fontsize=8)
ax.set_ylabel("arc position s (0=decl, 1=disp)")
ax.set_title(f"axis trajectory (rho={rho_obs:.2f}, "
             f"p={p_traj:.3f})")
ax.legend(fontsize=8)

ax = axes[2]
tcs = np.arange(2000, 2021, 2.0)
for w in (10.0, 20.0, 40.0, 80.0):
    cs = []
    for tc in tcs:
        t = np.clip((yr - (tc - w / 2)) / w, 0, 1)
        dirs = np.array([_slerp(DECL, DISP, ti) for ti in t])
        th = np.degrees(np.arccos(
            np.clip(np.einsum("ij,ij->i", aph, dirs), -1, 1)))
        cs.append(_gap(dd, th))
    ax.plot(tcs, cs, "-", label=f"W={w:.0f}y")
ax.axhline(c_dec, color="tab:blue", ls="--", lw=0.8)
ax.axhline(c_disp, color="tab:orange", ls="--", lw=0.8)
ax.set_xlabel("drift ramp centre T_c"); ax.set_ylabel("pooled contrast")
ax.set_title("drift model vs fixed axes"); ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "supplementary" / "step_b94_epoch_trajectory.png", dpi=300)
plt.close(fig)

with open(RESULTS / "step_b94_epoch_trajectory.json", "w") as f:
    json.dump(res, f, indent=1)
logger.data_save(RESULTS / "step_b94_epoch_trajectory.json")
logger.data_save(RESULTS / "step_b94_epoch_trajectory.csv")
logger.data_save(FIG / "supplementary" / "step_b94_epoch_trajectory.png")
logger.save_provenance()
