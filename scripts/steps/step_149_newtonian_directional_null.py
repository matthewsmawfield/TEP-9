"""step_149: Newtonian directional null (result b113).

step_063 showed that full planetary integrations reproduce the
catalogued orig->fut periapsis rotation per comet (rho = 0.999), and
step_042/063 regressed that rotation on a four-covariate encounter
budget (energy kick, closest approach, perihelion depth, inclination).
What has never been done is to apply the COMPLETE directional-analysis
pipeline -- the same battery the manuscript applies to the observed
discrepancy -- to the pure Newtonian rotation field itself.  This step
closes that gap.

For every channel v (per-comet scalar) the identical battery used on
the observed signal is run:

  T1  declared-cap contrast: Mann-Whitney (greater) inside vs outside
      the 60-deg cap about the TNO axis (49,-17) and the declared comet
      axis (34,-13); continuous Spearman vs theta.
  T2  shell morphology: median v in 0-30 / 30-45 / 45-60 / 60-75 deg
      shells about the TNO axis (plateau vs gradient).
  T3  free sky scan on the 5-deg 2232-axis grid (the step_037 look-
      elsewhere grid): best-contrast axis, its separation from the TNO
      and comet axes, and the fraction of trial axes beating the
      declared-axis p.
  T4  label-swap permutation null (20 000 draws): the cap median
      contrast with directions fixed -- the step_037 B2 statistic.
  T5  axis-coincidence null: 2 000 full permutation sky scans (15-deg
      coarse grid, MW z-score criterion).  For each direction-unlinked
      realization the best-contrast axis is recovered and its
      separation from the TNO axis recorded.  This prices the
      coincidence that the rotation field's own best axis lands near
      the independently derived resident axis.

Channels:
  sim  drot_sim -- the pure standard-dynamics rotation (N-body model)
  cat  drot_cat -- the catalogued rotation (observed)
  res  drot_cat - drot_sim -- leftover after subtracting the full
       conventional model (the review's decisive residual channel)
  kick |daa_sim| -- modelled energy kick magnitude
  rpk  drot_sim / |daa_sim| -- rotation per unit kick

  T6  arrival-geometry audit: in-cap fraction vs the cap solid-angle
      fraction, median theta -- whether the arrival directions alone
      predispose a cap result.
  T7  conventional-mediator audit: Spearman of theta against |daa_sim|,
      denc and drot_sim -- which encounter property carries the
      directional structure under standard dynamics.

Inputs:  results/step_b28_bidirectional_rotation.csv,
         data/raw/code/code_original.html
Outputs: results/step_b113_newtonian_null.json / .csv
         results/figures/step_b113_newtonian_null.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (
    DATA_RAW, RESULTS, lv, lb, sep, perih_dir, parse_code, tee_stdout)
logger = StepLogger("step_149_newtonian_directional_null")
tee_stdout(logger)
logger.header("Newtonian directional null -- full pipeline on the model field")

import csv
import json
import math
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr, rankdata, binomtest

SEED = 20260920
CAP = 60.0
N_PERM = 20000
N_AXIS_PERM = 2000
rng = np.random.default_rng(SEED)

TNO = lv(49.0, -17.0)          # resident cluster axis
CAX = lv(34.0, -13.0)          # declared comet transit axis
CAP_FRAC = (1.0 - math.cos(math.radians(CAP))) / 2.0   # 0.25 of the sky

# fine grid identical to step_037 B3 look-elsewhere
FINE = [lv(lam, b) for lam in np.arange(0, 360, 5)
        for b in np.arange(-75, 76, 5)]
# coarse grid identical to step_037 B4 split-half scan
COARSE = [lv(lam, b) for lam in np.arange(0, 360, 15)
          for b in np.arange(-60, 61, 15)]
FINE = np.array(FINE)          # (2232, 3)
COARSE = np.array(COARSE)      # (216, 3)


# ------------------------------------------------------------------
# data
# ------------------------------------------------------------------
orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))

rows = []
for r in csv.DictReader(open(RESULTS / "step_b28_bidirectional_rotation.csv")):
    k = r["desig"]
    if k not in orig:
        continue
    ro = orig[k]
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                   math.radians(ro["i"]))
    u = -po                                  # arrival (inbound) direction
    rows.append(dict(desig=k, cls=r["cls"], q=float(r["q"]),
                     i=float(r["i"]), u=u,
                     theta=float(r["theta"]),
                     cat=float(r["drot_cat"]), sim=float(r["drot_sim"]),
                     kick=abs(float(r["daa_sim"])),
                     denc=float(r["denc"])))

for r in rows:
    r["res"] = r["cat"] - r["sim"]
    r["rpk"] = r["sim"] / r["kick"] if r["kick"] > 0 else np.nan

logger.info(f"joined {len(rows)} comets (step_b28 x CODE original)")
subsets = {"matched": [r for r in rows if r["q"] < 3.1],
           "all_c1": rows}
logger.info(f"matched (q<3.1): {len(subsets['matched'])}; "
            f"all class-1: {len(subsets['all_c1'])}")


# ------------------------------------------------------------------
# directional battery
# ------------------------------------------------------------------
def thetas(U, ax):
    """Angular separation of every direction in U (n,3) from ax."""
    with np.errstate(all="ignore"):   # Accelerate BLAS emits spurious FP flags
        c = np.clip(U @ ax, -1.0, 1.0)
    return np.degrees(np.arccos(c))


def mw_greater(v, inc):
    if inc.sum() < 3 or (~inc).sum() < 3:
        return np.nan
    return float(mannwhitneyu(v[inc], v[~inc], alternative="greater").pvalue)


def channel_pipeline(U, v, label):
    ok = np.isfinite(v) & np.isfinite(U).all(axis=1)
    out = {"label": label, "n": int(ok.sum())}
    U, v = U[ok], v[ok]
    th_tno, th_cax = thetas(U, TNO), thetas(U, CAX)

    # T1 declared-cap contrast
    for tag, th in (("tno", th_tno), ("comet_axis", th_cax)):
        inc = th <= CAP
        rho, p = spearmanr(th, v)
        out[f"T1_{tag}"] = {
            "n_in": int(inc.sum()),
            "med_in": float(np.median(v[inc])),
            "med_out": float(np.median(v[~inc])),
            "p_mw_greater": mw_greater(v, inc),
            "rho": float(rho), "p_spearman": float(p)}

    # T2 shell morphology about the TNO axis
    shells = {}
    for lo, hi in ((0, 30), (30, 45), (45, 60), (60, 75), (75, 180)):
        m = (th_tno > lo) & (th_tno <= hi)
        shells[f"{lo}-{hi}"] = {"n": int(m.sum()),
                                "med": float(np.median(v[m])) if m.any() else None}
    out["T2_shells_tno"] = shells

    # T3 free sky scan (5-deg grid, MW-greater criterion)
    with np.errstate(all="ignore"):   # Accelerate BLAS emits spurious FP flags
        TH = np.degrees(np.arccos(np.clip(U @ FINE.T, -1, 1)))   # (n, 2232)
    best = {"p": 1.0, "ax": None}
    n_beat = 0
    p_decl = out["T1_tno"]["p_mw_greater"]
    n_eval = 0
    for j in range(FINE.shape[0]):
        inc = TH[:, j] <= CAP
        if inc.sum() < 3 or (~inc).sum() < 3:
            continue
        pj = mw_greater(v, inc)
        n_eval += 1
        if pj < best["p"]:
            best = {"p": pj, "ax": FINE[j]}
        if p_decl is not None and np.isfinite(p_decl) and pj <= p_decl:
            n_beat += 1
    out["T3_scan"] = {
        "n_grid": n_eval,
        "best_p": best["p"],
        "best_lb": list(lb(best["ax"])) if best["ax"] is not None else None,
        "sep_best_tno": float(sep(best["ax"], TNO)) if best["ax"] is not None else None,
        "sep_best_comet": float(sep(best["ax"], CAX)) if best["ax"] is not None else None,
        "look_elsewhere_frac": n_beat / n_eval if n_eval else None}

    # T4 label-swap cap null (B2 statistic on the TNO cap)
    inc = th_tno <= CAP
    obs = float(np.median(v[inc]) - np.median(v[~inc]))
    cnt = 0
    for _ in range(N_PERM):
        vs = rng.permutation(v)
        if np.median(vs[inc]) - np.median(vs[~inc]) >= obs:
            cnt += 1
    out["T4_perm"] = {"obs_meddiff": obs, "n_perm": N_PERM,
                      "p_perm": float((cnt + 1) / (N_PERM + 1))}
    return out, U, v, th_tno


def coarse_best(U, v):
    """Best-contrast axis on the 15-deg grid via MW z-score (vectorized).

    Returns (axis, z).  Ranks are computed once; per axis the U
    statistic is a masked rank sum."""
    n = len(v)
    with np.errstate(all="ignore"):   # Accelerate BLAS emits spurious FP flags
        TH = np.degrees(np.arccos(np.clip(U @ COARSE.T, -1, 1)))   # (n,216)
    inc = TH <= CAP                                            # (n,216)
    n_in = inc.sum(axis=0)
    valid = (n_in >= 3) & (n - n_in >= 3)
    ranks = rankdata(v)
    Ustat = (inc * ranks[:, None]).sum(axis=0) - n_in * (n_in + 1) / 2.0
    mu = n_in * (n - n_in) / 2.0
    sigma = np.sqrt(n_in * (n - n_in) * (n + 1) / 12.0)
    z = np.where(valid & (sigma > 0), (Ustat - mu - 0.5) / np.maximum(sigma, 1e-12), -np.inf)
    j = int(np.argmax(z))
    return COARSE[j], float(z[j])


def axis_coincidence(U, v, obs_sep):
    """T5: distribution of best-axis separation from TNO under the
    label-swap null, and of the scan-best z."""
    seps = np.empty(N_AXIS_PERM)
    zmax = np.empty(N_AXIS_PERM)
    for p in range(N_AXIS_PERM):
        ax, z = coarse_best(U, rng.permutation(v))
        seps[p] = sep(ax, TNO)
        zmax[p] = z
    return {"n_perm": N_AXIS_PERM,
            "obs_sep_tno": float(obs_sep),
            "p_sep_le_obs": float((np.sum(seps <= obs_sep) + 1) / (N_AXIS_PERM + 1)),
            "frac_within_30": float(np.mean(seps <= 30)),
            "frac_within_45": float(np.mean(seps <= 45)),
            "frac_within_60": float(np.mean(seps <= 60)),
            "median_sep": float(np.median(seps)),
            "seps": seps, "zmax": zmax}


# ------------------------------------------------------------------
# run
# ------------------------------------------------------------------
CHANNELS = ("sim", "cat", "res", "kick", "rpk")
res = {"seed": SEED, "cap_deg": CAP, "n_perm": N_PERM,
       "n_axis_perm": N_AXIS_PERM,
       "axes": {"tno": [49.0, -17.0], "comet": [34.0, -13.0]},
       "method": (
           "The complete directional battery applied to the pure "
           "Newtonian rotation field (drot_sim) alongside the catalogued "
           "field (drot_cat), the observed-minus-model residual (res), "
           "the energy-kick magnitude (kick) and rotation-per-kick (rpk). "
           "T5 prices the axis coincidence: the fraction of direction-"
           "unlinked rotation fields whose own best-contrast axis lands "
           "as close to the TNO axis as the model field's does.")}

perm_store = {}
for stag, sub in subsets.items():
    res[stag] = {}
    U = np.array([r["u"] for r in sub])
    for ch in CHANNELS:
        v = np.array([r[ch] for r in sub])
        out, Uo, vo, th = channel_pipeline(U, v, f"{stag}.{ch}")
        # T5 only where the channel produces a real best axis
        if out["T3_scan"]["sep_best_tno"] is not None:
            t5 = axis_coincidence(Uo, vo, out["T3_scan"]["sep_best_tno"])
            perm_store[f"{stag}.{ch}"] = t5
            out["T5_axis_coincidence"] = {k: v2 for k, v2 in t5.items()
                                          if k not in ("seps", "zmax")}
        res[stag][ch] = out
        t1, t3 = out["T1_tno"], out["T3_scan"]
        logger.info(
            f"{stag}.{ch:4s} n={out['n']:3d} | cap p={t1['p_mw_greater']} "
            f"rho={t1['rho']:+.3f} | perm p={out['T4_perm']['p_perm']:.4f} | "
            f"best axis {t3['best_lb']} p={t3['best_p']:.3g} "
            f"sep_TNO={t3['sep_best_tno']:.1f} | "
            f"coinc p={out.get('T5_axis_coincidence', {}).get('p_sep_le_obs')}")

# T6 arrival-geometry audit
for stag, sub in subsets.items():
    th = np.array([r["theta"] for r in sub])
    res[stag]["T6_arrival_geometry"] = {
        "n": len(sub), "n_in": int((th <= CAP).sum()),
        "frac_in": float(np.mean(th <= CAP)),
        "cap_sky_frac": CAP_FRAC,
        "p_binom": float(binomtest(int((th <= CAP).sum()), len(sub),
                                   CAP_FRAC).pvalue),
        "med_theta": float(np.median(th))}
    logger.info(f"{stag} arrival geometry: {int((th<=CAP).sum())}/{len(sub)} "
                f"in cap (sky frac {CAP_FRAC:.2f}), "
                f"med theta={np.median(th):.1f}")

# T7 conventional mediators (theta vs encounter properties)
for stag, sub in subsets.items():
    th = np.array([r["theta"] for r in sub])
    med = {}
    for key in ("kick", "denc", "sim"):
        vv = np.array([r[key] for r in sub])
        rho, p = spearmanr(th, vv)
        med[key] = {"rho": float(rho), "p": float(p)}
    res[stag]["T7_mediators"] = med
    logger.info(f"{stag} mediators vs theta: " +
                " ".join(f"{k} rho={m['rho']:+.3f} (p={m['p']:.3g})"
                         for k, m in med.items()))

# model-vs-catalogue axis agreement (the direct null comparison)
axagree = {}
for stag in subsets:
    a_sim = res[stag]["sim"]["T3_scan"]["best_lb"]
    a_cat = res[stag]["cat"]["T3_scan"]["best_lb"]
    if a_sim and a_cat:
        axagree[stag] = float(sep(lv(*a_sim), lv(*a_cat)))
res["axis_recovery_sim_vs_cat_sep"] = axagree
res["finding"] = (
    "If the modelled (pure Newtonian) rotation field reproduces the "
    "in-cap excess, the recovered axis and the shell morphology of the "
    "catalogued signal, the directional structure is a product of "
    "standard dynamics on this arrival geometry; the residual channel "
    "then bounds what is left for non-standard physics. T5 prices how "
    "unlikely the model axis's proximity to the independent TNO axis is "
    "under direction-unlinked fields.")

# ------------------------------------------------------------------
# outputs
# ------------------------------------------------------------------
csv_out = RESULTS / "step_b113_newtonian_null.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["desig", "cls", "q", "i", "l_arr", "b_arr", "theta",
                "drot_cat", "drot_sim", "res", "kick", "denc", "rpk"])
    for r in rows:
        l_, b_ = lb(r["u"])
        w.writerow([r["desig"], r["cls"], r["q"], r["i"],
                    f"{l_:.4f}", f"{b_:.4f}", f"{r['theta']:.4f}",
                    f"{r['cat']:.6f}", f"{r['sim']:.6f}", f"{r['res']:.6f}",
                    f"{r['kick']:.4f}", f"{r['denc']:.4f}",
                    f"{r['rpk']:.4f}" if np.isfinite(r["rpk"]) else ""])

json_out = RESULTS / "step_b113_newtonian_null.json"
with open(json_out, "w") as _fh:
    json.dump(res, _fh, indent=1, default=str)

# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sub = subsets["matched"]
okfig = np.array([np.isfinite(r["u"]).all() for r in sub])
U = np.array([r["u"] for r in sub])[okfig]
vs = np.array([r["sim"] for r in sub])[okfig]
th = np.array([r["theta"] for r in sub])[okfig]
inc = th <= CAP

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axes[0]
ax.scatter(th[~inc], vs[~inc], s=22, c="0.55", alpha=0.75,
           label=f"outside cap (n={int((~inc).sum())})")
ax.scatter(th[inc], vs[inc], s=26, c="crimson", alpha=0.85,
           label=f"inside cap (n={int(inc.sum())})")
ax.axvline(CAP, color="k", ls=":", lw=1)
ax.set_xlabel(r"arrival-direction separation from TNO axis, $\theta$ (deg)")
ax.set_ylabel(r"modelled rotation $d_{\rm sim}$ (deg)")
ax.set_title("Newtonian rotation field vs axis (matched)")
ax.legend(fontsize=8)

ax = axes[1]
with np.errstate(all="ignore"):   # Accelerate BLAS emits spurious FP flags
    THg = np.degrees(np.arccos(np.clip(U @ FINE.T, -1, 1)))
zmap = np.full(FINE.shape[0], np.nan)
ranks = rankdata(vs)
for j in range(FINE.shape[0]):
    ij = THg[:, j] <= CAP
    if ij.sum() < 3 or (~ij).sum() < 3:
        continue
    zmap[j] = -math.log10(max(mw_greater(vs, ij), 1e-16))
sc = ax.scatter([lb(a)[0] for a in FINE], [lb(a)[1] for a in FINE],
                c=zmap, cmap="viridis", s=14)
for ax_, mk, cl, lab in ((TNO, "*", "crimson", "TNO axis"),
                         (CAX, "s", "white", "comet axis")):
    l_, b_ = lb(ax_)
    ax.scatter([l_], [b_], marker=mk, c=cl, s=140, edgecolors="k",
               linewidths=0.6, label=lab, zorder=5)
ax.set_xlabel("ecliptic longitude (deg)"); ax.set_ylabel("ecliptic latitude (deg)")
ax.set_title(r"scan $-\log_{10}p$ on $d_{\rm sim}$ (matched)")
ax.legend(fontsize=8, loc="upper right")
fig.colorbar(sc, ax=ax, label=r"$-\log_{10} p$")

ax = axes[2]
key = "matched.sim"
if key in perm_store:
    seps_p = perm_store[key]["seps"]
    ax.hist(seps_p, bins=40, color="0.6", alpha=0.8,
            label="direction-unlinked null")
    obs = res["matched"]["sim"]["T3_scan"]["sep_best_tno"]
    ax.axvline(obs, color="crimson", lw=1.6,
               label=f"model best-axis sep = {obs:.1f} deg")
    ax.set_xlabel("best-axis separation from TNO axis (deg)")
    ax.set_ylabel("permutation realizations")
    ax.set_title("T5 axis-coincidence null")
    ax.legend(fontsize=8)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b113_newtonian_null.png", dpi=150)

print("wrote", json_out)
print("wrote", csv_out)
print("wrote", FIG / "step_b113_newtonian_null.png")
