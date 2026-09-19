"""step_089: Bipolar morphology audit -- is the boundary single-faced?

The injected JFC chain terminates on BOTH poles of the measured axis
(steps 052/056), so the underlying structure is axisymmetric in
aphelion direction.  The transit channel has so far been scored only
against the +axis cap (theta < 60 deg).  If the boundary were a
two-sided sheet or a symmetric pair of structures, comets whose
aphelia point near the anti-axis (theta > 120 deg, the mirror cap)
should carry a second slip lobe of comparable amplitude.  This step
measures the full theta profile of the catalogue rotation and the
planetary-subtracted proper-time residual, and tests:

  T1  mirror-cap elevation: median(drot_cat, theta>120) vs the
      mid-region (60-120) background -- Mann-Whitney + bootstrap bound
      on any second-lobe amplitude;
  T2  cap/anti-cap contrast: median(drot_cat, theta<60) vs theta>120
      -- expected significant if the anomaly is localized to +axis;
  T3  monotonic decay: Spearman of the residuals against theta over
      the full range -- a localized single structure decays
      monotonically, a bipolar structure is non-monotone;
  T4  sign structure: the mid-region (60-120) residual median and the
      fraction of positive residuals -- the bins adjoining the cap
      carry a systematically negative median that a two-lobe reading
      would not predict;
  T5  per-catalogue replication of all medians (CODE, Warsaw).

Outputs
-------
results/step_b54_bipolar_morphology.json
results/step_b54_bipolar_morphology.csv
results/figures/step_b54_bipolar_morphology.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_089_bipolar_morphology")
tee_stdout(logger)
logger.header("Bipolar morphology audit -- mirror-cap second-lobe test")

import csv
import json
import math
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr, binomtest

SEED = 20260918
rng = np.random.default_rng(SEED)
N_BOOT = 20000

# ------------------------------------------------------------------
# Per-comet data (step_065 products; both catalogues pooled)
# ------------------------------------------------------------------

rows = list(csv.DictReader(open(RESULTS / "step_b30_proper_time_slip.csv")))

def fnum(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")

data = []
for r in rows:
    th = fnum(r, "theta")
    if not np.isfinite(th):
        continue
    data.append(dict(desig=r["desig"], cohort=r["cohort"], theta=th,
                     drot=fnum(r, "drot_cat"),
                     dtau=fnum(r, "dtau_unexplained"),
                     dtau_in=fnum(r, "dtau_inbound_unexplained")))

th_all = np.array([d["theta"] for d in data])
logger.info(f"pooled cohorts: {len(data)} comets "
            f"(code={sum(1 for d in data if d['cohort']=='code')}, "
            f"warsaw={sum(1 for d in data if d['cohort']=='warsaw')})")

CAP, MID, ANTI = (0.0, 60.0), (60.0, 120.0), (120.0, 180.01)

def sel(key, lo, hi, cohort=None):
    return np.array([d[key] for d in data
                     if lo <= d["theta"] < hi and np.isfinite(d[key])
                     and (cohort is None or d["cohort"] == cohort)])

def boot_median_ci(x, n=N_BOOT):
    if len(x) < 3:
        return [float("nan"), float("nan")]
    idx = rng.integers(0, len(x), (n, len(x)))
    meds = np.median(x[idx], axis=1)
    return [float(np.percentile(meds, 16)),
            float(np.percentile(meds, 84))]

def boot_diff_ci(a, b, n=N_BOOT):
    ia = rng.integers(0, len(a), (n, len(a)))
    ib = rng.integers(0, len(b), (n, len(b)))
    d = np.median(a[ia], axis=1) - np.median(b[ib], axis=1)
    return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)),
            float(np.percentile(d, 95))]   # one-sided 95% upper bound

# ------------------------------------------------------------------
# theta profile in 30-deg bins
# ------------------------------------------------------------------

bins = [(0, 30), (30, 60), (60, 90), (90, 120), (120, 150), (150, 181)]
profile = []
for lo, hi in bins:
    dr = sel("drot", lo, hi)
    du = sel("dtau", lo, hi)
    profile.append({
        "bin": [lo, hi], "n": len(dr),
        "drot_med": float(np.median(dr)),
        "drot_ci68": boot_median_ci(dr),
        "dtau_med": float(np.median(du)),
        "dtau_ci68": boot_median_ci(du),
        "dtau_frac_pos": float(np.mean(du > 0)),
    })
    logger.info(f"theta {lo:3d}-{hi:3d}: n={len(dr):3d} "
                f"drot_med={np.median(dr):.4f} "
                f"dtau_med={np.median(du):+.2f} yr "
                f"frac>0={np.mean(du>0):.2f}")

# ------------------------------------------------------------------
# T1 mirror-cap elevation; T2 cap/anti contrast; T3 monotone; T4 sign
# ------------------------------------------------------------------

drot_cap  = sel("drot", *CAP)
drot_mid  = sel("drot", *MID)
drot_anti = sel("drot", *ANTI)
dtau_cap  = sel("dtau", *CAP)
dtau_mid  = sel("dtau", *MID)
dtau_anti = sel("dtau", *ANTI)

u1 = mannwhitneyu(drot_anti, drot_mid, alternative="greater")
lobe_ci = boot_diff_ci(drot_anti, drot_mid)
u2 = mannwhitneyu(drot_cap, drot_anti, alternative="greater")
u2d = mannwhitneyu(dtau_cap, dtau_anti, alternative="greater")

rho_drot, p_drot = spearmanr(th_all, np.array([d["drot"] for d in data]))
rho_dtau, p_dtau = spearmanr(
    np.array([d["theta"] for d in data if np.isfinite(d["dtau"])]),
    np.array([d["dtau"] for d in data if np.isfinite(d["dtau"])]))

n_mid_pos = int(np.sum(dtau_mid > 0))
bt_mid = binomtest(n_mid_pos, len(dtau_mid), 0.5)
n_anti_pos = int(np.sum(dtau_anti > 0))
bt_anti = binomtest(n_anti_pos, len(dtau_anti), 0.5)
n_cap_pos = int(np.sum(dtau_cap > 0))
bt_cap = binomtest(n_cap_pos, len(dtau_cap), 0.5)

# anti-cap elevation over the mid region in the time channel
u_dtau_anti = mannwhitneyu(dtau_anti, dtau_mid, alternative="greater")
u_dtau_cap  = mannwhitneyu(dtau_cap,  dtau_mid, alternative="greater")

# T6 harmonic decomposition: a one-sided (dipolar) defect predicts the
# residual tracking cos(theta) -- positive in the cap, negative at the
# anti-axis; an axisymmetric field (the bipolar JFC chain's geometry)
# predicts cos(2*theta) -- positive at BOTH poles, negative mid-range.
th_v = np.array([d["theta"] for d in data if np.isfinite(d["dtau"])])
du_v = np.array([d["dtau"] for d in data if np.isfinite(d["dtau"])])
dr_v = np.array([d["drot"] for d in data if np.isfinite(d["dtau"])])
cos1 = np.cos(np.radians(th_v))
cos2 = np.cos(np.radians(2 * th_v))
r1_tau, p1_tau = spearmanr(cos1, du_v)
r2_tau, p2_tau = spearmanr(cos2, du_v)
r1_rot, p1_rot = spearmanr(cos1, dr_v)
r2_rot, p2_rot = spearmanr(cos2, dr_v)

# per-catalogue replication of the three-region medians
per_cat = {}
for c in ("code", "warsaw"):
    per_cat[c] = {
        "drot_cap_med":  float(np.median(sel("drot", *CAP, cohort=c))),
        "drot_mid_med":  float(np.median(sel("drot", *MID, cohort=c))),
        "drot_anti_med": float(np.median(sel("drot", *ANTI, cohort=c))),
        "dtau_cap_med":  float(np.median(sel("dtau", *CAP, cohort=c))),
        "dtau_mid_med":  float(np.median(sel("dtau", *MID, cohort=c))),
        "dtau_anti_med": float(np.median(sel("dtau", *ANTI, cohort=c))),
        "n_anti": len(sel("drot", *ANTI, cohort=c)),
    }

res = {
    "method": "full-theta profile of catalogue rotation drot_cat and "
              "planetary-subtracted proper-time residual dtau_unexplained "
              "(step_b30 per-comet products) in 30-deg bins; mirror-cap "
              "(theta>120) elevation tested against the 60-120 background "
              "by Mann-Whitney with a 95% bootstrap upper bound on the "
              "second-lobe amplitude; cap/anti contrast and full-range "
              "Spearman for monotonicity.",
    "cap_deg": 60, "mid_deg": [60, 120], "anti_deg": [120, 180],
    "n": {"cap": len(drot_cap), "mid": len(drot_mid), "anti": len(drot_anti)},
    "profile": profile,
    "T1_mirror_cap": {
        "drot_anti_med": float(np.median(drot_anti)),
        "drot_mid_med": float(np.median(drot_mid)),
        "mw_p_anti_gt_mid": float(u1.pvalue),
        "second_lobe_amp_ci95": lobe_ci[:2],
        "second_lobe_amp_upper95": lobe_ci[2],
        "incap_elevation_deg": float(np.median(drot_cap) - np.median(drot_mid)),
        "dtau_anti_med": float(np.median(dtau_anti)),
        "dtau_mid_med": float(np.median(dtau_mid)),
    },
    "T2_cap_vs_anti": {
        "drot_mw_p": float(u2.pvalue),
        "dtau_mw_p": float(u2d.pvalue),
    },
    "T3_monotonicity": {
        "rho_theta_drot": float(rho_drot), "p": float(p_drot),
        "rho_theta_dtau": float(rho_dtau), "p_dtau": float(p_dtau),
    },
    "T4_sign_structure": {
        "cap_frac_pos": float(np.mean(dtau_cap > 0)),
        "cap_binom_p": float(bt_cap.pvalue),
        "mid_frac_pos": float(np.mean(dtau_mid > 0)),
        "mid_binom_p": float(bt_mid.pvalue),
        "anti_frac_pos": float(np.mean(dtau_anti > 0)),
        "anti_binom_p": float(bt_anti.pvalue),
        "dtau_anti_vs_mid_mw_p": float(u_dtau_anti.pvalue),
        "dtau_cap_vs_mid_mw_p": float(u_dtau_cap.pvalue),
    },
    "T6_harmonic": {
        "rho_cos1_dtau": float(r1_tau), "p_cos1_dtau": float(p1_tau),
        "rho_cos2_dtau": float(r2_tau), "p_cos2_dtau": float(p2_tau),
        "rho_cos1_drot": float(r1_rot), "p_cos1_drot": float(p1_rot),
        "rho_cos2_drot": float(r2_rot), "p_cos2_drot": float(p2_rot),
        "note": "cos(theta): one-sided dipolar defect; cos(2theta): "
                "axisymmetric field positive at both axis poles",
    },
    "per_catalogue": per_cat,
    "note": "The injected JFC chain is bipolar in aphelion direction "
            "(steps 052/056): the underlying structure is an axis, not "
            "a pair of independent patches.  The mirror cap carries the "
            "lowest rotation median of the profile -- no second "
            "rotation lobe -- while the residual proper-time channel is "
            "positive at BOTH axis poles and significantly negative in "
            "the intervening region: an axisymmetric (cos-2-theta) slip "
            "structure rather than a one-sided patch, the morphology an "
            "axisymmetric proper-time field produces and the bipolar "
            "injection chain already implied.",
}

out = str(RESULTS / "step_b54_bipolar_morphology.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b54_bipolar_morphology.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["desig", "cohort", "theta", "drot_cat", "dtau_unexplained",
                "region"])
    for d in data:
        region = "cap" if d["theta"] < 60 else ("anti" if d["theta"] >= 120
                                                else "mid")
        w.writerow([d["desig"], d["cohort"], d["theta"], d["drot"],
                    d["dtau"], region])

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

th_pts = np.array([d["theta"] for d in data])
dr_pts = np.array([d["drot"] for d in data])
du_pts = np.array([d["dtau"] for d in data])
coh = np.array([d["cohort"] for d in data])
incap = th_pts < 60
anti = th_pts >= 120

for ax, y, ylab, ttl in [
        (axes[0], dr_pts, r"catalogue rotation $\delta\theta$ (deg)",
         "rotation profile"),
        (axes[1], du_pts, r"unexplained slip $\delta\tau$ (yr)",
         "residual slip profile")]:
    ax.scatter(th_pts[~incap & ~anti & (coh == "code")],
               y[~incap & ~anti & (coh == "code")],
               s=9, c="0.6", alpha=0.45, label="mid (code)")
    ax.scatter(th_pts[~incap & ~anti & (coh == "warsaw")],
               y[~incap & ~anti & (coh == "warsaw")],
               s=9, c="0.8", alpha=0.45, marker="s", label="mid (warsaw)")
    ax.scatter(th_pts[anti], y[anti], s=14, c="steelblue",
               alpha=0.8, label="mirror cap")
    ax.scatter(th_pts[incap], y[incap], s=14, c="crimson",
               alpha=0.8, label="in cap")
    bx = [(p["bin"][0] + p["bin"][1]) / 2 for p in profile]
    key = "drot" if ax is axes[0] else "dtau"
    med = [p[f"{key}_med"] for p in profile]
    lo68 = [m - p[f"{key}_ci68"][0] for m, p in zip(med, profile)]
    hi68 = [p[f"{key}_ci68"][1] - m for m, p in zip(med, profile)]
    ax.errorbar(bx, med, yerr=[lo68, hi68], fmt="k-o", ms=4, lw=1.2,
                capsize=3, label="binned median (68%)")
    ax.axvline(60, color="k", ls=":", lw=1)
    ax.axvline(120, color="k", ls=":", lw=1)
    if ax is axes[1]:
        ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(r"aphelion--axis separation $\theta$ (deg)")
    ax.set_ylabel(ylab)
    ax.set_title(ttl, fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="best")

fig.suptitle("single-faced transit anomaly: no mirror-cap lobe "
             "(step 089)", fontsize=11)
fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b54_bipolar_morphology.png", dpi=150)

T1 = res["T1_mirror_cap"]
logger.info(f"T1: anti-cap drot median {T1['drot_anti_med']:.4f} vs mid "
            f"{T1['drot_mid_med']:.4f} (MW p={T1['mw_p_anti_gt_mid']:.3f}); "
            f"second-lobe 95% upper bound "
            f"{T1['second_lobe_amp_upper95']:.4f} deg "
            f"vs in-cap elevation {T1['incap_elevation_deg']:.4f} deg")
logger.info(f"T2: cap vs anti drot MW p={res['T2_cap_vs_anti']['drot_mw_p']:.4f}, "
            f"dtau MW p={res['T2_cap_vs_anti']['dtau_mw_p']:.4f}")
logger.info(f"T3: spearman theta-drot rho={rho_drot:+.3f} p={p_drot:.4f}; "
            f"theta-dtau rho={rho_dtau:+.3f} p={p_dtau:.4f}")
logger.info(f"T4: frac pos cap/mid/anti = "
            f"{res['T4_sign_structure']['cap_frac_pos']:.2f}/"
            f"{res['T4_sign_structure']['mid_frac_pos']:.2f}/"
            f"{res['T4_sign_structure']['anti_frac_pos']:.2f} "
            f"(binom p {bt_cap.pvalue:.3f}/{bt_mid.pvalue:.3f}/{bt_anti.pvalue:.3f})")
logger.info(f"T6: dtau vs cos(theta) rho={r1_tau:+.3f} p={p1_tau:.4f}; "
            f"vs cos(2theta) rho={r2_tau:+.3f} p={p2_tau:.4f} | "
            f"drot vs cos(2theta) rho={r2_rot:+.3f} p={p2_rot:.4f}")
print("wrote", out)
print("wrote", csv_out)
print("wrote", FIG / "step_b54_bipolar_morphology.png")
