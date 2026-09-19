"""step_065: Operational time-equivalent of the transit rotation.

Every Phase-3 step expresses the comet anomaly in angular units: the
orig->fut periapsis-direction rotation d_of and its residual after the
measured encounter budget is regressed out (step_063).  Under the
domain-boundary reading that residual is a phase error -- a proper-time
offset accumulated during the boundary crossing that a standard-clock
reconstruction mis-assigns as an in-plane rotation.  This step converts
the measured rotation into an operational time-equivalent; it is not a
derived proper-time measurement until a field/clock response model is
specified.

Definition
----------
At the 250 AU barycentric reference sphere a periapsis-direction
rotation d_th (radians) corresponds to the comet's sky position having
been swept through the same angle.  The comet sweeps sky angle at the
rate Omega_b = h / r_b^2, where h = sqrt(mu q (1+e)) is the specific
angular momentum of the boundary orbit and r_b = 250 AU.  The operational
time-equivalent is therefore

    dtau = d_th * r_b^2 / h        [yr]

-- the time an unperturbed comet on the same boundary orbit needs to
subtend the unexplained rotation.  No lapse-field amplitude is assumed.

Two quantities are reported per comet:
  dtau_total   : the full orig->fut rotation expressed in years
  dtau_unexplained : the residual of dtau after regression on the
                     measured encounter budget (energy kick, minimum
                     planetary approach, perihelion depth, inclination)

The fractional slip dtau/(t_back + t_fwd) expresses the offset relative
to the comet's measured transit duration inside the 250 AU sphere -- a
dimensionless lapse contrast.

Inputs
------
results/step_b28_bidirectional_rotation.csv   (CODE cohort, step_063)
results/step_b29_warsaw_bidirectional.csv     (Warsaw cohort, step_064)

Outputs
-------
results/step_b30_proper_time_slip.json
results/step_b30_proper_time_slip.csv
results/figures/step_b30_proper_time_slip.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_065_proper_time_slip")
tee_stdout(logger)
logger.header("Operational time-equivalent (transit rotation in time units)")

import csv
import json
import math
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

MU    = 4 * math.pi ** 2          # AU^3 yr^-2
R_B   = 250.0                     # barycentric reference sphere (AU)
CAP   = 60.0                      # pre-declared cap (deg)
SEED  = 20260918
NBOOT = 20000

rng = np.random.default_rng(SEED)

# ------------------------------------------------------------------
# Load the bidirectional cohorts produced by steps 063/064
# ------------------------------------------------------------------

def load_csv(path, cohort):
    rows = []
    for r in csv.DictReader(open(path)):
        try:
            rec = dict(
                desig=r["desig"], cohort=cohort,
                q=float(r["q"]), i=float(r["i"]),
                theta=float(r["theta"]),
                drot=float(r["drot_sim"]),
                drot_cat=float(r["drot_cat"]),
                daa_sim=float(r["daa_sim"]),
                denc=float(r["denc"]),
                aa_back=float(r["aa_back"]),
                t_back=abs(float(r["t_back"])),
                t_fwd=abs(float(r["t_fwd"])))
            # inbound-leg channel where present (Warsaw, step_064)
            for k in ("ddir_sim", "ddir_cat"):
                if k in r and r[k] not in ("", "nan"):
                    rec[k] = float(r[k])
            rows.append(rec)
        except (ValueError, KeyError):
            continue
    return rows

code   = load_csv(RESULTS / "step_b28_bidirectional_rotation.csv", "code")
warsaw = load_csv(RESULTS / "step_b29_warsaw_bidirectional.csv", "warsaw")
logger.info(f"cohorts: CODE n={len(code)}, Warsaw n={len(warsaw)}")

# ------------------------------------------------------------------
# Per-comet proper-time conversion
# ------------------------------------------------------------------

def enrich(rows):
    for r in rows:
        a = 1e6 / r["aa_back"] if r["aa_back"] != 0 else float("inf")
        e = 1.0 - r["q"] / a
        h = math.sqrt(MU * r["q"] * (1.0 + e))          # AU^2/yr
        om_b = h / R_B ** 2                              # rad/yr at boundary
        r["h"] = h
        r["om_b_deg_yr"] = math.degrees(om_b)
        r["dtau_total"] = math.radians(r["drot"]) / om_b
        r["dtau_cat"]   = math.radians(r["drot_cat"]) / om_b
        if "ddir_sim" in r:
            r["dtau_inbound"] = math.radians(r["ddir_sim"]) / om_b
            r["dtau_inbound_cat"] = math.radians(r["ddir_cat"]) / om_b
        r["t_transit"]  = r["t_back"] + r["t_fwd"]
        r["frac_slip"]  = r["dtau_total"] / r["t_transit"]
    return rows

enrich(code); enrich(warsaw)
allrows = code + warsaw

# ------------------------------------------------------------------
# Statistics
# ------------------------------------------------------------------

def lin_resid(y, X_cols):
    X = np.column_stack([np.ones(len(y))] + X_cols)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    with np.errstate(all="ignore"):
        return y - X @ coef

def boot_med_ci(x, n=NBOOT):
    idx = rng.integers(0, len(x), (n, len(x)))
    meds = np.median(x[idx], axis=1)
    return float(np.percentile(meds, 16)), float(np.percentile(meds, 84))

def run(sub):
    th  = np.array([r["theta"]     for r in sub])
    dt  = np.array([r["dtau_total"] for r in sub])
    K   = np.abs(np.array([r["daa_sim"] for r in sub]))
    D   = np.array([r["denc"] for r in sub])
    Q   = np.array([r["q"]    for r in sub])
    I   = np.array([r["i"]    for r in sub])
    Tt  = np.array([r["t_transit"] for r in sub])
    inc = th < CAP
    out = {"n": len(sub), "n_in": int(inc.sum()),
           "med_dtau_total_in":  float(np.median(dt[inc])) if inc.any() else None,
           "med_dtau_total_out": float(np.median(dt[~inc])) if (~inc).any() else None,
           "med_transit_yr": float(np.median(Tt))}
    # raw cap contrast in time units
    if inc.any() and (~inc).any():
        u = mannwhitneyu(dt[inc], dt[~inc], alternative="greater")
        out["dtau_in_gt_out_p"] = float(u.pvalue)
        rho, p = spearmanr(th, dt)
        out["dtau_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
        # residual in log space (mirrors step_063 instrument)
        resid_log = lin_resid(np.log(dt),
            [np.log10(K + 1.0), np.log10(D), Q, I])
        u = mannwhitneyu(resid_log[inc], resid_log[~inc], alternative="greater")
        out["dtau_resid_log_cap"] = {"p": float(u.pvalue),
            "med_resid_in": float(np.median(resid_log[inc])),
            "med_resid_out": float(np.median(resid_log[~inc]))}
        rho, p = spearmanr(th, resid_log)
        out["dtau_resid_log_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
        # residual in linear space -> the implied offset in years
        resid_lin = lin_resid(dt,
            [np.log10(K + 1.0), np.log10(D), Q, I])
        u = mannwhitneyu(resid_lin[inc], resid_lin[~inc], alternative="greater")
        lo, hi = boot_med_ci(resid_lin[inc])
        out["dtau_unexplained"] = {
            "med_in":  float(np.median(resid_lin[inc])),
            "ci68_in": [lo, hi],
            "med_out": float(np.median(resid_lin[~inc])),
            "p_in_gt_out": float(u.pvalue)}
        rho, p = spearmanr(th, resid_lin)
        out["dtau_unexplained_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
        for r, rr in zip(sub, resid_lin):
            r["dtau_unexplained"] = float(rr)
        # fractional lapse contrast
        fs = dt / Tt
        resid_fs = lin_resid(np.log(fs),
            [np.log10(K + 1.0), np.log10(D), Q, I])
        u = mannwhitneyu(resid_fs[inc], resid_fs[~inc], alternative="greater")
        out["frac_slip"] = {"med_in": float(np.median(fs[inc])),
            "med_out": float(np.median(fs[~inc])),
            "resid_cap_p": float(u.pvalue)}
        # inbound-leg channel (Warsaw cohort)
        if all("dtau_inbound" in r for r in sub):
            di = np.array([r["dtau_inbound"] for r in sub])
            resid_in = lin_resid(di,
                [np.log10(K + 1.0), np.log10(D), Q, I])
            u = mannwhitneyu(resid_in[inc], resid_in[~inc], alternative="greater")
            lo, hi = boot_med_ci(resid_in[inc])
            out["dtau_inbound_unexplained"] = {
                "med_in":  float(np.median(resid_in[inc])),
                "ci68_in": [lo, hi],
                "med_out": float(np.median(resid_in[~inc])),
                "p_in_gt_out": float(u.pvalue),
                "med_total_in": float(np.median(di[inc])),
                "med_total_out": float(np.median(di[~inc]))}
            for r, rr in zip(sub, resid_in):
                r["dtau_inbound_unexplained"] = float(rr)
    return out

res = {
    "method": "operational time-equivalent dt_equiv = d_th * r_b^2 / h; "
              "boundary rotation converted to the time an unperturbed "
              "comet on the same orbit needs to sweep the angle.  "
              "Residuals regressed on the measured encounter budget "
              "(|kick|, min planet approach, q, i).",
    "definition": "dtau is the time-equivalent that a standard-clock "
                  "reconstruction would associate with the measured phase "
                  "error; it is not a derived proper-time measurement "
                  "without an explicit field/clock response model; no "
                  "model amplitude assumed",
    "r_b_AU": R_B, "cap_deg": CAP, "seed": SEED, "nboot": NBOOT,
    "code":   {"n": len(code),
               "matched": run([r for r in code if r["q"] < 3.1]),
               "all":     run(code)},
    "warsaw": {"n": len(warsaw),
               "matched": run([r for r in warsaw if r["q"] < 3.1]),
               "all":     run(warsaw)},
    "pooled": {"n": len(allrows),
               "matched": run([r for r in allrows if r["q"] < 3.1]),
               "all":     run(allrows)},
}

out = str(RESULTS / "step_b30_proper_time_slip.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b30_proper_time_slip.csv")
fields = ["desig", "cohort", "q", "i", "theta", "drot", "drot_cat",
          "ddir_sim", "ddir_cat",
          "daa_sim", "denc", "aa_back", "t_back", "t_fwd", "t_transit",
          "h", "om_b_deg_yr", "dtau_total", "dtau_cat", "dtau_inbound",
          "dtau_inbound_cat", "frac_slip",
          "dtau_unexplained", "dtau_inbound_unexplained"]
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    w.writeheader(); w.writerows(allrows)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

ax = axes[0]
for cohort, mk in [("code", "o"), ("warsaw", "s")]:
    sub = [r for r in allrows if r["cohort"] == cohort]
    th = np.array([r["theta"] for r in sub])
    dt = np.array([r["dtau_total"] for r in sub])
    inc = th < CAP
    ax.scatter(th[~inc], dt[~inc], s=20, marker=mk, facecolors="none",
               edgecolors="0.55", label=f"{cohort} outside ($n={int((~inc).sum())}$)")
    ax.scatter(th[inc], dt[inc], s=24, marker=mk, c="crimson",
               label=f"{cohort} inside ($n={int(inc.sum())}$)")
ax.axvline(CAP, color="k", ls=":", lw=1)
ax.set_xlabel(r"aphelion--axis separation $\theta$ (deg)")
ax.set_ylabel(r"$\delta\tau_{\rm total}$ (yr)")
ax.set_yscale("log")
ax.legend(frameon=False, fontsize=8)
ax.set_title("implied proper-time offset vs axis distance", fontsize=10)

ax = axes[1]
sub = [r for r in allrows if "dtau_unexplained" in r]
resid = np.array([r["dtau_unexplained"] for r in sub])
th = np.array([r["theta"] for r in sub])
inc = th < CAP
bins = np.linspace(np.percentile(resid, 1), np.percentile(resid, 99), 30)
ax.hist(resid[~inc], bins=bins, color="0.55", alpha=0.7, density=True,
        label=f"outside ($n={int((~inc).sum())}$)")
ax.hist(resid[inc], bins=bins, color="crimson", alpha=0.6, density=True,
        label=f"inside ($n={int(inc.sum())}$)")
ax.axvline(0, color="k", ls=":", lw=1)
ax.set_xlabel(r"unexplained proper-time offset $\delta\tau$ (yr)")
ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=8)
ax.set_title("residual after encounter-budget regression", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b30_proper_time_slip.png", dpi=150)

for tag, dd in [("code", res["code"]), ("warsaw", res["warsaw"]), ("pooled", res["pooled"])]:
    d = dd["all"]
    logger.info(f"{tag} (all): n={d['n']} in={d['n_in']} | "
                f"dtau_total in/out = {d['med_dtau_total_in']:.2f}/{d['med_dtau_total_out']:.2f} yr | "
                f"unexplained in = {d['dtau_unexplained']['med_in']:.2f} yr "
                f"ci68 {d['dtau_unexplained']['ci68_in']} p={d['dtau_unexplained']['p_in_gt_out']:.4f}")
    if "dtau_inbound_unexplained" in d:
        di = d["dtau_inbound_unexplained"]
        logger.info(f"{tag} inbound: med total in/out = {di['med_total_in']:.2f}/"
                    f"{di['med_total_out']:.2f} yr | unexplained in = {di['med_in']:.2f} yr "
                    f"ci68 {di['ci68_in']} p={di['p_in_gt_out']:.4f}")
print("wrote", out)
print("wrote", csv_out)
print("wrote", FIG / "step_b30_proper_time_slip.png")
