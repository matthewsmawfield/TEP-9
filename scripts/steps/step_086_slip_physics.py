#!/usr/bin/env python3
"""Step 086 -- slip-physics discriminator: where in the orbit does
the anomaly act?

The inner-shell profile (step 076) shows the in-cap rotation excess
is already fully developed at the innermost probed shell of 8 AU --
inside the region where the comet is observed and the reconstruction
residual is generated.  That leaves a discriminating question the
radial profile alone cannot answer: is the slip a temporal offset
imprinted during the observed arc, a temporal offset imprinted at a
distant boundary crossing, or a fixed angular displacement?

For near-parabolic comets the three readings predict different
power-law scalings of the measured rotation with perihelion
distance q (h ~ q^{1/2}):

  angular systematic   d_theta ~ q^0      (a fixed apsidal error:
                                          zonal catalogue error or
                                          directional fit bias)
  boundary time slip   d_theta = dtau*h/r_b^2 ~ q^{+0.5}
                       (fixed proper-time offset at a ~250 AU
                                          crossing; more angular
                                          momentum h sweeps more
                                          angle per unit time slip)
  arc-epoch time slip  d_theta = omega_dot_arc * dtau ~ q^{-1.5}
                       (fixed time offset absorbed by the fit at
                                          the observed-arc epoch;
                                          the perihelion sweep rate
                                          dominates and scales as
                                          q^{-3/2})

Equivalently, on the boundary-normalized time slip
dtau = d_theta * r_b^2 / h the predictions are
  angular    slope(dtau ~ q) = -0.5
  boundary   slope            =  0
  arc-epoch  slope            = -2.0

One dataset, one slope measurement on each of two normalizations,
three separated model points: (0, -0.5), (+0.5, 0), (-1.5, -2).
The out-of-cap cohort supplies the control -- a sky- or
catalogue-wide systematic would carry its scaling into both caps;
a boundary mechanism acts only on the crossing population.

Data: per-comet catalogue rotations drot_cat and boundary-
normalized slips dtau_cat from step_b30_proper_time_slip.csv
(CODE n=131, Warsaw n=98; theta vs the transit axis).

Outputs
-------
results/step_b51_slip_physics.json
results/step_b51_slip_physics.csv   (per-comet inputs used)
results/figures/step_b51_slip_physics.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_086_slip_physics")
tee_stdout(logger)
logger.header("Slip-physics discriminator: angular vs boundary vs arc-epoch")

import csv
import json
import numpy as np
from scipy.stats import theilslopes, spearmanr

rng = np.random.default_rng(20260918)
N_BOOT = 20000
CAP = 60.0

MODELS = {
    "angular":   {"s_theta": 0.0,  "s_tau": -0.5,
                  "note": "fixed angular displacement (zonal/fit "
                          "systematic): d_theta ~ q^0"},
    "boundary":  {"s_theta": 0.5,  "s_tau": 0.0,
                  "note": "fixed proper-time slip at the ~250 AU "
                          "boundary: d_theta = dtau*h/r^2 ~ q^+0.5"},
    "arc_epoch": {"s_theta": -1.5, "s_tau": -2.0,
                  "note": "fixed time slip absorbed at the observed-"
                          "arc epoch: d_theta ~ perihelion sweep "
                          "rate ~ q^-1.5"},
}

# ------------------------------------------------------------------
# Load per-comet products
# ------------------------------------------------------------------

rows = []
with open(RESULTS / "step_b30_proper_time_slip.csv") as f:
    for r in csv.DictReader(f):
        try:
            q = float(r["q"]); th = float(r["theta"])
            dcat = float(r["drot_cat"]); dsim = float(r["drot"])
            tcat = float(r["dtau_cat"])
            dun = float(r["dtau_unexplained"])
        except (ValueError, KeyError):
            continue
        rows.append(dict(desig=r["desig"], cohort=r["cohort"], q=q,
                         theta=th, drot_cat=dcat, drot_sim=dsim,
                         dtau_cat=tcat, dtau_unexplained=dun,
                         dres=dcat - dsim))

logger.info(f"loaded {len(rows)} comets "
            f"(code {sum(r['cohort']=='code' for r in rows)}, "
            f"warsaw {sum(r['cohort']=='warsaw' for r in rows)})")

# ------------------------------------------------------------------
# Slope machinery
# ------------------------------------------------------------------

def boot_slope(x, y, n_boot=N_BOOT):
    """Theil-Sen slope of y~x with a bootstrap 68% interval."""
    x, y = np.asarray(x), np.asarray(y)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 6:
        return float("nan"), float("nan"), float("nan"), len(x)
    s0 = theilslopes(y, x)[0]
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    sx = x[idx]; sy = y[idx]
    # vectorized Theil-Sen is heavy; resample and use polyfit for CI
    sl = np.array([np.polyfit(sx[j], sy[j], 1)[0] for j in range(n_boot)])
    lo, hi = np.percentile(sl, [16, 84])
    return float(s0), float(lo), float(hi), int(len(x))

def slope_block(sub, xkey, ykey):
    lx = np.log10([r[xkey] for r in sub])
    ly = np.log10([r[ykey] for r in sub])
    s, lo, hi, nn = boot_slope(lx, ly)
    rho, p = spearmanr(lx, ly) if nn >= 6 else (float("nan"), float("nan"))
    return {"n": nn, "slope": s, "ci68": [lo, hi],
            "spearman_rho": float(rho), "p_2sided": float(p)}

def analyse(sub):
    out = {}
    out["s_theta"] = slope_block(sub, "q", "drot_cat")
    out["s_tau"] = slope_block(sub, "q", "dtau_cat")
    # residual channels: the planetary baseline carries its own
    # q-dependence (~q^-0.6, identical in both caps -- deeper plungers
    # rotate more under encounters), so the mechanism test must run on
    # the unexplained residual, not the raw rotation
    pos = [r for r in sub if r["dres"] > 0]
    out["s_theta_resid"] = slope_block(pos, "q", "dres")
    out["s_theta_resid"]["n_pos_frac"] = len(pos) / len(sub) if sub else float("nan")
    pos_t = [r for r in sub if r["dtau_unexplained"] > 0]
    out["s_tau_resid"] = slope_block(pos_t, "q", "dtau_unexplained")
    out["s_tau_resid"]["n_pos_frac"] = len(pos_t) / len(sub) if sub else float("nan")
    # signed Spearman on the residual slip -- no log, keeps negatives
    if len(sub) >= 6:
        rho, p = spearmanr([r["q"] for r in sub],
                           [r["dtau_unexplained"] for r in sub])
        out["spearman_q_vs_dtau_resid"] = {"rho": float(rho), "p_2sided": float(p)}
    # distance of the measured (s_theta, s_tau) point from each model,
    # on both the raw and residual channels
    for tag, (sk1, sk2) in {"raw": ("s_theta", "s_tau"),
                            "resid": ("s_theta_resid", "s_tau_resid")}.items():
        st, sa = out[sk1]["slope"], out[sk2]["slope"]
        for m, mm in MODELS.items():
            out[f"dist_{tag}_{m}"] = float(
                np.hypot(st - mm["s_theta"], sa - mm["s_tau"])) \
                if np.isfinite(st) and np.isfinite(sa) else float("nan")
    return out

inc = [r for r in rows if r["theta"] < CAP]
outc = [r for r in rows if r["theta"] >= CAP]
res = {"method": "power-law scaling of the catalogue rotation "
       "drot_cat and of the boundary-normalized slip dtau_cat with "
       "perihelion distance q; the three slip mechanisms predict "
       "separated (s_theta, s_tau) model points.",
       "cap_deg": CAP,
       "models": MODELS,
       "cohorts": {}}

for name, sub in (("pooled_in", inc), ("pooled_out", outc),
                  ("code_in", [r for r in inc if r["cohort"] == "code"]),
                  ("code_out", [r for r in outc if r["cohort"] == "code"]),
                  ("warsaw_in", [r for r in inc if r["cohort"] == "warsaw"]),
                  ("warsaw_out", [r for r in outc if r["cohort"] == "warsaw"])):
    res["cohorts"][name] = analyse(sub)
    a = res["cohorts"][name]
    logger.info(f"{name}: n={a['s_theta']['n']} "
                f"s_theta={a['s_theta']['slope']:+.2f} "
                f"s_tau={a['s_tau']['slope']:+.2f} | resid "
                f"s_theta={a['s_theta_resid']['slope']:+.2f} "
                f"s_tau={a['s_tau_resid']['slope']:+.2f} "
                f"(pos {a['s_theta_resid'].get('n_pos_frac', float('nan')):.2f})")

# residual channel: median unexplained slip vs q in 3 q-bins (in-cap),
# a distribution-free companion to the slope fits
qbins = [(0, 1.0), (1.0, 2.0), (2.0, 99.0)]
res["unexplained_vs_q_incap"] = {}
for lo, hi in qbins:
    sub = [r for r in inc if lo <= r["q"] < hi]
    if len(sub) >= 3:
        med = float(np.median([r["dtau_unexplained"] for r in sub]))
        res["unexplained_vs_q_incap"][f"q_{lo}_{hi}"] = {
            "n": len(sub), "median_dtau_unexplained": med}

# ------------------------------------------------------------------
# Write
# ------------------------------------------------------------------

out = str(RESULTS / "step_b51_slip_physics.json")
json.dump(res, open(out, "w"), indent=1, default=float)
csv_out = str(RESULTS / "step_b51_slip_physics.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("wrote", out)
print("wrote", csv_out)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))

ax = axes[0]
for sub, c, lab in ((inc, "crimson", "in-cap"), (outc, "0.6", "out-of-cap")):
    ax.scatter([r["q"] for r in sub], [r["drot_cat"] for r in sub],
               s=14, c=c, alpha=0.6, label=lab)
ax.set_xscale("log"); ax.set_yscale("log")
a = res["cohorts"]["pooled_in"]
qline = np.logspace(np.log10(0.3), np.log10(4), 50)
med_q = np.median([r["q"] for r in inc]); med_d = np.median([r["drot_cat"] for r in inc])
ax.plot(qline, med_d * (qline / med_q) ** a["s_theta"]["slope"], "r--",
        lw=1.5, label=f"in-cap slope {a['s_theta']['slope']:+.2f}")
for m, sty in (("boundary", ":"), ("angular", "-."), ("arc_epoch", "-")):
    ax.plot(qline, med_d * (qline / med_q) ** MODELS[m]["s_theta"],
            "k", ls=sty, lw=0.9, alpha=0.6, label=f"{m} ({MODELS[m]['s_theta']:+.1f})")
ax.set_xlabel("perihelion distance q (AU)")
ax.set_ylabel("catalogue rotation $\\delta\\theta$ (deg)")
ax.legend(frameon=False, fontsize=7)
ax.set_title("angular channel vs q", fontsize=10)

ax = axes[1]
for sub, c, lab in ((inc, "crimson", "in-cap"), (outc, "0.6", "out-of-cap")):
    v = np.array([r["dtau_cat"] for r in sub])
    ok = v > 0
    ax.scatter(np.array([r["q"] for r in sub])[ok], v[ok],
               s=14, c=c, alpha=0.6, label=lab)
ax.set_xscale("log"); ax.set_yscale("log")
a = res["cohorts"]["pooled_in"]
med_t = np.median([r["dtau_cat"] for r in inc if r["dtau_cat"] > 0])
ax.plot(qline, med_t * (qline / med_q) ** a["s_tau"]["slope"], "r--",
        lw=1.5, label=f"in-cap slope {a['s_tau']['slope']:+.2f}")
for m, sty in (("boundary", ":"), ("angular", "-."), ("arc_epoch", "-")):
    ax.plot(qline, med_t * (qline / med_q) ** MODELS[m]["s_tau"],
            "k", ls=sty, lw=0.9, alpha=0.6, label=f"{m} ({MODELS[m]['s_tau']:+.1f})")
ax.set_xlabel("perihelion distance q (AU)")
ax.set_ylabel("boundary-normalized slip $\\delta\\tau$ (yr)")
ax.legend(frameon=False, fontsize=7)
ax.set_title("time channel vs q", fontsize=10)

ax = axes[2]
for m, (mk, mm) in enumerate(MODELS.items()):
    ax.scatter(mm["s_theta"], mm["s_tau"], marker="*", s=220, c="k",
               zorder=3, label=f"{mk} model" if m == 0 else None)
    ax.annotate(mk, (mm["s_theta"], mm["s_tau"]),
                xytext=(6, -12), textcoords="offset points", fontsize=9)
for name, c in (("pooled_in", "crimson"), ("code_in", "darkred"),
                ("warsaw_in", "teal"), ("pooled_out", "0.55")):
    a = res["cohorts"][name]
    st, sa = a["s_theta"]["slope"], a["s_tau"]["slope"]
    if not (np.isfinite(st) and np.isfinite(sa)):
        continue
    xe = [[max(st - a["s_theta"]["ci68"][0], 0.001)],
          [max(a["s_theta"]["ci68"][1] - st, 0.001)]]
    ye = [[max(sa - a["s_tau"]["ci68"][0], 0.001)],
          [max(a["s_tau"]["ci68"][1] - sa, 0.001)]]
    ax.errorbar(st, sa, xerr=xe, yerr=ye,
                fmt="o", ms=7, c=c, capsize=3, label=name)
    # residual-channel point, open marker
    rt, ra = a["s_theta_resid"]["slope"], a["s_tau_resid"]["slope"]
    if np.isfinite(rt) and np.isfinite(ra):
        ax.scatter(rt, ra, marker="o", s=60, facecolors="none",
                   edgecolors=c, linewidths=1.5)
ax.scatter([], [], marker="o", s=60, facecolors="none",
           edgecolors="k", label="residual channel")
ax.set_xlabel("$s_\\theta$: slope of $\\log\\delta\\theta$ vs $\\log q$")
ax.set_ylabel("$s_\\tau$: slope of $\\log\\delta\\tau$ vs $\\log q$")
ax.legend(frameon=False, fontsize=8)
ax.set_title("measured vs predicted scaling", fontsize=10)
ax.grid(alpha=0.2)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b51_slip_physics.png", dpi=150)
print(f"wrote {FIG / 'step_b51_slip_physics.png'}")
