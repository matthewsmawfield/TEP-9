"""step_068: Wall vs distributed field -- the structure of the implied
proper-time slip.

Step_065 expressed the transit anomaly in time units (dtau).  Two
physically distinct field geometries produce the same cap-localized
excess and must be separated on the measured scaling:

  * a localized wall or thin boundary, crossed once: the phase slip is
    applied at the crossing instant, so dtau is independent of the
    time the comet spends inside the 250 AU sphere (log-log slope
    b ~= 0 of dtau vs t_transit);

  * a distributed lapse gradient filling the interior: the slip
    accumulates at a steady rate, so dtau ~ t_transit (b ~= 1).

Both readings are TEP-compatible (Section 6.8 keeps both open); they
are separated here by the measured exponent, by the variance reduction
of the fractional slip, and by the inbound/outbound localization of
the Warsaw channel.

Additional coherence checks reported:

  * signed omega change dw = w_fut - w_orig for the CODE cohort --
    whether the in-cap rotation has a preferred handedness;
  * within-Warsaw correlation of the inbound slip with the full
    orig->fut slip, in cap vs out -- same-crossing-event coupling;
  * in-cap flatness of dtau vs axis distance (wall plateau vs the
    1/b^2 profile a point mass must write).

Inputs
------
results/step_b30_proper_time_slip.csv   (per-comet dtau, step_065)
data/raw/code/code_original.html, code_future.html  (signed omega)

Outputs
-------
results/step_b33_wall_vs_field.json
results/figures/step_b33_wall_vs_field.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, parse_code, tee_stdout
logger = StepLogger("step_068_wall_vs_field")
tee_stdout(logger)
logger.header("Wall vs distributed field -- slip structure")

import csv
import json
import math
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr, binomtest

CAP   = 60.0
SEED  = 20260918
NBOOT = 20000
rng   = np.random.default_rng(SEED)

rows = list(csv.DictReader(open(RESULTS / "step_b30_proper_time_slip.csv")))
for r in rows:
    for k in ("theta", "t_transit", "dtau_total", "dtau_unexplained",
              "frac_slip", "q", "i"):
        r[k] = float(r[k]) if r[k] not in ("", "nan") else float("nan")
    for k in ("dtau_inbound", "dtau_inbound_unexplained"):
        r[k] = float(r[k]) if r.get(k) not in ("", "nan", None) else float("nan")

logger.info(f"comets: n={len(rows)} "
            f"(code {sum(r['cohort']=='code' for r in rows)}, "
            f"warsaw {sum(r['cohort']=='warsaw' for r in rows)})")

# ------------------------------------------------------------------
# Signed omega change for the CODE cohort (catalogue legs)
# ------------------------------------------------------------------

orig = parse_code(DATA_RAW / "code" / "code_original.html")
fut  = parse_code(DATA_RAW / "code" / "code_future.html")
n_signed = 0
for r in rows:
    if r["cohort"] != "code":
        continue
    o, f = orig.get(r["desig"]), fut.get(r["desig"])
    if o and f:
        dw = (f["w"] - o["w"] + 540.0) % 360.0 - 180.0
        r["dw_signed"] = dw
        n_signed += 1
logger.info(f"signed dw recovered for {n_signed} CODE comets")

# ------------------------------------------------------------------
# Scaling test: does the unexplained slip grow with time inside?
# A wall applies the slip once -> residual uncorrelated with T.
# A distributed field accumulates eta*T -> positive correlation.
# ------------------------------------------------------------------

def slope_boot(x, y, n=NBOOT):
    """OLS slope of y~x with bootstrap 68% interval."""
    x, y = np.asarray(x), np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 6:
        return None
    b_ols = float(np.polyfit(x, y, 1)[0])
    idx = rng.integers(0, len(x), (n, len(x)))
    bs = np.array([np.polyfit(x[i], y[i], 1)[0] for i in idx])
    return {"b": b_ols,
            "ci68": [float(np.percentile(bs, 16)),
                     float(np.percentile(bs, 84))],
            "ci95": [float(np.percentile(bs, 2.5)),
                     float(np.percentile(bs, 97.5))],
            "n": int(len(x))}

def scaling(sub):
    """residual-slip scaling + raw-total context."""
    th = np.array([r["theta"] for r in sub])
    T  = np.array([r["t_transit"] for r in sub])
    dt = np.array([r["dtau_total"] for r in sub])
    ru = np.array([r["dtau_unexplained"] for r in sub])
    K  = np.abs(np.array([float(r["daa_sim"]) for r in sub]))
    D  = np.array([float(r["denc"]) for r in sub])
    Q  = np.array([r["q"] for r in sub]); I = np.array([r["i"] for r in sub])
    X  = np.column_stack([np.ones(len(sub)), np.log10(K + 1.0),
                          np.log10(D), Q, I])
    coef, *_ = np.linalg.lstsq(X, np.log(dt), rcond=None)
    with np.errstate(all="ignore"):
        resid_log = np.log(dt) - X @ coef
    out = {}
    for tag, mask in [("in", th < CAP), ("out", th >= CAP)]:
        o = {}
        if mask.sum() >= 6:
            rho, p = spearmanr(T[mask], ru[mask])
            o["resid_vs_T"] = {"rho": float(rho), "p_2sided": float(p)}
            rho, p = spearmanr(np.log(T[mask]), resid_log[mask])
            o["residlog_vs_logT"] = {"rho": float(rho), "p_2sided": float(p)}
            o["rate_yr_per_yr"] = slope_boot(T[mask], ru[mask])
            o["raw_loglog"] = slope_boot(np.log(T[mask]),
                                         np.log(dt[mask]))
            o["n"] = int(mask.sum())
        out[tag] = o if o else None
    return out

def var_reduction(sub):
    """log-space scatter of dtau vs frac_slip inside the cap."""
    s = [r for r in sub if r["theta"] < CAP]
    a = np.log([r["dtau_total"] for r in s])
    b = np.log([r["frac_slip"] for r in s])
    return {"sig_log_dtau": float(np.std(a)),
            "sig_log_frac": float(np.std(b)),
            "n": len(s)}

# ------------------------------------------------------------------
# Localization: inbound share of the slip (Warsaw)
# ------------------------------------------------------------------

def inbound_share(sub):
    s = [r for r in sub if np.isfinite(r["dtau_inbound"])]
    if len(s) < 10:
        return None
    th = np.array([r["theta"] for r in s])
    share = np.array([r["dtau_inbound"] / r["dtau_total"] for r in s])
    inc = th < CAP
    out = {"n": len(s), "n_in": int(inc.sum()),
           "med_share_in":  float(np.median(share[inc]))  if inc.any() else None,
           "med_share_out": float(np.median(share[~inc])) if (~inc).any() else None}
    if inc.any() and (~inc).any():
        u = mannwhitneyu(share[inc], share[~inc], alternative="greater")
        out["p_share_in_gt_out"] = float(u.pvalue)
        rho, p = spearmanr(np.array([r["dtau_inbound_unexplained"] for r in s]),
                           np.array([r["dtau_unexplained"] for r in s]))
        out["slip_coupling"] = {"rho": float(rho), "p_2sided": float(p)}
        rho, p = spearmanr(th, share)
        out["share_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
    return out

# ------------------------------------------------------------------
# Coherence checks
# ------------------------------------------------------------------

def coherence(sub):
    out = {}
    th = np.array([r["theta"] for r in sub])
    inc = th < CAP
    # in-cap flatness of the slip (wall plateau vs 1/b^2)
    dt = np.array([r["dtau_total"] for r in sub])
    if inc.sum() >= 6:
        rho, p = spearmanr(th[inc], dt[inc])
        out["dtau_vs_theta_incap"] = {"rho": float(rho),
                                      "p_2sided": float(p)}
    # signed omega handedness (CODE only)
    s = [r for r in sub if "dw_signed" in r]
    if len(s) >= 10:
        ths = np.array([r["theta"] for r in s])
        dw  = np.array([r["dw_signed"] for r in s])
        incs = ths < CAP
        for tag, m in [("in", incs), ("out", ~incs)]:
            if m.sum() >= 5:
                pos = int((dw[m] > 0).sum())
                out[f"dw_sign_{tag}"] = {
                    "n": int(m.sum()), "n_pos": pos,
                    "frac_pos": pos / int(m.sum()),
                    "p_binom": float(binomtest(pos, int(m.sum()), 0.5).pvalue)}
        if incs.any() and (~incs).any():
            u = mannwhitneyu(np.abs(dw[incs]), np.abs(dw[~incs]),
                             alternative="greater")
            out["abs_dw_in_gt_out_p"] = float(u.pvalue)
    return out

code   = [r for r in rows if r["cohort"] == "code"]
warsaw = [r for r in rows if r["cohort"] == "warsaw"]
pooled = rows
m_code   = [r for r in code   if r["q"] < 3.1]
m_warsaw = [r for r in warsaw if r["q"] < 3.1]
m_pooled = [r for r in pooled if r["q"] < 3.1]

res = {
    "method": "log-log slope b of dtau vs transit time separates a "
              "single-crossing wall (b~0) from a distributed lapse "
              "gradient (b~1); localization and coherence checks on "
              "the same per-comet slips",
    "cap_deg": CAP, "seed": SEED, "nboot": NBOOT,
    "t_transit_range_yr": [
        float(min(r["t_transit"] for r in rows)),
        float(max(r["t_transit"] for r in rows))],
    "scaling": {
        "code_all":       scaling(code),
        "code_matched":   scaling(m_code),
        "warsaw_all":     scaling(warsaw),
        "pooled_all":     scaling(pooled),
        "pooled_matched": scaling(m_pooled),
    },
    "variance_reduction": {
        "code_incap":   var_reduction(code),
        "pooled_incap": var_reduction(pooled),
        "pooled_matched_incap": var_reduction(m_pooled),
    },
    "warsaw_inbound_share": inbound_share(warsaw),
    "coherence": {
        "code":   coherence(code),
        "warsaw": coherence(warsaw),
        "pooled": coherence(pooled),
    },
}

out = str(RESULTS / "step_b33_wall_vs_field.json")
json.dump(res, open(out, "w"), indent=1, default=float)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

ax = axes[0]
for co, mk in [("code", "o"), ("warsaw", "s")]:
    s = [r for r in rows if r["cohort"] == co]
    th = np.array([r["theta"] for r in s]); inc = th < CAP
    tt = np.array([r["t_transit"] for r in s])
    dt = np.array([r["dtau_total"] for r in s])
    ax.scatter(tt[~inc], dt[~inc], s=20, marker=mk, facecolors="none",
               edgecolors="0.55", label=f"{co} outside")
    ax.scatter(tt[inc], dt[inc], s=24, marker=mk, c="crimson",
               label=f"{co} inside")
sc = res["scaling"]["pooled_matched"]["in"]
if sc and sc.get("raw_loglog"):
    xs = np.linspace(math.log(500), math.log(760), 50)
    inc_m = np.array([r["theta"] for r in m_pooled]) < CAP
    b0 = sc["raw_loglog"]["b"]
    b1 = np.mean(np.log(np.array([r["dtau_total"] for r in m_pooled])[inc_m])) \
         - b0 * np.mean(np.log(np.array([r["t_transit"] for r in m_pooled])[inc_m]))
    ax.plot(np.exp(xs), np.exp(b1 + b0 * xs), "r--", lw=1,
            label=f"in-cap slope b={b0:.2f}")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"transit time inside 250 AU $T$ (yr)")
ax.set_ylabel(r"$\delta\tau_{\rm total}$ (yr)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("slip vs transit time -- wall ($b{=}0$) vs field ($b{=}1$)",
             fontsize=10)

ax = axes[1]
s = [r for r in warsaw if np.isfinite(r["dtau_inbound"])]
th = np.array([r["theta"] for r in s]); inc = th < CAP
share = np.array([r["dtau_inbound"] / r["dtau_total"] for r in s])
bins = np.linspace(0, 1.5, 25)
ax.hist(share[~inc], bins=bins, color="0.55", alpha=0.7, density=True,
        label=f"outside ($n={int((~inc).sum())}$)")
ax.hist(share[inc], bins=bins, color="crimson", alpha=0.6, density=True,
        label=f"inside ($n={int(inc.sum())}$)")
ax.axvline(1.0, color="k", ls=":", lw=1)
ax.set_xlabel(r"inbound share of slip $\delta\tau_{\rm in}/\delta\tau_{\rm total}$")
ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Warsaw: crossing localization", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b33_wall_vs_field.png", dpi=150)

for tag, d in res["scaling"].items():
    if d.get("in"):
        i, o = d["in"], d.get("out") or {}
        logger.info(f"{tag}: in-cap resid-vs-T rho="
                    f"{i['resid_vs_T']['rho']:.3f} "
                    f"p={i['resid_vs_T']['p_2sided']:.3f} | "
                    f"out rho={o.get('resid_vs_T',{}).get('rho',float('nan')):.3f} "
                    f"p={o.get('resid_vs_T',{}).get('p_2sided',float('nan')):.3f}")
ws = res["warsaw_inbound_share"]
if ws:
    logger.info(f"warsaw inbound share in/out = "
                f"{ws['med_share_in']:.2f}/{ws['med_share_out']:.2f} "
                f"p={ws['p_share_in_gt_out']:.4f} | slip coupling "
                f"rho={ws['slip_coupling']['rho']:.3f} "
                f"p={ws['slip_coupling']['p_2sided']:.4f}")
print("wrote", out)
print("wrote", FIG / "step_b33_wall_vs_field.png")
