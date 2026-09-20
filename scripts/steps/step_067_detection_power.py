"""step_067: Prospective detection power -- how many future detached-TNO
discoveries separate the domain-boundary prediction from the
footprint-only null?

step_054 issued the dated falsifiable predictions: the intrinsic-axis
model expects ~59 per cent of future detached discoveries inside the
pre-declared 60 deg cap (mean varpi 49 deg), against a footprint-only
baseline of ~42 per cent (mean 10.5 deg).  This step turns those rates
into a power statement by Monte-Carlo: for each candidate future sample
size N, draw realizations of the discovery count under both models and
measure how decisively a future observation separates them.

Channels
--------
C1 cap fraction   : in-cap count ~ Binomial(N, p) with p_wall = 0.59
                    (mixture prediction, step_b13/b19) vs p_fp = 0.42
                    (footprint-convolved baseline).  The discriminator is
                    the one-sided binomial test of the observed count
                    against p_fp; power = fraction of wall realizations
                    rejecting at alpha.
C2 mean direction : per-object varpi drawn from the measured mixture
                    (clustered fraction f = 0.43, von Mises
                    sigma = 34.4 deg about 49 deg + footprint remainder)
                    vs the footprint-only distribution centred on
                    10.5 deg.  Discriminator: circular mean resultant
                    length R; power estimated against the footprint
                    null's R distribution.
C3 comet channel  : fresh class-1-like comet cohort; in-cap fraction
                    0.352 of the sample (CODE-measured); rate of
                    rotation > 0.2 deg = 0.579 in-cap vs 0.229 out-cap.
                    Power of the in/out rate contrast vs cohort size.

Inputs (all measured in earlier steps, no new assumptions)
----------------------------------------------------------
results/step_b19_prospective_predictions.json   (rates, caps, baselines)
results/step_b13_patch_scale.json               (mixture f, sigma)
results/step_b15_evidence_summary.json          (context)

Outputs
-------
results/step_b32_detection_power.json
results/figures/supplementary/step_b32_detection_power.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_067_detection_power")
tee_stdout(logger)
logger.header("Prospective detection power (Monte-Carlo)")

import json
import math
import numpy as np
from scipy.stats import binomtest, norm, vonmises

SEED  = 20260918
NMC   = 40000
rng   = np.random.default_rng(SEED)

pred  = json.load(open(RESULTS / "step_b19_prospective_predictions.json"))
patch = json.load(open(RESULTS / "step_b13_patch_scale.json"))

P_WALL   = float(pred["P1_next_detached"]["predicted_frac_in_cap"])   # 0.59
P_FP     = float(pred["P1_next_detached"]["footprint_only_baseline"]) # 0.42
MU_WALL  = float(pred["axis_deg"])                                    # 49 deg
MU_FP    = float(pred["P2_direction"]["footprint_only_baseline_deg"]) # 10.5 deg
CAP      = float(pred["cap_deg"])                                     # 60 deg

# mixture parameters (step_b13): clustered fraction f, von Mises sigma
def find_mixture(d):
    """Locate f and sigma in the patch-scale result."""
    if isinstance(d, dict):
        if "f_clustered" in d and "sigma_deg" in d:
            return float(d["f_clustered"]), float(d["sigma_deg"])
        for v in d.values():
            r = find_mixture(v)
            if r: return r
    return None

mix = find_mixture(patch)
F_CLUS, SIG_CLUS = mix if mix else (0.43, 34.4)
KAPPA = 1.0 / math.radians(SIG_CLUS) ** 2      # von Mises kappa

logger.info(f"rates: wall p={P_WALL}, footprint p={P_FP}; "
            f"mixture f={F_CLUS}, sigma={SIG_CLUS} deg")

N_GRID = [20, 30, 50, 75, 100, 150, 200, 300, 500]
ALPHAS = [0.05, 0.01, 0.001]

# ------------------------------------------------------------------
# C1: cap-fraction channel
# ------------------------------------------------------------------

c1 = {}
for N in N_GRID:
    k_crit = {a: int(np.quantile(rng.binomial(N, P_FP, 200000), 1 - a)) + 1
              for a in ALPHAS}
    # exact critical counts via binomial survival function
    from scipy.stats import binom
    k_crit = {a: int(binom.isf(a, N, P_FP)) for a in ALPHAS}
    pw = {a: float(1 - binom.cdf(k - 1, N, P_WALL)) for a, k in k_crit.items()}
    # expected sigma separation: difference in proportions / joint sigma
    sd = math.sqrt(P_WALL * (1 - P_WALL) / N + P_FP * (1 - P_FP) / N)
    z = (P_WALL - P_FP) / sd
    c1[N] = {"k_crit": {str(a): k for a, k in k_crit.items()},
             "power": {str(a): pw[a] for a in ALPHAS},
             "z_sep": float(z)}

def N_for(power_target, alpha, chan):
    for N in N_GRID:
        if chan[N]["power"][str(alpha)] >= power_target:
            return N
    return None


def N_exact_binom(power_target, alpha):
    """Smallest n whose exact one-sided binomial cap-fraction test
    reaches power_target -- no N_GRID rounding."""
    for n in range(1, 600):
        kc = int(binom.isf(alpha, n, P_FP))
        if float(binom.sf(kc - 1, n, P_WALL)) >= power_target:
            return n
    return None

# ------------------------------------------------------------------
# C2: mean-direction channel
# ------------------------------------------------------------------

def draw_varpi_wall(N):
    """Mixture: fraction f clustered about MU_WALL (von Mises sigma),
    remainder footprint-distributed about MU_FP."""
    n_c = rng.binomial(N, F_CLUS)
    out = np.concatenate([
        rng.vonmises(math.radians(MU_WALL), KAPPA, n_c),
        rng.vonmises(math.radians(MU_FP), 0.5, N - n_c)])
    return out

def draw_varpi_fp(N):
    """Footprint-only: broad distribution about MU_FP.  The footprint
    itself is broad (its realized R ~ the measured value); modelled as
    a weak von Mises + uniform mix consistent with the observed
    footprint scatter."""
    return rng.vonmises(math.radians(MU_FP), 0.5, N)

def R_of(x):
    return float(np.abs(np.exp(1j * x).mean()))

c2 = {}
for N in N_GRID:
    R_fp   = np.array([R_of(draw_varpi_fp(N))   for _ in range(NMC // 4)])
    R_wall = np.array([R_of(draw_varpi_wall(N)) for _ in range(NMC // 4)])
    crit = {a: float(np.quantile(R_fp, 1 - a)) for a in ALPHAS}
    pw = {a: float(np.mean(R_wall > c)) for a, c in crit.items()}
    # direction offset: probability the wall sample's circular mean lands
    # inside the CI lower edge (17 deg) vs footprint mean 10.5
    mu_w = np.array([math.degrees(np.angle(np.exp(1j * draw_varpi_wall(N)).mean()))
                     for _ in range(2000)])
    c2[N] = {"R_crit": {str(a): c for a, c in crit.items()},
             "power": {str(a): pw[a] for a in ALPHAS},
             "med_mean_wall_deg": float(np.median(mu_w % 360))}

# ------------------------------------------------------------------
# C3: comet-channel power (fresh cohort)
# ------------------------------------------------------------------

FR_CAP   = float(pred["P4_comet_channel"]["frac_cap_of_sample"])      # 0.352
R_IN     = float(pred["P4_comet_channel"]["rate_gt_0p2_in"])          # 0.579
R_OUT    = float(pred["P4_comet_channel"]["rate_gt_0p2_out"])         # 0.229

c3 = {}
for M in N_GRID:
    n_in  = rng.binomial(M, FR_CAP,  NMC // 4)
    n_out = M - n_in
    x_in  = rng.binomial(np.maximum(n_in, 1),  R_IN,  (NMC // 4))
    x_out = rng.binomial(np.maximum(n_out, 1), R_OUT, (NMC // 4))
    # one-sided normal approx of the two-proportion z test
    p_pool = (x_in + x_out) / np.maximum(n_in + n_out, 1)
    with np.errstate(all="ignore"):
        z = (x_in / np.maximum(n_in, 1) - x_out / np.maximum(n_out, 1)) / \
            np.sqrt(p_pool * (1 - p_pool) * (1 / np.maximum(n_in, 1) +
                                             1 / np.maximum(n_out, 1)))
    pw = {a: float(np.mean(z > norm.ppf(1 - a))) for a in ALPHAS}
    c3[M] = {"power": {str(a): pw[a] for a in ALPHAS},
             "exp_n_in": float(FR_CAP * M)}

res = {
    "method": "Monte-Carlo discrimination power vs future sample size; "
              "wall and footprint rates from step_b19 (measured mixture "
              "and footprint-convolved baseline); NMC realizations per "
              "grid point; seeded RNG.",
    "seed": SEED, "nmc": NMC,
    "inputs": {"p_wall": P_WALL, "p_fp": P_FP,
               "mu_wall_deg": MU_WALL, "mu_fp_deg": MU_FP,
               "cap_deg": CAP, "f_clustered": F_CLUS,
               "sigma_deg": SIG_CLUS,
               "comet_frac_in_cap": FR_CAP,
               "comet_rate_in": R_IN, "comet_rate_out": R_OUT},
    "C1_cap_fraction": c1,
    "C2_mean_direction": c2,
    "C3_comet_channel": c3,
    "N_needed": {
        "C1_95pct_at_1pct":  N_for(0.95, 0.01, c1),
        "C1_95pct_at_0p1pct": N_for(0.95, 0.001, c1),
        "C2_95pct_at_1pct":  N_for(0.95, 0.01, c2),
        "C3_95pct_at_5pct":  N_for(0.95, 0.05, c3),
        "C3_95pct_at_1pct":  N_for(0.95, 0.01, c3),
    },
    "N_exact_C1_binomial": {
        "C1_95pct_at_5pct":   N_exact_binom(0.95, 0.05),
        "C1_95pct_at_1pct":   N_exact_binom(0.95, 0.01),
        "C1_95pct_at_0p1pct": N_exact_binom(0.95, 0.001),
        "note": "smallest n giving >=95 per cent power on the exact "
                "one-sided binomial cap-fraction test; the N_needed "
                "grid values above round up to the N_GRID mesh"},
}

out = str(RESULTS / "step_b32_detection_power.json")
json.dump(res, open(out, "w"), indent=1, default=float)

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4))

ax = axes[0]
for a, c in [(0.05, "0.5"), (0.01, "crimson"), (0.001, "navy")]:
    ax.plot(N_GRID, [c1[N]["power"][str(a)] for N in N_GRID],
            marker="o", ms=4, color=c, label=fr"$\alpha={a}$")
ax.axhline(0.95, color="k", ls=":", lw=1)
ax.set_xlabel("future detached-TNO discoveries $N$")
ax.set_ylabel("power (reject footprint null)")
ax.set_ylim(0, 1.02)
ax.legend(frameon=False, fontsize=8)
ax.set_title("C1: cap-fraction channel", fontsize=10)

ax = axes[1]
for a, c in [(0.05, "0.5"), (0.01, "crimson"), (0.001, "navy")]:
    ax.plot(N_GRID, [c2[N]["power"][str(a)] for N in N_GRID],
            marker="o", ms=4, color=c, label=fr"$\alpha={a}$")
ax.axhline(0.95, color="k", ls=":", lw=1)
ax.set_xlabel("future detached-TNO discoveries $N$")
ax.set_ylabel("power (mean-direction $R$)")
ax.set_ylim(0, 1.02)
ax.legend(frameon=False, fontsize=8)
ax.set_title("C2: mean-direction channel", fontsize=10)

ax = axes[2]
for a, c in [(0.05, "0.5"), (0.01, "crimson"), (0.001, "navy")]:
    ax.plot(N_GRID, [c3[N]["power"][str(a)] for N in N_GRID],
            marker="o", ms=4, color=c, label=fr"$\alpha={a}$")
ax.axhline(0.95, color="k", ls=":", lw=1)
ax.set_xlabel("fresh class-1-like comet cohort $M$")
ax.set_ylabel("power (in/out rate contrast)")
ax.set_ylim(0, 1.02)
ax.legend(frameon=False, fontsize=8)
ax.set_title("C3: comet channel", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b32_detection_power.png", dpi=300)

for N in N_GRID:
    logger.info(f"N={N}: C1 z={c1[N]['z_sep']:.2f} "
                f"power@1%={c1[N]['power']['0.01']:.2f} | "
                f"C2 power@1%={c2[N]['power']['0.01']:.2f} | "
                f"C3 power@1%={c3[N]['power']['0.01']:.2f}")
logger.info(f"N needed: {res['N_needed']}")
logger.data_save(out)
logger.data_save(FIG / "supplementary" / "step_b32_detection_power.png")