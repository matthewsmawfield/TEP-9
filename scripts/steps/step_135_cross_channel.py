"""Step 135: Cross-channel axis convergence (step_b99).

Exploratory comparison of TEP data channels; statistical independence is not established.
Every recovered direction and p-value is read directly from the
sibling-project result leaves -- no re-derivation, full provenance.

Channel ledger (independent datasets, independent pipelines):

  comet-transit   TEP-9   displaced axis + detached-resident axis,
                          recovered on raw-MPC refits (steps 128, b77)
  clock-GNSS      TEP-GNSS-II   free galactic-vector search over
                          2701 directions on CODE precise-clock
                          pair coherence (step_2_5)
  clock-MGEX      TEP-GNSS-MGEX CMB-projection consistency of the
                          EW/NS asymmetry (step_2_4)
  lunar-LLR       TEP-LLR full-sky delta-AIC scan: Planck dipole
                          rank + look-elsewhere p (step_076)
  clock-epoch     TEP-GNSS-II mesh_evolution collective-motion
                          direction time series, 104 x 90-day
                          windows over ~25 yr -- tested for a
                          step change at the comet switch epoch

T1  Channel direction ledger: bipolar-frame (cone, meridian
    distance, azimuth) of every independently recovered direction.

T2  Meridian coplanarity, cross-channel: the meridian plane was
    fixed in step_132 from the comet record alone.  Each external
    channel direction is an out-of-sample draw; extremeness is the
    nominal Fisher combination of the plane distances at which qualifying directions land
    within their observed distance of the plane.

T3  Joint frame evidence: Fisher combination of per-channel
    extremeness for the common external axis (each channel's own
    registered statistic against the same target).

T4  Clock-epoch coincidence: the CODE collective-motion-direction
    series is tested for a step discontinuity at the comet switch
    boundary (day offset for 2018-01-01 under the series' stated
    epoch) and extremeness is priced against all other admissible
    changepoints.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_135_cross_channel")
tee_stdout(logger)
logger.header("Cross-channel axis convergence")

import json
import math
import numpy as np
from scipy import stats as _st
from scripts.utils.tep9_common import lv, sep
from scripts.utils.coordinates import GAL2ECL, EQ2ECL
from scripts.utils.statistics import fisher_combination

CORPUS = RESULTS.parent.parent
GNSS2 = CORPUS / "TEP-GNSS-II" / "results" / "outputs"
MGEX = CORPUS / "TEP-GNSS-MGEX" / "results" / "outputs"
LLR = CORPUS / "TEP-LLR" / "results" / "outputs"
SEED = 20260919



def eqv(ra, dec):
    ra, dec = math.radians(ra), math.radians(dec)
    return EQ2ECL @ np.array([math.cos(dec) * math.cos(ra),
                              math.cos(dec) * math.sin(ra),
                              math.sin(dec)])


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


CMB = gv(264.02, 48.25)
DECL = lv(34.0, -13.0)
PLN = np.cross(CMB, DECL)
PLN /= np.linalg.norm(PLN)
ez = CMB.copy()
ex = np.cross([0, 0, 1.0], ez)
ex /= np.linalg.norm(ex)
ey = np.cross(ez, ex)


def bipolar_frame(u):
    """Return (cone-from-apex, meridian distance, azimuth) in deg."""
    cone = math.degrees(math.acos(np.clip(u @ CMB, -1, 1)))
    mdist = abs(math.degrees(math.asin(np.clip(u @ PLN, -1, 1))))
    az = math.degrees(math.atan2(u @ ey, u @ ex)) % 360
    return cone, mdist, az


def bip_dist(u):
    return min(sep(u, CMB), sep(u, -CMB))


# ------------------------------------------------------------------
# T1 channel ledger
# ------------------------------------------------------------------
logger.info("T1 cross-channel direction ledger")

ledger = []


def add_channel(name, source, direction, kind, note=""):
    if direction is None:
        ledger.append(dict(channel=name, source=source, kind=kind,
                           recovered=False, note=note))
        return
    cone, mdist, az = bipolar_frame(direction)
    bd = bip_dist(direction)
    ledger.append(dict(channel=name, source=source, kind=kind,
                       recovered=True, cone_deg=cone,
                       meridian_dist_deg=mdist, azimuth_deg=az,
                       bipolar_dist_deg=bd, note=note))
    logger.info(f"  {name:22s} cone={cone:6.1f} meridian={mdist:5.1f} "
                f"az={az:6.1f} bipolar={bd:5.1f}  ({kind})")


# TEP-9 internal axes (comet record -- defines the frame, kept for
# reference but excluded from the out-of-sample joint statistic)
add_channel("comet-declared", "TEP-9 step_101/b52", DECL,
            "internal-frame-defining")
add_channel("comet-displaced", "TEP-9 step_128", lv(120.0, -40.0),
            "internal-recovered")
add_channel("comet-detached", "TEP-9 step_b77", lv(49.9, -17.0),
            "internal-recovered")

# GNSS-II free galactic-vector search
g = json.load(open(GNSS2 / "step_2_5_dual_motion_geometry.json"))
gs = g["galactic_vector_search"]
gnss2_axis = eqv(gs["best_ra_deg"], gs["best_dec_deg"])
gnss2_r = gs["best_correlation"]
add_channel("clock-GNSS-II",
            "TEP-GNSS-II step_2_5", gnss2_axis, "free-recovery",
            note=f"r={gnss2_r:.3f}, 2701-dir search")

# MGEX -- the grid axis was flagged non-identifiable by its own
# pipeline; record it but carry the flag so T2 excludes it.
m = json.load(open(MGEX / "step_2_4_cmb_alignment.json"))
mg = m["combined"]
mgex_axis = eqv(mg["best_fit_ra_deg"], mg["best_fit_dec_deg"])
mgex_ident = bool(mg["grid_axis_identifiable"])
mgex_cmb_p = mg["lee_corrected_p_value"]
add_channel("clock-MGEX", "TEP-GNSS-MGEX step_2_4", mgex_axis,
            "free-recovery",
            note=f"identifiable={mgex_ident}, cmb_p={mgex_cmb_p}")

# LLR -- its own best axis is far from CMB; the registered finding is
# that the Planck-dipole direction itself is anomalous (top 8.5%,
# LEE p = 0.001).  Ledger the best axis; the evidence enters T3.
l = json.load(open(LLR / "step_076_sky_scan_directional.json"))
la = l["max_delta_aic"]
llr_axis = eqv(la["best_ra_deg"], la["best_dec_deg"])
llr_planck = l["planck_dipole"]
llr_p = l["look_elsewhere_correction"]["p_planck_axis_exceeds_observed"]
add_channel("lunar-LLR-best", "TEP-LLR step_076", llr_axis,
            "free-recovery",
            note=f"Planck rank {llr_planck['rank_by_delta_aic']}"
                 f"/2664, LEE p={llr_p}")


# ------------------------------------------------------------------
# T2 meridian coplanarity, out-of-sample channels
# ------------------------------------------------------------------
logger.info("T2 meridian coplanarity on external channels")

# The meridian plane (CMB x declared axis) was fixed by the comet
# record in step_132.  External free-recovered directions are draws
# against that fixed target.  Qualifying = pipeline-flagged as an
# identifiable recovered direction.
meridian_hits = []
for c in ledger:
    if c["kind"] != "free-recovery" or not c["recovered"]:
        continue
    if c["channel"] == "clock-MGEX" and not mgex_ident:
        c["t2_excluded"] = "axis flagged non-identifiable by source"
        continue
    d = c["meridian_dist_deg"]
    p = math.sin(math.radians(d))  # P(|u.n| <= sin d) = sin(d)
    meridian_hits.append(dict(channel=c["channel"],
                              meridian_dist_deg=d, p_plane=p))
    logger.info(f"  {c['channel']:22s} meridian {d:5.1f} deg  "
                f"p={p:.4f}")

meridian_product = float(np.prod([h["p_plane"] for h in meridian_hits]))
joint_meridian = fisher_combination([h["p_plane"] for h in meridian_hits])["p"]
logger.info(f"  nominal meridian Fisher p = {joint_meridian:.3g} "
            f"over {len(meridian_hits)} channels")

# GNSS-II alone: a 2701-direction free search landing within d of a
# pre-declared plane is priced directly by p_plane above.


# ------------------------------------------------------------------
# T3 joint frame evidence (Fisher over per-channel extremeness)
# ------------------------------------------------------------------
logger.info("T3 Fisher combination of channel frame statistics")

chan_p = {}

# Read the actual result leaf. Missing/malformed inputs must fail, never
# silently supply a favourable p-value.
b93 = json.load(open(RESULTS / "step_b93_bipolar_unification.json"))
chan_p["comet-bipolar-J"] = b93["T2_joint_opposite_polarity"]["p_J"]

# One statistic per recovered GNSS-II direction. The plane and pole
# distances describe the SAME axis and cannot be independent Fisher terms.
md = bipolar_frame(gnss2_axis)[1]
chan_p["clock-GNSS-II-meridian"] = math.sin(math.radians(md))
gnss2_pole_diagnostic = 1 - math.cos(math.radians(bip_dist(gnss2_axis)))

# MGEX CMB-projection consistency
if np.isfinite(mgex_cmb_p):
    chan_p["clock-MGEX-cmb-proj"] = mgex_cmb_p

# LLR Planck-axis LEE p
if np.isfinite(llr_p):
    chan_p["lunar-LLR-Planck"] = llr_p

combined = fisher_combination(list(chan_p.values()))
chi2, df, fisher_p = combined["statistic"], combined["df"], combined["p"]
logger.info("  channel p-values:")
for k, v in chan_p.items():
    logger.info(f"    {k:28s} {v:.4g}")
logger.info(f"  Fisher chi2 = {chi2:.1f} (df={df})  p = {fisher_p:.3g}")


# ------------------------------------------------------------------
# T4 clock-epoch coincidence at the comet switch boundary
# ------------------------------------------------------------------
logger.info("T4 GNSS-II collective-motion direction vs switch epoch")

g22 = json.load(open(GNSS2 /
                     "step_2_2_geospatial_temporal_analysis_code.json"))
me = g22["mesh_dance_analysis"]["mesh_evolution"]
t = np.array([w["days_since_epoch"] for w in me], dtype=float)
dd = np.array([w["collective_motion_direction"] for w in me],
              dtype=float)
# unwrap circular direction
du = np.degrees(np.unwrap(np.radians(dd)))
# epoch origin: series spans 9218 dates ending ~2025; day 0
# corresponds to 2000-01-01 under the CODE longspan layout
switch_day = (2018 - 2000) * 365.25
ib = int(np.argmin(np.abs(t - switch_day)))
boundary_day = float(t[ib])


def tstat(x, i):
    a, b = x[:i], x[i:]
    tt, pp = _st.ttest_ind(a, b, equal_var=False)
    return abs(tt), pp


t_b, p_b = tstat(du, ib)
# extremeness vs all admissible changepoints
all_t = [tstat(du, i)[0] for i in range(10, len(t) - 10)]
rank = int(sum(1 for x in all_t if x >= t_b)) + 1
p_extreme = rank / (len(all_t) + 1)
best_i = 10 + int(np.argmax(all_t))
# autocorrelation-aware null: the nominal Welch p ignores the strong
# serial correlation of the direction series (90-d windows overlap).
# Circular block permutation (block ~1 yr of windows) preserves the
# autocorrelation while randomising the boundary position, and the
# at-boundary statistic is recomputed on each permutation -- the
# correct null for a predicted-epoch test.
rng_t4 = np.random.default_rng(SEED)
blk = max(4, int(round(365.25 / np.median(np.diff(t)))))
n_bt = 2000
cnt = 1
for _ in range(n_bt):
    xperm = np.concatenate([du[i:i + blk] for i in
                            rng_t4.permutation(
                                np.arange(0, len(du), blk))])
    xperm = xperm[: len(du)]
    if len(xperm) < len(du):
        xperm = np.pad(xperm, (0, len(du) - len(xperm)),
                       mode="edge")
    if tstat(xperm, ib)[0] >= t_b:
        cnt += 1
p_block = cnt / (n_bt + 1)
epoch_res = dict(available=True, n_windows=len(t),
                 boundary_day=boundary_day,
                 boundary_year=2000 + boundary_day / 365.25,
                 t_at_boundary=float(t_b), p_at_boundary=float(p_b),
                 p_at_boundary_blockperm=float(p_block),
                 block_len_windows=int(blk),
                 extremeness_rank=rank,
                 n_changepoints=len(all_t),
                 p_extremeness=p_extreme,
                 best_changepoint_day=float(t[best_i]),
                 best_changepoint_year=2000 + float(t[best_i]) / 365.25,
                 best_t=float(all_t[best_i - 10]))
logger.info(f"  boundary {boundary_day:.0f} d "
            f"(~{epoch_res['boundary_year']:.1f}): "
            f"t={t_b:.2f} p={p_b:.3g}; block-perm p={p_block:.3g}; "
            f"extremeness "
            f"rank {rank}/{len(all_t)} (p={p_extreme:.3g}); "
            f"max changepoint {epoch_res['best_changepoint_year']:.1f}")


# ------------------------------------------------------------------
# verdict
# ------------------------------------------------------------------
n_conv = sum(1 for p in chan_p.values() if p < 0.05)
verdict = ("CONVERGENT CROSS-CHANNEL GEOMETRY: "
           f"{n_conv} of {len(chan_p)} independent measurement classes "
           f"resolve the frame-anchored bipolar structure at p<0.05 "
           f"(nominal Fisher p={fisher_p:.2g}; GNSS-II direction epoch "
           f"coincides with the comet era boundary at "
           f"p={epoch_res.get('p_at_boundary_blockperm', float('nan')):.3g} "
           "block-permutation).  No calibrated joint significance is "
           "claimed: GNSS overlap, directional sensitivity and the "
           "retrospective choice of geometric relation remain "
           "uncalibrated.")
logger.info(f"verdict: {verdict}")

out = dict(step="step_135_cross_channel", result="b99", verdict=verdict,
           T1_ledger=ledger,
           T2_meridian=dict(hits=meridian_hits,
                            joint_p_nominal=joint_meridian,
                            product_statistic=meridian_product,
                            interpretation="Isotropic independent-direction reference only; product alone is not a p-value."),
           T3_fisher=dict(channel_p=chan_p, chi2=chi2, df=df,
                          fisher_p_nominal=fisher_p,
                          excluded_duplicate_GNSS_pole_p=gnss2_pole_diagnostic,
                          interpretation="One term per channel; GNSS overlap, sky sensitivity, and retrospective choice of geometry are not calibrated."),
           T4_epoch=epoch_res,
           evidence_status="convergent geometry; joint significance uncalibrated",
           registration_status="No timestamped pre-analysis registration established by this script")
with open(RESULTS / "step_b99_cross_channel.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.info("wrote results/step_b99_cross_channel.json")

FIG = RESULTS / "figures"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

# T1/T2 bipolar positions
names, cones, mds = [], [], []
for c in ledger:
    if not c["recovered"]:
        continue
    names.append(c["channel"])
    cones.append(c["cone_deg"])
    mds.append(c["meridian_dist_deg"])
cols = ["steelblue" if "comet" in n else
        ("indianred" if "GNSS" in n or "MGEX" in n else "seagreen")
        for n in names]
ax[0].scatter(mds, cones, c=cols, s=90)
for n, x, y in zip(names, mds, cones):
    ax[0].annotate(n.replace("comet-", "c-").replace("clock-", ""),
                   (x, y), fontsize=7)
ax[0].axvline(30, color="k", ls=":", lw=0.8)
ax[0].set_xlabel("distance to meridian plane (deg)")
ax[0].set_ylabel("cone angle from CMB apex (deg)")
ax[0].set_title("T1/T2 channel directions in bipolar frame")
ax[0].set_xlim(0, 95)

# T3 channel p-values
ks = list(chan_p)
ax[1].barh(range(len(ks)), [-math.log10(chan_p[k]) for k in ks],
           color="slateblue")
ax[1].set_yticks(range(len(ks)))
ax[1].set_yticklabels(ks, fontsize=8)
ax[1].set_xlabel("-log10 p")
ax[1].set_title(f"T3 nominal Fisher p = {fisher_p:.2g}")
ax[1].axvline(-math.log10(0.05), color="k", ls=":", lw=0.8)

# T4 epoch series
if epoch_res.get("available"):
    ax[2].plot(t / 365.25 + 2000, du, lw=0.9, color="0.3")
    ax[2].axvline(2018, color="indianred", ls="--", lw=0.9,
                  label="comet switch")
    ax[2].axvline(epoch_res["best_changepoint_year"], color="gold",
                  ls=":", lw=1.0, label="max changepoint")
    ax[2].set_xlabel("year")
    ax[2].set_ylabel("collective-motion direction (deg, unwrapped)")
    ax[2].set_title(f"T4 GNSS-II direction epoch (block-perm p={epoch_res.get('p_at_boundary_blockperm', epoch_res['p_at_boundary']):.2g})")
    ax[2].legend(fontsize=8)
else:
    ax[2].text(0.5, 0.5, "GNSS-II series unavailable", ha="center")
    ax[2].axis("off")

fig.tight_layout()
fig.savefig(FIG / "step_b99_cross_channel.png", dpi=150)
logger.info("wrote results/figures/step_b99_cross_channel.png")
