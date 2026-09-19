"""Step 139: Retrospective geometry audit (step_b103).

The input directions were already inspected in steps 131-135. Freezing them
now does not turn those data into a prospective test. This audit reports
geometric tail areas, retains each ISO's three orbit directions as one
rotation block, and excludes the resident direction that motivated the
transit sector. None of these tail areas is a discovery significance or
posterior probability for TEP. Future data require a timestamped protocol.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))

import json
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import combine_pvalues, mannwhitneyu

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (
    DATA_RAW, RESULTS, lv, sep)
from scripts.utils.coordinates import GAL2ECL, EQ2ECL
from scripts.utils.statistics import monte_carlo_tail, nearest_cone_probability

logger = StepLogger("step_139_frozen_geometry")

CORPUS = RESULTS.parent.parent
GNSS2 = CORPUS / "TEP-GNSS-II" / "results" / "outputs"
MGEX = CORPUS / "TEP-GNSS-MGEX" / "results" / "outputs"
LLR = CORPUS / "TEP-LLR" / "results" / "outputs"

SEED = 20260919
N_SHUF = 20000

# ------------------------------------------------------------------
# frozen frame
# ------------------------------------------------------------------



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


CMB = gv(264.02, 48.25)          # frozen: corpus-declared frame
DECL = lv(34.0, -13.0)           # frozen: declared pre-2018 axis
DISP = lv(120.0, -40.0)          # recovered on post-2017 record
RES = lv(49.9, -17.0)            # held-out: resident detached axis

# meridian plane from frozen inputs only
PLN = np.cross(CMB, DECL)
PLN /= np.linalg.norm(PLN)
# bipolar frame basis (azimuth about the CMB axis)
ez = CMB.copy()
ex = np.cross([0, 0, 1.0], ez)
ex /= np.linalg.norm(ex)
ey = np.cross(ez, ex)


def bipolar(u):
    """(cone from apex, meridian distance to frozen plane, azimuth)."""
    cone = math.degrees(math.acos(np.clip(u @ CMB, -1, 1)))
    md = abs(math.degrees(math.asin(np.clip(u @ PLN, -1, 1))))
    az = math.degrees(math.atan2(u @ ey, u @ ex)) % 360.0
    return cone, md, az


def az_shuffle(u, n, rng):
    """Rotate u about the CMB axis by uniform azimuths."""
    cone, _, _ = bipolar(u)
    th = math.radians(cone)
    ph = rng.uniform(0, 2 * math.pi, n)
    return (math.cos(th) * ez[None, :]
            + math.sin(th) * (np.cos(ph)[:, None] * ex[None, :]
                              + np.sin(ph)[:, None] * ey[None, :]))


# ------------------------------------------------------------------
# held-out direction ledger
# ------------------------------------------------------------------
dirs = []

# post-2017 recovered displaced axis (held-out relative to the
# frozen plane; semi-held-out relative to the mirror prediction)
dirs.append(dict(channel="comet-displaced", u=DISP,
                 tier="post2017-recovered",
                 note="free-scan axis, independent-refit confirmed"))

# resident detached axis (TNO population -- never entered the fit)
dirs.append(dict(channel="comet-detached-resident", u=RES,
                 tier="training-related",
                 note="resident population motivated the transit sector; exclude from combined score"))

# ISO leg axes from the harmonic-bipolar ledger (b97): reconstruct
# unit vectors from the stored bipolar theta/phi coordinates.
try:
    b97 = json.load(open(RESULTS / "step_b97_harmonic_bipolar.json"))
    for leg in b97["T3_iso_meridian"]["legs"]:
        th = math.radians(leg["theta"])
        ph = math.radians(leg["phi"])
        u = (math.cos(th) * ez
             + math.sin(th) * (math.cos(ph) * ex + math.sin(ph) * ey))
        dirs.append(dict(channel=f"ISO-{leg['des']}-{leg['leg']}",
                         u=u / np.linalg.norm(u), tier="retrospective",
                         note=f"ISO {leg['des']} {leg['leg']} axis"))
    logger.info(f"  loaded {sum(1 for d in dirs if 'ISO' in d['channel'])} "
                f"ISO leg axes from step_b97")
except FileNotFoundError as exc:
    logger.warning(f"  ISO ledger unavailable: {exc}")

# external channel axes (same readers as step_135)
try:
    g = json.load(open(GNSS2 / "step_2_5_dual_motion_geometry.json"))
    gs = g["galactic_vector_search"]
    u = eqv(gs["best_ra_deg"], gs["best_dec_deg"])
    dirs.append(dict(channel="clock-GNSS-II", u=u / np.linalg.norm(u),
                     tier="retrospective",
                     note=f"2701-dir free search, r={gs['best_correlation']:.3f}"))
except FileNotFoundError as exc:
    logger.warning(f"  GNSS-II axis unavailable: {exc}")
try:
    m = json.load(open(MGEX / "step_2_4_cmb_alignment.json"))
    mg = m["combined"]
    u = eqv(mg["best_fit_ra_deg"], mg["best_fit_dec_deg"])
    dirs.append(dict(channel="clock-MGEX", u=u / np.linalg.norm(u),
                     tier="non-identifiable",
                     note=f"identifiable={mg['grid_axis_identifiable']}"))
except FileNotFoundError as exc:
    logger.warning(f"  MGEX axis unavailable: {exc}")
try:
    l = json.load(open(LLR / "step_076_sky_scan_directional.json"))
    la = l["max_delta_aic"]
    u = eqv(la["best_ra_deg"], la["best_dec_deg"])
    dirs.append(dict(channel="lunar-LLR", u=u / np.linalg.norm(u),
                     tier="retrospective", note="free sky-scan axis"))
except FileNotFoundError as exc:
    logger.warning(f"  LLR axis unavailable: {exc}")

for d in dirs:
    cone, md, az = bipolar(d["u"])
    d.update(cone_deg=cone, meridian_dist_deg=md, azimuth_deg=az)
    logger.info(f"  {d['channel']:26s} [{d['tier']:>18s}] "
                f"cone={cone:6.1f} merid={md:6.1f} az={az:6.1f}")

res = {"step": "step_139_frozen_geometry", "result": "b103",
       "description": "Retrospective geometric sensitivity audit with orbit-block rotations",
       "evidence_status": "exploratory; no prospective confirmation",
       "frozen_inputs": {
           "cmb_apex_gal_lb": [264.02, 48.25],
           "declared_axis_ecl_lb": [34.0, -13.0],
           "hypothesis": "one frame-anchored bipolar lapse field"},
       "n_shuffle": N_SHUF, "seed": SEED,
       "direction_ledger": [
           {k: (round(v, 4) if isinstance(v, float) else v)
            for k, v in d.items() if k != "u"} for d in dirs]}

rng = np.random.default_rng(SEED)

# ------------------------------------------------------------------
# T1 mirror prediction -- parameter-free
# ------------------------------------------------------------------
logger.info("T1 mirror prediction: equatorial reflection of the "
            "declared axis about the CMB axis")
mirror = DECL - 2.0 * (DECL @ CMB) * CMB
mirror /= np.linalg.norm(mirror)
sep_obs = sep(mirror, DISP)
# null: azimuth shuffle of the recovered displaced axis about the
# CMB axis -- preserves its recovered cone angle
sh = az_shuffle(DISP, N_SHUF, rng)
seps = np.degrees(np.arccos(np.clip(sh @ mirror, -1, 1)))
p_mirr = monte_carlo_tail(seps, sep_obs, "less")
res["T1_mirror_prediction"] = {
    "predicted_mirror_cone_deg": bipolar(mirror)[0],
    "predicted_mirror_az_deg": bipolar(mirror)[2],
    "recovered_displaced_cone_deg": bipolar(DISP)[0],
    "recovered_displaced_az_deg": bipolar(DISP)[2],
    "sep_deg": float(sep_obs),
    "p_azimuth_shuffle": p_mirr,
    "p_fixed_target": float((1 - math.cos(math.radians(sep_obs)))
                            / 2.0),
    "note": ("Retrospective mirror comparison: parameter-free after choosing "
             "the mirror rule, but the rule was examined on these directions. "
             "Neither target-selection history nor survey sensitivity is calibrated.")}

logger.info(f"  predicted mirror at cone {res['T1_mirror_prediction']['predicted_mirror_cone_deg']:.1f} "
            f"az {res['T1_mirror_prediction']['predicted_mirror_az_deg']:.1f}; "
            f"recovered displaced {sep_obs:.2f} deg away; "
            f"p={p_mirr:.4f}")

# ------------------------------------------------------------------
# T2 frozen meridian -- all recovered directions as held-out draws
# ------------------------------------------------------------------
logger.info("T2 frozen meridian {CMB x declared} on recovered "
            "directions")
scored = [d for d in dirs if d["tier"] not in {"training-related", "non-identifiable"}]
t2 = []
for d in scored:
    md = d["meridian_dist_deg"]
    p_plane = math.sin(math.radians(md))  # P(|u.n| <= sin d)
    t2.append(dict(channel=d["channel"], tier=d["tier"],
                   meridian_dist_deg=float(md),
                   p_plane=float(p_plane)))
    logger.info(f"  {d['channel']:26s} meridian {md:6.2f}  "
                f"p={p_plane:.4f}")
# A single random rotation is applied to all directions of the SAME ISO.
# Independent rotations would destroy the orbital constraints and inflate
# the effective sample size. Rotations preserve every within-orbit dot product.
pvec = np.array([t["p_plane"] for t in t2])
chi2_obs = float(-2.0 * np.log(np.clip(pvec, 1e-300, 1)).sum())
groups = {}
for d in scored:
    group = "-".join(d["channel"].split("-")[:2]) if d["channel"].startswith("ISO-") else d["channel"]
    groups.setdefault(group, []).append(d["u"])
null_score = np.zeros(N_SHUF)
for members in groups.values():
    ph = rng.uniform(0, 2*np.pi, N_SHUF)
    for u in members:
        parallel = (u @ CMB)*CMB
        rotated = (parallel + np.cos(ph)[:, None]*(u-parallel)
                   + np.sin(ph)[:, None]*np.cross(CMB, u))
        p_null = np.abs(rotated @ PLN)
        null_score += -2*np.log(np.clip(p_null, 1e-300, 1))
p_joint = monte_carlo_tail(null_score, chi2_obs)
res["T2_frozen_meridian"] = {
    "plane_definition": "normal = CMB x declared (frozen inputs only)",
    "per_direction": t2,
    "joint_score": chi2_obs, "n_rotation_blocks": len(groups),
    "n_shuffle": N_SHUF,
    "interpretation": "Retrospective conditional-azimuth reference, preserving ISO orbit geometry; no global selection calibration.",
    "joint_p_azimuth_shuffle": float(p_joint),
    "excluded": [{"channel": d["channel"], "note": d["note"]}
                 for d in dirs if d["tier"] in {"training-related", "non-identifiable"}]}
logger.info(f"  joint meridian chi2={chi2_obs:.1f} on {len(t2)} "
            f"directions, shuffle p={p_joint:.4f}")

# ------------------------------------------------------------------
# T3 frozen cone band -- folded cone vs nearest comet lobe
# ------------------------------------------------------------------
logger.info("T3 frozen cone band on held-out directions")
lobes = {"declared": min(bipolar(DECL)[0], 180-bipolar(DECL)[0]),
         "displaced": min(bipolar(DISP)[0], 180-bipolar(DISP)[0])}
t3 = []
for d in scored:
    if d["channel"] == "comet-displaced":
        continue  # the displaced direction defines a target cone
    folded = min(d["cone_deg"], 180.0 - d["cone_deg"])
    near = min(lobes.values(), key=lambda x: abs(x - folded))
    delta = abs(folded - near)
    p_cone = nearest_cone_probability(folded, list(lobes.values()))
    t3.append(dict(channel=d["channel"], tier=d["tier"],
                   folded_cone_deg=float(folded),
                   nearest_lobe_deg=float(near),
                   delta_deg=float(delta),
                   p_cone=float(p_cone)))
    logger.info(f"  {d['channel']:26s} folded cone {folded:6.1f} "
                f"(lobe {near:.1f}, delta {delta:5.1f}) p={p_cone:.3f}")
res["T3_frozen_cone"] = {
    "comet_lobe_cones_deg": lobes,
    "per_direction": t3,
    "note": ("Retrospective folded cone angle, scoring the union of both lobe bands; isotropic null "
             "per direction -- polar directions (GNSS-II, LLR) are "
             "honestly off-band")}

# ------------------------------------------------------------------
# T4 provenance partition + joint predictive score
# ------------------------------------------------------------------
tiers = {}
for d in dirs:
    tiers.setdefault(d["tier"], []).append(d["channel"])
res["T4_provenance"] = {
    "partition": tiers,
    "prospectively_held_out_channels": [],
    "note": ("All directions were available during geometry exploration. The "
             "resident sample informed the axis; the displaced direction defines "
             "one cone; three legs of each ISO share one orbit. A future test "
             "must freeze hypotheses, cuts, nulls and outcomes before data inspection.")}
res["verdict"] = "RETROSPECTIVE GEOMETRY ONLY; prospective confirmation not established"
logger.info("verdict: " + res["verdict"])

out = RESULTS / "step_b103_frozen_geometry.json"
json.dump(res, open(out, "w"), indent=1, default=float)
logger.info(f"wrote {out}")

# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axes[0]
for d in dirs:
    col = {"retrospective": "tab:blue",
           "training-related": "grey",
           "post2017-recovered": "tab:red",
           "non-identifiable": "grey"}[d["tier"]]
    ax.scatter(d["azimuth_deg"], d["cone_deg"], c=col, s=40,
               zorder=3)
    ax.annotate(d["channel"].replace("comet-", "").replace(
        "clock-", "").replace("lunar-", ""),
        (d["azimuth_deg"], d["cone_deg"]),
        textcoords="offset points", xytext=(5, 4), fontsize=7)
for nm, u, mk in (("declared", DECL, "k^"), ("mirror-pred", mirror,
                                              "g*"),
                  ("displaced", DISP, "rs"), ("resident", RES, "bD"),
                  ("CMB apex", CMB, "mo")):
    cone, md, az = bipolar(u)
    ax.scatter(az, cone, marker=mk[-1], c=mk[0], s=90 if mk == "g*" else 60,
               zorder=4, label=nm)
ax.set_xlabel("azimuth about CMB axis (deg)")
ax.set_ylabel("cone angle from apex (deg)")
ax.legend(fontsize=7, loc="upper right")
ax.set_title("retrospective direction ledger")

ax = axes[1]
lab = [t["channel"].replace("comet-", "") for t in t2]
ax.barh(range(len(t2)), [t["meridian_dist_deg"] for t in t2],
        color=["tab:red" if t["tier"] == "post2017-recovered"
               else "tab:blue" for t in t2])
ax.set_yticks(range(len(t2)))
ax.set_yticklabels(lab, fontsize=7)
ax.invert_yaxis()
ax.set_xlabel("distance to frozen meridian plane (deg)")
ax.set_title(f"orbit-block meridian reference (p={p_joint:.3g})")
fig.tight_layout()
fig.savefig(FIG / "step_b103_frozen_geometry.png", dpi=150)
logger.info("wrote results/figures/step_b103_frozen_geometry.png")
print(json.dumps(res["T1_mirror_prediction"], indent=1))
