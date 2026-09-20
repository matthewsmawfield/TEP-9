#!/usr/bin/env python3
"""step_126: Ten-channel global cross-survey and multi-lineage synthesis.

This step aggregates ten overlapping observational summaries across six
distinct catalogue lineages and multiple dynamical populations, all testing
convergence on the declared outer solar system proper-time domain boundary:

  Resident population channels (identical per-axis in-cap Rayleigh
  machinery on each survey's own per-object elements):
    C1: SBDB secure detached TNOs (a>150, q>30, cc<=3, N=44) -- varpi clustering
    C2: DES Y6 detached TNOs (Bernardinelli et al. 2022, N=16) -- independent survey/fit
    C3: MPCORB census detached TNOs (MPC lineage, N=63) -- independent catalogue
    C4: MPCORB distant belt mean-plane lean azimuth (80-400 AU non-resonant)

  Cometary transit & spatial channels:
    C5: CODE comets orig->fut rotation without energy exchange (N=131, class-1)
    C6: Warsaw comets inbound-leg rotation discrepancy (N=98, near-parabolic)
    C7: One-apparition comets uncertainty-normalized rotation (N=30)
    C8: MPC CometEls census near-parabolic aphelion dipole (fragment-merged census, as step_108)
    C9: SBDB post-2017 prospective comets aphelion dipole (N=287, 2018-2026)

  Inflow / Outflow bipolar channel:
    C10: Jupiter-family comets (JFCs) antipodal perihelion concentration (N=378)

A 20,000-draw random-axis comparison evaluates a conditional sky-area rank: the fraction of random sky directions whose joint Fisher statistic
S = -2 sum ln p reaches or exceeds the value measured at the declared boundary axis.

Inputs
------
data/raw/sbdb/sbdb_outer_ss.json
data/raw/des/y6_des_tnos_color.fits
data/raw/mpc/MPCORB.DAT.gz
results/step_b28_bidirectional_rotation.csv
results/step_b30_proper_time_slip.csv
results/step_b72_cometels_dipole.csv
results/step_b88_lineage_aphelion.csv
results/step_b55_mean_plane.json
results/step_b17_injection_chain.json

Outputs
-------
results/step_b90_global_synthesis.json
results/step_b90_global_synthesis.csv
results/figures/step_b90_global_synthesis.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.statistics import monte_carlo_p
from scripts.utils.tep9_common import DATA_RAW, RESULTS, lv, lb, sep, tee_stdout
logger = StepLogger("step_126_global_cross_survey_synthesis")
tee_stdout(logger)
logger.header("Ten-channel global cross-survey and multi-lineage synthesis")

import csv
import json
import math
import numpy as np
from scipy.stats import norm, rankdata, chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20261031
NAX = 20000
CAP = 60.0
rng = np.random.default_rng(SEED)

AX_TNO = lv(49.0, -17.0)          # detached-TNO cluster axis
AX_COM = lv(34.0, -13.0)          # declared comet transit axis
ANTI_AX = lv(229.0, 17.0)         # antipodal axis

# ------------------------------------------------------------------
# 1. Load Resident Channels (C1: SBDB, C2: DES, C3: MPCORB, C4: Mean-plane)
# ------------------------------------------------------------------

# C1: SBDB Detached TNOs
def fnum(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")

def perih_dir(om, Om, inc):
    co, so = math.cos(om), math.sin(om)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(inc), math.sin(inc)
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci, so * si])

_d = json.load(open(DATA_RAW / "sbdb" / "sbdb_outer_ss.json"))
sbdb = [dict(zip(_d["fields"], rec)) for rec in _d["data"]]
tnos_sbdb = []
for r in sbdb:
    a, q, cc = fnum(r, "a"), fnum(r, "q"), fnum(r, "condition_code")
    if not (a > 150 and q > 30 and 0 <= cc <= 3):
        continue
    Om, w, i = fnum(r, "om"), fnum(r, "w"), fnum(r, "i")
    if np.isfinite(Om) and np.isfinite(w) and np.isfinite(i):
        tnos_sbdb.append(dict(
            dir=perih_dir(math.radians(w), math.radians(Om), math.radians(i)),
            varpi=math.radians((Om + w) % 360.0)))
logger.info(f"C1 SBDB secure detached TNOs: {len(tnos_sbdb)}")

# C2: DES Detached TNOs -- per-object elements from the Y6 FITS,
# matched to the resident-ledger cohort by MPC designation so the
# channel is evaluated with the same per-axis in-cap Rayleigh
# machinery as C1 rather than a direction-proximity score.
des_res = json.load(open(RESULTS / "step_b77_resident_ledger.json"))
des_entry = next(s for s in des_res["per_survey"] if s["survey"] == "DES")
p_des_obs = float(des_entry["p_vs_baseline"])
logger.info(f"C2 DES detached TNOs (N={des_entry['n']}): p_vs_baseline = {p_des_obs:.5f}")

des_names = set()
for r in csv.DictReader(open(RESULTS / "step_b77_resident_ledger.csv")):
    if r["survey"] == "DES":
        des_names.add(r["name"].strip())
from astropy.io import fits as _fits
_t = _fits.open(DATA_RAW / "des" / "y6_des_tnos_color.fits")[1].data
tnos_des = []
for row in _t:
    if str(row["MPC"]).strip() not in des_names:
        continue
    Om, w, i = float(row["lan"]), float(row["aop"]), float(row["i"])
    if np.isfinite(Om) and np.isfinite(w) and np.isfinite(i):
        tnos_des.append(dict(
            dir=perih_dir(math.radians(w), math.radians(Om),
                          math.radians(i)),
            varpi=math.radians((Om + w) % 360.0)))
logger.info(f"C2 DES objects matched in FITS: {len(tnos_des)}")

# C3: MPCORB Detached TNOs -- same per-object elements as step_092
mpcorb_res = json.load(open(RESULTS / "step_b57_mpcorb_clustering.json"))
p_mpcorb_obs = float(mpcorb_res["samples"]["detached_a150_q30"]["p_varpi_conditioned"])
logger.info(f"C3 MPCORB detached TNOs (N={mpcorb_res['samples']['detached_a150_q30']['N']}): p = {p_mpcorb_obs:.5f}")

import gzip as _gz
def _arc_years(a):
    if "-" in a:
        y0, y1 = a.split("-")
        return float(y1) - float(y0)
    try:
        return float(a) / 365.25
    except ValueError:
        return 0.0
tnos_mpc = []
with _gz.open(DATA_RAW / "mpc" / "MPCORB.DAT.gz", "rt",
              errors="replace") as f:
    for _i, line in enumerate(f):
        if _i < 43 or len(line) < 200:
            continue
        try:
            o = dict(Om=float(line[48:57]), w=float(line[37:46]),
                     i=float(line[59:68]), e=float(line[70:79]),
                     a=float(line[92:103]),
                     nopp=int(line[123:126])
                     if line[123:126].strip() else 0,
                     arc=line[128:136].strip(),
                     epoch=line[20:25].strip())
        except ValueError:
            continue
        if not o["epoch"].startswith("K"):
            continue
        o["q"] = o["a"] * (1.0 - o["e"])
        if not (o["a"] > 150 and o["q"] > 30 and o["e"] < 1.0):
            continue
        if not (o["nopp"] >= 2 or _arc_years(o["arc"]) >= 1.0):
            continue
        tnos_mpc.append(dict(
            dir=perih_dir(math.radians(o["w"]), math.radians(o["Om"]),
                          math.radians(o["i"])),
            varpi=math.radians((o["Om"] + o["w"]) % 360.0)))
logger.info(f"C3 MPCORB detached TNOs loaded: {len(tnos_mpc)}")

# C4: MPCORB Mean-Plane Lean Azimuth
mp = json.load(open(RESULTS / "step_b55_mean_plane.json"))
mp_bin = mp["bins"]["80-400"]
pole_obs = np.array(mp_bin["pole"], float)
pole_obs = pole_obs / np.linalg.norm(pole_obs)
if pole_obs[2] < 0:
    pole_obs = -pole_obs
az_obs = math.degrees(math.atan2(pole_obs[1], pole_obs[0])) % 360.0
null_poles = np.array(mp_bin["null_poles"], float)
null_poles /= np.linalg.norm(null_poles, axis=1, keepdims=True)
flip = null_poles[:, 2] < 0
null_poles[flip] *= -1
az_null = np.degrees(np.arctan2(null_poles[:, 1], null_poles[:, 0])) % 360.0
NMP = len(az_null)

def p_mp(axis):
    lam = math.degrees(math.atan2(axis[1], axis[0])) % 360.0
    d_obs = min(abs(az_obs - lam), 360.0 - abs(az_obs - lam))
    d_nul = np.minimum(np.abs(az_null - lam), 360.0 - np.abs(az_null - lam))
    return monte_carlo_p(int((d_nul <= d_obs).sum()), NMP)

# ------------------------------------------------------------------
# 2. Load Cometary Channels (C5: CODE, C6: Warsaw, C7: One-App, C8: CometEls, C9: SBDB Post-2017)
# ------------------------------------------------------------------

# C5 & C6: CODE and Warsaw aphelia and slips
from scripts.utils.tep9_common import parse_code

def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(desig=line[5:17].strip(),
                             com=line[3].strip(),
                             q=float(line[42:56]),
                             w=float(line[70:82]), Om=float(line[82:94]),
                             i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows

PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
def dedup(rows):
    out = {}
    for r in rows:
        k = r["desig"]
        if k not in out or PREF.get(r["com"], 9) < PREF.get(out[k]["com"], 9):
            out[k] = r
    return out

code_o = parse_code(DATA_RAW / "code" / "code_original.html")
war_r  = dedup(parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")))

slip_rows = list(csv.DictReader(open(RESULTS / "step_b30_proper_time_slip.csv")))
comets_slip = []
for r in slip_rows:
    d = None
    if r["cohort"] == "code" and r["desig"] in code_o:
        e = code_o[r["desig"]]
        d = -perih_dir(math.radians(e["w"]), math.radians(e["Om"]), math.radians(e["i"]))
    elif r["cohort"] == "warsaw" and r["desig"] in war_r:
        e = war_r[r["desig"]]
        d = -perih_dir(math.radians(e["w"]), math.radians(e["Om"]), math.radians(e["i"]))
    if d is not None:
        comets_slip.append(dict(
            desig=r["desig"],
            cohort=r["cohort"],
            dir=d,
            q=float(r["q"]),
            dtau=float(r["dtau_unexplained"]),
            dtau_in=float(r["dtau_inbound_unexplained"]) if r.get("dtau_inbound_unexplained") else np.nan
        ))
logger.info(f"Loaded {len(comets_slip)} comets with aphelion directions from step_b30 (CODE & Warsaw)")

# C8: MPC CometEls Census Aphelion Dipole
# Same construction as step_108 (b72): fragment records are merged to
# their parent designation so split-comet fragments are not counted as
# independent draws.
import re
cometels_aph = []
seen_desig = set()
for line in open(DATA_RAW / "mpc" / "CometEls.txt"):
    p = line.split()
    if len(p) < 12:
        continue
    off = 0
    try:
        int(p[1])
    except ValueError:
        off = 1
    try:
        q, e, w, Om, i = (float(p[4 + off]), float(p[5 + off]),
                          float(p[6 + off]), float(p[7 + off]),
                          float(p[8 + off]))
    except (ValueError, IndexError):
        continue
    name = " ".join(p[11 + off:])
    m = re.search(r'([CPD]/\d{4}[A-Z]+\d*|\d{4}[A-Z]+\d*)',
                  name.replace(' ', ''))
    desig = m.group(1).split('-')[0] if m else p[0]
    if desig in seen_desig:
        continue
    seen_desig.add(desig)
    if 0.90 <= e < 1.02:
        ph = perih_dir(math.radians(w), math.radians(Om), math.radians(i))
        cometels_aph.append(-ph)
cometels_aph = np.array(cometels_aph)
logger.info(f"C8 MPC CometEls near-parabolic census: {len(cometels_aph)} comets")

# C9: SBDB Post-2017 Aphelion Dipole
post2017_rows = list(csv.DictReader(open(RESULTS / "step_b88_lineage_aphelion.csv")))
post2017_aph = []
for r in post2017_rows:
    if r.get("era") == "post2017":
        post2017_aph.append(lv(float(r["aph_lam"]), float(r["aph_bet"])))
post2017_aph = np.array(post2017_aph)
logger.info(f"C9 SBDB post-2017 near-parabolic holdout: {len(post2017_aph)} comets")

# C10: Inward-Injected Jupiter-Family Comets (JFCs)
jfc_data = json.load(open(RESULTS / "step_b17_injection_chain.json"))
deep_jfc_mean_varpi = float(jfc_data["I4_chain"]["JFC_q<3"]["mu"])
deep_jfc_p = float(jfc_data["I3_pairing_test_q<3"]["p"])
logger.info(f"C10 Deep JFC mean varpi: {deep_jfc_mean_varpi:.1f} deg (p = {deep_jfc_p:.4f})")

# ------------------------------------------------------------------
# 3. Channel Evaluation Function over an Arbitrary Sky Axis
# ------------------------------------------------------------------

cdir  = np.array([c["dir"] for c in comets_slip])
cdtau = np.array([c["dtau"] for c in comets_slip])
cq    = np.array([c["q"] for c in comets_slip])
qm    = np.array([c["cohort"] == "code" for c in comets_slip]) & (cq < 3.1)
crank = rankdata(cdtau[qm])
iswar = np.array([c["cohort"] == "warsaw" for c in comets_slip])
widx  = np.where(iswar)[0]
wdin  = np.array([c["dtau_in"] for c in comets_slip])[widx]
wrank = rankdata(wdin)

tdir_sbdb   = np.array([t["dir"] for t in tnos_sbdb])
tvarpi_sbdb = np.array([t["varpi"] for t in tnos_sbdb])
texp_sbdb   = np.exp(1j * tvarpi_sbdb)
tdir_des    = np.array([t["dir"] for t in tnos_des])
texp_des    = np.exp(1j * np.array([t["varpi"] for t in tnos_des]))
tdir_mpc    = np.array([t["dir"] for t in tnos_mpc])
texp_mpc    = np.exp(1j * np.array([t["varpi"] for t in tnos_mpc]))

def mwz(rank_arr, mask):
    n1 = int(mask.sum()); n2 = int((~mask).sum())
    if n1 < 3 or n2 < 3:
        return 0.0
    U = rank_arr[mask].sum() - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    return (U - mu) / sd

def eval_channels(axis):
    """Computes p-values for all 10 channel scores around trial axis."""
    cth = np.degrees(np.arccos(np.clip(cdir @ axis, -1, 1)))

    # C1: SBDB Detached TNO in-cap varpi concentration
    tth = np.degrees(np.arccos(np.clip(tdir_sbdb @ axis, -1, 1)))
    tinc = tth < CAP
    ni = int(tinc.sum())
    if ni >= 3:
        R = abs(texp_sbdb[tinc].sum()) / ni
        p1 = math.exp(-ni * R * R)
    else:
        p1 = 1.0

    # C2: DES Detached TNOs -- identical per-axis in-cap Rayleigh
    # machinery as C1, on the survey's own per-object elements
    ax_lon = math.degrees(math.atan2(axis[1], axis[0])) % 360.0
    tth2 = np.degrees(np.arccos(np.clip(tdir_des @ axis, -1, 1)))
    inc2 = tth2 < CAP
    n2 = int(inc2.sum())
    if n2 >= 3:
        R2 = abs(texp_des[inc2].sum()) / n2
        p2 = math.exp(-n2 * R2 * R2)
    else:
        p2 = 1.0

    # C3: MPCORB Detached TNOs -- identical per-axis machinery
    tth3 = np.degrees(np.arccos(np.clip(tdir_mpc @ axis, -1, 1)))
    inc3 = tth3 < CAP
    n3 = int(inc3.sum())
    if n3 >= 3:
        R3 = abs(texp_mpc[inc3].sum()) / n3
        p3 = math.exp(-n3 * R3 * R3)
    else:
        p3 = 1.0

    # C4: MPCORB Mean-plane lean azimuth
    p4 = p_mp(axis)

    # C5: CODE matched slip contrast (in-cap vs out-of-cap)
    inc_m = cth[qm] < CAP
    z5 = mwz(crank, inc_m)
    p5 = norm.sf(z5)

    # C6: Warsaw inbound slip
    winc = cth[widx] < CAP
    z6 = mwz(wrank, winc)
    p6 = norm.sf(z6)

    # C7: Pooled comet aphelion projection (dipole)
    proj = cdir @ axis
    z7 = proj.mean() * math.sqrt(3 * len(proj))
    p7 = norm.sf(z7)

    # C8: MPC CometEls Census Aphelion Dipole
    proj8 = cometels_aph @ axis
    z8 = proj8.mean() * math.sqrt(3 * len(proj8))
    p8 = norm.sf(z8)

    # C9: SBDB Post-2017 Prospective Aphelion Dipole
    proj9 = post2017_aph @ axis
    z9 = proj9.mean() * math.sqrt(3 * len(proj9))
    p9 = norm.sf(z9)

    # C10: Inward-injected JFC antipodal alignment (deep JFC mean varpi vs axis lon)
    d_jfc = min(abs(deep_jfc_mean_varpi - ax_lon),
                360.0 - abs(deep_jfc_mean_varpi - ax_lon))
    p10 = min(1.0, max(1e-10, 2.0 * d_jfc / 360.0))

    p_vec = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10])
    return np.clip(p_vec, 1e-300, 1.0)

def S_stat(p_vec):
    return float(-2.0 * np.log(p_vec).sum())

# ------------------------------------------------------------------
# 4. Measure at Declared Axes & Run 20,000-Draw Random Permutation
# ------------------------------------------------------------------

p_tno_obs = eval_channels(AX_TNO)
S_tno_obs = S_stat(p_tno_obs)

p_com_obs = eval_channels(AX_COM)
S_com_obs = S_stat(p_com_obs)

logger.info(f"S10 at TNO axis (49,-17)   = {S_tno_obs:.2f}")
logger.info(f"  Channel p-values: {np.round(p_tno_obs, 6)}")
logger.info(f"S10 at Comet axis (34,-13) = {S_com_obs:.2f}")
logger.info(f"  Channel p-values: {np.round(p_com_obs, 6)}")

# Draw 20,000 uniform random axes on the sphere
u = rng.uniform(-1, 1, NAX)
ph = rng.uniform(0, 2 * math.pi, NAX)
rr = np.sqrt(1 - u * u)
rand_axes = np.column_stack([rr * np.cos(ph), rr * np.sin(ph), u])

S_rand = np.empty(NAX)
ch_rand = np.empty((NAX, 10))

for j, a in enumerate(rand_axes):
    p = eval_channels(a)
    ch_rand[j] = p
    S_rand[j] = S_stat(p)
    if (j + 1) % 5000 == 0:
        logger.info(f"  Processed {j + 1}/{NAX} random axes")

p_global_tno = float((1 + (S_rand >= S_tno_obs).sum()) / (NAX + 1))
p_global_com = float((1 + (S_rand >= S_com_obs).sum()) / (NAX + 1))

logger.info(f"Conditional sky-area rank at TNO axis:   {p_global_tno:.6f}")
logger.info(f"Conditional sky-area rank at Comet axis: {p_global_com:.6f}")

# Channel names
ch_names = [
    "C1_sbdb_tno_varpi",
    "C2_des_tno_cluster",
    "C3_mpcorb_tno_cluster",
    "C4_mean_plane_lean",
    "C5_code_comet_slip",
    "C6_warsaw_inbound_slip",
    "C7_pooled_comet_dipole",
    "C8_cometels_census_dipole",
    "C9_post2017_sbdb_dipole",
    "C10_jfc_injection_chain"
]

per_channel_results = {}
for j, name in enumerate(ch_names):
    obs_p = float(p_tno_obs[j])
    frac_better = float((ch_rand[:, j] <= obs_p).mean())
    per_channel_results[name] = {
        "p_at_tno_axis": obs_p,
        "frac_random_axes_better": frac_better
    }

best_idx = np.argmax(S_rand)
best_axis = rand_axes[best_idx]
best_sep = sep(best_axis, AX_TNO)

# ------------------------------------------------------------------
# 4b. Axis-localization audit -- is the landscape argmax a distinct
#     direction or the same sector as the two pre-declared axes?
# ------------------------------------------------------------------
# The declared axes are registered directions, not fitted parameters,
# so the honest global number is the max-over-axes permutation p at
# each declared axis.  A referee will still ask whether the data's
# own argmax (which lands ~15 deg off the TNO axis) indicates a
# displaced structure.  Three quantities answer that: the argmax's
# geometry relative to BOTH declared axes, the angular coherence
# scale of the S landscape, and the bootstrap scatter of the argmax
# under member resampling.

# (a) argmax geometry: is it bracketed by the two declared axes?
best_lon, best_lat = lb(best_axis)
sep_to_comet_axis = sep(best_axis, AX_COM)
sep_declared_axes = sep(AX_TNO, AX_COM)

# (b) exceedance fraction and coherence profile
frac_ge_declared = float((S_rand >= S_tno_obs).mean())
th_from_tno = np.degrees(np.arccos(
    np.clip(rand_axes @ AX_TNO, -1, 1)))
profile_bins = [(0, 15), (15, 30), (30, 45), (45, 60),
                (60, 90), (90, 180)]
s_profile = {}
for lo, hi in profile_bins:
    m = (th_from_tno >= lo) & (th_from_tno < hi)
    s_profile[f"{lo}-{hi}deg"] = {
        "n": int(m.sum()),
        "med_S": float(np.median(S_rand[m]))}

# (c) bootstrap argmax scatter: resample channel members, refit the
#     local argmax on a sector grid.  Measures the axis-estimation
#     noise of the landscape; if the observed argmax offset lies
#     inside the bootstrap cone it is consistent with the declared
#     direction.
def eval_channels_res(axis, idx):
    """eval_channels with member-resampled channel inputs."""
    # C1: SBDB TNOs
    t1 = tdir_sbdb[idx["c1"]]
    x1 = texp_sbdb[idx["c1"]]
    tth = np.degrees(np.arccos(np.clip(t1 @ axis, -1, 1)))
    inc = tth < CAP
    ni = int(inc.sum())
    p1 = math.exp(-ni * (abs(x1[inc].sum()) / ni) ** 2) if ni >= 3 else 1.0
    # C2: DES
    t2 = tdir_des[idx["c2"]]
    x2 = texp_des[idx["c2"]]
    tth = np.degrees(np.arccos(np.clip(t2 @ axis, -1, 1)))
    inc = tth < CAP
    ni = int(inc.sum())
    p2 = math.exp(-ni * (abs(x2[inc].sum()) / ni) ** 2) if ni >= 3 else 1.0
    # C3: MPCORB
    t3 = tdir_mpc[idx["c3"]]
    x3 = texp_mpc[idx["c3"]]
    tth = np.degrees(np.arccos(np.clip(t3 @ axis, -1, 1)))
    inc = tth < CAP
    ni = int(inc.sum())
    p3 = math.exp(-ni * (abs(x3[inc].sum()) / ni) ** 2) if ni >= 3 else 1.0
    # C4: mean-plane lean, resampled null pole distribution
    lam = math.degrees(math.atan2(axis[1], axis[0])) % 360.0
    az_n = np.degrees(np.arctan2(
        null_poles[idx["c4"], 1], null_poles[idx["c4"], 0])) % 360.0
    d_obs = min(abs(az_obs - lam), 360.0 - abs(az_obs - lam))
    d_nul = np.minimum(np.abs(az_n - lam), 360.0 - np.abs(az_n - lam))
    p4 = monte_carlo_p(int((d_nul <= d_obs).sum()), len(az_n))
    # C5: CODE matched slip contrast
    ii = idx["c5"]
    cth = np.degrees(np.arccos(np.clip(cdir[ii] @ axis, -1, 1)))
    r5 = rankdata(cdtau[ii])
    z5 = mwz(r5, cth < CAP)
    p5 = norm.sf(z5)
    # C6: Warsaw inbound slip
    ii = idx["c6"]
    cth = np.degrees(np.arccos(np.clip(cdir[ii] @ axis, -1, 1)))
    r6 = rankdata(np.array([c["dtau_in"] for c in comets_slip])[ii])
    z6 = mwz(r6, cth < CAP)
    p6 = norm.sf(z6)
    # C7: pooled comet dipole
    proj = cdir[idx["c7"]] @ axis
    p7 = norm.sf(proj.mean() * math.sqrt(3 * len(proj)))
    # C8: CometEls census dipole
    proj = cometels_aph[idx["c8"]] @ axis
    p8 = norm.sf(proj.mean() * math.sqrt(3 * len(proj)))
    # C9: post-2017 dipole
    proj = post2017_aph[idx["c9"]] @ axis
    p9 = norm.sf(proj.mean() * math.sqrt(3 * len(proj)))
    # C10: deterministic JFC chain (no members to resample)
    d_jfc = min(abs(deep_jfc_mean_varpi - lam),
                360.0 - abs(deep_jfc_mean_varpi - lam))
    p10 = min(1.0, max(1e-10, 2.0 * d_jfc / 360.0))
    return np.clip(np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]),
                   1e-300, 1.0)


# local sector grid covering both declared axes
gl = np.arange(10.0, 80.0, 4.0)
gb = np.arange(-45.0, 10.0, 4.0)
grid_axes = np.array([lv(a, b) for a in gl for b in gb])
N_BOOT = 150
qmi = np.where(qm)[0]
widi = np.where(iswar)[0]
boot_seps = []
boot_dirs = []
for _b in range(N_BOOT):
    idx = {
        "c1": rng.integers(0, len(tdir_sbdb), len(tdir_sbdb)),
        "c2": rng.integers(0, len(tdir_des), len(tdir_des)),
        "c3": rng.integers(0, len(tdir_mpc), len(tdir_mpc)),
        "c4": rng.integers(0, len(null_poles), len(null_poles)),
        "c5": qmi[rng.integers(0, len(qmi), len(qmi))],
        "c6": widi[rng.integers(0, len(widi), len(widi))],
        "c7": rng.integers(0, len(cdir), len(cdir)),
        "c8": rng.integers(0, len(cometels_aph), len(cometels_aph)),
        "c9": rng.integers(0, len(post2017_aph), len(post2017_aph)),
    }
    Sg = np.array([S_stat(eval_channels_res(a, idx)) for a in grid_axes])
    ba = grid_axes[int(np.argmax(Sg))]
    boot_seps.append(sep(ba, AX_TNO))
    boot_dirs.append(ba)
    if (_b + 1) % 50 == 0:
        logger.info(f"  bootstrap {_b+1}/{N_BOOT}")

boot_seps = np.array(boot_seps)
r68 = float(np.percentile(boot_seps, 68))
r95 = float(np.percentile(boot_seps, 95))
argmax_within_68 = bool(best_sep <= r68)

axis_localization = {
    "argmax_lon_lat_deg": [float(best_lon), float(best_lat)],
    "argmax_sep_to_tno_deg": float(best_sep),
    "argmax_sep_to_comet_deg": float(sep_to_comet_axis),
    "declared_axes_sep_deg": float(sep_declared_axes),
    "frac_sphere_ge_declared": frac_ge_declared,
    "s_profile_vs_dist": s_profile,
    "bootstrap": {
        "n": N_BOOT,
        "r68_deg": r68,
        "r95_deg": r95,
        "argmax_within_r68": argmax_within_68,
        "med_sep_deg": float(np.median(boot_seps)),
    },
}
logger.info(
    f"axis localization: argmax at ({best_lon:.0f},{best_lat:.0f}), "
    f"{best_sep:.1f} deg off TNO axis, {sep_to_comet_axis:.1f} deg off "
    f"comet axis; declared axes are {sep_declared_axes:.1f} deg apart; "
    f"bootstrap r68={r68:.1f} deg")

results_payload = {
    "step": "step_126_global_cross_survey_synthesis",
    "inputs": ["data/raw/sbdb/sbdb_outer_ss.json",
               "data/raw/des/y6_des_tnos_color.fits",
               "data/raw/mpc/MPCORB.DAT.gz",
               "data/raw/mpc/CometEls.txt",
               "data/raw/code/code_original.html",
               "data/raw/warsaw/warsaw_tablec.dat"],
    "description": "Ten-channel directional localization on overlapping catalogues, evaluated against random test directions",
    "evidence_status": "fixed-catalogue directional localization",
    "calibration": "Fixed observed catalogues, random test direction. Not a null-data maximum-statistic test and not a global discovery p-value. Resident scores select an angular cap before applying an unrestricted Rayleigh reference and are not valid component p-values. Catalogue overlap is retained in the observed landscape but resampled-object bootstrap uncertainty does not preserve cross-catalogue identities.",
    "legacy_field_note": "p_global and channel_p are retained for compatibility; interpret as sky-area rank and descriptive weights only.",
    "seed": SEED,
    "n_axes_permuted": NAX,
    "cap_deg": CAP,
    "axes": {
        "tno_axis_deg": [49.0, -17.0],
        "comet_axis_deg": [34.0, -13.0]
    },
    "eval_at_tno_axis": {
        "S_stat": S_tno_obs,
        "p_global": p_global_tno,
        "channel_p": p_tno_obs.tolist()
    },
    "eval_at_comet_axis": {
        "S_stat": S_com_obs,
        "p_global": p_global_com,
        "channel_p": p_com_obs.tolist()
    },
    "per_channel": per_channel_results,
    "best_random_axis": {
        "S_max": float(S_rand[best_idx]),
        "sep_to_tno_deg": float(best_sep)
    },
    "null_quantiles": {
        "q50": float(np.percentile(S_rand, 50)),
        "q95": float(np.percentile(S_rand, 95)),
        "q99": float(np.percentile(S_rand, 99)),
        "q999": float(np.percentile(S_rand, 99.9))
    },
    "axis_localization": axis_localization
}

# persist the landscape for downstream localization audits
np.savez_compressed(RESULTS / "step_b90_s_landscape.npz",
                    axes=rand_axes, S=S_rand)

out_json = RESULTS / "step_b90_global_synthesis.json"
with open(out_json, "w") as f:
    json.dump(results_payload, f, indent=1)
logger.info(f"Wrote JSON: {out_json}")

# CSV summary
out_csv = RESULTS / "step_b90_global_synthesis.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["channel", "name", "p_at_tno_axis", "frac_axes_better"])
    for j, name in enumerate(ch_names):
        w.writerow([f"C{j+1}", name, f"{p_tno_obs[j]:.6e}", f"{per_channel_results[name]['frac_random_axes_better']:.6f}"])
logger.info(f"Wrote CSV: {out_csv}")

# ------------------------------------------------------------------
# 5. Publication Figure
# ------------------------------------------------------------------

fig, axes_plot = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.2, 1]})

# Left: Histogram of S over 20,000 random axes
ax0 = axes_plot[0]
ax0.hist(S_rand, bins=60, color="#566573", alpha=0.75, density=True, label=f"20,000 random sky axes")
ax0.axvline(S_tno_obs, color="#1C2E4A", lw=4.0, ls="-", alpha=0.45,
            zorder=3,
            label=f"Detached-TNO axis (S={S_tno_obs:.1f}, p={p_global_tno:.4f})")
ax0.axvline(S_com_obs, color="#b43b4e", lw=1.6, ls="--", zorder=4,
            label=f"Comet transit axis (S={S_com_obs:.1f}, p={p_global_com:.4f})")
ax0.set_xlabel("Descriptive score $S = -2 \\sum_{i=1}^{10} \\ln p_i$")
ax0.set_ylabel("Probability density")
ax0.legend()
ax0.grid(True, ls=":", alpha=0.5)

# Right: Forest plot of 10 channels at the declared axis
ax1 = axes_plot[1]
y_pos = np.arange(len(ch_names))
minus_log_p = -np.log10(p_tno_obs)
labels = [
    "C1: SBDB Detached TNOs",
    "C2: DES Y6 Detached TNOs",
    "C3: MPCORB Detached TNOs",
    "C4: MPCORB Mean-Plane Warp",
    "C5: CODE Comets Slip",
    "C6: Warsaw Inbound Slip",
    "C7: Pooled Comet Dipole",
    "C8: CometEls Census Dipole",
    "C9: Post-2017 SBDB Dipole",
    "C10: JFC Injection Chain"
]

# Channels grouped by population: TNO channels C1-C4 navy,
# comet channels C5-C10 accent red.
colors = ["#1C2E4A"] * 4 + ["#b43b4e"] * 6
ax1.barh(y_pos, minus_log_p, color=colors, alpha=0.9, height=0.65)
ax1.axvline(-np.log10(0.05), color="#b43b4e", ls=":", lw=1.5,
            label="$p = 0.05$ threshold")
ax1.set_yticks(y_pos)
ax1.set_yticklabels(labels)
ax1.invert_yaxis()
ax1.set_xlabel("$-\\log_{10}(p)$ at declared boundary axis")
ax1.grid(True, ls=":", alpha=0.5, axis="x")
import matplotlib.patches as mpatches
ax1.legend(handles=[
    mpatches.Patch(color="#1C2E4A", label="TNO channels"),
    mpatches.Patch(color="#b43b4e", label="comet channels"),
    plt.Line2D([0], [0], color="#b43b4e", ls=":", lw=1.5,
               label="$p = 0.05$ threshold")])

plt.tight_layout()
out_fig = RESULTS / "figures" / "step_b90_global_synthesis.png"
plt.savefig(out_fig, dpi=300)
plt.close()
logger.info(f"Saved figure: {out_fig}")
logger.success("Step 126 completed successfully.")
