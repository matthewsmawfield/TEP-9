#!/usr/bin/env python3
"""Step 090 -- bias-free mean-plane warp of the distant Kuiper belt
(Siraj, Chyba & Tremaine 2025 estimator) and the direction of the warp
====================================================================

Siraj, Chyba & Tremaine (2025, MNRAS Letters slaf091; arXiv:2508.14156)
introduced a likelihood mean-plane estimator that is, by construction,
independent of the survey footprint: conditioning on each object's
observed position r-hat removes the on-sky selection function w(r)
from the likelihood, leaving only the intrinsic orbital-plane term.
They report a warp of the distant Kuiper-belt mean plane relative to
the invariable plane in the 80-200 AU (2.5 sigma) and 80-400 AU
(2.7 sigma) semimajor-axis bins -- i0 ~ 15 deg toward node Omega0 ~
120 deg -- and show in their own n-body tests that NEITHER the
Brown-Batygin (2021) Planet Nine NOR their own Siraj+25 Planet-X
realization can generate it; a third, distinct perturber ('Planet Y',
~0.06-1 M_earth, a = 100-200 AU, i >~ 10 deg) would be required.

This step (i) re-implements the exact estimator on real MPCORB
elements -- an orbit-fit lineage independent of JPL SBDB -- (ii)
reproduces the warp in their semimajor-axis bins on two samples
(their published 46-object Table-1 list, and a pipeline-native
non-resonant cut), and (iii) asks the question their paper leaves
open: WHERE does the disk lean?  The mean-pole offset azimuth is
measured against the detached-TNO boundary axis (lam=49 deg) already
recovered by the resident and transit channels, and against the
published perturber planes.

Likelihood (their eq. 16)
-------------------------
    log L = gamma * m_hat . sum_i J_hat_i
            - sum_i log I0( gamma * sin(theta_i) ) + const

with m_hat the trial mean-plane pole, J_hat_i = r_i x v_i / |r_i x v_i|
the angular-momentum direction, and theta_i the angle between m_hat
and the observed position r_hat_i.  Delta log L is profiled over
gamma at the invariable plane (i_inv = 1.6 deg, Om_inv = 107 deg);
significance = sqrt(2 Delta log L) (their convention).  The false-
alarm probability is measured by Monte-Carlo: synthetic disks
distributed symmetrically about the invariable plane, drawn at the
OBSERVED position vectors (the von-Mises-on-the-tangent-circle
generative null of SCT25, conditioned on the real on-sky geometry so
the C_i denominators match the data exactly).

Positions are propagated from MPCORB osculating elements
(a, e, i, Om, w, M at the file epoch, K2669 = 2026-06-09 for >99 per
cent of rows) by two-body Kepler solution -- heliocentric ecliptic;
the heliocentric-vs-barycentric offset is far below the estimator's
resolution and is recorded as a caveat.

Non-resonant cut: objects within +/-1.5 per cent in semimajor axis of
a Neptune mean-motion resonance (p+q):p with p+q <= 12 are excluded --
an analytic proxy for the SBDynT 90-clone filter used by SCT25, run
at +/-1.0 and +/-2.0 per cent as sensitivity checks.

Outputs
-------
results/step_b55_mean_plane.json
results/figures/step_b55_mean_plane.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_090_mean_plane")
tee_stdout(logger)
logger.header("Bias-free mean-plane warp (SCT25 estimator, MPCORB lineage)")

import gzip
import json
import math
import numpy as np
from scipy.optimize import minimize
from scipy.special import i0

rng = np.random.default_rng(20260918)

# ------------------------------------------------------------------
# MPCORB fixed-width parser
# ------------------------------------------------------------------

MPCORB = DATA_RAW / "mpc" / "MPCORB.DAT.gz"

MONTH_C = "123456789ABC"
DAY_C = "123456789ABCDEFGHIJKLMNOPQRSTUV"

def unpack_epoch(s):
    """MPC packed epoch -> (year, month, day).  I=1800s,J=1900s,K=2000s."""
    cen = {"I": 1800, "J": 1900, "K": 2000}.get(s[0])
    if cen is None:
        return None
    return (cen + int(s[1:3]), MONTH_C.index(s[3]) + 1,
            DAY_C.index(s[4]) + 1)

def load_mpcorb():
    """All MPCORB rows with usable elements -> list of dicts."""
    objs = []
    with gzip.open(MPCORB, "rt", errors="replace") as f:
        for i, line in enumerate(f):
            if i < 43 or len(line) < 200:
                continue
            try:
                o = dict(
                    desig=line[0:7].strip(),
                    epoch=line[20:25].strip(),
                    M=float(line[26:35]), w=float(line[37:46]),
                    Om=float(line[48:57]), i=float(line[59:68]),
                    e=float(line[70:79]), a=float(line[92:103]),
                    nopp=int(line[123:126]) if line[123:126].strip() else 0,
                    arc=line[128:136].strip(),
                    rms=float(line[137:141]) if line[137:141].strip() else float("nan"),
                    name=line[166:194].strip())
            except ValueError:
                continue
            ep = unpack_epoch(o["epoch"])
            if ep is None:
                continue
            o["epoch_ymd"] = ep
            o["q"] = o["a"] * (1.0 - o["e"])
            objs.append(o)
    return objs

# ------------------------------------------------------------------
# Kepler propagation -> position / angular-momentum unit vectors
# ------------------------------------------------------------------

MU = 4 * math.pi ** 2          # AU^3 yr^-2

def kepler_E(M, e, tol=1e-12):
    E = M + e * math.sin(M)
    for _ in range(60):
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E

def state_vec(a, e, i, Om, w, M):
    """Osculating elements (deg; M rad internally) -> r, v (AU, AU/yr)."""
    E = kepler_E(M, e)
    cE, sE = math.cos(E), math.sin(E)
    x_, y_ = a * (cE - e), a * math.sqrt(1 - e * e) * sE
    fac = math.sqrt(MU * a) / (a * (1 - e * cE))
    vx_, vy_ = -fac * sE, fac * math.sqrt(1 - e * e) * cE
    cO, sO = math.cos(Om), math.sin(Om)
    co, so = math.cos(w), math.sin(w)
    ci, si = math.cos(i), math.sin(i)
    R = np.array([[cO * co - sO * so * ci, -cO * so - sO * co * ci, sO * si],
                  [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
                  [so * si, co * si, ci]])
    return R @ np.array([x_, y_, 0.0]), R @ np.array([vx_, vy_, 0.0])

def rhat_jhat(o):
    r, v = state_vec(o["a"], o["e"], math.radians(o["i"]),
                     math.radians(o["Om"]), math.radians(o["w"]),
                     math.radians(o["M"]))
    rhat = r / np.linalg.norm(r)
    j = np.cross(r, v)
    return rhat, j / np.linalg.norm(j)

# ------------------------------------------------------------------
# Neptune MMR exclusion (analytic proxy for the SCT25 clone filter)
# ------------------------------------------------------------------

A_NEP = 30.07

def mmr_centres(a_lo=30.0, a_hi=1000.0, max_sum=12):
    """Principal Neptune MMRs (p+q):p with p+q <= max_sum -- the
    resonances wide enough to trap appreciable populations."""
    out = []
    for p in range(1, max_sum):
        for q in range(1, max_sum - p + 1):
            ar = A_NEP * ((p + q) / p) ** (2.0 / 3.0)
            if a_lo < ar < a_hi:
                out.append(((p + q), p, ar))
    return sorted(out, key=lambda t: t[2])

MMRS = mmr_centres()

def near_mmr(a, frac):
    return any(abs(a - ar) / ar < frac for _, _, ar in MMRS)

# ------------------------------------------------------------------
# Published non-resonant sample (SCT25 Table 1, 80-400 AU bin)
# ------------------------------------------------------------------

SCT25_TABLE1 = ["2016 UP273", "2013 GZ136", "2008 JO41", "2003 YQ179",
                "1999 CF119", "2015 TV361", "2013 RK109", "2013 RJ158",
                "2012 HW87", "2020 KV11", "2011 BR163", "2014 YD50",
                "2010 ER65", "2000 OM67", "2013 VU71", "2015 GY55",
                "2015 RC279", "1999 RZ215", "2014 BE70", "2014 VM43",
                "2015 DW224", "2017 BZ236", "2019 SS149", "2013 RE124",
                "2014 MJ70", "2014 KZ101", "2015 RB279", "2014 JW80",
                "2013 JO64", "2003 HB57", "2015 SO20", "2007 VJ305",
                "2015 OC193", "2017 CL54", "2001 FP185", "2015 RY245",
                "2016 SG58", "2023 KQ14", "2012 VP113", "2015 DY248",
                "2015 GT50", "2013 SL102", "2014 TU115", "2013 FL28",
                "2013 RF98", "2014 SX403"]

def prov_desig(name):
    """Readable-name field -> trailing provisional designation."""
    return name.split(")")[-1].strip()

# ------------------------------------------------------------------
# SCT25 likelihood estimator
# ------------------------------------------------------------------

def pole(l_deg, b_deg):
    l, b = math.radians(l_deg), math.radians(b_deg)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l), math.sin(b)])

def pole_to_iOm(m):
    """Pole unit vector -> (i0, Omega0) of the corresponding plane."""
    i0 = math.degrees(math.acos(np.clip(m[2], -1, 1)))
    Om0 = math.degrees(math.atan2(m[0], -m[1])) % 360.0
    return i0, Om0

def logL(params, rh, jh, jsum):
    l, b, lg = params
    m = pole(l, b)
    g = math.exp(float(lg))
    sin_th = np.sqrt(np.clip(1.0 - (rh @ m) ** 2, 0.0, 1.0))
    return g * float(m @ jsum) - float(np.log(i0(g * sin_th)).sum())

def _starts(rh, jh, full):
    """Starting points: the mean-J direction first, plus a coarse
    grid when the full search is requested."""
    jm = jh.mean(axis=0); jm /= np.linalg.norm(jm)
    l0 = math.degrees(math.atan2(jm[1], jm[0])) % 360
    b0 = math.degrees(math.asin(np.clip(jm[2], -1, 1)))
    st = [[l0, b0, math.log(30.0)], [l0, -b0, math.log(30.0)]]
    if full:
        for l_ in np.arange(0, 360, 60):
            for b_ in (-45, 0, 45):
                for lg_ in (math.log(10.0), math.log(100.0)):
                    st.append([l_, b_, lg_])
    return st

def max_logL(rh, jh, full=True):
    """Maximize logL over (l, b, log gamma)."""
    jsum = jh.sum(axis=0)
    best = None
    for s0 in _starts(rh, jh, full):
        res = minimize(lambda p: -logL(p, rh, jh, jsum), s0,
                       method="Nelder-Mead",
                       options=dict(maxiter=1500, xatol=1e-5,
                                    fatol=1e-7))
        if best is None or res.fun < best.fun:
            best = res
    l, b, lg = best.x
    return pole(l, b), math.exp(float(lg)), -best.fun

def prof_logL_at(m_fixed, rh, jh):
    """Profile likelihood with m pinned to a fixed direction."""
    jsum = jh.sum(axis=0)
    sin_th = np.sqrt(np.clip(1.0 - (rh @ m_fixed) ** 2, 0.0, 1.0))
    lt = np.log(i0(np.linspace(0.0, 400.0, 4001)[None, :]
                   * sin_th[:, None])).sum(axis=0)
    gs = np.linspace(0.0, 400.0, 4001)
    ll = gs * float(m_fixed @ jsum) - lt
    return float(ll.max())

# ------------------------------------------------------------------
# Unwarped-disk Monte Carlo (SCT25 FAP construction)
# ------------------------------------------------------------------

def synth_disk(rh_obs, gamma, m_sym):
    """Conditional unwarped-disk draw: at each OBSERVED position r-hat,
    draw J-hat on the tangent unit circle with the von Mises weight
    exp(gamma * J.m_sym) -- the SCT25 generative null, conditioned on
    the real position geometry so the C_i denominators match."""
    n = len(rh_obs)
    jh = np.empty((n, 3))
    for k in range(n):
        r_ = rh_obs[k]
        while True:
            u = np.cross(r_, rng.normal(size=3))
            nu = np.linalg.norm(u)
            if nu > 1e-9:
                break
        u /= nu
        v = np.cross(r_, u)
        phi = rng.uniform(0, 2 * np.pi)
        j_ = u * math.cos(phi) + v * math.sin(phi)
        if rng.uniform() < math.exp(gamma * (float(j_ @ m_sym) - 1.0)):
            jh[k] = j_
        else:
            k_retry = 0
            while True:
                phi = rng.uniform(0, 2 * np.pi)
                j_ = u * math.cos(phi) + v * math.sin(phi)
                if rng.uniform() < math.exp(
                        gamma * (float(j_ @ m_sym) - 1.0)):
                    jh[k] = j_
                    break
                k_retry += 1
                if k_retry > 200:
                    jh[k] = j_
                    break
    return rh_obs, jh

def delta_logL(rh, jh, m_null, full=True):
    m_, g_, lmax = max_logL(rh, jh, full=full)
    return lmax - prof_logL_at(m_null, rh, jh), m_, g_

# ------------------------------------------------------------------
# Load catalogue and build samples
# ------------------------------------------------------------------

logger.info("loading MPCORB ...")
objs = load_mpcorb()
logger.info(f"  {len(objs)} orbits")

# propagate only objects that can enter any bin or the validation
# sample (a in 30-450 AU, q > 28) -- the rest cannot contribute
cand = [o for o in objs
        if 30.0 < o["a"] <= 450.0 and o["q"] > 28.0]
name_set = set(SCT25_TABLE1)
for o in objs:
    if prov_desig(o["name"]) in name_set and o not in cand:
        cand.append(o)
logger.info(f"  {len(cand)} candidates in 30-450 AU, q>28")
for o in cand:
    o["rhat"], o["jhat"] = rhat_jhat(o)

# invariable plane (SCT25 convention)
M_INV = np.array([math.sin(math.radians(1.6)) * math.sin(math.radians(107.0)),
                  -math.sin(math.radians(1.6)) * math.cos(math.radians(107.0)),
                  math.cos(math.radians(1.6))])
TNO_AXIS = np.array([math.cos(math.radians(-17.0)) * math.cos(math.radians(49.0)),
                     math.cos(math.radians(-17.0)) * math.sin(math.radians(49.0)),
                     math.sin(math.radians(-17.0))])

def arc_years(o):
    a = o["arc"]
    if "-" in a:                      # multi-opposition "YYYY-YYYY"
        try:
            y0, y1 = a.split("-")
            return float(y1) - float(y0)
        except ValueError:
            return 0.0
    try:                              # single-opposition arc in days
        return float(a) / 365.25
    except ValueError:
        return 0.0

def select(lo, hi, mmr_frac=0.015):
    """Pipeline-native non-resonant cut: a in (lo,hi], q>30,
    multi-opposition or >=1 yr arc, outside all MMR bands."""
    return [o for o in cand
            if lo < o["a"] <= hi and o["q"] > 30.0
            and (o["nopp"] >= 2 or arc_years(o) >= 1.0)
            and not near_mmr(o["a"], mmr_frac)]

def measure(sample, label, m_null=M_INV, n_mc=400):
    rh = np.array([o["rhat"] for o in sample])
    jh = np.array([o["jhat"] for o in sample])
    if len(sample) < 8:
        logger.warning(f"{label}: only {len(sample)} objects -- skipped")
        return None
    m_, g_, lmax = max_logL(rh, jh)
    lnull = prof_logL_at(m_null, rh, jh)
    dl = lmax - lnull
    i0, Om0 = pole_to_iOm(m_)
    # FAP Monte Carlo on the unwarped-disk null (SCT25 construction),
    # conditioned on the observed position geometry.  The fitted pole
    # azimuth of every null draw is kept so the lean direction can be
    # tested against the lambda=49 deg boundary axis.
    dls = np.empty(n_mc)
    azs = np.empty(n_mc)
    null_poles = np.empty((n_mc, 3))
    for k in range(n_mc):
        rs, js = synth_disk(rh, 20.0, m_null)
        dls[k], m_s, _ = delta_logL(rs, js, m_null, full=False)
        if m_s[2] < 0:
            m_s = -m_s
        null_poles[k] = m_s
        azs[k] = math.degrees(math.atan2(m_s[1], m_s[0])) % 360.0
    fap = float((int((dls >= dl).sum()) + 1) / (n_mc + 1))
    m_obs = m_ if m_[2] >= 0 else -m_
    az_obs = math.degrees(math.atan2(m_obs[1], m_obs[0])) % 360.0
    d_obs = min(abs(az_obs - 49.0), 360.0 - abs(az_obs - 49.0))
    d_nul = np.minimum(np.abs(azs - 49.0), 360.0 - np.abs(azs - 49.0))
    out = dict(n=len(sample), i0_deg=round(i0, 2), Om0_deg=round(Om0, 1),
               gamma=round(g_, 1), delta_logL=round(float(dl), 3),
               sigma=round(float(math.sqrt(2 * max(dl, 0.0))), 2),
               fap_mc=round(fap, 4), n_mc=n_mc,
               pole=list(np.round(m_, 4)),
               null_poles=np.round(null_poles, 4).tolist(),
               lean_azimuth_deg=round(az_obs, 1),
               lean_az_vs_49deg=round(d_obs, 1),
               p_lean_toward_49=round(float(
                   (int((d_nul <= d_obs).sum()) + 1)
                   / (n_mc + 1)), 4))
    logger.info(f"{label}: N={len(sample)}  i0={i0:.1f}  Om0={Om0:.0f}  "
                f"gamma={g_:.0f}  dlogL={dl:.2f}  "
                f"sigma={out['sigma']:.2f}  FAP={fap:.4f}")
    return out

# ------------------------------------------------------------------
# A. published-sample validation (SCT25 Table 1 -> MPCORB)
# ------------------------------------------------------------------

logger.subheader("A. SCT25 Table-1 published-sample validation")
name_map = {}
for o in cand:
    pd = prov_desig(o["name"])
    if pd:
        name_map[pd] = o
pub = [name_map[d] for d in SCT25_TABLE1 if d in name_map]
missing = [d for d in SCT25_TABLE1 if d not in name_map]
logger.info(f"matched {len(pub)}/{len(SCT25_TABLE1)} published objects; "
            f"missing: {missing}")
res_pub = measure(pub, "SCT25 Table-1 (80-400 AU)")

# ------------------------------------------------------------------
# B. pipeline-native non-resonant bins
# ------------------------------------------------------------------

logger.subheader("B. Pipeline-native bins (MMR-excluded)")
BINS = [(50, 80), (80, 200), (200, 400), (80, 400)]
res_bins = {}
for lo, hi in BINS:
    s = select(lo, hi)
    res_bins[f"{lo}-{hi}"] = measure(s, f"a {lo}-{hi} AU")

# sensitivity to the MMR half-width
sens = {}
for fr in (0.010, 0.020):
    s = select(80, 400, mmr_frac=fr)
    sens[f"mmr_{fr}"] = measure(s, f"80-400 AU, MMR cut {fr}", n_mc=200)

# ------------------------------------------------------------------
# C. warp direction vs the boundary axis
# ------------------------------------------------------------------

logger.subheader("C. Direction of the warp")
def sep_deg(u, v):
    return math.degrees(math.acos(np.clip(float(u @ v), -1, 1)))

def az_of(m):
    return math.degrees(math.atan2(m[1], m[0])) % 360.0

dir_tests = {}
for key, r in list(res_bins.items()) + ([("Table1", res_pub)] if res_pub else []):
    if r is None:
        continue
    m = np.array(r["pole"])
    if m[2] < 0:
        m = -m
    lean = r["lean_azimuth_deg"]
    seps = dict(
        pole_vs_tno_axis_deg=round(sep_deg(m, TNO_AXIS), 1),
        lean_azimuth_deg=lean,
        lean_vs_axis_lon49_deg=r["lean_az_vs_49deg"],
        p_lean_toward_49=r["p_lean_toward_49"],
        node_Om0_vs_BB21_Om9_96d9_deg=round(
            min(abs(r["Om0_deg"] - 96.9),
                360 - abs(r["Om0_deg"] - 96.9)), 1))
    dir_tests[key] = seps
    logger.info(f"{key}: lean az {lean:.0f} deg, "
                f"|az-49|={r['lean_az_vs_49deg']} deg "
                f"(p={r['p_lean_toward_49']}), "
                f"Om0-P9node {seps['node_Om0_vs_BB21_Om9_96d9_deg']} deg")

# leave-one-out stability of the fitted lean azimuth -- is the
# direction carried by a handful of objects?
jack = {}
for key in ("80-400", "200-400"):
    lo, hi = (int(v) for v in key.split("-"))
    s = select(lo, hi)
    if len(s) < 12:
        continue
    azs_j = np.empty(len(s))
    for k in range(len(s)):
        sub = s[:k] + s[k + 1:]
        rh = np.array([o["rhat"] for o in sub])
        jh = np.array([o["jhat"] for o in sub])
        m_, g_, _ = max_logL(rh, jh, full=False)
        if m_[2] < 0:
            m_ = -m_
        azs_j[k] = az_of(m_)
    ph = np.deg2rad(azs_j)
    R = float(abs(np.exp(1j * ph).mean()))
    csd = math.degrees(math.sqrt(-2.0 * math.log(max(R, 1e-12))))
    d49 = np.minimum(np.abs(azs_j - 49.0), 360.0 - np.abs(azs_j - 49.0))
    jack[key] = dict(n=len(s), circ_sd_deg=round(csd, 1),
                     med_az_deg=round(float(np.median(azs_j)), 1),
                     frac_within_30deg_of_49=round(float((d49 <= 30).mean()), 3))
    logger.info(f"{key} jackknife: median az "
                f"{jack[key]['med_az_deg']:.0f} deg, circ sd {csd:.1f} deg, "
                f"{jack[key]['frac_within_30deg_of_49']:.2f} within 30 deg "
                f"of 49")

# ------------------------------------------------------------------
# D. internal cross-check: footprint-aware mean pole (step-018 style)
# ------------------------------------------------------------------

logger.subheader("D. Mean-pole cross-check (same bins)")
def pole_vec(om, i):
    n = np.stack([np.sin(i) * np.sin(om), -np.sin(i) * np.cos(om),
                  np.cos(i)], axis=-1)
    n[n[..., 2] < 0] *= -1
    return n

def mean_pole_az(om, i):
    n = pole_vec(np.asarray(om), np.asarray(i))
    m = n.mean(axis=0); m /= np.linalg.norm(m)
    return math.degrees(math.acos(np.clip(m[2], -1, 1))), az_of(m)

xc = {}
for lo, hi in BINS:
    s = select(lo, hi)
    if len(s) < 8:
        continue
    off, az = mean_pole_az(np.deg2rad([o["Om"] for o in s]),
                           np.deg2rad([o["i"] for o in s]))
    xc[f"{lo}-{hi}"] = dict(n=len(s), naive_pole_tilt_deg=round(off, 2),
                            naive_pole_az_deg=round(az, 1))
    logger.info(f"{lo}-{hi} AU: naive mean-pole tilt {off:.2f} deg, "
                f"azimuth {az:.0f} deg")

# ------------------------------------------------------------------
# write
# ------------------------------------------------------------------

out = dict(meta=dict(
    method="SCT25 likelihood mean-plane estimator on MPCORB osculating "
           "elements (independent orbit-fit lineage; M at epoch K2669 "
           "= 2026-06-09 -> Kepler-propagated r-hat, J-hat)",
    likelihood="logL = gamma*m.sum(J) - sum log I0(gamma*sin theta_i); "
               "Delta logL profiled over gamma at the invariable pole",
    invariable=dict(i_inv_deg=1.6, Om_inv_deg=107.0),
    mmr_exclusion="+/-1.5% in a around Neptune (p+q):p MMRs, p+q<=12 "
                  "(analytic proxy for the SCT25 90-clone SBDynT filter); "
                  "sensitivity at 1.0/2.0%",
    quality="nopp>=2 or arc>=1 yr; q>30 AU in the native bins",
    fap_mc="unwarped disks about the invariable plane drawn at the "
           "observed r-hat positions, J-hat on the tangent circle with "
           "von Mises weight exp(gamma*J.m_inv), gamma=20 (SCT25 "
           "generative null conditioned on the real on-sky geometry); "
           "per-draw fitted-pole azimuths test the lean toward "
           "lambda=49 deg",
    caveat="heliocentric osculating elements (MPCORB) vs barycentric "
           "(SCT25/SBDB); epochs not perfectly synchronous for a small "
           "minority of rows",
    reference="Siraj, Chyba & Tremaine 2025, MNRAS Lett. slaf091 "
              "(arXiv:2508.14156)"),
    published_sample=dict(n_matched=len(pub), missing=missing,
                          result=res_pub,
                          sct25_published=dict(i0_deg=15.0, Om0_deg=120.0,
                                               sigma=2.74)),
    bins=res_bins, mmr_sensitivity=sens,
    warp_direction=dir_tests, lean_jackknife=jack,
    mean_pole_crosscheck=xc)

RES = RESULTS
json.dump(out, open(RES / "step_b55_mean_plane.json", "w"), indent=1,
          default=float)
print("wrote", RES / "step_b55_mean_plane.json")

# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))

ax = axes[0]
labels = [k for k in res_bins if res_bins[k] is not None]
i0s = [res_bins[k]["i0_deg"] for k in labels]
sgs = [res_bins[k]["sigma"] for k in labels]
ax.bar(range(len(labels)), i0s, color="steelblue", edgecolor="k", lw=0.4)
for j, s in enumerate(sgs):
    ax.text(j, i0s[j] + 0.3, f"{s:.1f}$\\sigma$", ha="center", fontsize=8)
ax.axhline(1.6, color="r", ls="--", lw=1, label="invariable $i_{inv}=1.6°$")
if res_pub:
    ax.axhline(res_pub["i0_deg"], color="purple", ls=":", lw=1,
               label=f"Table-1 fit $i_0$={res_pub['i0_deg']}°")
ax.set_xticks(range(len(labels)))
ax.set_xticklabels([f"{k}\n(N={res_bins[k]['n']})" for k in labels],
                   fontsize=8)
ax.set(ylabel="best-fit mean-plane tilt $i_0$ (deg)",
       title="warp amplitude per bin\n(SCT25 likelihood, MPCORB)")
ax.legend(fontsize=7)

ax = axes[1]
for k in labels:
    m = np.array(res_bins[k]["pole"])
    if m[2] < 0:
        m = -m
    ax.plot(az_of(m), math.degrees(math.acos(np.clip(m[2], -1, 1))),
            "o", ms=8, label=f"{k} AU")
ax.plot(49.0, 73.0, "*", ms=14, color="crimson",
        label="TNO axis poleward offset")
ax.plot(az_of(M_INV), math.degrees(math.acos(np.clip(M_INV[2], -1, 1))),
        "x", ms=9, color="k", label="invariable pole")
ax.set(xlabel="pole azimuth (deg)", ylabel="pole offset (deg)",
       title="warp pole direction\n(lean azimuth vs 49° axis)")
ax.set_xlim(0, 360); ax.legend(fontsize=7)

ax = axes[2]
for k in labels:
    ax.bar(k, res_bins[k]["fap_mc"], color="teal", edgecolor="k", lw=0.4)
ax.axhline(0.05, color="r", ls="--", lw=1, label="5% FAP")
ax.set(ylabel="false-alarm probability", title="unwarped-disk Monte Carlo\n(SCT25 construction)")
ax.legend(fontsize=7)

fig.tight_layout()
FIG = RES / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b55_mean_plane.png", dpi=150)
print(f"wrote {FIG / 'step_b55_mean_plane.png'}")
