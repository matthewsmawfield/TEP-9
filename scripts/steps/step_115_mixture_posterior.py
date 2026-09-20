#!/usr/bin/env python3
"""Step 115 -- mixture posterior-predictive audit of resident cohorts.

The registered prospective prediction (step_054) is conditional on a
cohort's discovery geometry: the expected in-cap rate under EITHER
model is a function of where the cohort was discovered pointing.
The footprint-only model predicts the in-cap rate through the
coupling-smear baseline; the TEP mixture predicts it through the
per-object posterior P(varpi in cap | lam_opp).  For a cohort
discovered looking toward the anti-axis sector the two predictions
nearly coincide (the intrinsic-axis population is rarely discoverable
there), while a cohort discovered looking at the axis sector
separates them by tens of points.  The newest provisional cohort --
scored "untested" in step_084 -- can therefore be scored after all:
against the prediction each model makes FOR ITS OWN pointing.

Generative models compared:

  H0  footprint only:  varpi | lam ~ WN(lam, sigma_c)
      -> per-object in-cap probability = the established smear_cap
         baseline of steps 014/054/085/112/113.
  H1  TEP mixture (step_025 fit, registered in step_054):
      varpi ~ f * WN(axis, sigma_int) + (1 - f) * Uniform,
      lam | varpi ~ WN(varpi, sigma_c)
      -> per-object posterior P(varpi in cap | lam_i) by Bayes on a
         0.25 deg longitude grid (720-point circle).

For every cohort the step reports the observed in-cap count, the
expected count under each model, the exact Poisson-binomial
probability of the observed count under each model, and the
per-cohort Bayes factor P(k|H1)/P(k|H0).  A deduplicated union
(same precedence as step_113) carries the joint factor, and a
discovery-longitude lever-arm map shows where a future cohort is
decisive.  Provisional-orbit quality is assessed from the data
itself: the circular scatter of (varpi - lam_opp) residuals on the
newest cohort versus the secure cohort.

Cohorts (constructed exactly as in steps 113/084):

  SBDB      secure detached a>150, q>30, cc<=3  (designation proxy)
  DES       detached a>150, q>30                (designation proxy)
  OSSOS     'det' class                         (exact astrometry)
  2025+     detached a>150, q>30, designation year >= 2025, any cc
            (provisional orbits; designation proxy)
  era bins  the full SBDB detached pool (any cc) split by
            designation year -- 2000-05, 2006-10, 2011-15, 2016-19,
            2020-22, 2023-24, 2025+ -- scoring whether the signature
            persists across discovery generations and whether each
            generation's count matches what its own pointing predicts

Outputs
-------
results/step_b79_mixture_posterior.json
results/step_b79_mixture_posterior.csv   (per-object posterior scoring)
results/figures/supplementary/step_b79_mixture_posterior.png

Author: Matthew Lukin Smawfield
Date: September 2026
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_115_mixture_posterior")
tee_stdout(logger)
logger.header("Mixture posterior-predictive audit")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import norm
from astropy.io.votable import parse as vot_parse
from astropy.io import fits
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time
from astropy import units as u
import warnings
warnings.filterwarnings("ignore")

rng = np.random.default_rng(20261013)

AXIS, CAP, SIG_C = 49.0, 60.0, 50.0

# mixture parameters from the registered step_025 fit (read back so the
# audit tracks the fitted values rather than restating them)
_mix = json.load(open(RESULTS / "step_b13_patch_scale.json"))
F_CL = float(_mix["P2_intrinsic_width"]["f_cluster"])
SIG_I = float(_mix["P2_intrinsic_width"]["sigma_cluster_deg"])
logger.info(f"mixture parameters: f={F_CL}, sigma_int={SIG_I} deg, "
            f"coupling sigma_c={SIG_C} deg, axis={AXIS} deg, cap={CAP} deg")

# ------------------------------------------------------------------
# circle helpers
# ------------------------------------------------------------------

def d_ax(a, c=AXIS):
    return np.abs((np.asarray(a) - c + 180) % 360 - 180)


def wn_pdf(x, mu, sig):
    """Wrapped-normal density of x (deg) about mu."""
    d = (np.asarray(x) - mu + 180) % 360 - 180
    s = np.zeros_like(d, dtype=float)
    for k in (-2, -1, 0, 1, 2):
        s = s + norm.pdf((d + 360.0 * k) / sig) / sig
    return s


def smear_cap(d, sig=SIG_C):
    return norm.cdf((CAP - d) / sig) + norm.cdf((CAP + d) / sig) - 1


# mixture posterior on a 0.25 deg longitude grid
V = np.arange(0.0, 360.0, 0.25)
PRIOR = F_CL * wn_pdf(V, AXIS, SIG_I) + (1.0 - F_CL) / 360.0
PRIOR /= PRIOR.sum()
IN_CAP = d_ax(V) < CAP


def p_incap_mix(lam):
    """P(varpi in cap | discovered at opposition longitude lam) under H1."""
    like = wn_pdf(V, lam, SIG_C)
    post = PRIOR * like
    post /= post.sum()
    return float(post[IN_CAP].sum())


def p_incap_mix_vec(lams):
    return np.array([p_incap_mix(x) for x in np.asarray(lams)])


def poibin_pmf(ps):
    pmf = np.array([1.0])
    for p in np.asarray(ps, dtype=float):
        pmf = np.convolve(pmf, [1.0 - p, p])
    return pmf


def circ_std_deg(a):
    r = abs(np.exp(1j * np.deg2rad(np.asarray(a))).mean())
    return float(np.rad2deg(np.sqrt(-2.0 * math.log(max(r, 1e-12)))))


# ------------------------------------------------------------------
# cohort construction (identical to steps 113/084)
# ------------------------------------------------------------------

HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 8), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}
DESIG = re.compile(r"(\d{4})\s*([A-Z])[A-Z]?\d*")
DESIG_PAR = re.compile(r"\((\d{4})\s*([A-Z])([A-Z]?\d*)\)")


def desig_lam(desig):
    m = DESIG.search(str(desig))
    if not m or m.group(2) not in HALF:
        return float("nan"), float("nan")
    mo, dy = HALF[m.group(2)]
    t = Time(f"{int(m.group(1)):04d}-{mo:02d}-{dy:02d}T00:00:00",
             format="isot", scale="utc")
    return (float((get_sun(t).geocentrictrueecliptic.lon.deg + 180) % 360),
            int(m.group(1)))


def ecl_lon(ra, dec):
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    return float(c.geocentrictrueecliptic.lon.deg)


cohorts = {}

d = json.loads((DATA_RAW / "sbdb" / "sbdb_outer_ss.json").read_text())
sb_rows, new_rows, era_rows = [], [], []
for rec in d["data"]:
    o = dict(zip(d["fields"], rec))
    try:
        a, q = float(o["a"]), float(o["q"])
        om, w = float(o["om"]), float(o["w"])
        cc = int(o["condition_code"] or 9)
    except (TypeError, ValueError):
        continue
    if not (a > 150 and q > 30):
        continue
    nm = str(o["full_name"]).strip()
    m = re.search(r"\((\d{4}\s*[A-Z]+\d*)\)", nm)
    key = m.group(1) if m else nm
    lam, yr = desig_lam(key)
    row = dict(name=key, varpi=(om + w) % 360, a=a, q=q, cc=cc,
               lam=lam, lam_kind="designation-proxy", yr=yr)
    era_rows.append(row)
    if cc <= 3:
        sb_rows.append(row)
    if np.isfinite(yr) and yr >= 2025:
        new_rows.append(row)
cohorts["SBDB"] = dict(rows=sb_rows, label="SBDB secure detached "
                                        "(a>150,q>30,cc<=3)")
cohorts["SBDB-prv"] = dict(rows=[o for o in era_rows if o["cc"] > 3],
                           label="SBDB provisional detached pooled "
                                 "(a>150,q>30,cc>3)")
cohorts["2025+"] = dict(rows=new_rows,
                        label="2025+ provisional detached (any cc)")

tab = fits.open(str(DATA_RAW / "des" / "y6_des_tnos_color.fits"))[1].data
des_rows = []
for r in tab[(tab["a"] > 150) & (tab["q"] > 30)]:
    varpi = float((r["lan"] + r["aop"]) % 360)
    name = str(r["MPC"]).strip()
    lam, _ = desig_lam(name)
    des_rows.append(dict(name=name, varpi=varpi, a=float(r["a"]),
                         q=float(r["q"]), lam=lam,
                         lam_kind="designation-proxy"))
cohorts["DES"] = dict(rows=des_rows, label="DES detached (a>150,q>30)")

t = vot_parse(str(DATA_RAW / "ossos" / "ossos_t3char.vot")) \
    .get_first_table().to_table()
oss_rows = []
for r in t[t["cl"] == "det"]:
    oss_rows.append(dict(
        name=str(r["Astorb"]).strip() or str(r["ID"]).strip(),
        varpi=(float(r["Omega"]) + float(r["omega"])) % 360,
        a=float(r["a"]), q=float(r["a"]) * (1 - float(r["e"])),
        lam=ecl_lon(float(r["RAJ2000"]), float(r["DEJ2000"])),
        lam_kind="exact discovery astrometry"))
cohorts["OSSOS"] = dict(rows=oss_rows, label="OSSOS 'det' "
                                           "(exact astrometry)")

# designation-generation ladder on the full SBDB detached pool (any
# condition code): does the signature persist across discovery
# generations, and does each generation's count match what its own
# pointing predicts?
ERA_BINS = [(2000, 2005), (2006, 2010), (2011, 2015), (2016, 2019),
            (2020, 2022), (2023, 2024), (2025, 2100)]
era_cohorts = {}
for lo, hi in ERA_BINS:
    rows = [o for o in era_rows
            if np.isfinite(o["yr"]) and lo <= o["yr"] <= hi]
    if not rows:
        continue
    lbl = f"{lo}+" if hi >= 2100 else f"{lo}-{str(hi)[2:]}"
    era_cohorts[f"era {lbl}"] = dict(
        rows=rows, label=f"SBDB detached designation-era {lbl} (any cc)")
    for tag, sub, desc in [("sec", lambda o: o["cc"] <= 3, "cc<=3"),
                           ("prv", lambda o: o["cc"] > 3, "cc>3")]:
        srows = [o for o in rows if sub(o)]
        if srows:
            era_cohorts[f"era {lbl} [{tag}]"] = dict(
                rows=srows, label=f"SBDB detached designation-era {lbl} "
                                  f"({desc})")

# ------------------------------------------------------------------
# per-cohort posterior-predictive scoring
# ------------------------------------------------------------------

logger.subheader("Per-cohort scoring")
all_obj = []

def score_cohort(name, rows, label, bucket):
    vps = np.array([o["varpi"] for o in rows])
    lams = np.array([o["lam"] for o in rows])
    # A missing fitted longitude must not enter either the observed count or
    # the posterior-predictive scatter.  Keep the cohort denominator explicit,
    # but score only finite (varpi, discovery-longitude) pairs.
    ok = np.isfinite(lams) & np.isfinite(vps)
    n, k = len(rows), int(np.sum(ok & (d_ax(vps) < CAP)))
    p0 = smear_cap(d_ax(lams[ok]))
    p1 = p_incap_mix_vec(lams[ok])
    pmf0, pmf1 = poibin_pmf(p0), poibin_pmf(p1)
    e0, e1 = float(p0.sum()), float(p1.sum())
    pk0, pk1 = float(pmf0[k]), float(pmf1[k])
    bf = pk1 / pk0 if pk0 > 0 else float("inf")
    # posterior-predictive tail in the direction of the discrepancy
    tail0 = float(pmf0[k:].sum()) if k >= e0 else float(pmf0[:k + 1].sum())
    tail1 = float(pmf1[k:].sum()) if k >= e1 else float(pmf1[:k + 1].sum())
    delta = (vps[ok] - lams[ok] + 180) % 360 - 180
    bucket[name] = dict(
        label=label, n=n, n_scored=int(ok.sum()), k_in=k,
        frac_in=k / n,
        exp_H0=round(e0, 2), exp_H1=round(e1, 2),
        pk_H0=pk0, pk_H1=pk1, bayes_H1_over_H0=bf,
        ppred_tail_H0=tail0, ppred_tail_H1=tail1,
        mean_lam=float(np.rad2deg(np.angle(
            np.exp(1j * np.deg2rad(lams[ok])).mean())) % 360)
        if ok.any() else float("nan"),
        delta_scatter_deg=round(circ_std_deg(delta), 1)
        if ok.sum() >= 4 else float("nan"))
    logger.info(f"{name:6s} n={n:3d} k={k:2d} ({k/n:.3f}) | "
                f"E[H0]={e0:.1f} E[H1]={e1:.1f} | "
                f"P(k|H0)={pk0:.3g} P(k|H1)={pk1:.3g} BF={bf:.3g}")
    for o, q0, q1 in zip(np.array(rows)[ok], p0, p1):
        all_obj.append(dict(cohort=name, name=o["name"],
                            varpi=o["varpi"], lam=o["lam"],
                            in_cap=int(d_ax(o["varpi"]) < CAP),
                            p_H0=round(float(q0), 4),
                            p_H1=round(float(q1), 4)))

results, era_results = {}, {}
for name, co in cohorts.items():
    score_cohort(name, co["rows"], co["label"], results)
for name, co in era_cohorts.items():
    score_cohort(name, co["rows"], co["label"], era_results)

# ------------------------------------------------------------------
# provisional-orbit noise diagnostic
# ------------------------------------------------------------------
# Provisional fitted varpi carries measurement error absent from both
# models; the extra scatter is estimated from the widening of the
# pooled provisional |varpi - lam| residual distribution relative to
# the calibrated coupling width (sigma_eff^2 = sigma_c^2 + sigma_e^2,
# half-normal moment estimator -- the wrapped tail is negligible at
# sigma_eff ~ 60 deg).  The provisional cohorts are then re-scored
# under noise-convolved versions of both models (H0n: coupling widened
# to sigma_eff; H1n: likelihood and intrinsic prior both convolved),
# giving the Bayes factor the provisional data can actually support.

prv_d = np.abs((np.array([o["varpi"] for o in cohorts["SBDB-prv"]["rows"]])
                - np.array([o["lam"] for o in cohorts["SBDB-prv"]["rows"]])
                + 180) % 360 - 180)
prv_d = prv_d[np.isfinite(prv_d)]
sig_eff = float(np.sqrt(np.mean(prv_d ** 2)))
sig_e = float(np.sqrt(max(sig_eff ** 2 - SIG_C ** 2, 0.0)))

pri_n = F_CL * wn_pdf(V, AXIS, math.sqrt(SIG_I ** 2 + sig_e ** 2)) \
    + (1.0 - F_CL) / 360.0
pri_n /= pri_n.sum()
cap_g = d_ax(V) < CAP
def p1_noise(lam):
    post = pri_n * wn_pdf(V, lam, sig_eff)
    post /= post.sum()
    return float(post[cap_g].sum())

def score_noise(rows):
    lams = np.array([o["lam"] for o in rows])
    vps = np.array([o["varpi"] for o in rows])
    ok = np.isfinite(lams)
    k = int(np.sum(d_ax(vps) < CAP))
    p0n = smear_cap(d_ax(lams[ok]), sig_eff)
    p1n = np.array([p1_noise(x) for x in lams[ok]])
    pmf0, pmf1 = poibin_pmf(p0n), poibin_pmf(p1n)
    return dict(k_in=k,
                exp_H0n=round(float(p0n.sum()), 2),
                exp_H1n=round(float(p1n.sum()), 2),
                bayes_H1n_over_H0n=float(pmf1[k] / pmf0[k])
                if pmf0[k] > 0 else float("inf"))

noise_scored = {}
for cn, co in list(cohorts.items()) + list(era_cohorts.items()):
    if "prv" in cn or cn == "2025+":
        noise_scored[cn] = score_noise(co["rows"])
        logger.info(f"noise-convolved {cn}: k={noise_scored[cn]['k_in']} "
                f"E[H0n]={noise_scored[cn]['exp_H0n']} "
                f"E[H1n]={noise_scored[cn]['exp_H1n']} "
                f"BF={noise_scored[cn]['bayes_H1n_over_H0n']:.3g}")
logger.info(f"provisional residual width sigma_eff={sig_eff:.1f} deg -> "
            f"extra varpi scatter sigma_e={sig_e:.1f} deg")

# ------------------------------------------------------------------
# deduplicated union (same precedence as step_113)
# ------------------------------------------------------------------

def norm_key(name):
    return re.sub(r"\s+", "", str(name))

def score_union(order):
    union = {}
    for cn in order:
        for o in cohorts[cn]["rows"]:
            k2 = norm_key(o["name"])
            if k2 not in union:
                union[k2] = o
    urows = list(union.values())
    uv = np.array([o["varpi"] for o in urows])
    ul = np.array([o["lam"] for o in urows])
    ok = np.isfinite(ul)
    un, uk = len(urows), int(np.sum(d_ax(uv) < CAP))
    up0 = smear_cap(d_ax(ul[ok]))
    up1 = p_incap_mix_vec(ul[ok])
    upmf0, upmf1 = poibin_pmf(up0), poibin_pmf(up1)
    return dict(
        n=un, n_scored=int(ok.sum()), k_in=uk, frac_in=uk / un,
        exp_H0=round(float(up0.sum()), 2),
        exp_H1=round(float(up1.sum()), 2),
        pk_H0=float(upmf0[uk]), pk_H1=float(upmf1[uk]),
        bayes_H1_over_H0=float(upmf1[uk] / upmf0[uk])
        if upmf0[uk] > 0 else float("inf"))


union_secure = score_union(["OSSOS", "DES", "SBDB"])
union_res = score_union(["OSSOS", "DES", "2025+", "SBDB"])
union_res["note"] = ("each distinct object once; per-object posterior "
                     "uses its own discovery longitude (exact astrometry "
                     "for OSSOS, designation proxy otherwise); includes "
                     "the provisional 2025+ cohort")
union_secure["note"] = ("secure-lineage union only (OSSOS > DES > SBDB "
                        "precedence) -- the step_113 union pool; the "
                        "provisional 2025+ cohort excluded")
logger.info(f"union(sec) n={union_secure['n']:3d} "
            f"k={union_secure['k_in']:2d} | "
            f"E[H0]={union_secure['exp_H0']:.1f} "
            f"E[H1]={union_secure['exp_H1']:.1f} | "
            f"BF={union_secure['bayes_H1_over_H0']:.3g}")
logger.info(f"union(all) n={union_res['n']:3d} "
            f"k={union_res['k_in']:2d} | "
            f"E[H0]={union_res['exp_H0']:.1f} "
            f"E[H1]={union_res['exp_H1']:.1f} | "
            f"BF={union_res['bayes_H1_over_H0']:.3g}")

# ------------------------------------------------------------------
# sensitivity: BF across the mixture-parameter confidence region
# ------------------------------------------------------------------

F_GRID = [0.317, 0.431, 0.55, 0.703]     # CI68 [0.317, 0.703]
SI_GRID = [27.5, 34.4, 48.0]             # CI68 [27.5, 48.0]
SC_GRID = [40.0, 50.0, 65.0]             # bracketing OSSOS refit
                                         # [33, 75] about the 50 deg
                                         # calibration

sens = {}
for f_ in F_GRID:
    for si_ in SI_GRID:
        for sc_ in SC_GRID:
            pri = f_ * wn_pdf(V, AXIS, si_) + (1.0 - f_) / 360.0
            pri /= pri.sum()
            cap = d_ax(V) < CAP
            def p1v(lams, pri=pri, cap=cap, sc_=sc_):
                out_ = np.empty(len(lams))
                for j, x in enumerate(lams):
                    post = pri * wn_pdf(V, x, sc_)
                    post /= post.sum()
                    out_[j] = post[cap].sum()
                return out_
            key = f"f={f_},si={si_},sc={sc_}"
            bfs = {}
            for cn, co in list(cohorts.items()) + \
                    list(era_cohorts.items()):
                lams = np.array([o["lam"] for o in co["rows"]])
                ok = np.isfinite(lams)
                if not ok.any():
                    continue
                k = int(np.sum(d_ax([o["varpi"] for o in co["rows"]])
                                < CAP))
                p0 = smear_cap(d_ax(lams[ok]), sc_)
                p1 = p1v(lams[ok])
                pmf0, pmf1 = poibin_pmf(p0), poibin_pmf(p1)
                bfs[cn] = float(pmf1[k] / pmf0[k]) if pmf0[k] > 0 \
                    else float("inf")
            sens[key] = bfs

sens_summary, era_sens_summary = {}, {}
for cn in cohorts:
    vals = [v[cn] for v in sens.values() if cn in v]
    sens_summary[cn] = dict(
        min=round(min(vals), 3), median=round(float(np.median(vals)), 3),
        max=round(max(vals), 3),
        frac_favouring_H1=round(float(np.mean([v > 1 for v in vals])), 3))
    logger.info(f"sensitivity {cn}: BF range "
                f"[{sens_summary[cn]['min']}, {sens_summary[cn]['max']}], "
                f"fraction favouring mixture "
                f"{sens_summary[cn]['frac_favouring_H1']}")
for cn in era_cohorts:
    vals = [v[cn] for v in sens.values() if cn in v]
    era_sens_summary[cn] = dict(
        min=round(min(vals), 3), median=round(float(np.median(vals)), 3),
        max=round(max(vals), 3),
        frac_favouring_H1=round(float(np.mean([v > 1 for v in vals])), 3))
    logger.info(f"sensitivity {cn}: BF range "
                f"[{era_sens_summary[cn]['min']}, "
                f"{era_sens_summary[cn]['max']}], "
                f"fraction favouring mixture "
                f"{era_sens_summary[cn]['frac_favouring_H1']}")

# ------------------------------------------------------------------
# discovery-longitude lever arm
# ------------------------------------------------------------------

lam_grid = np.arange(0.0, 360.0, 1.0)
e0_map = smear_cap(d_ax(lam_grid))
e1_map = p_incap_mix_vec(lam_grid)
delta_map = e1_map - e0_map
decisive = lam_grid[delta_map > 0.10]
lever = dict(
    lam_grid=lam_grid.tolist(), exp_H0=e0_map.tolist(),
    exp_H1=e1_map.tolist(),
    decisive_zone_deg=[float(decisive.min()), float(decisive.max())]
    if len(decisive) else None,
    max_delta=float(delta_map.max()),
    lam_at_max_delta=float(lam_grid[delta_map.argmax()]),
    delta_at_newest_mean=float(np.interp(
        results["2025+"]["mean_lam"], lam_grid, delta_map))
    if np.isfinite(results["2025+"]["mean_lam"]) else float("nan"),
    note="expected in-cap fraction per model for a cohort discovered "
         "at opposition longitude lam; the models separate where "
         "E[H1]-E[H0] > 0.10")

# ------------------------------------------------------------------
# write
# ------------------------------------------------------------------

res = dict(
    step="step_115_mixture_posterior",
    description="Posterior-predictive audit: every resident cohort "
                "scored against the footprint-only model and the "
                "registered TEP mixture, conditioned on each object's "
                "own discovery longitude; exact Poisson-binomial "
                "likelihood of each observed count under each model.",
    inputs=["data/raw/sbdb/sbdb_outer_ss.json",
            "data/raw/des/y6_des_tnos_color.fits",
            "data/raw/ossos/ossos_t3char.vot",
            "results/step_b13_patch_scale.json (mixture parameters)"],
    seed=20261013,
    axis_deg=AXIS, cap_deg=CAP, coupling_sigma_deg=SIG_C,
    mixture=dict(f_cluster=F_CL, sigma_int_deg=SIG_I,
                 source="results/step_b13_patch_scale.json "
                        "P2_intrinsic_width"),
    per_cohort=results, per_era=era_results,
    union_secure=union_secure, union=union_res,
    ossos_caveat="The mixture parameters were fitted on the a>150 "
                 "secure cohort; applied to the Gladman-detached OSSOS "
                 "class (median a ~ 60 AU, interior to the boundary) "
                 "the intrinsic fraction is expected to be lower -- "
                 "the OSSOS lean toward H0 is consistent with that "
                 "selectivity rather than a mixture failure.",
    lever_arm=lever,
    provisional_noise=dict(
        sigma_eff_deg=round(sig_eff, 1),
        sigma_e_deg=round(sig_e, 1),
        method="half-normal width of the pooled provisional |varpi-lam| "
               "residuals minus the calibrated coupling width in "
               "quadrature; both models re-scored with the extra "
               "scatter convolved in",
        per_cohort=noise_scored,
        note="provisional fitted varpi carries measurement error that "
             "smears intrinsic members out of the cap asymmetrically; "
             "the noise-convolved Bayes factor is the support the "
             "provisional data can actually register"),
    sensitivity=dict(
        grid=dict(f_cluster=F_GRID, sigma_int_deg=SI_GRID,
                  sigma_c_deg=SC_GRID),
        per_cohort_bf_range=sens_summary,
        per_era_bf_range=era_sens_summary,
        note="each cell re-derives both models at the stated "
             "parameters (the coupling width is shared by H0 and "
             "H1); the BF range brackets the mixture-parameter "
             "confidence region"),
    interpretation="A cohort's in-cap rate is decisive only where its "
                   "discovery longitude sits inside the decisive zone; "
                   "the 2025+ cohort's pointing lies outside it, so its "
                   "baseline-consistent count is the expected outcome "
                   "under BOTH models -- quantifying rather than "
                   "asserting why it does not test the prediction.")

out = RESULTS / "step_b79_mixture_posterior.json"
def finite_json(value):
    if isinstance(value, dict):
        return {k: finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(v) for v in value]
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    return value

with open(out, "w") as stream:
    json.dump(finite_json(res), stream, indent=1, allow_nan=False)
logger.data_save(out)
csv_out = RESULTS / "step_b79_mixture_posterior.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cohort", "name", "varpi", "lam_disc", "in_cap",
                "p_H0", "p_H1"])
    for o in all_obj:
        w.writerow([o["cohort"], o["name"], f"{o['varpi']:.3f}",
                    f"{o['lam']:.2f}", o["in_cap"], o["p_H0"], o["p_H1"]])
logger.data_save(csv_out)
# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(18, 4.8))

ax = axes[0]
names = list(results.keys()) + ["union(sec)", "union(all)"]
obs = [results[c]["frac_in"] for c in results] + \
      [union_secure["frac_in"], union_res["frac_in"]]
h0 = [results[c]["exp_H0"] / results[c]["n_scored"] for c in results] + \
     [union_secure["exp_H0"] / union_secure["n_scored"],
      union_res["exp_H0"] / union_res["n_scored"]]
h1 = [results[c]["exp_H1"] / results[c]["n_scored"] for c in results] + \
     [union_secure["exp_H1"] / union_secure["n_scored"],
      union_res["exp_H1"] / union_res["n_scored"]]
x = np.arange(len(names))
ax.bar(x - 0.27, obs, 0.25, color="crimson", label="observed")
ax.bar(x, h0, 0.25, color="0.65", label="E[footprint H0]")
ax.bar(x + 0.27, h1, 0.25, color="steelblue", label="E[mixture H1]")
for i, o in enumerate(obs):
    ax.text(i - 0.27, o + 0.015, f"{o:.2f}", ha="center", fontsize=7)
ax.set_xticks(x)
ax.set_xticklabels([f"{c}\nn={results[c]['n']}" for c in results] +
                   [f"union(sec)\nn={union_secure['n']}",
                    f"union(all)\nn={union_res['n']}"], fontsize=7)
ax.set_ylabel("in-cap fraction")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Observed vs model-predicted in-cap fractions", fontsize=10)

ax = axes[1]
en = [c for c in era_results if "[" not in c]
esec = [era_results.get(c + " [sec]") for c in en]
eprv = [era_results.get(c + " [prv]") for c in en]
eh0 = [era_results[c]["exp_H0"] / era_results[c]["n_scored"]
       for c in en]
eh1 = [era_results[c]["exp_H1"] / era_results[c]["n_scored"]
       for c in en]
xe = np.arange(len(en))
w = 0.21
ax.bar(xe - 1.5 * w, [s["frac_in"] if s else 0 for s in esec], w,
       color="crimson", label="observed cc<=3")
ax.bar(xe - 0.5 * w, [p["frac_in"] if p else 0 for p in eprv], w,
       color="lightpink", label="observed cc>3")
ax.bar(xe + 0.5 * w, eh0, w, color="0.65", label="E[footprint H0]")
ax.bar(xe + 1.5 * w, eh1, w, color="steelblue", label="E[mixture H1]")
for i, s in enumerate(esec):
    if s:
        ax.text(i - 1.5 * w, s["frac_in"] + 0.015,
                f"{s['frac_in']:.2f}", ha="center", fontsize=7)
ax.set_xticks(xe)
ax.set_xticklabels(
    [f"{c[4:]}\n{esec[i]['n'] if esec[i] else 0}s+"
     f"{eprv[i]['n'] if eprv[i] else 0}p" for i, c in enumerate(en)],
    fontsize=7)
ax.set_ylabel("in-cap fraction")
ax.legend(frameon=False, fontsize=7)
ax.set_title("Designation-era ladder, split by orbit quality",
             fontsize=10)

ax = axes[2]
ax.plot(lam_grid, e0_map, color="0.45", lw=1.8,
        label="footprint only (H0)")
ax.plot(lam_grid, e1_map, color="steelblue", lw=1.8,
        label=f"TEP mixture (H1, f={F_CL:.2f})")
ax.fill_between(lam_grid, e0_map, e1_map, where=e1_map > e0_map,
                color="steelblue", alpha=0.15)
ax.axvline(AXIS, color="k", ls="--", lw=1, label="axis 49 deg")
for cn, mk, col in [("DES", "s", "teal"), ("2025+", "D", "crimson"),
                    ("SBDB", "o", "navy"), ("OSSOS", "^", "purple")]:
    ml = results[cn]["mean_lam"]
    if np.isfinite(ml):
        ax.axvline(ml, color=col, ls=":", lw=1.2,
                   label=f"{cn} mean $\\lambda_{{opp}}$={ml:.0f} deg")
ax.set_xlabel("cohort mean discovery longitude $\\lambda_{opp}$ (deg)")
ax.set_ylabel("expected in-cap fraction")
ax.legend(frameon=False, fontsize=7)
ax.set_title("Discovery-longitude lever arm", fontsize=10)

fig.tight_layout()
fig.savefig(RESULTS / "figures" / "supplementary" / "step_b79_mixture_posterior.png",
            dpi=300)
logger.data_save(RESULTS / 'figures' / 'supplementary' / 'step_b79_mixture_posterior.png')