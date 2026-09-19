#!/usr/bin/env python3
"""Step 124 -- lineage-independent aphelion dipole (SBDB eras).

The prospective analysis (step_117) found that the post-2017 SBDB
cohort fails the registered three-leg residual test, and the
factorial audit (steps_118-120) diagnosed a unipolar,
cohort-private residual field pointing ~(120, -40) that is
invisible to the CODE machinery at the same sky.  The residual
channel is reconstruction machinery; the aphelion-direction
channel is not -- each comet's aphelion unit vector is a function
of the published osculating elements (w, Om, i) alone, with no
three-leg integration and no catalogue boundary solution.  The
residual failure therefore does not automatically doom the
spatial datum: if the newly designated population leans toward
the declared axis in its aphelia, the sky anomaly exists in an
independently discovered, independently fitted population and
the caveat narrows to the reconstruction channel.

This step runs the identical tide-aware dipole test of step_061/
step_108 on the SBDB record split by designation era:

  post-2017   C/2018+ (the step_117 holdout cohort definition:
              0.95<=e<1.5, arc>=30 d, nobs>=20, q>=0.1)
  pre-2018    C/<2018 (the step_120 cohort, same cuts)

T1  d_par = <aph . axis> under the tide-aware null (preserve
    |b_gal|, isotropize galactic longitude) on the cap axis, the
    detached-TNO axis, the anti-axis, the galactic anticenter,
    and the displaced (120,-40) axis the residual audit found.
T2  ecliptic-footprint null (preserve |b_ecl|) -- the survey
    footprint scans the ecliptic; does the lean survive a null
    that keeps the ecliptic band occupancy?
T3  quality strata -- bound vs hyperbolic, perihelion-spanning
    vs one-sided arcs (f_pre = (tp - first_obs)/(arc window)),
    short/long arcs, small/large nobs, condition code: the same
    cells where step_119/122 localized the residual systematic.
T4  census e-bin gradient (0.90<=e<1.02) per era for CometEls
    comparability.
T5  free spherical-mean dipole per cohort: does the population
    self-organize toward the declared axis, the displaced axis,
    or neither?
T6  era combination: Fisher p across the two SBDB eras and the
    comparison against the CODE/Warsaw/CometEls values already
    measured.

Inputs
------
data/raw/sbdb/sbdb_comets_all.json
results/step_b30_proper_time_slip.csv      (training membership)
results/step_b84_pre2018_sbdb.csv          (f_pre arc-fraction)

Outputs
-------
results/step_b88_lineage_aphelion.json
results/step_b88_lineage_aphelion.csv
results/figures/step_b88_lineage_aphelion.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_124_lineage_aphelion")
tee_stdout(logger)
logger.header("Lineage-independent aphelion dipole (SBDB eras)")

import csv
import json
import math
import re
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scripts.utils.coordinates import angular_separation
from scipy.stats import chi2, kstest
from scripts.utils.statistics import fisher_combination

SEED = 20261030
N_MC = 20000
rng = np.random.default_rng(SEED)



def perih_dir(om, Om, inc):
    co, so = np.cos(om), np.sin(om)
    cO, sO, ci, si = (np.cos(Om), np.sin(Om),
                      np.cos(inc), np.sin(inc))
    return np.array([cO * co - sO * so * ci,
                     sO * co + cO * so * ci, so * si])


def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b) * math.cos(l),
                     math.cos(b) * math.sin(l), math.sin(b)])


def gv(l, b):
    l, b = math.radians(l), math.radians(b)
    return GAL2ECL @ np.array([math.cos(b) * math.cos(l),
                               math.cos(b) * math.sin(l),
                               math.sin(b)])


def gal_of(v):
    g = ECL2GAL @ v
    return (math.degrees(math.atan2(g[1], g[0])) % 360,
            math.degrees(math.asin(np.clip(g[2], -1, 1))))


def ecl_lb(v):
    return (math.degrees(math.atan2(v[1], v[0])) % 360,
            math.degrees(math.asin(np.clip(v[2], -1, 1))))


def sep_deg(a, b):
    return angular_separation(a, b)


AXES = {
    "cap_34_-13": lv(34.0, -13.0),
    "detached_50_-17": lv(49.9, -17.0),
    "anti_214_13": lv(214.0, 13.0),
    "gal_antictr": gv(180, 0),
    "displaced_120_-40": lv(120.0, -40.0),
}


def dipole_null(aphs, axis, frame="gal", n_mc=N_MC):
    """MC null for <aph.axis>: preserve |b| in the given frame,
    isotropize the frame longitude, random hemisphere sign."""
    n = len(aphs)
    obs = float(np.mean([np.dot(a, axis) for a in aphs]))
    if frame == "gal":
        bs = np.array([abs(gal_of(a)[1]) for a in aphs])
    else:
        bs = np.array([abs(ecl_lb(a)[1]) for a in aphs])
    cnt = 0
    for _ in range(n_mc):
        lrand = rng.uniform(0, 2 * np.pi, n)
        bpick = bs[rng.integers(0, n, n)]
        sgn = rng.choice([-1.0, 1.0], n)
        bfr = np.deg2rad(bpick * sgn)
        if frame == "gal":
            vg = np.stack([np.cos(bfr) * np.cos(lrand),
                           np.cos(bfr) * np.sin(lrand),
                           np.sin(bfr)], axis=-1)
            with np.errstate(all="ignore"):
                ve = vg @ GAL2ECL.T
        else:
            ve = np.stack([np.cos(bfr) * np.cos(lrand),
                           np.cos(bfr) * np.sin(lrand),
                           np.sin(bfr)], axis=-1)
        with np.errstate(all="ignore"):
            stat = float(np.mean(ve @ axis))
        if stat >= obs:
            cnt += 1
    return obs, (cnt + 1) / (n_mc + 1)


def dipole_block(recs, label):
    out = {"label": label, "n": len(recs)}
    aphs = np.array([r["aph"] for r in recs])
    t1 = {}
    for name, ax in AXES.items():
        obs, p = dipole_null(aphs, ax, "gal")
        t1[name] = {"d_par": round(obs, 4), "p_tide_null": float(p)}
    out["T1_tide_aware"] = t1
    obs_e, p_e = dipole_null(aphs, AXES["cap_34_-13"], "ecl")
    out["T2_ecliptic_null_cap"] = {"d_par": round(obs_e, 4),
                                   "p_ecl_null": float(p_e)}
    logger.metric(
        f"dipole[{label}]",
        f"n={len(recs)} d(cap34)={t1['cap_34_-13']['d_par']:+.4f} "
        f"p={t1['cap_34_-13']['p_tide_null']:.3g} "
        f"d(det50)={t1['detached_50_-17']['d_par']:+.4f} "
        f"d(displ)={t1['displaced_120_-40']['d_par']:+.4f}")
    return out


def free_dipole(aphs):
    A = np.array(aphs)
    mean_v = A.mean(axis=0)
    R = float(np.linalg.norm(mean_v))
    dv = mean_v / R if R > 0 else np.array([1, 0, 0])
    lam = float(np.degrees(np.arctan2(dv[1], dv[0])) % 360)
    bet = float(np.degrees(np.arcsin(np.clip(dv[2], -1, 1))))
    return dict(R=round(R, 4), lam=round(lam, 2), bet=round(bet, 2),
                sep_cap=round(sep_deg(dv, AXES["cap_34_-13"]), 2),
                sep_det=round(sep_deg(dv, AXES["detached_50_-17"]), 2),
                sep_displ=round(
                    sep_deg(dv, AXES["displaced_120_-40"]), 2))


# ------------------------------------------------------------------ parse
sbdb = json.load(open(DATA_RAW / "sbdb" / "sbdb_comets_all.json"))
fields = sbdb["fields"]
rows = [dict(zip(fields, rec)) for rec in sbdb["data"]]
logger.info(f"SBDB comet table: {len(rows)} rows")


def fnum(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")


def desig_year(name):
    m = re.match(r"\s*[CP]/(\d{4})", name)
    return int(m.group(1)) if m else None


def is_fragment(name):
    return bool(re.match(r"\s*C/\d{4}\s+\S+-\w", name))


# f_pre (arc fraction before perihelion) from the step_120 CSV
fpre = {}
try:
    for r in csv.DictReader(
            open(RESULTS / "step_b84_pre2018_sbdb.csv")):
        try:
            fpre[r["desig"].strip()] = float(r["f_pre"])
        except (KeyError, ValueError):
            pass
except FileNotFoundError:
    logger.info("step_b84 CSV absent; f_pre strata skipped")

# training-set membership
def norm(d):
    d = d.strip().upper().replace(' ', '')
    m = re.search(r'([CPD]/\d{4}[A-Z]+\d*|\d{4}[A-Z]+\d*)', d)
    return m.group(1).split('-')[0] if m else d


train = {norm(r['desig'])
         for r in csv.DictReader(
             open(RESULTS / "step_b30_proper_time_slip.csv"))}

recs = []
for r in rows:
    name = str(r["full_name"]).strip()
    if not name.startswith("C/"):
        continue
    yr = desig_year(name)
    if yr is None or is_fragment(name):
        continue
    e, q = fnum(r, "e"), fnum(r, "q")
    arc, nobs = fnum(r, "data_arc"), fnum(r, "n_obs_used")
    w, Om, inc = fnum(r, "w"), fnum(r, "om"), fnum(r, "i")
    if not all(np.isfinite(v) for v in (e, q, w, Om, inc)):
        continue
    aph = -perih_dir(math.radians(w), math.radians(Om),
                     math.radians(inc))
    recs.append(dict(name=name, yr=yr, q=q, e=e, i=inc,
                     arc=arc, nobs=nobs,
                     cc=fnum(r, "condition_code"),
                     f_pre=fpre.get(name, float("nan")),
                     in_train=norm(name) in train,
                     aph=aph))
logger.info(f"parsed {len(recs)} C/ comets with finite elements")


def era_cohort(post):
    """Mirror the step_117/120 near-parabolic definitions."""
    out = []
    for r in recs:
        if (r["yr"] >= 2018) != post:
            continue
        if not (0.95 <= r["e"] < 1.5):
            continue
        if not (np.isfinite(r["q"]) and r["q"] >= 0.1):
            continue
        if not (np.isfinite(r["arc"]) and r["arc"] >= 30
                and np.isfinite(r["nobs"]) and r["nobs"] >= 20):
            continue
        out.append(r)
    return out


post17 = era_cohort(True)
pre18 = era_cohort(False)
logger.info(f"near-parabolic cohorts: post-2017 n={len(post17)}, "
            f"pre-2018 n={len(pre18)}")

res = {}
res["post2017_near_parabolic"] = dipole_block(post17, "SBDB post-2017")
res["pre2018_near_parabolic"] = dipole_block(pre18, "SBDB pre-2018")

# deep-plunger subsets (primary cohort restriction of step_117)
for tag, co in (("post2017", post17), ("pre2018", pre18)):
    sub = [r for r in co if r["q"] < 3.1]
    if len(sub) >= 15:
        res[f"{tag}_deep_plunger_q_lt_3p1"] = dipole_block(
            sub, f"{tag} q<3.1")

# ------------------------------------------------------------------ T3 strata
strata = {}
for tag, co in (("post2017", post17), ("pre2018", pre18)):
    st = {}
    cuts = {
        "bound_e_lt_1": [r for r in co if r["e"] < 1.0],
        "hyperbolic_e_ge_1": [r for r in co if r["e"] >= 1.0],
        "arc_lt_365d": [r for r in co if r["arc"] < 365],
        "arc_ge_365d": [r for r in co if r["arc"] >= 365],
        "nobs_lt_300": [r for r in co if r["nobs"] < 300],
        "nobs_ge_300": [r for r in co if r["nobs"] >= 300],
        "cc_le_2": [r for r in co
                    if np.isfinite(r["cc"]) and r["cc"] <= 2],
        "span_arc": [r for r in co
                     if np.isfinite(r["f_pre"])
                     and 0.15 < r["f_pre"] < 0.85],
        "one_sided_arc": [r for r in co
                          if np.isfinite(r["f_pre"])
                          and (r["f_pre"] <= 0.15
                               or r["f_pre"] >= 0.85)],
        "in_training": [r for r in co if r["in_train"]],
        "not_in_training": [r for r in co if not r["in_train"]],
    }
    for sname, sub in cuts.items():
        if len(sub) < 15:
            st[sname] = {"n": len(sub)}
            continue
        obs, p = dipole_null([r["aph"] for r in sub],
                             AXES["cap_34_-13"], "gal")
        st[sname] = {"n": len(sub), "d_par": round(obs, 4),
                     "p_tide_null": float(p)}
    strata[tag] = st
    for sname, s in st.items():
        if "d_par" in s:
            logger.metric(f"stratum[{tag}:{sname}]",
                          f"n={s['n']} d_par={s['d_par']:+.4f} "
                          f"p={s['p_tide_null']:.3g}")
res["T3_strata"] = strata

# ------------------------------------------------------------------ T4 e-bins
t4 = {}
for tag, co in (("post2017", post17), ("pre2018", pre18)):
    rows_t = []
    for lo, hi in ((0.90, 0.95), (0.95, 0.985), (0.985, 1.0),
                   (1.0, 1.02)):
        sub = [r for r in recs
               if ((r["yr"] >= 2018) == (tag == "post2017"))
               and lo <= r["e"] < hi and r["q"] >= 0.1]
        if len(sub) < 10:
            rows_t.append(dict(e_bin=f"[{lo},{hi})", n=len(sub)))
            continue
        obs, p = dipole_null([r["aph"] for r in sub],
                             AXES["cap_34_-13"], "gal")
        rows_t.append(dict(e_bin=f"[{lo},{hi})", n=len(sub),
                           d_par=round(obs, 4),
                           p_tide_null=float(p)))
    t4[tag] = rows_t
res["T4_e_gradient"] = t4

# ------------------------------------------------------------------ T5 free dipoles
res["T5_free_dipole"] = {
    "post2017": free_dipole([r["aph"] for r in post17]),
    "pre2018": free_dipole([r["aph"] for r in pre18]),
}
for tag in ("post2017", "pre2018"):
    f = res["T5_free_dipole"][tag]
    logger.metric(f"free_dipole[{tag}]",
                  f"({f['lam']:.1f},{f['bet']:.1f}) R={f['R']:.3f} "
                  f"sep_cap={f['sep_cap']:.0f} "
                  f"sep_displ={f['sep_displ']:.0f}")

# ------------------------------------------------------------------ T6 era Fisher
p_post = res["post2017_near_parabolic"]["T1_tide_aware"][
    "cap_34_-13"]["p_tide_null"]
p_pre = res["pre2018_near_parabolic"]["T1_tide_aware"][
    "cap_34_-13"]["p_tide_null"]
combined = fisher_combination([p_post, p_pre])
fisher, p_comb = combined["statistic"], combined["p"]
res["T6_era_combination"] = dict(
    p_post2017=p_post, p_pre2018=p_pre,
    fisher_stat=fisher, df=combined["df"], p_combined=p_comb,
    mc_draws_per_test=N_MC, mc_p_floor=1/(N_MC+1),
    interpretation=combined["interpretation"],
    limitation="Shared selection effects are not removed by disjoint eras. Zero Monte Carlo exceedances indicate finite resolution.")
logger.metric("era_fisher", f"X2={fisher:.2f} p={p_comb:.3g}")

# ------------------------------------------------------------------ figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.4))

# panel 1: d_par by cohort + strata
labels, vals, cols = [], [], []
for tag, lab in (("post2017_near_parabolic", "post-2017"),
                 ("pre2018_near_parabolic", "pre-2018")):
    d = res[tag]["T1_tide_aware"]["cap_34_-13"]
    labels.append(lab); vals.append(d["d_par"])
    cols.append("crimson" if tag.startswith("post") else "tab:blue")
for tag in ("post2017", "pre2018"):
    for sname in ("bound_e_lt_1", "hyperbolic_e_ge_1",
                  "span_arc", "one_sided_arc"):
        s = res["T3_strata"][tag].get(sname, {})
        if "d_par" in s:
            labels.append(f"{tag[:4]}\n{sname.split('_')[0]}")
            vals.append(s["d_par"])
            cols.append("salmon" if tag == "post2017"
                        else "cornflowerblue")
ax[0].bar(range(len(vals)), vals, color=cols)
ax[0].set_xticks(range(len(vals)))
ax[0].set_xticklabels(labels, fontsize=6, rotation=45, ha="right")
ax[0].axhline(0, color="k", lw=0.8)
ax[0].set(ylabel=r"$d_\parallel$ toward cap axis",
          title="aphelion dipole by era and stratum")

# panel 2: axis specificity for both cohorts
xpos = np.arange(len(AXES))
w_ = 0.36
for j, tag in enumerate(("post2017_near_parabolic",
                         "pre2018_near_parabolic")):
    ds = [res[tag]["T1_tide_aware"][k]["d_par"] for k in AXES]
    ax[1].bar(xpos + (j - 0.5) * w_, ds, w_,
              color="crimson" if j == 0 else "tab:blue",
              label="post-2017" if j == 0 else "pre-2018")
ax[1].set_xticks(xpos)
ax[1].set_xticklabels([k.replace("_", "\n") for k in AXES],
                      fontsize=6.5)
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set(ylabel=r"$d_\parallel$", title="axis specificity")
ax[1].legend(fontsize=8)

# panel 3: e-gradient
for j, tag in enumerate(("post2017", "pre2018")):
    rows_t = t4[tag]
    xs = np.arange(len(rows_t))
    ys = [t.get("d_par", 0) for t in rows_t]
    ax[2].bar(xs + (j - 0.5) * w_, ys, w_,
              color="crimson" if j == 0 else "tab:blue",
              label=tag)
ax[2].set_xticks(np.arange(4))
ax[2].set_xticklabels([t["e_bin"] for t in t4["post2017"]],
                      fontsize=7)
ax[2].axhline(0, color="k", lw=0.8)
ax[2].set(xlabel="eccentricity bin", ylabel=r"$d_\parallel$",
          title="eccentricity gradient")
ax[2].legend(fontsize=8)
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b88_lineage_aphelion.png", dpi=150)
plt.close(fig)
logger.data_save(FIG / "step_b88_lineage_aphelion.png")

# ------------------------------------------------------------------ output
d_post = res["post2017_near_parabolic"]["T1_tide_aware"]["cap_34_-13"]
d_pre = res["pre2018_near_parabolic"]["T1_tide_aware"]["cap_34_-13"]
d_post_d = res["post2017_near_parabolic"]["T1_tide_aware"][
    "detached_50_-17"]
d_pre_d = res["pre2018_near_parabolic"]["T1_tide_aware"][
    "detached_50_-17"]
verdict = (
    f"On the reconstruction-free aphelion channel, the post-2017 "
    f"SBDB holdout cohort (n={len(post17)}) leans toward the cap "
    f"axis at d_par={d_post['d_par']:+.3f} "
    f"(tide-aware p={d_post['p_tide_null']:.3g}) and toward the "
    f"detached-TNO axis at {d_post_d['d_par']:+.3f} "
    f"(p={d_post_d['p_tide_null']:.3g}); the pre-2018 cohort "
    f"(n={len(pre18)}) leans at {d_pre['d_par']:+.3f} "
    f"(p={d_pre['p_tide_null']:.3g}) / "
    f"{d_pre_d['d_par']:+.3f} (p={d_pre_d['p_tide_null']:.3g}).  "
    f"Combined era Fisher p={p_comb:.3g}.  This channel uses "
    f"published osculating elements only -- no three-leg "
    f"reconstruction -- so it adjudicates whether the spatial "
    f"anomaly lives in the discovered population or only in the "
    f"Poznan-lineage fit record.")

res_full = dict(
    step="step_124_lineage_aphelion",
    description=("Aphelion-dipole channel on the SBDB record "
                 "split by designation era: the reconstruction-"
                 "independent spatial datum applied to the same "
                 "cohorts where the residual channel was tested "
                 "(steps_117-120).  Tide-aware and ecliptic-"
                 "footprint nulls, quality strata matching the "
                 "diagnosed residual systematic, free-dipole "
                 "recovery, and era combination."),
    inputs=["data/raw/sbdb/sbdb_comets_all.json",
            "results/step_b30_proper_time_slip.csv",
            "results/step_b84_pre2018_sbdb.csv"],
    seed=SEED, n_mc=N_MC,
    n_post2017=len(post17), n_pre2018=len(pre18),
    results=res,
    verdict=verdict,
    caveats=[
        "Aphelion directions derive from SBDB osculating "
        "elements; the discovered population itself may carry "
        "survey footprint structure, which the tide-aware and "
        "ecliptic-band nulls only partially absorb (step_109 "
        "bounds the residual footprint at the few-percent level "
        "on the same channel).",
        "The post-2017 and pre-2018 cohorts share no objects by "
        "construction (designation-era split), so their dipoles "
        "are compositionally independent; each is still a "
        "discovered-population sample, not an injected recovery.",
        "A flat post-2017 dipole would not isolate which side "
        "fails -- the anomaly could be confined to the older "
        "population (survival bias: comets seen repeatedly "
        "before 2018) rather than to a fitting lineage."])

out = RESULTS / "step_b88_lineage_aphelion.json"
json.dump(res_full, open(out, "w"), indent=1, default=float)
logger.data_save(out)

csv_out = RESULTS / "step_b88_lineage_aphelion.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "yr", "era", "q", "e", "arc", "nobs",
                "cc", "f_pre", "in_training",
                "aph_lam", "aph_bet"])
    for r in post17 + pre18:
        lam, bet = ecl_lb(r["aph"])
        w.writerow([r["name"], r["yr"],
                    "post2017" if r["yr"] >= 2018 else "pre2018",
                    r["q"], r["e"], r["arc"], r["nobs"], r["cc"],
                    r["f_pre"], int(r["in_train"]),
                    round(lam, 3), round(bet, 3)])
logger.data_save(csv_out)
logger.info("verdict: " + verdict)
logger.info("lineage aphelion dipole complete")
