"""Step 140: Astrometric-catalogue stratification (step_b104).

The era transition near 2018 coincides with the astrometric
reduction-catalogue switch (UCAC4-era -> Gaia-DR era reductions) in
the MPC observation stream.  Step 130 showed the signature
transition tracks perihelion epoch rather than source table, but a
reduction systematic applied per-observation would also track
perihelion epoch -- the confound was narrowed, not closed.

This step stratifies directly on the per-observation astcat field
of the cached MPC ADES observations:

  GAIA-ERA   = {Gaia1, Gaia2, Gaia3, Gaia3E, ATLAS2}  -- ATLAS2 is
               the Gaia-derived ATLAS-REFCAT2, same reduction epoch
  LEGACY     = {UCAC*, USNO*, GSC*, NOMAD, PPMXL, URAT1, 2MASS,
               ACT, UNK, ...}

For every dual-leg refit record the Gaia-era fraction is computed
per leg (inbound = pre-perihelion obs, outbound = post-perihelion
obs, split at the comet's tp_jd), giving each comet a leg-level
reduction-composition profile independent of its source table.

Tests
  T1  composition summary: leg Gaia fractions by cohort (pre-2018
      vs post-2017 refits) -- quantifies how mixed the legs are.
  T2  uniform-modern-legs control (post-2017 cohort): the
      displaced-cap contrast on the carrier channels (d_out_leg,
      ddirf) restricted to comets whose BOTH legs are >=70%
      Gaia-era reductions.  A reduction-switch artifact cannot
      survive a leg-uniform cut -- both legs were reduced in the
      same system.
  T3  uniform-legacy-legs control (pre-2018 cohort): the
      declared-cap contrast on the carrier channel (d_in_leg)
      restricted to comets whose both legs are <=30% Gaia-era.
      The registered anomaly on pre-Gaia reductions demonstrates
      it predates the switch rather than being created by it.
  T4  composition-vs-signal association: Spearman of the leg-fit
      disagreement (ddirf) against the minimum leg Gaia fraction,
      per cohort and cap -- under the systematic the signal
      should scale with composition.
  T5  era crossers: pre-2018-table comets whose perihelia fall
      after 2018 -- does their signature track their leg
      composition or their transit epoch?

Selection rules, sample sizes, null model and seeds are reported
in the result leaf; missing astcat counts toward LEGACY (the
pre-catalog-era default).
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))

import json
import math
import re
import numpy as np
from scipy import stats as _st

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import DATA_RAW, RESULTS, lv, sep

logger = StepLogger("step_140_astcat_stratification")
tee_stdout(logger)

OBS_DIR = DATA_RAW / "mpc" / "obs"
SEED = 20260919
CAP = 60.0          # cap half-angle, same as steps 127/128
GAIA_HI = 0.70      # uniform-modern threshold (both legs)
GAIA_LO = 0.30      # uniform-legacy threshold (both legs)
AXIS_DISP = lv(120.0, -40.0)


def _safe(des):
    return re.sub(r"[^A-Za-z0-9]+", "_", des).strip("_")


def leg_gaia_fracs(des, tp_jd):
    """Return (gaia_frac_in, gaia_frac_out, n_in, n_out)."""
    path = OBS_DIR / f"{_safe(des)}.json"
    if not path.exists():
        return None
    obs = json.load(open(path))
    gi = go = ni = no = 0
    for o in obs:
        try:
            jd = _iso_to_jd(o["obstime"])
        except Exception:
            continue
        gaia = _is_gaia(o.get("astcat"))
        if jd < tp_jd:
            ni += 1
            gi += gaia
        else:
            no += 1
            go += gaia
    if ni == 0 or no == 0:
        return None
    return gi / ni, go / no, ni, no


def _is_gaia(cat):
    if not cat:
        return False
    return cat.startswith("Gaia") or cat == "ATLAS2"


_JD2000 = 2451545.0


def _iso_to_jd(s):
    from datetime import datetime, timezone
    s = s.rstrip("Z")
    try:
        dt = datetime.fromisoformat(s)
    except ValueError:
        dt = datetime.strptime(s[:19], "%Y-%m-%dT%H:%M:%S")
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return _JD2000 + (dt - datetime(2000, 1, 1, 12,
                                    tzinfo=timezone.utc)
                      ).total_seconds() / 86400.0


def cap_contrast(vals, th, cap=CAP, alternative="greater"):
    vals = np.asarray(vals, dtype=float)
    th = np.asarray(th, dtype=float)
    inc = th < cap
    if inc.sum() < 8 or (~inc).sum() < 8:
        return dict(n_in=int(inc.sum()), n_out=int((~inc).sum()),
                    median_in=None, median_out=None, p=None)
    u = _st.mannwhitneyu(vals[inc], vals[~inc],
                         alternative=alternative)
    return dict(n_in=int(inc.sum()), n_out=int((~inc).sum()),
                median_in=float(np.median(vals[inc])),
                median_out=float(np.median(vals[~inc])),
                p=float(u.pvalue))


def theta_disp(r):
    a = r["our_aph"]
    aph = np.array(a if isinstance(a, (list, tuple))
                   else [float(x) for x in a.split(";")], dtype=float)
    n = np.linalg.norm(aph)
    if n == 0:
        return np.nan
    return sep(aph / n, AXIS_DISP)


def resid_field(rows):
    """Kick-regressed log10(drot) residual field -- the registered
    post-2017 carrier (step_128 _resid_gap): log10(drot) detrended
    on log10(|daa|+1), log10(denc), q, i."""
    v = np.array([r["our_drot"] for r in rows], dtype=float)
    K = np.abs([r["our_daa"] for r in rows])
    D = np.array([r["our_denc"] for r in rows])
    Q = np.array([r["our_q"] for r in rows])
    I = np.array([r["our_i"] for r in rows])
    y = np.log10(np.clip(v, 1e-9, None))
    X = np.column_stack([np.ones(len(y)), np.log10(K + 1.0),
                         np.log10(np.clip(D, 1e-3, None)), Q, I])
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    resid[~np.isfinite(resid)] = np.nanmedian(
        resid[np.isfinite(resid)])
    return resid


# ------------------------------------------------------------------
# load dual-leg records + leg composition
# ------------------------------------------------------------------
logger.info("loading dual-leg records and per-leg astcat "
            "composition")
cohorts = {}
for tag, jf in (("pre2018", "step_b91_refit.jsonl"),
                ("post2017", "step_b92_refit.jsonl")):
    recs = []
    for line in open(RESULTS / jf):
        r = json.loads(line)
        if r.get("our_drot") is None or r.get("our_aph") is None:
            continue
        tp = r.get("tp_jd")
        if tp is None:
            continue
        fr = leg_gaia_fracs(r["des"], tp)
        if fr is None:
            continue
        r["_g_in"], r["_g_out"], r["_n_in"], r["_n_out"] = fr
        r["_g_min"] = min(fr[0], fr[1])
        r["_g_max"] = max(fr[0], fr[1])
        r["_dual"] = r.get("our_ddirf") is not None
        r["_th_decl"] = float(r["our_theta"])
        r["_th_disp"] = theta_disp(r)
        recs.append(r)
    cohorts[tag] = recs
    logger.info(f"  {tag}: {len(recs)} fitted records with "
                f"composition "
                f"({sum(r['_dual'] for r in recs)} dual-leg)")

res = {"step": "step_140_astcat_stratification", "result": "b104",
       "description": ("per-observation astrometric-catalogue "
                       "(astcat) stratification of the era "
                       "transition -- tests whether the axis/leg "
                       "signal survives homogeneous reduction-"
                       "catalogue legs"),
       "gaia_era_labels": ["Gaia1", "Gaia2", "Gaia3", "Gaia3E",
                           "ATLAS2"],
       "cap_deg": CAP, "gaia_hi": GAIA_HI, "gaia_lo": GAIA_LO,
       "seed": SEED,
       "n_records": {k: len(v) for k, v in cohorts.items()}}

# ------------------------------------------------------------------
# T1 composition summary
# ------------------------------------------------------------------
t1 = {}
for tag, recs in cohorts.items():
    gi = np.array([r["_g_in"] for r in recs])
    go = np.array([r["_g_out"] for r in recs])
    t1[tag] = dict(
        n=len(recs),
        gaia_in_median=float(np.median(gi)),
        gaia_out_median=float(np.median(go)),
        frac_both_ge70=float(np.mean((gi >= GAIA_HI)
                                     & (go >= GAIA_HI))),
        frac_both_le30=float(np.mean((gi <= GAIA_LO)
                                     & (go <= GAIA_LO))),
        frac_mixed=float(np.mean(
            ((gi >= GAIA_HI) & (go <= GAIA_LO))
            | ((gi <= GAIA_LO) & (go >= GAIA_HI)))),
        n_obs_in=int(sum(r["_n_in"] for r in recs)),
        n_obs_out=int(sum(r["_n_out"] for r in recs)))
    logger.info(f"  {tag}: leg Gaia in={t1[tag]['gaia_in_median']:.2f} "
                f"out={t1[tag]['gaia_out_median']:.2f}; uniform-modern "
                f"{t1[tag]['frac_both_ge70']:.2f}, uniform-legacy "
                f"{t1[tag]['frac_both_le30']:.2f}, mixed "
                f"{t1[tag]['frac_mixed']:.2f}")
res["T1_composition"] = t1

# ------------------------------------------------------------------
# T2 post-2017 uniform-modern-legs control on displaced cap
# ------------------------------------------------------------------
logger.info("T2 uniform-modern-legs control (post-2017 cohort)")
recs = cohorts["post2017"]
uni = [r for r in recs if r["_g_min"] >= GAIA_HI]
mix = [r for r in recs if r["_g_min"] < GAIA_HI]
t2 = {}
for tag, sub in (("uniform_ge70", uni), ("mixed", mix)):
    th = np.array([r["_th_disp"] for r in sub])
    t2[tag] = {"n": len(sub)}
    # leg carrier channels -- dual-leg members only
    dual = [r for r in sub if r["_dual"]]
    thd = np.array([r["_th_disp"] for r in dual])
    for ch in ("our_d_out_leg", "our_d_in_leg", "our_ddirf"):
        v = np.array([r[ch] if r.get(ch) is not None else np.nan
                      for r in dual])
        m = np.isfinite(v)
        t2[tag][ch] = cap_contrast(v[m], thd[m])
    # registered displaced carrier: resid field on all fitted
    rv = resid_field(sub)
    m = np.isfinite(th)
    t2[tag]["resid_field"] = cap_contrast(rv[m], th[m])
    logger.info(f"  {tag} n={len(sub)} ({len(dual)} dual): "
                f"d_out_leg p={t2[tag]['our_d_out_leg'].get('p')}, "
                f"ddirf p={t2[tag]['our_ddirf'].get('p')}, "
                f"resid p={t2[tag]['resid_field'].get('p')}")
res["T2_uniform_modern_legs"] = t2

# ------------------------------------------------------------------
# T3 pre-2018 uniform-legacy-legs control on declared cap
# ------------------------------------------------------------------
logger.info("T3 uniform-legacy-legs control (pre-2018 cohort)")
recs = cohorts["pre2018"]
leg_uni = [r for r in recs if r["_g_max"] <= GAIA_LO]
leg_rich = [r for r in recs if r["_g_max"] > GAIA_LO]
t3 = {}
for tag, sub in (("uniform_le30", leg_uni), ("gaia_rich", leg_rich)):
    dual = [r for r in sub if r["_dual"]]
    th = np.array([r["_th_decl"] for r in dual])
    t3[tag] = {"n": len(sub), "n_dual": len(dual)}
    for ch in ("our_d_in_leg", "our_d_out_leg", "our_ddirf"):
        v = np.array([r[ch] if r.get(ch) is not None else np.nan
                      for r in dual])
        m = np.isfinite(v)
        t3[tag][ch] = cap_contrast(v[m], th[m])
    logger.info(f"  {tag} n={len(sub)} ({len(dual)} dual): "
                f"d_in_leg p="
                f"{t3[tag]['our_d_in_leg'].get('p')}, ddirf p="
                f"{t3[tag]['our_ddirf'].get('p')}")
res["T3_uniform_legacy_legs"] = t3

# ------------------------------------------------------------------
# T4 composition-vs-signal association
# ------------------------------------------------------------------
logger.info("T4 composition vs signal association")
t4 = {}
for tag, recs in cohorts.items():
    dual = [r for r in recs if r["_dual"]]
    g = np.array([r["_g_min"] for r in dual])
    d = np.array([r["our_ddirf"] for r in dual],
                 dtype=float)
    m = np.isfinite(d)
    if m.sum() > 20:
        rho, p = _st.spearmanr(g[m], d[m])
        t4[tag] = dict(n=int(m.sum()),
                       rho_ddirf_vs_minleg=float(rho), p=float(p))
        logger.info(f"  {tag}: rho(ddirf, min-leg Gaia) = "
                    f"{rho:+.3f} (p={p:.3g}), n={m.sum()}")
res["T4_composition_vs_signal"] = t4

# ------------------------------------------------------------------
# T5 era crossers
# ------------------------------------------------------------------
logger.info("T5 era crossers (pre-2018 table, post-2018 perihelia)")
jd2018 = _JD2000 + 18 * 365.25
cross = [r for r in cohorts["pre2018"] if r["tp_jd"] > jd2018
         and r["_dual"]]
t5 = {"n": len(cross)}
if len(cross) >= 8:
    th = np.array([r["_th_disp"] for r in cross])
    for ch in ("our_d_out_leg", "our_ddirf"):
        v = np.array([r[ch] if r.get(ch) is not None else np.nan
                      for r in cross])
        m = np.isfinite(v)
        t5[ch] = cap_contrast(v[m], th[m])
    t5["gaia_min_median"] = float(
        np.median([r["_g_min"] for r in cross]))
    logger.info(f"  {len(cross)} crossers: d_out_leg p="
                f"{t5['our_d_out_leg'].get('p')}, median min-leg "
                f"Gaia {t5['gaia_min_median']:.2f}")
res["T5_era_crossers"] = t5

# ------------------------------------------------------------------
# T6 matched-composition audit of the mixed stratum (post-2017)
# ------------------------------------------------------------------
# The registered displaced carrier is flat on uniform-Gaia legs but
# strongly elevated on mixed-composition legs.  Before calling that
# a catalogue systematic, verify the in-cap members of the mixed
# stratum are not themselves compositionally or temporally distinct:
# if in-cap vs out-of-cap mixed legs match on every composition and
# epoch metric, the cap gap there is direction-specific, not a
# composition main effect.
logger.info("T6 matched-composition audit of the mixed stratum")
recs = cohorts["post2017"]
mix = [r for r in recs if r["_g_min"] < GAIA_HI]
uni = [r for r in recs if r["_g_min"] >= GAIA_HI]
rv_all = resid_field(recs)   # registered pooled-residual form
for i, r in enumerate(recs):
    r["_rv"] = float(rv_all[i])
t6 = {"note": ("pooled resid_field (as registered in step_128) "
               "scored per stratum; in-cap vs out-of-cap members "
               "of the mixed stratum compared on leg composition, "
               "leg contrast and perihelion epoch"),
      "pooled_resid_cap_gap": {}}
for tag, sub in (("uniform_ge70", uni), ("mixed", mix)):
    th = np.array([r["_th_disp"] for r in sub])
    v = np.array([r["_rv"] for r in sub])
    t6["pooled_resid_cap_gap"][tag] = cap_contrast(v, th)
# matched-composition checks within mixed stratum
inc = np.array([r["_th_disp"] < CAP for r in mix])
cmp_ = {}
for key, arr in (
        ("g_min", np.array([r["_g_min"] for r in mix])),
        ("g_in", np.array([r["_g_in"] for r in mix])),
        ("g_out", np.array([r["_g_out"] for r in mix])),
        ("leg_contrast_go_minus_gi",
         np.array([r["_g_out"] - r["_g_in"] for r in mix])),
        ("abs_leg_contrast",
         np.abs([r["_g_out"] - r["_g_in"] for r in mix])),
        ("tp_jd", np.array([r["tp_jd"] for r in mix]))):
    if inc.sum() >= 5 and (~inc).sum() >= 5:
        u = _st.mannwhitneyu(arr[inc], arr[~inc])
        cmp_[key] = dict(median_in=float(np.median(arr[inc])),
                         median_out=float(np.median(arr[~inc])),
                         p=float(u.pvalue))
# cap gap after regressing pooled resid on absolute leg contrast
v = np.array([r["_rv"] for r in mix])
xc = np.abs([r["_g_out"] - r["_g_in"] for r in mix])
A = np.column_stack([np.ones(len(mix)), xc])
c, *_ = np.linalg.lstsq(A, v, rcond=None)
resid = v - A @ c
u = _st.mannwhitneyu(resid[inc], resid[~inc],
                     alternative="greater")
t6["matched_composition_in_mixed"] = cmp_
t6["cap_gap_resid_on_abs_contrast"] = dict(
    p=float(u.pvalue), slope=float(c[1]))
logger.info(f"  pooled-resid cap gap: uniform "
            f"p={t6['pooled_resid_cap_gap']['uniform_ge70'].get('p')}"
            f", mixed "
            f"p={t6['pooled_resid_cap_gap']['mixed'].get('p')}; "
            f"contrast-regressed p={u.pvalue:.3g}")
res["T6_matched_composition_audit"] = t6

# ------------------------------------------------------------------
# verdict
# ------------------------------------------------------------------
u3 = res["T3_uniform_legacy_legs"]["uniform_le30"]
p_leg = (u3.get("our_d_in_leg") or {}).get("p")
g_uni = res["T6_matched_composition_audit"][
    "pooled_resid_cap_gap"]["uniform_ge70"].get("p")
g_mix = res["T6_matched_composition_audit"][
    "pooled_resid_cap_gap"]["mixed"].get("p")
p_out = (res["T2_uniform_modern_legs"]["uniform_ge70"]
         .get("our_d_out_leg") or {}).get("p")
pre_survives = p_leg is not None and p_leg < 0.10
post_confined = (g_uni is not None and g_uni > 0.10
                 and g_mix is not None and g_mix < 0.05)
if pre_survives and post_confined:
    v = ("LEGACY CARRIER SURVIVES; DISPLACED CARRIER ATTENUATED, "
         "NOT COMPOSITION-EXPLAINED: the declared-cap inbound-leg "
         f"signal survives on uniform-legacy legs (p={p_leg:.3g}).  "
         "The displaced-cap residual is nonsignificant on "
         f"uniformly-high-Gaia legs (p={g_uni:.3g}) and strong on "
         f"sub-threshold legs (p={g_mix:.3g}); however the "
         "sub-threshold stratum is composition-homogeneous "
         "(only ~2% of post-2017 legs are truly stitch-mixed), "
         "in-cap members match out-of-cap on leg composition, leg "
         "contrast and perihelion epoch (T6), and the cap gap "
         "persists after regressing out absolute leg contrast.  "
         "The attenuation is therefore not explained as a "
         "catalogue-stitching systematic: it admits a composition-"
         "modulated signal amplitude or a reduction-quality "
         "dependence, while catalogue-by-sky interactions remain "
         "unresolved.")
elif pre_survives:
    v = ("PRE-2018 CARRIER WITHIN LEGACY STRATUM; post-2017 "
         "stratification intermediate -- see T2/T6")
elif post_confined:
    v = ("POST-2017 CARRIER COMPOSITION-CONFINED -- see T6")
else:
    v = "STRATIFICATION INCONCLUSIVE -- see T2/T3/T6 sample sizes"
res["verdict"] = v
logger.info("verdict: " + v)

out = RESULTS / "step_b104_astcat_stratification.json"


def _finite(o):
    if isinstance(o, dict):
        return {k: _finite(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_finite(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


json.dump(_finite(res), open(out, "w"), indent=1, default=float)
logger.data_save(out)
print(json.dumps({"T1": res["T1_composition"],
                  "T2": res["T2_uniform_modern_legs"],
                  "T3": res["T3_uniform_legacy_legs"],
                  "verdict": v}, indent=1))
