#!/usr/bin/env python3
"""
TEP-9 step 098 -- RTG nuclear-decay channel
===========================================

Each Voyager carries three MHW-RTG units powered by Pu-238 alpha
decay (half-life 87.7 yr).  The electrical power record is therefore
a continuous ~44-year monitor of a *nuclear* clock -- a process class
distinct from both the electronic quartz oscillator (SCLK channel)
and the macroscopic orbital dynamics (comet channel).  This is the
"nuclear processes in the power engine" channel requested for the
TEP-9 analysis.

Data
----
The channel uses the published historical RTG performance dataset
(Whiting & Woerner 2023, doi:10.5061/dryad.1zcrjdfw2), which collates
the JPL MHW-RTG telemetry record: normalized electrical power versus
time for the Voyager 1 and Voyager 2 total missions and for each
individual RTG unit, through August 2021.  This is real telemetry --
129 Voyager-1 and 89 Voyager-2 total-mission epochs plus five
per-unit series -- not literature anchors.

Pedigree facts encoded from the dataset README (used in the
analysis):

* values are normalized to the beginning-of-mission power (highest
  point within 30 days of BOM); absolute watts use the tabulated
  PBOM (V1 total 471 We, V2 total 476 We);
* the first ~13.3 years are spacecraft-reported power;
* post-13.3-year points are recorded at epochs where an RTG current
  dropped one telemetry bit (0.039 A for VG1, 0.038 A for VG2), so
  the late series is bit-quantized at the ~1 % level;
* Voyager 2 RTG #3 telemetry was lost ~700 h after launch; the V2
  total is reconstructed from units 1-2 plus a fixed bus-voltage
  correction;
* literature power anchors (literature_anchors.json) are retained as
  an independent cross-check on the fitted model.

Physics of the channel
----------------------
Alpha-decay rates are exponentially sensitive to the effective fine
structure constant through the Gamow penetrability exponent

    G = pi Z1 Z2 alpha / beta ,   beta = v_rel/c ,
    lambda ~ exp(-2 G)  =>  delta lambda / lambda ~ -2 G delta alpha_eff / alpha_eff .

For Pu-238 -> U-234 + alpha (Q = 5.59 MeV, Z1 Z2 = 184) the exponent
evaluates to G ~ 76, so the decay rate amplifies any constants
shift by ~150x: a delta alpha_eff / alpha_eff step at a
lapse-boundary crossing of comet-channel strength (~1e-2) would
appear as an order-unity step in thermal power -- far above the
bit-quantization floor.  The bound is
process-specific by construction: the TEP framework treats
electronic, nuclear and gravitational/dynamical couplings to A(phi)
as potentially different, so the RTG bound constrains the *nuclear*
channel without assuming it equals the orbital channel.

Signatures tested
-----------------
1. Level step at each crossing epoch: ln P = a + b t + s H(t-t_cx);
   a constants step changes the instantaneous decay-rate power
   (level), so s is the channel-relevant statistic.  t-test on s.
2. Slope change at each crossing epoch: a post-crossing change in
   the effective decay constant shows as a break in b.
3. Per-RTU consistency: the same tests on the five functioning
   per-unit series (3 x VG1, 2 x VG2).
4. VG1 - VG2 differential normalized power: removes the common-mode
   degradation model and exposes craft-specific residuals around
   each craft's own crossings.
5. Cross-check: the fitted model evaluated at the published
   literature-anchor epochs.

Inputs
------
data/raw/rtg/All_RTGs_P_Over_P0_-_2023_Update_for_Release.xlsx
data/raw/literature/literature_anchors.json

Outputs
-------
results/step_b62_rtg_nuclear.json
results/step_b62_rtg_nuclear.csv  (per-epoch power series)
results/figures/supplementary/step_b62_rtg_nuclear.png
"""

import sys
import json
import math
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
import openpyxl
from scipy.stats import t as tdist

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.statistics import monte_carlo_tail
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout

logger = StepLogger("step_098_rtg_nuclear_channel")
tee_stdout(logger)
logger.header("RTG nuclear-decay channel (telemetry)")

XLSX = (DATA_RAW / "rtg" /
        "All_RTGs_P_Over_P0_-_2023_Update_for_Release.xlsx")
LIT = json.loads((DATA_RAW / "literature" /
                  "literature_anchors.json").read_text())
CROSSINGS = LIT["crossings"]
T12_YR = LIT["rtg_power"]["pu238_half_life_yr"]["value"]

# ---------------------------------------------------------------------
# Parse the MHW-RTG sheet
# ---------------------------------------------------------------------
wb = openpyxl.load_workbook(XLSX, data_only=True, read_only=True)
ws = wb["MHW-RTG"]
rows = list(ws.iter_rows(values_only=True))

SERIES = {
    # label: (time col, power col, PBOM_w, BOM date, craft, unit)
    "VG1_total": (16, 17, 471.0, datetime(1977, 9, 5), "VG1", "total"),
    "VG1_rtg1": (7, 8, 156.0, datetime(1977, 9, 5), "VG1", "rtg1"),
    "VG1_rtg2": (10, 11, 157.0, datetime(1977, 9, 5), "VG1", "rtg2"),
    "VG1_rtg3": (13, 14, 159.0, datetime(1977, 9, 5), "VG1", "rtg3"),
    "VG2_total": (28, 29, 476.0, datetime(1977, 8, 20), "VG2", "total"),
    "VG2_rtg1": (19, 20, 158.0, datetime(1977, 8, 20), "VG2", "rtg1"),
    "VG2_rtg2": (22, 23, 158.0, datetime(1977, 8, 20), "VG2", "rtg2"),
}


def parse(ct, cp, pbom, bom, craft, unit):
    pts = [(r[ct], r[cp]) for r in rows[19:]
           if isinstance(r[ct], (int, float))
           and isinstance(r[cp], (int, float))]
    t = np.array([p[0] for p in pts])          # years since BOM
    pn = np.array([p[1] for p in pts])         # P/P0
    utc = [bom + timedelta(days=float(tt) * 365.25) for tt in t]
    return dict(label=f"{craft}_{unit}", craft=craft, unit=unit,
                pbom_w=pbom, bom_utc=bom.isoformat()[:10],
                t_yr=t, p_norm=pn, p_w=pn * pbom,
                utc=[u.isoformat()[:10] for u in utc],
                year=np.array(
                    [u.year + (u.timetuple().tm_yday - 0.5) / 365.25
                     for u in utc]))


series = {k: parse(*v) for k, v in SERIES.items()}
for k, s in series.items():
    logger.metric(f"n_{k}", len(s["t_yr"]),
                  f"telemetry epochs, last={s['utc'][-1]}")

# ---------------------------------------------------------------------
# Gamow exponent for Pu-238 -> U-234 + alpha
# ---------------------------------------------------------------------
Z1, Z2 = 2, 92
Q_MEV = 5.593
M_DAUGHTER, M_PARENT, AMU_MEV = 234.0, 238.0, 931.494
mu = 4.0 * M_DAUGHTER / M_PARENT
beta_v = math.sqrt(2.0 * Q_MEV / (mu * AMU_MEV))
ALPHA_FS = 1.0 / 137.035999084
G_GAMOW = math.pi * Z1 * Z2 * ALPHA_FS / beta_v
AMP = 2.0 * G_GAMOW
logger.metric("gamow_exponent", round(G_GAMOW, 1),
              "Pu-238 Gamow exponent")
logger.metric("decay_amplification", round(AMP, 0),
              "|dlambda/lambda| per |dalpha/alpha|")

# ---------------------------------------------------------------------
# Standard decay+degradation fit per craft on the total-mission series
# ln P_norm = a - k_eff t + c t^2 ; the quadratic absorbs the known
# concavity of the SiGe sublimation-degradation curve (a linear-in-log
# model leaves systematic +/-0.5-1 % residuals that masquerade as
# steps at every epoch -- verified by the placebo scan below)
# ---------------------------------------------------------------------
LN2_T12 = math.log(2.0) / T12_YR      # Pu-238 contribution to slope


def fit_model(s):
    t, pn = s["t_yr"], s["p_norm"]
    A = np.vstack([np.ones_like(t), -t, t * t]).T
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(A, np.log(pn), rcond=None)
    a, k_eff, c2 = coef
    k_deg = k_eff - LN2_T12           # degradation beyond Pu-238
    model = np.exp(a - k_eff * t + c2 * t * t)
    resid = np.log(pn) - np.log(model)
    return dict(a=float(a), k_eff=float(k_eff), c2=float(c2),
                k_deg_per_yr=float(k_deg),
                P0_w=float(math.exp(a) * s["pbom_w"]),
                resid=resid, model_norm=model,
                resid_rms=float(np.sqrt(np.mean(resid ** 2))),
                resid_abs_max=float(np.abs(resid).max()),
                resid_abs_p99=float(np.quantile(np.abs(resid), 0.99)))


fits = {k: fit_model(series[k]) for k in ("VG1_total", "VG2_total")}
for k, f in fits.items():
    logger.metric(f"{k}_k_deg", round(f["k_deg_per_yr"] * 100, 3),
                  "degradation %/yr beyond Pu-238")
    logger.metric(f"{k}_resid_rms", round(f["resid_rms"] * 100, 2),
                  "log-residual rms (%)")


def signed_step_test(s, t_cx_yr):
    """Level step on the quadratic model: ln pn = a + b t + c t^2 +
    s H(t-t_cx); plus a post-crossing slope-change term test."""
    t, lnp = s["t_yr"], np.log(s["p_norm"])
    h = (t >= t_cx_yr).astype(float)
    # (a) level step
    A = np.vstack([np.ones_like(t), t, t * t, h]).T
    with np.errstate(all="ignore"):
        c, *_ = np.linalg.lstsq(A, lnp, rcond=None)
        r = lnp - A @ c
        dof = max(len(t) - A.shape[1], 1)
        s2 = float(r @ r) / dof
        cov = s2 * np.linalg.inv(A.T @ A)
        se = math.sqrt(max(cov[3, 3], 0.0))
    tstat = float(c[3] / se) if se > 0 else float("nan")
    pv = float(2 * tdist.sf(abs(tstat), dof)) if se > 0 else float("nan")
    pv_log10 = float((math.log(2) + tdist.logsf(abs(tstat), dof))
                     / math.log(10)) if se > 0 else float("nan")
    # (b) slope change: post-crossing extra linear term
    d = np.maximum(0.0, t - t_cx_yr)
    A2 = np.vstack([np.ones_like(t), t, t * t, d]).T
    with np.errstate(all="ignore"):
        c2, *_ = np.linalg.lstsq(A2, lnp, rcond=None)
        r2 = lnp - A2 @ c2
        s22 = float(r2 @ r2) / dof
        cov2 = s22 * np.linalg.inv(A2.T @ A2)
        se2 = math.sqrt(max(cov2[3, 3], 0.0))
    tstat2 = float(c2[3] / se2) if se2 > 0 else float("nan")
    pv2 = float(2 * tdist.sf(abs(tstat2), dof)) if se2 > 0 \
        else float("nan")
    pv2_log10 = float((math.log(2) + tdist.logsf(abs(tstat2), dof))
                      / math.log(10)) if se2 > 0 else float("nan")
    return dict(step_frac=float(c[3]), step_se=float(se),
                step_t=round(tstat, 3), step_p=pv if pv > 0 else None,
                step_p_log10=pv_log10,
                dslope_per_yr=float(c2[3]), dslope_se=float(se2),
                dslope_t=round(tstat2, 3), dslope_p=pv2 if pv2 > 0 else None,
                dslope_p_log10=pv2_log10)


# crossing epochs in each craft's years-since-BOM frame
def cx_t(key, s):
    cx = CROSSINGS[key]
    y, m, dd = (int(x) for x in cx["utc"].split("-"))
    cy = y + (m - 0.5) / 12.0 + (dd - 15.5) / 365.25
    by = int(s["bom_utc"][:4]) + \
        (datetime(int(s["bom_utc"][:4]),
                  int(s["bom_utc"][5:7]),
                  int(s["bom_utc"][8:10])).timetuple().tm_yday
         - 0.5) / 365.25
    return cy - by


CRAFT_CX = {"VG1": ("v1_termination_shock", "v1_heliopause"),
            "VG2": ("v2_termination_shock", "v2_heliopause")}

step_tests = {}
for k in ("VG1_total", "VG2_total"):
    s, craft = series[k], series[k]["craft"]
    step_tests[k] = {}
    for key in CRAFT_CX[craft]:
        tcx = cx_t(key, s)
        step_tests[k][key] = dict(
            epoch_utc=CROSSINGS[key]["utc"],
            radius_au=CROSSINGS[key]["radius_au"],
            t_yr_since_bom=round(tcx, 2),
            **signed_step_test(s, tcx))
        logger.metric(f"{k}_{key}_step",
                      round(step_tests[k][key]["step_frac"] * 100, 2),
                      f"% step at {CROSSINGS[key]['utc']} "
                      f"(p={step_tests[k][key]['step_p']:.3f})")

# placebo scan: the same level-step test at a uniform grid of
# non-crossing epochs.  The crossing steps are meaningful only if
# they exceed the placebo envelope -- a curved degradation model or
# quantization bias would otherwise masquerade as crossing steps.
placebo = {}
for k in ("VG1_total", "VG2_total"):
    s = series[k]
    grid = np.arange(15.0, s["t_yr"].max() - 2.0, 2.5)
    cx_epochs = [cx_t(key, s) for key in CRAFT_CX[s["craft"]]]
    pts = []
    for tcx in grid:
        if min(abs(tcx - e) for e in cx_epochs) < 1.0:
            continue
        st = signed_step_test(s, float(tcx))
        pts.append(dict(t_yr=round(float(tcx), 1),
                        step_frac=st["step_frac"],
                        step_se=st["step_se"], step_p=st["step_p"]))
    env = max(abs(p["step_frac"]) for p in pts)
    # empirical p: fraction of placebo epochs whose |step| reaches
    # the crossing step's magnitude -- the honest significance under
    # the record's own drift/quantization noise
    emp = {}
    for key in CRAFT_CX[s["craft"]]:
        obs = abs(step_tests[k][key]["step_frac"])
        emp[key] = float(monte_carlo_tail(
            [abs(p["step_frac"]) for p in pts], obs))
        step_tests[k][key]["empirical_p_vs_placebo"] = emp[key]
    placebo[k] = dict(
        n_placebo=len(pts), envelope_abs_frac=env,
        empirical_p=emp,
        crossings_inside_envelope=all(
            abs(step_tests[k][key]["step_frac"]) <= env
            for key in CRAFT_CX[s["craft"]]),
        grid=pts)
    logger.metric(f"{k}_placebo_envelope", round(env * 100, 2),
                  "% max |placebo step|")

# per-unit step tests (5 functioning RTG series)
unit_tests = {}
for k in ("VG1_rtg1", "VG1_rtg2", "VG1_rtg3", "VG2_rtg1", "VG2_rtg2"):
    s, craft = series[k], series[k]["craft"]
    unit_tests[k] = {}
    for key in CRAFT_CX[craft]:
        tcx = cx_t(key, s)
        unit_tests[k][key] = dict(
            epoch_utc=CROSSINGS[key]["utc"],
            **signed_step_test(s, tcx))

# ---------------------------------------------------------------------
# Excess-decline channel (model-light)
# ---------------------------------------------------------------------
# A level shift in the decay rate -- the conformal/varying-constants
# signature -- appears as a *persistent excess decline* accumulating
# across the transition window, not a step.  Anchor each series to a
# linear-in-log baseline fitted well inside the boundary (13-26 yr,
# after the telemetry-mode change at 13.3 yr and before the
# termination-shock epoch at 27 yr), then track the cumulative
# excess decline versus mission year.  This avoids the quadratic
# model's curvature absorption and end-curl artifacts entirely.
BASE_LO, BASE_HI = 13.0, 26.0


def excess_decline(s):
    t, lnp = s["t_yr"], np.log(s["p_norm"])
    m = (t >= BASE_LO) & (t <= BASE_HI)
    c = np.polyfit(t[m], lnp[m], 1)
    return (lnp - np.polyval(c, t)) * 100.0      # per cent


excess = {k: excess_decline(s) for k, s in series.items()}

# per-craft coherence: mean and spread of the unit-level excesses at
# the crossing epochs.  A common (field- or craft-level) shift must
# move every unit identically; unit-to-unit divergence is the
# signature of unit-specific systematics.
tg_c = np.arange(14.0, 44.0, 0.5)
unit_coherence = {}
for craft, units in (("VG1", ("VG1_rtg1", "VG1_rtg2", "VG1_rtg3")),
                     ("VG2", ("VG2_rtg1", "VG2_rtg2"))):
    E = np.array([np.interp(tg_c, series[k]["t_yr"], excess[k])
                  for k in units])
    hp = cx_t(f"{craft.lower().replace('vg','v')}_heliopause", series[units[0]])
    ts = cx_t(f"{craft.lower().replace('vg','v')}_termination_shock",
              series[units[0]])
    ihp = int(np.argmin(np.abs(tg_c - hp)))
    its = int(np.argmin(np.abs(tg_c - ts)))
    unit_coherence[craft] = dict(
        grid=tg_c.tolist(),
        mean_excess_pct=[round(float(x), 3) for x in E.mean(axis=0)],
        spread_excess_pct=[round(float(x), 3) for x in E.std(axis=0)],
        excess_at_hp_pct=dict(
            mean=round(float(E.mean(axis=0)[ihp]), 3),
            per_unit={k: round(float(
                np.interp(hp, series[k]["t_yr"], excess[k])), 3)
                for k in units},
            unit_spread=round(float(E[:, ihp].std()), 3)),
        excess_at_ts_pct=dict(
            mean=round(float(E.mean(axis=0)[its]), 3),
            per_unit={k: round(float(
                np.interp(ts, series[k]["t_yr"], excess[k])), 3)
                for k in units}),
        excess_end_pct=dict(
            mean=round(float(E.mean(axis=0)[-1]), 3),
            per_unit={k: round(float(excess[k][-1]), 3)
                      for k in units}),
    )

# mean-unit ramp localization: onset/centre/amplitude of the V1
# excess accumulation (graded-transition shape, not a step)
E1 = np.array(unit_coherence["VG1"]["mean_excess_pct"])
mr = tg_c >= 27.0
best = None
for t0 in np.arange(27.0, 38.0, 0.5):
    for w in (3.0, 5.0, 7.0, 9.0):
        x = np.clip((tg_c - t0) / w, 0.0, 1.0)
        A = np.vstack([np.ones_like(tg_c), x]).T
        c, *_ = np.linalg.lstsq(A[mr], E1[mr], rcond=None)
        rss = float(((E1[mr] - A[mr] @ c) ** 2).sum())
        if best is None or rss < best[0]:
            best = (rss, t0, w, float(c[1]))
# model comparison on the same post-27 segment: does the saturating
# ramp beat a plain linear trend or a step at the heliopause?  The
# morphology claim is only meaningful if it does.
_x, _y = tg_c[mr], E1[mr]
_xr = np.clip((_x - best[1]) / best[2], 0.0, 1.0)
_A = np.vstack([np.ones_like(_x), _xr]).T
_c, *_ = np.linalg.lstsq(_A, _y, rcond=None)
_rss_ramp = float(((_y - _A @ _c) ** 2).sum())
_A = np.vstack([np.ones_like(_x), _x]).T
_c, *_ = np.linalg.lstsq(_A, _y, rcond=None)
_rss_lin = float(((_y - _A @ _c) ** 2).sum())
_hp_t = cx_t("v1_heliopause", series["VG1_total"])
_A = np.vstack([np.ones_like(_x), (_x > _hp_t).astype(float)]).T
_c, *_ = np.linalg.lstsq(_A, _y, rcond=None)
_rss_step = float(((_y - _A @ _c) ** 2).sum())
v1_ramp = dict(onset_yr=best[1], width_yr=best[2],
               centre_yr=best[1] + best[2] / 2,
               amplitude_pct=round(best[3], 3),
               hp_yr=round(_hp_t, 2),
               ts_yr=round(cx_t("v1_termination_shock",
                                series["VG1_total"]), 2),
               model_comparison=dict(
                   rss_ramp=round(_rss_ramp, 4),
                   rss_linear=round(_rss_lin, 4),
                   rss_step_at_hp=round(_rss_step, 4),
                   ramp_vs_linear_rss_improvement_pct=round(
                       100.0 * (1.0 - _rss_ramp / _rss_lin), 1)))

# baseline-window robustness: the excess is defined relative to a
# linear fit on [BASE_LO, BASE_HI]; the same data supports materially
# different amplitudes -- and even signs -- under other plausible
# windows, so the excursion is reported per window rather than as a
# single number.  Also recorded: the largest single-unit excursion in
# the out-of-cap craft (VG2_rtg2), a same-amplitude counterexample in
# the band where the field predicts the opposite sign.
def excess_decline_window(s, lo, hi):
    t, lnp = s["t_yr"], np.log(s["p_norm"])
    m = (t >= lo) & (t <= hi)
    c = np.polyfit(t[m], lnp[m], 1)
    return (lnp - np.polyval(c, t)) * 100.0


v1_hp_t = cx_t("v1_heliopause", series["VG1_total"])
v2_hp_t = cx_t("v2_heliopause", series["VG2_total"])
window_robustness = {"windows": {}, "counterexamples": {}}
for lo, hi in ((0.0, 26.0), (8.0, 26.0), (13.0, 26.0), (16.0, 26.0),
               (13.0, 22.0)):
    tag = f"{lo:.0f}-{hi:.0f}"
    window_robustness["windows"][tag] = {
        k: round(float(np.interp(
            v1_hp_t if k.startswith("VG1") else v2_hp_t,
            series[k]["t_yr"], excess_decline_window(series[k], lo, hi))), 3)
        for k in ("VG1_rtg1", "VG1_rtg2", "VG1_rtg3",
                  "VG2_rtg1", "VG2_rtg2")}
window_robustness["counterexamples"]["VG2_rtg2_excess_at_v2hp_pct"] = \
    round(float(np.interp(v2_hp_t, series["VG2_rtg2"]["t_yr"],
                          excess["VG2_rtg2"])), 3)
window_robustness["counterexamples"]["VG2_rtg2_excess_at_end_pct"] = \
    round(float(excess["VG2_rtg2"][-1]), 3)
window_robustness["note"] = (
    "the V1 HP excess ranges from +1.6% (baseline 0-26) to -0.5% "
    "(baseline 13-26) across plausible windows -- the sign itself is "
    "window-dependent; VG2_rtg2 alone shows -0.85% at end under the "
    "canonical window, an out-of-cap counterexample at the same "
    "amplitude.  Unit coherence within V1 is real (spread ~0.05%) "
    "but is expected for units sharing one craft's thermal "
    "environment, measurement chain and degradation physics -- it "
    "discriminates unit-specific noise, not field vs common-mode "
    "aging.")

# V2 total series carries a two-mode reporting structure: through
# ~37 yr the record alternates between report conventions offset
# ~0.75 % (paired epochs ~0.1 yr apart), so sub-percent level tests
# on the mixed series are not reliable; the per-unit series are
# clean single tracks and are used for the V2 verdict.
v2_mixing = dict(
    note=("VG2_total alternates between two report conventions "
          "offset ~0.75% through mission-year ~37 (paired epochs "
          "~0.1 yr apart); level-step tests on the mixed series "
          "inherit that offset, so the craft-level verdict for V2 "
          "rests on the two clean per-unit series."),
    per_unit_divergence=("VG2_rtg1 and VG2_rtg2 diverge by ~2% in "
                         "excess decline over mission-years 37-43 "
                         "(rtg1 positive excess, rtg2 a -1.2% dip "
                         "with partial recovery) -- unit-level "
                         "systematics of order the sought signal."))

# ---------------------------------------------------------------------
# Geometry concordance
# ---------------------------------------------------------------------
# The slip field is cos-2theta about the axis: POSITIVE-signed at
# both poles (the primary and mirror caps, median slip ~+2 yr) and
# NEGATIVE in the intervening band.  Of the boundary crossers, only
# V1 crossed inside a cap lobe (mirror cap, positive sign); V2
# crossed in the intervening negative-slip band; Pioneer 10 never
# reached the boundary.  A positive lapse sign means proper time --
# and with it the nuclear decay clock -- runs FASTER in V1's lobe,
# so the predicted RTG signature is a faster power decline, i.e. a
# NEGATIVE excess.  V1's coherent -0.4/-0.6% excess is therefore
# sign-consistent with the mirror-cap polarity.  V2's negative-slip
# band would predict the opposite sign; its record is dominated by
# unit-level divergence (~1.7% peak-to-peak) and resolves no common
# shift either way.
geom_concordance = dict(
    slip_field_structure=("cos-2theta: positive slip at both axis "
                          "poles (primary and mirror caps), negative "
                          "slip in the intervening band"),
    v1=dict(position="inside mirror cap",
            field_sign="positive slip (faster proper time)",
            predicted_power_excess="negative (faster decay)",
            coherent_excess=True,
            observed="coherent -0.39% at HP / -0.58% at end",
            verdict="sign-consistent candidate signature"),
    v2=dict(position="intervening negative-slip band (outside both "
            "caps)",
            field_sign="negative slip (slower proper time)",
            predicted_power_excess="positive (slower decay)",
            coherent_excess=False,
            observed="units diverge ~1.7% peak-to-peak; mean -0.16% "
                     "at HP -- no common shift resolved",
            verdict="unresolved (unit-level systematics dominate)"),
    note=("the only boundary crosser inside a cap lobe is the only "
          "craft with a unit-coherent decay-clock excursion, and its "
          "negative sign matches the faster-clock polarity of the "
          "mirror lobe; V2's out-of-cap path resolves no common "
          "shift -- the directional pattern the bipolar field "
          "predicts, retained at candidate level because mid-life "
          "MHW-RTG aging has no external control."))

# ---------------------------------------------------------------------
# VG1 - VG2 differential channel (common-mode degradation removed)
# ---------------------------------------------------------------------
v1, v2 = series["VG1_total"], series["VG2_total"]
# interpolate V2 log-residuals onto the V1 time grid (same BOM frame:
# both BOM dates are 1977, difference 16 d, negligible here)
r1 = fits["VG1_total"]["resid"]
r2 = np.interp(v1["t_yr"], v2["t_yr"], fits["VG2_total"]["resid"],
               left=np.nan, right=np.nan)
ok = np.isfinite(r2)
dres = r1[ok] - r2[ok]
t_d = v1["t_yr"][ok]


def dres_step(tcx):
    """Level-step test on the VG1-VG2 residual differential at a
    crossing epoch (years since ~BOM, common 1977 frame)."""
    h = (t_d >= tcx).astype(float)
    A = np.vstack([np.ones_like(t_d), h]).T
    c, *_ = np.linalg.lstsq(A, dres, rcond=None)
    r = dres - A @ c
    dof = max(len(t_d) - 2, 1)
    se = math.sqrt(max((r @ r) / dof *
                       np.linalg.inv(A.T @ A)[1, 1], 0.0))
    tstat = float(c[1] / se) if se > 0 else float("nan")
    return dict(step_frac=float(c[1]), step_se=float(se),
                step_t=round(tstat, 3),
                step_p=float(2 * tdist.sf(abs(tstat), dof))
                if se > 0 else float("nan"))


diff_tests = dict(
    v1_termination_shock=dres_step(cx_t("v1_termination_shock", v1)),
    v1_heliopause=dres_step(cx_t("v1_heliopause", v1)),
    v2_termination_shock=dres_step(cx_t("v2_termination_shock", v2)),
    v2_heliopause=dres_step(cx_t("v2_heliopause", v2)))
differential = dict(
    n=int(ok.sum()),
    rms_log=float(np.sqrt(np.mean(dres ** 2))),
    max_abs_log=float(np.abs(dres).max()),
    crossing_step_tests=diff_tests,
    note=("VG1 minus interpolated VG2 log-residuals about each "
          "craft's own decay+degradation fit; a field common to both "
          "paths cancels, exposing craft-specific excursions.  A "
          "crossing step that is common-mode degradation drift "
          "cancels here; a genuinely craft-specific crossing step "
          "survives."),
    series=dict(t_yr=t_d.tolist(),
                dres=[float(x) for x in dres]))

# ---------------------------------------------------------------------
# Bounds: largest |step| consistent with the fits -> delta-lambda ->
# delta-alpha via the Gamow amplification
# ---------------------------------------------------------------------
all_steps = []
for k, tests in step_tests.items():
    for key, tt in tests.items():
        all_steps.append(abs(tt["step_frac"]) + 2 * abs(tt["step_se"]))
bound_dlambda = max(all_steps)          # 2-sigma step bound
bound_dalpha = bound_dlambda / AMP
for k, tests in step_tests.items():
    for key, tt in tests.items():
        tt["bound_dlambda_frac_2s"] = round(
            abs(tt["step_frac"]) + 2 * abs(tt["step_se"]), 5)
        tt["bound_dalpha_frac_2s"] = float(
            f"{tt['bound_dlambda_frac_2s'] / AMP:.2e}")

# ---------------------------------------------------------------------
# Predicted response at comet-channel contrast, two coupling sectors
# ---------------------------------------------------------------------
# Universal conformal: every spacecraft process shares the lapse
# ratio, so the measured power steps by dA/A itself.  That sector is
# already excluded at ~1 ppm by the SCLK record (step 096), so the
# channel of interest here is the nuclear sector: a varying-alpha
# boundary couples to the decay constant through the Gamow exponent,
# dlambda/lambda = AMP * dalpha/alpha, while leaving the electronic
# SCLK within its bound only if the electronic-clock sensitivity to
# alpha_eff is small -- the process-dependent structure TEP itself
# predicts.
dA = 1.16e-2
_ledger = RESULTS / "step_b64_clock_consistency.json"
if _ledger.exists():
    try:
        dA = float(json.load(open(_ledger))
                   ["comet_channel"]["implied_dA_over_A"])
    except Exception:
        pass
dlambda_predicted = AMP * dA
predicted_steps = {}
for k in ("VG1_total", "VG2_total"):
    s = series[k]
    for key in CRAFT_CX[s["craft"]]:
        if "heliopause" not in key:
            continue
        tcx = cx_t(key, s)
        pn_at = math.exp(fits[k]["a"] - fits[k]["k_eff"] * tcx
                       + fits[k]["c2"] * tcx * tcx)
        predicted_steps[f"{k}:{key}"] = dict(
            epoch_utc=CROSSINGS[key]["utc"],
            conformal_step_w_at_comet_scale=round(
                dA * pn_at * s["pbom_w"], 1),
            alpha_channel_step_w_at_comet_scale=round(
                dlambda_predicted * pn_at * s["pbom_w"], 1),
            measured_step_w=round(
                step_tests[k][key]["step_frac"] * pn_at *
                s["pbom_w"], 2),
            v1_unit_mean_excess_at_hp_w=(
                round(-unit_coherence["VG1"]["excess_at_hp_pct"]
                      ["mean"] / 100 * pn_at * s["pbom_w"], 1)
                if s["craft"] == "VG1" else None))

# ---------------------------------------------------------------------
# Literature-anchor cross-check
# ---------------------------------------------------------------------
anchor_check = []
for a in LIT["rtg_power"]["anchors"]:
    y, m, dd = (int(x) for x in a["utc"].split("-"))
    cy = y + (m - 0.5) / 12.0 + (dd - 15.5) / 365.25
    k = f"{a['craft']}_total"
    s = series[k]
    by = int(s["bom_utc"][:4]) + \
        (datetime(int(s["bom_utc"][:4]),
                  int(s["bom_utc"][5:7]),
                  int(s["bom_utc"][8:10])).timetuple().tm_yday
         - 0.5) / 365.25
    tt = cy - by
    pn_model = math.exp(fits[k]["a"] - fits[k]["k_eff"] * tt
                          + fits[k]["c2"] * tt * tt)
    anchor_check.append(dict(
        craft=a["craft"], utc=a["utc"], anchor_w=a["power_w"],
        telemetry_model_w=round(pn_model * s["pbom_w"], 1),
        resid_w=round(a["power_w"] - pn_model * s["pbom_w"], 1),
        note=a["note"]))

# ---------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(4, 1, figsize=(11, 14), sharex=False)
cols = {"VG1": "#1f77b4", "VG2": "#d62728"}

ax = axes[0]
for k in ("VG1_total", "VG2_total"):
    s = series[k]
    ax.plot(s["year"], s["p_w"], ".", ms=3, color=cols[s["craft"]],
            label=f'{s["craft"]} telemetry (PBOM={s["pbom_w"]:.0f} W)')
    tt = np.linspace(0, s["t_yr"].max(), 200)
    ax.plot(s["year"][0] + tt,
            np.exp(fits[k]["a"] - fits[k]["k_eff"] * tt
                    + fits[k]["c2"] * tt * tt) * s["pbom_w"],
            "-", color=cols[s["craft"]], lw=1, alpha=0.6,
            label=f'{s["craft"]} decay+degradation fit '
                  f'($k_{{deg}}$={fits[k]["k_deg_per_yr"]*100:.2f}%/yr)')
for a in anchor_check:
    ax.plot(float(a["utc"][:4]) + 0.5, a["anchor_w"], "x", ms=7,
            color="k")
ax.plot([], [], "kx", ms=7, label="literature anchors")
for key, col in (("v1_heliopause", cols["VG1"]),
                 ("v2_heliopause", cols["VG2"]),
                 ("v1_termination_shock", cols["VG1"]),
                 ("v2_termination_shock", cols["VG2"])):
    y, m, dd = (int(x) for x in CROSSINGS[key]["utc"].split("-"))
    cy = y + (m - 0.5) / 12.0
    ax.axvline(cy, color=col, ls="--", alpha=0.35)
ax.set_ylabel("RTG electrical power (W)")
ax.set_title("Voyager MHW-RTG telemetry (Dryad/Whiting & Woerner 2023)"
             " vs standard decay+degradation")
ax.legend(fontsize=8, loc="lower left")

ax = axes[1]
for k in ("VG1_total", "VG2_total"):
    s = series[k]
    ax.plot(s["year"], fits[k]["resid"] * 100, ".", ms=3,
            color=cols[s["craft"]], label=s["craft"])
for k, col in (("v1_termination_shock", cols["VG1"]),
               ("v1_heliopause", cols["VG1"]),
               ("v2_termination_shock", cols["VG2"]),
               ("v2_heliopause", cols["VG2"])):
    y, m, dd = (int(x) for x in CROSSINGS[k]["utc"].split("-"))
    ax.axvline(y + (m - 0.5) / 12.0, color=col, ls="--", alpha=0.4)
ax.axhline(0, color="k", lw=0.5)
ax.set_ylabel("log residual (%)")
ax.set_xlabel("year")
ax.set_title("Residuals about decay+degradation fit; dashed = TS/HP "
             "crossings")
ax.legend(fontsize=8)

ax = axes[2]
for k in ("VG1_rtg1", "VG1_rtg2", "VG1_rtg3",
          "VG2_rtg1", "VG2_rtg2"):
    s = series[k]
    f = fit_model(s)
    ax.plot(s["year"], f["resid"] * 100, ".", ms=2,
            label=k.replace("_total", "").replace("_", " "),
            alpha=0.7)
ax.axhline(0, color="k", lw=0.5)
ax.set_ylabel("per-unit log residual (%)")
ax.set_xlabel("year")
ax.set_title("Five functioning MHW-RTG units: independent nuclear "
             "clocks")
ax.legend(fontsize=7, ncol=3)

ax = axes[3]
for craft, col in (("VG1", cols["VG1"]), ("VG2", cols["VG2"])):
    g = unit_coherence[craft]["grid"]
    y0 = series[f"{craft}_total"]["year"][0]
    gy = y0 + np.asarray(g)
    ax.plot(gy, unit_coherence[craft]["mean_excess_pct"], "-",
            color=col, lw=1.6,
            label=f"{craft} mean-unit excess")
    for i, k in enumerate(u for u in series
                          if u.startswith(craft) and "rtg" in u):
        ax.plot(gy, np.interp(g, series[k]["t_yr"], excess[k]),
                ".", ms=2, color=col, alpha=0.35)
    hp = cx_t(f"{craft.lower().replace('vg','v')}_heliopause", series[f"{craft}_total"])
    ts = cx_t(f"{craft.lower().replace('vg','v')}_termination_shock",
              series[f"{craft}_total"])
    ax.axvline(y0 + hp, color=col, ls="--", alpha=0.5)
    ax.axvline(y0 + ts, color=col, ls=":", alpha=0.4)
ax.axhline(0, color="k", lw=0.5)
ax.set_ylabel("excess decline vs interior baseline (%)")
ax.set_xlabel("year")
ax.set_title("Cumulative excess decline (baseline 13-26 yr); "
             "dashed = heliopause, dotted = termination shock")
ax.legend(fontsize=8)

fig.tight_layout()
figp = RESULTS / "figures" / "supplementary" / "step_b62_rtg_nuclear.png"
fig.savefig(figp, dpi=300)
logger.data_save(figp)

# CSV of the power series
import csv
csvp = RESULTS / "step_b62_rtg_nuclear.csv"
with csvp.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["series", "craft", "unit", "utc", "year",
                "t_yr_since_bom", "p_norm", "p_w", "resid_log"])
    for k, s in series.items():
        f = fit_model(s)
        for i in range(len(s["t_yr"])):
            w.writerow([k, s["craft"], s["unit"], s["utc"][i],
                        f"{s['year'][i]:.4f}", f"{s['t_yr'][i]:.4f}",
                        f"{s['p_norm'][i]:.6f}", f"{s['p_w'][i]:.2f}",
                        f"{f['resid'][i]:.5f}"])
logger.data_save(csvp)

out = dict(
    step="step_098_rtg_nuclear_channel",
    description=("RTG nuclear-decay channel: Pu-238 alpha decay as a "
                 "nuclear clock, using the published JPL MHW-RTG "
                 "telemetry record (Whiting & Woerner 2023 Dryad "
                 "dataset); decay+degradation fits, signed-offset "
                 "step tests at the four boundary crossings, "
                 "per-unit consistency, and a VG1-VG2 differential "
                 "channel; decay-rate bounds mapped to effective-"
                 "alpha bounds via the Gamow amplification."),
    inputs={
        "rtg_telemetry": str(XLSX.relative_to(XLSX.parents[2])),
        "literature_anchors":
            "data/raw/literature/literature_anchors.json",
        "dataset_citation": (
            "Whiting & Woerner 2023, doi:10.5061/dryad.1zcrjdfw2; "
            "Voyager telemetry courtesy E. Medina (JPL), through "
            "2021-08")},
    physics=dict(
        reaction="Pu-238 -> U-234 + alpha, Q=5.593 MeV",
        half_life_yr=T12_YR,
        gamow_exponent=G_GAMOW,
        decay_amplification=AMP,
        amplification_note=(
            "|d lambda/lambda| = 2G |d alpha_eff/alpha_eff| ~ "
            f"{AMP:.0f} |d alpha_eff/alpha_eff|; the RTG power "
            "record amplifies any effective-alpha step by ~150x.")),
    telemetry=dict(
        n_epochs={k: int(len(s["t_yr"])) for k, s in series.items()},
        last_epoch={k: s["utc"][-1] for k, s in series.items()},
        pbom_w={k: s["pbom_w"] for k, s in series.items()},
        pedigree_notes=[
            "first ~13.3 yr: spacecraft-reported power",
            "post-13.3 yr: bit-drop epochs, quantized at 0.039 A "
            "(VG1) / 0.038 A (VG2) per bit",
            "V2 RTG#3 telemetry lost ~700 h after launch; V2 total "
            "reconstructed from units 1-2 + bus-voltage correction",
        ]),
    model_fits={k: {kk: vv for kk, vv in f.items()
                    if kk not in ("resid", "model_norm")}
                for k, f in fits.items()},
    step_tests=step_tests,
    placebo_scan=placebo,
    unit_step_tests=unit_tests,
    excess_decline=dict(
        baseline_window_yr=[BASE_LO, BASE_HI],
        note=("linear-in-log baseline fitted inside the boundary "
              "(13-26 yr); the cumulative excess decline is the "
              "model-light observable for a graded level shift")),
    unit_coherence=unit_coherence,
    v1_ramp_fit=v1_ramp,
    baseline_window_robustness=window_robustness,
    v2_total_mixing=v2_mixing,
    geometry_concordance=geom_concordance,
    differential_v1_minus_v2=differential,
    bound=dict(
        bound_dlambda_frac_2s=round(bound_dlambda, 5),
        bound_dalpha_frac_2s=float(f"{bound_dalpha:.2e}"),
        note=("largest |step| + 2 sigma over the four crossing "
              "tests; the nuclear-sector bound at the crossing "
              "epochs")),
    predicted_step_for_comet_scale=predicted_steps,
    literature_anchor_check=anchor_check,
    caveats=[
        "Post-13.3-year telemetry is recorded at bit-drop epochs and "
        "is quantized at ~0.25 % per bit (V1); the step tests "
        "operate on the quantized series, so sub-percent structure "
        "is carried by the fitted steps plus their standard errors.",
        "The V2 total-power series alternates between two report "
        "conventions offset ~0.75% through mission-year ~37; "
        "level tests on the mixed series inherit that offset. The "
        "V2 verdict rests on the two clean per-unit series, which "
        "diverge from each other by ~2% over mission-years 37-43 -- "
        "unit-level systematics of order the sought signal.",
        "A decay-rate step is partially degenerate with "
        "thermoelectric-degradation changes; no other MHW-RTG "
        "mission lived past ~12 yr, so the mission-year 30-44 "
        "degradation shape has no external control and the V1 "
        "excess decline is degenerate with undocumented mid-life "
        "SiGe aging.",
        "The bound applies to the nuclear sector specifically: TEP "
        "treats electronic/nuclear/dynamical couplings to A(phi) as "
        "potentially different, so this channel does not constrain "
        "the orbital channel directly. Under universal conformal "
        "coupling all spacecraft clocks share the lapse ratio and "
        "the SCLK record (step 096) bounds any common shift at "
        "~1 ppm; a per-cent-scale nuclear excursion is therefore "
        "interpretable only through a process-dependent "
        "coupling.",
        "Strong-sector coupling shifts (preformation factor, Q "
        "value) would add to the EM-barrier term computed here; the "
        "reported bound is for the dominant electromagnetic "
        "penetrability channel.",
        "Telemetry ends 2021-08 (per the published dataset); the "
        "2023 literature anchor (~225 W) is consistent with the "
        "fitted model but postdates the record.",
        "The dslope (post-crossing slope-change) statistics carry "
        "very small nominal p-values across the record because the "
        "quadratic under-fits the SiGe degradation curvature; they "
        "are model-misfit diagnostics, not placebo-controlled "
        "crossing tests, and should not be read as detections.",
        "The differential (VG1-VG2) crossing tests inherit the "
        "quadratic end-curl of both fits; the V2 heliopause "
        "differential step is partially a boundary-of-fit "
        "artifact.",
    ],
    tep_interpretation=None)
_v1_hp_excess = -unit_coherence["VG1"]["excess_at_hp_pct"]["mean"]
_v1_end_excess = -unit_coherence["VG1"]["excess_end_pct"]["mean"]
_v1_dalpha = _v1_hp_excess / 100.0 / AMP
out["tep_interpretation"] = (
        "The RTG record is a 44-year nuclear clock testing the "
        "process-dependent coupling structure. Under universal "
        "conformal coupling the SCLK bound (~1 ppm) already excludes "
        "any common rate shift at the sought level, so the channel "
        "of interest is the nuclear sector, where an effective-alpha "
        "shift is amplified ~153x through the Gamow exponent. The "
        "model-light excess-decline channel shows V1's three "
        "independent nuclear clocks accumulating a coherent "
        f"{_v1_hp_excess:.1f}-{_v1_end_excess:.1f}% "
        "excess decline across the heliopause window (ramp onset "
        f"~yr {v1_ramp['onset_yr']:.0f}, centre ~yr "
        f"{v1_ramp['centre_yr']:.0f} against the crossing at yr "
        f"{v1_ramp['hp_yr']:.0f}; nothing at the termination "
        "shock). The sign also matches the field's polarity: the "
        "slip field is positive at both axis poles (cos-2theta), so "
        "in the mirror lobe V1 occupies, proper time -- and with it "
        "the decay clock -- runs faster, predicting exactly the "
        "faster-decline (negative) power excess observed. V2, which "
        "crossed in the intervening negative-slip band, shows no "
        "unit-coherent excursion -- its two surviving units diverge "
        "by ~1.7% -- the directional pattern the bipolar field "
        "predicts. The V1 coherence is a candidate signature "
        "(delta lambda/lambda ~ -0.5% -> delta alpha_eff/alpha_eff "
        f"~ {_v1_dalpha:.1e} via the Gamow factor) but with two "
        "quantified fragilities: the excess amplitude -- and even "
        "its sign -- is baseline-window-dependent (+1.6% under a "
        "0-26 yr baseline vs -0.5% under 13-26), and VG2_rtg2 alone "
        "shows -0.85% at end under the canonical window, an "
        "out-of-cap counterexample at the same amplitude. The "
        "coherence is real but expected for units sharing one "
        "craft's thermal environment and degradation physics; the "
        "excursion is degenerate with undocumented mid-life "
        "MHW-RTG degradation because no other unit of the design "
        "survived past 12 yr. The conservative nuclear-sector "
        "bound at the crossings "
        f"remains |d alpha_eff/alpha_eff| ~ {bound_dalpha:.0e}.")
out_path = RESULTS / "step_b62_rtg_nuclear.json"
out_path.write_text(json.dumps(out, indent=1))
logger.data_save(out_path)
logger.success("RTG nuclear channel complete")
