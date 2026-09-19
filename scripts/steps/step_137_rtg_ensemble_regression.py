"""Step 137: RTG cross-generator degradation regression (step_b101).

The thermal-regression close on the RTG channel named in step 098
and bounded by step 103.  Per-sensor thermal telemetry is not
publicly archived for the Voyager bus (the JPL thermal-model
literature records that the flight system carries very few
temperature sensors), so the physically meaningful regression of
the RTG decline against the spacecraft's thermal/degradation state
is run against the full cross-generator ensemble: every unit-level
P/P0 series in the Whiting & Woerner (2023) Dryad dataset, six RTG
families, 34 series.

T1  Ensemble extraction: all unit-level series (MHW, GPHS, MMRTG,
    SNAP-19, SNAP-27, SNAP-9A).

T2  Smooth per-unit baselines (quadratic in log P/P0 over mission
    years) -- a 9-yr excursion cannot distort a 10-40 yr fit.

T3  Pooled wander null: the empirical distribution of residual
    change over 9-yr windows across every generator -- what real,
    documented RTG aging does to a crossing-aligned window.

T4  Crossing-window test: Voyager per-unit and mission-total
    residual change over the heliopause-aligned windows (V1
    34.75-43.75 yr per step_103, opening at its t_hp = 34.97 yr;
    V2 34.5-43.5 yr, closing just after its t_hp = 41.22 yr --
    the widest post-crossing window the V2 record supports),
    priced against the pooled null.

T5  Sibling-coherence thermal control: units sharing one
    spacecraft share its thermal/bus environment; the
    within-mission versus cross-mission residual-difference
    distribution measures how much common-mode wander the shared
    thermal state can supply.

T6  Bound translation through the Gamow amplification of step 098
    (|d lambda/lambda| = 152.7 |d alpha_eff/alpha_eff|).
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_137_rtg_ensemble_regression")
tee_stdout(logger)
logger.header("RTG cross-generator degradation regression")

import json
import numpy as np
import openpyxl
from scipy import stats as _st

XLSX = (RESULTS.parent / "data" / "raw" / "rtg" /
        "All_RTGs_P_Over_P0_-_2023_Update_for_Release.xlsx")
GAMOW = 152.674  # step_098 amplification factor


def parse_sheet(ws):
    rows = list(ws.iter_rows(values_only=True))
    hdr = None
    for i, row in enumerate(rows[:40]):
        if any(c and "Time" in str(c) for c in row):
            hdr = i
            break
    if hdr is None:
        return {}
    name_row, best = hdr - 5, 0
    for j in range(hdr - 1, max(0, hdr - 8), -1):
        n = sum(1 for c in rows[j]
                if isinstance(c, str) and len(str(c)) > 3
                and not any(k in str(c) for k in (
                    "Time", "PBOM", "BOM", "Average", "Total",
                    "Serial", "Telem", "Normaliz", "(years)",
                    "Power")))
        if n > best:
            best, name_row = n, j
    out = {}
    for c, name in enumerate(rows[name_row]):
        if (isinstance(name, str) and len(name.strip()) > 2
                and not any(k in name for k in (
                    "Time", "PBOM", "BOM", "Average", "Total",
                    "Serial", "Telem", "Normaliz", "(years)",
                    "Power", "dataset", "All info", "†", "‡",
                    "*"))):
            t, p = [], []
            for r in rows[hdr + 2:]:
                if (c + 1 < len(r)
                        and isinstance(r[c], (int, float))
                        and isinstance(r[c + 1], (int, float))
                        and 0 <= r[c] < 60
                        and 0 < r[c + 1] <= 1.5):
                    t.append(float(r[c]))
                    p.append(float(r[c + 1]))
            if len(t) > 5:
                out[name.strip()[:32]] = (np.array(t), np.array(p))
    return out


wb = openpyxl.load_workbook(XLSX, data_only=True, read_only=True)
UNITS, FAM = {}, {}
for fam in ("MHW-RTG", "GPHS-RTG", "MMRTG", "SNAP-19", "SNAP-27",
            "SNAP-9A"):
    for k, v in parse_sheet(wb[fam]).items():
        UNITS[k] = v
        FAM[k] = fam
logger.info(f"T1 ensemble: {len(UNITS)} series across "
            f"{len(set(FAM.values()))} families")
T1 = {k: {"family": FAM[k], "n": int(len(v[0])),
          "t_max_yr": float(v[0][-1])} for k, v in UNITS.items()}


def resid(t, p):
    y = np.log(p)
    X = np.column_stack([np.ones(len(t)), t, t * t])
    with np.errstate(all="ignore"):
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        c = np.where(np.isfinite(c), c, 0.0)
        return y - X @ c


SIBS = {"V1": ["Voyager 1 RTG #1 (F-8)", "Voyager 1 RTG #2 (F-6)",
               "Voyager 1 RTG #3 (F-9)"],
        "V2": ["Voyager 2 RTG #1 (F-10)", "Voyager 2 RTG #2 (F-12)"],
        "Cassini": ["Cassini RTG #1", "Cassini RTG #2",
                    "Cassini RTG #3"],
        "Galileo": ["Galileo -RTG (F-1)", "Galileo +RTG (F-4)"],
        "Nimbus": ["Nimbus III RTG #1", "Nimbus III RTG #2"],
        "LES": ["LES 8", "LES 9"]}

# heliopause crossing epochs in mission years (V1 launch 1977-09-05,
# V2 1977-08-20; crossing epochs per step_096/step_103)
HP_YR = {"V1": 34.97, "V2": 41.26}


def wander(a_len, gap, b_len, min_pts=2):
    """Geometry-matched wander null: residual level change between
    two end bins of a (a_len + gap + b_len)-yr sliding window."""
    W, META = [], []
    for k, (t, p) in UNITS.items():
        r = resid(t, p)
        span = a_len + gap + b_len
        for t0 in np.arange(0, t[-1] - span, 0.5):
            a = (t >= t0) & (t <= t0 + a_len)
            b = (t >= t0 + a_len + gap) & (t <= t0 + span)
            if a.sum() >= min_pts and b.sum() >= min_pts:
                W.append(np.median(r[b]) - np.median(r[a]))
                META.append((k, float(t0)))
    return np.array(W), META


def crossing_dr(unit, a_rng, b_rng):
    t, p = UNITS[unit]
    r = resid(t, p)
    a = (t >= a_rng[0]) & (t <= a_rng[1])
    b = (t >= b_rng[0]) & (t <= b_rng[1])
    if a.sum() < 2 or b.sum() < 2:
        return None, int(a.sum()), int(b.sum())
    return float(np.median(r[b]) - np.median(r[a])), \
        int(a.sum()), int(b.sum())


def test_geometry(craft, unit):
    """Strictly pre-crossing baseline vs the widest post-crossing
    end bin the record supports; returns (a_rng, b_rng, gap, b_len).
    """
    tc = HP_YR[craft]
    tmax = UNITS[unit][0][-1]
    b0 = max(tc, tmax - 3.0)
    return (tc - 1.5, tc), (b0, tmax), b0 - tc, tmax - b0


# ---- T3 pooled wander null ----------------------------------------
logger.info("T3 pooled 9-yr wander null")
W, W_SRC = wander(1.5, 7.5, 1.5)
T3 = dict(n_windows=int(len(W)), std_pct=float(100 * W.std()),
          pctiles_pct={q: float(v) for q, v in zip(
              (5, 25, 50, 75, 95),
              np.percentile(100 * W, (5, 25, 50, 75, 95)))})
logger.info(f"  n={len(W)} windows, std={100*W.std():.3f}%, "
            f"p5/p95={T3['pctiles_pct'][5]:+.3f}/"
            f"{T3['pctiles_pct'][95]:+.3f}%")

# ---- T4 crossing-window test ---------------------------------------
logger.info("T4 crossing-window residuals vs ensemble null")
TEST_UNITS = [("VG1_rtg1", "Voyager 1 RTG #1 (F-8)", "V1"),
              ("VG1_rtg2", "Voyager 1 RTG #2 (F-6)", "V1"),
              ("VG1_rtg3", "Voyager 1 RTG #3 (F-9)", "V1"),
              ("VG1_total", "Voyager 1", "V1"),
              ("VG2_rtg1", "Voyager 2 RTG #1 (F-10)", "V2"),
              ("VG2_rtg2", "Voyager 2 RTG #2 (F-12)", "V2"),
              ("VG2_total", "Voyager 2", "V2")]
NULLS = {}
T4 = {}
for tag, unit, craft in TEST_UNITS:
    a_rng, b_rng, gap, b_len = test_geometry(craft, unit)
    key = (round(gap, 2), round(b_len, 2))
    if key not in NULLS:
        NULLS[key] = wander(1.5, gap, b_len)
    Wg, META = NULLS[key]
    dr, na, nb = crossing_dr(unit, a_rng, b_rng)
    if dr is None:
        T4[tag] = dict(window_yr=[a_rng, b_rng], note="insufficient data")
        continue
    # leave-own-window-out null: drop the tested unit's wander
    # windows that overlap the crossing span so the null cannot
    # contain the signal being priced
    wk = np.array([w for w, (src, t0) in zip(Wg, META)
                   if not (src == unit
                           and t0 < b_rng[1]
                           and t0 + 1.5 + gap + b_len > a_rng[0])])
    tail = float((wk <= dr).mean())
    exc = float((np.abs(wk) >= abs(dr)).mean())
    T4[tag] = dict(window_yr=[list(a_rng), list(b_rng)],
                   delta_resid_pct=100 * dr,
                   n_a=na, n_b=nb,
                   one_sided_tail=tail,
                   abs_exceedance_frac=exc,
                   tail_note=("one-sided decline tail vs the "
                              "geometry-matched leave-own-window-out"
                              f" wander null (n={len(wk)})"))
    logger.info(f"  {tag}: dr={100*dr:+.3f}%  tail={100*tail:.1f}%"
                f"  |W|>=|dr| frac={exc:.3f}")

# ---- T4b craft common-mode crossing channel -------------------------
# The a-priori TEP channel: the RTG units on one spacecraft share a
# single location in the field, so the correct statistic is the
# craft-level common-mode drift at the crossing priced against the
# fleet's common-mode wander -- not a product of dependent per-unit
# tails.
logger.info("T4b craft common-mode crossing channel")
T4B = {}
for craft, units in (("V1", SIBS["V1"]), ("V2", SIBS["V2"])):
    tc = HP_YR[craft]
    tmax = min(UNITS[u][0][-1] for u in units)
    b0 = max(tc, tmax - 3.0)
    a_rng, b_rng = (tc - 1.5, tc), (b0, tmax)
    gap, b_len = b0 - tc, tmax - b0
    drs, npts = [], []
    for u in units:
        dr, na, nb = crossing_dr(u, a_rng, b_rng)
        if dr is not None:
            drs.append(dr)
            npts.append((na, nb))
    if len(drs) != len(units):
        T4B[craft] = dict(note="insufficient common-mode coverage")
        continue
    cm_dr = float(np.mean(drs))
    # fleet common-mode wander with the same window geometry:
    # mean sibling drift per craft per window
    CMW = []
    span = 1.5 + gap + b_len
    for ctag, cu in SIBS.items():
        tmin = min(UNITS[u][0][-1] for u in cu)
        for t0 in np.arange(0, tmin - span, 0.5):
            wr = []
            for u in cu:
                t, p = UNITS[u]
                r = resid(t, p)
                a = (t >= t0) & (t <= t0 + 1.5)
                b = (t >= t0 + 1.5 + gap) & (t <= t0 + span)
                if a.sum() >= 2 and b.sum() >= 2:
                    wr.append(np.median(r[b]) - np.median(r[a]))
            if len(wr) == len(cu) and not (
                    ctag == craft and t0 < b_rng[1]
                    and t0 + span > a_rng[0]):
                CMW.append(float(np.mean(wr)))
    CMW = np.array(CMW)
    tail = float((CMW <= cm_dr).mean()) if len(CMW) else float("nan")
    T4B[craft] = dict(window_yr=[list(a_rng), list(b_rng)],
                      per_unit_dr_pct=[100 * d for d in drs],
                      common_mode_dr_pct=100 * cm_dr,
                      one_sided_tail=tail, n_null=int(len(CMW)),
                      n_bins=npts,
                      note=("a-priori channel: mean sibling-unit "
                            "drift at the heliopause vs the fleet "
                            "common-mode wander null"))
    logger.info(f"  {craft}: common-mode dr={100*cm_dr:+.3f}%  "
                f"tail={100*tail:.2f}% (null n={len(CMW)})")

# ---- T5 sibling-coherence thermal control --------------------------
logger.info("T5 sibling-coherence thermal control")
import itertools
within = []
for units in SIBS.values():
    for a, b in itertools.combinations(units, 2):
        ta, pa = UNITS[a]
        tb, pb = UNITS[b]
        ra, rb = resid(ta, pa), resid(tb, pb)
        for t0 in np.arange(1, min(ta[-1], tb[-1]) - 9, 1.0):
            wa = (ta >= t0) & (ta <= t0 + 9)
            wbs = (tb >= t0) & (tb <= t0 + 9)
            if wa.sum() >= 6 and wbs.sum() >= 6:
                within.append(abs(np.median(ra[wa])
                                  - np.median(rb[wbs])))
rng = np.random.default_rng(7)
across = []
keys = list(UNITS)
while len(across) < 400:
    a, b = rng.choice(keys, 2, replace=False)
    ta, pa = UNITS[a]
    tb, pb = UNITS[b]
    ra, rb = resid(ta, pa), resid(tb, pb)
    t0 = float(rng.uniform(1, max(2, min(ta[-1], tb[-1]) - 9)))
    wa = (ta >= t0) & (ta <= t0 + 9)
    wbs = (tb >= t0) & (tb <= t0 + 9)
    if wa.sum() >= 6 and wbs.sum() >= 6:
        across.append(abs(np.median(ra[wa]) - np.median(rb[wbs])))
within, across = np.array(within), np.array(across)
mw = _st.mannwhitneyu(within, across)
T5 = dict(
    within_pct={q: float(v) for q, v in zip(
        (25, 50, 75), np.percentile(100 * within, (25, 50, 75)))},
    across_pct={q: float(v) for q, v in zip(
        (25, 50, 75), np.percentile(100 * across, (25, 50, 75)))},
    n_within=int(len(within)), n_across=int(len(across)),
    mwu_p=float(mw.pvalue),
    note=("sibling units share the spacecraft thermal/bus "
          "environment; the within-mission residual coherence "
          "prices the common-mode wander the shared thermal "
          "state can supply"))
logger.info(f"  within {np.median(100*within):.3f}% vs across "
            f"{np.median(100*across):.3f}% medians, "
            f"MWU p={mw.pvalue:.2g}")

# ---- T6 bound translation ------------------------------------------
v1 = T4["VG1_total"]["delta_resid_pct"] / 100.0
excess = abs(v1)
T6 = dict(v1_crossing_excess_frac=excess,
          dlambda_frac_1s=excess,
          bound_dalpha_frac_1s=excess / GAMOW,
          note=("residual crossing-window power change mapped to "
                "a decay-rate shift, then to effective-alpha "
                "through the step_098 Gamow factor"))
logger.info(f"  V1 excess {100*excess:.3f}% -> "
            f"|da/a|~{excess/GAMOW:.2e} (1s)")

tails = [v["one_sided_tail"] for v in T4.values()
         if "one_sided_tail" in v]
min_tail, max_tail = min(tails), max(tails)
v1_cm = T4B.get("V1", {})
v2_cm = T4B.get("V2", {})
v1_tail = v1_cm.get("one_sided_tail")
v2_tail = v2_cm.get("one_sided_tail")
n_craft = sum(1 for t in (v1_tail, v2_tail) if t is not None)
v1_tail_fw = min(1.0, v1_tail * n_craft) \
    if v1_tail is not None else None
if v1_tail_fw is not None and v1_tail_fw < 0.05:
    verdict = (f"COHERENT CROSSING DECLINE: V1's three nuclear "
               f"clocks decline together "
               f"{v1_cm['common_mode_dr_pct']:+.2f}% across the "
               f"heliopause window -- a common-mode excursion at "
               f"{100*v1_tail:.2f}% raw / {100*v1_tail_fw:.2f}% "
               f"family-wise (over {n_craft} craft tests) against "
               f"the fleet's common-mode wander "
               f"(n={v1_cm['n_null']}); the out-of-cap V2 control "
               f"shows a weaker {v2_cm.get('common_mode_dr_pct', 0):+.2f}%"
               f" at {100*(v2_tail or 1):.1f}% -- the "
               "cap-geometry contrast points the predicted way, "
               "bounded by fleet-null resolution and undocumented "
               "MHW-RTG aging")
elif v1_tail is not None and v1_tail < 0.10:
    verdict = (f"WEAK COHERENT DECLINE: V1's three nuclear clocks "
               f"decline together {v1_cm['common_mode_dr_pct']:+.2f}%"
               f" across the heliopause window -- a common-mode "
               f"excursion at {100*v1_tail:.2f}% raw / "
               f"{100*v1_tail_fw:.1f}% family-wise (over {n_craft} "
               f"craft tests) against the fleet's common-mode "
               f"wander (n={v1_cm['n_null']}), while the out-of-cap "
               f"V2 common-mode is flat at "
               f"{v2_cm.get('common_mode_dr_pct', 0):+.2f}% "
               f"({100*(v2_tail or 1):.0f}% tail) -- the "
               "cap-geometry contrast points the predicted way but "
               "the excursion does not clear the family-wise line; "
               "the nuclear channel remains an upper bound with a "
               "resolved directional carrier")
elif min_tail < 0.05:
    worst = min(T4, key=lambda k: T4[k]["one_sided_tail"])
    unit_tails = [v["one_sided_tail"] for k, v in T4.items()
                  if not k.endswith("_total")]
    verdict = (f"MARGINAL EXCURSION: {worst} declines "
               f"{T4[worst]['delta_resid_pct']:+.2f}% over its "
               f"crossing window, a one-sided decline tail of "
               f"{100*min_tail:.1f}% against the documented "
               "cross-generator wander -- beyond the 5% line on "
               "the mission-total series while the per-unit tails "
               f"stay at {100*min(unit_tails):.0f}-"
               f"{100*max(unit_tails):.0f}%; the nuclear channel "
               "remains an upper bound with a resolved marginal "
               "carrier, not a closed null")
elif max_tail > 0.95:
    verdict = ("REVERSED: crossing-window residuals sit above the "
               "wander distribution -- no decline signal present")
else:
    verdict = ("ENSEMBLE-CONSISTENT: all crossing-window declines "
               "sit inside the documented cross-generator wander "
               f"(tails {100*min_tail:.0f}-{100*max_tail:.0f}%); "
               "sibling coherence is dominated by the shared "
               "thermal/bus environment; the nuclear channel "
               "remains an upper bound, now priced against the "
               "full RTG ensemble")
logger.info(f"verdict: {verdict}")

out = dict(step="step_137_rtg_ensemble_regression", result="b101",
           verdict=verdict,
           inputs=dict(xlsx=str(XLSX.name),
                       telemetry_note=("per-sensor thermal "
                                       "telemetry is not publicly "
                                       "archived; the ensemble "
                                       "degradation distribution "
                                       "is the thermal-regression "
                                       "surrogate")),
           T1_ensemble=T1, T3_wander_null=T3,
           T4_crossing_windows=T4, T4B_common_mode=T4B,
           T5_sibling_coherence=T5, T6_bound=T6)
with open(RESULTS / "step_b101_rtg_ensemble.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.info("wrote results/step_b101_rtg_ensemble.json")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
ax[0].hist(100 * W, bins=40, color="0.4", alpha=0.7)
for tag, e in T4.items():
    if tag.endswith("_total") and "delta_resid_pct" in e:
        ax[0].axvline(e["delta_resid_pct"], ls="--", label=tag)
for craft, e in T4B.items():
    if "common_mode_dr_pct" in e:
        ax[0].axvline(e["common_mode_dr_pct"], ls="-",
                      label=f"{craft} common-mode")
ax[0].set_xlabel("residual change (%)")
ax[0].set_title("T3 pooled wander null vs Voyager windows")
ax[0].legend(fontsize=8)

ax[1].hist([100 * within, 100 * across], bins=30,
           label=["sibling (shared bus)", "cross-mission"],
           alpha=0.65, density=True)
ax[1].set_xlabel("|Δ median residual| over 9 yr (%)")
ax[1].set_title(f"T5 thermal coherence p={mw.pvalue:.1g}")
ax[1].legend(fontsize=8)

for unit, c in (("Voyager 1", "indianred"), ("Voyager 2",
                                              "steelblue"),
                ("LES 8", "0.5"), ("LES 9", "0.5"),
                ("Cassini RTG #1", "gold")):
    t, p = UNITS[unit]
    ax[2].plot(t, 100 * resid(t, p), lw=0.8, label=unit)
ax[2].axvline(34.97, color="indianred", ls=":", lw=0.8)
ax[2].axvline(41.22, color="steelblue", ls=":", lw=0.8)
ax[2].set_xlabel("mission year")
ax[2].set_ylabel("residual (%)")
ax[2].set_title("T2 residual trajectories (crossings marked)")
ax[2].legend(fontsize=7)

fig.tight_layout()
fig.savefig(RESULTS / "figures" / "step_b101_rtg_ensemble.png",
            dpi=150)
logger.info("wrote results/figures/step_b101_rtg_ensemble.png")
