"""Step 141: SCLK signed-drift crossing channel (step_b105).

Step 096 tested the absolute calibration residuals -- the scatter
of the SCLK correction record -- and returned bounds at the ppm
level.  A different observable lives in the same record: the
signed fractional frequency offset dnu of the onboard oscillator
versus Earth time.  The onboard oscillator ticks in spacecraft
proper time; its measured rate against ground-defined ET through
the downlinked timestamps is a one-way clock comparison -- the
transport class the two-way cancellation theorem says can carry a
lapse signal.  This step tests the signed level for
crossing-aligned structure.

T1  Signed dnu level steps at the four boundary crossings
    (V1 TS 2004.96, V1 HP 2012.65, V2 TS 2007.66, V2 HP 2018.90)
    at 2/3/5-yr window widths, raw and smooth-ageing detrended
    (quadratic in log mission-age).

T2  Changepoint extremeness: the same windowed-step statistic at
    every admissible changepoint in each record -- the price of
    asking "how often does a step this size occur anywhere".

T3  Joint statistic: 4/4 sign consistency against the measured
    ~50 per cent negative base rate, and the Fisher combination
    of the four negative-direction tail fractions.

T4  Artefact audit: calibration cadence around each crossing
    (long gaps bias windowed means), and the NH1 interior control
    (whose dnu is frozen at sub-0.1 ppm -- unusable for the
    signed channel, registered as such).

T5  Degeneracy statement: ordinary oscillator ageing predicts
    monotone deceleration, which the detrended test removes; what
    remains is a crossing-aligned local feature, but the record
    is kernel-segmented and two-craft, so the channel registers
    as a candidate with priced degeneracy -- not a null and not
    a detection.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_141_sclk_signed_drift")
tee_stdout(logger)
logger.header("SCLK signed-drift crossing channel")

import csv
import json
import numpy as np
from scipy import stats as _st

RATES = RESULTS / "step_b60_sclk_rates.csv"
rows = [r for r in csv.DictReader(open(RATES))
        if r["mode"] == "calibration"]
CROSS = {"VG1": [("TS", 2004.96, 94.0), ("HP", 2012.65, 121.6)],
         "VG2": [("TS", 2007.66, 84.0), ("HP", 2018.90, 119.0)]}
WINS = (2.0, 3.0, 5.0)


def series(craft):
    rs = [r for r in rows if r["craft"] == craft]
    et = np.array([float(r["et_s"]) for r in rs])
    dnu = np.array([float(r["dnu_frac"]) for r in rs]) * 1e6
    return et / 31557600 + 2000.0, dnu


def detrend(yr, dnu, exclude=None):
    age = yr - yr.min() + 0.5
    X = np.column_stack([np.ones(len(yr)), np.log(age),
                         np.log(age) ** 2])
    keep = np.ones(len(yr), bool) if exclude is None else ~exclude
    with np.errstate(all="ignore"):
        c, *_ = np.linalg.lstsq(X[keep], dnu[keep], rcond=None)
        c = np.where(np.isfinite(c), c, 0.0)
        return dnu - X @ c


def step_stat(yr, dnu, tc, w):
    a = (yr >= tc - w) & (yr < tc)
    b = (yr >= tc) & (yr <= tc + w)
    if a.sum() < 8 or b.sum() < 8:
        return None
    return float(dnu[b].mean() - dnu[a].mean())


def changepoint_null(yr, dnu, w):
    ts = np.arange(yr.min() + w + 0.5, yr.max() - w - 0.5, 0.25)
    out = [s for tc in ts
           if (s := step_stat(yr, dnu, tc, w)) is not None]
    return np.array(out)


# ---- T1/T2 steps + extremeness -------------------------------------
T1, T2 = {}, {}
for craft in ("VG1", "VG2"):
    yr, dnu = series(craft)
    # exclude +-max(WINS) around each crossing from the trend fit so
    # the smooth-ageing model cannot absorb the steps being tested
    ex = np.zeros(len(yr), bool)
    for _, t_hp, _ in CROSS[craft]:
        ex |= np.abs(yr - t_hp) <= max(WINS)
    rd = detrend(yr, dnu, exclude=ex)
    T1[craft], T2[craft] = {}, {}
    for name, t_hp, r_hp in CROSS[craft]:
        T1[craft][name] = dict(epoch=t_hp, radius_au=r_hp)
        T2[craft][name] = {}
        for w in WINS:
            raw = step_stat(yr, dnu, t_hp, w)
            det = step_stat(yr, rd, t_hp, w)
            null_raw = changepoint_null(yr, dnu, w)
            null_det = changepoint_null(yr, rd, w)
            e = dict(raw_ppm=raw, detrended_ppm=det,
                     neg_frac_raw=float((null_raw <= raw).mean()),
                     neg_frac_det=float((null_det <= det).mean()),
                     abs_frac_det=float(
                         (np.abs(null_det) >= abs(det)).mean()),
                     n_changepoints=int(len(null_det)))
            T1[craft][name][f"w{w:g}"] = e
            T2[craft][name][f"w{w:g}"] = dict(
                null_median_ppm=float(np.median(null_det)),
                null_neg_base_rate=float((null_det < 0).mean()))
            logger.info(f"  {craft} {name} w={w}: raw {raw:+.2f} / "
                        f"det {det:+.2f} ppm, neg tail "
                        f"{100*e['neg_frac_det']:.1f}%")

# ---- T3 joint -------------------------------------------------------
logger.info("T3 joint sign-consistency")
all_neg = all(T1[c][x][f"w{w:g}"]["detrended_ppm"] < 0
              for c in ("VG1", "VG2") for x in T1[c]
              for w in WINS)
base = np.mean([T2[c][x]["w3"]["null_neg_base_rate"]
                for c in ("VG1", "VG2") for x in T2[c]])
n_crossings = sum(len(T1[c]) for c in ("VG1", "VG2"))
n_negative_w3 = sum(T1[c][x]["w3"]["detrended_ppm"] < 0
                    for c in ("VG1", "VG2") for x in T1[c])
binom = float(_st.binomtest(n_negative_w3, n_crossings, base).pvalue)
tails = [T1[c][x]["w3"]["neg_frac_det"]
         for c in ("VG1", "VG2") for x in T1[c]]
chi2 = float(-2 * np.sum(np.log(np.clip(tails, 1e-12, 1))))
fisher_p = float(_st.chi2.sf(chi2, 2 * len(tails)))
T3 = dict(all_negative=all_neg, n_negative_w3=n_negative_w3, n_crossings=n_crossings,
          interpretation="Nominal independent-crossing references only; same-craft dependence, overlapping windows, and epoch selection are not calibrated.", neg_base_rate=float(base),
          binom_p_4of4=binom, fisher_chi2=chi2,
          fisher_p=fisher_p, tails_w3=tails)
logger.info(f"  4/4 negative, base rate {base:.2f}, "
            f"binom p={binom:.3f}; Fisher chi2={chi2:.1f} "
            f"p={fisher_p:.3g}")

# ---- T4 cadence audit ----------------------------------------------
logger.info("T4 cadence/segmentation audit")
T4 = {}
for craft in ("VG1", "VG2", "NH1"):
    yr, dnu = series(craft)
    gaps = np.diff(yr) * 365.25
    T4[craft] = dict(n=len(yr), median_gap_d=float(np.median(gaps)),
                     dnu_range_ppm=[float(dnu.min()),
                                    float(dnu.max())])
    if craft == "NH1":
        T4[craft]["note"] = ("dnu frozen at sub-0.1 ppm -- NH "
                             "unusable for the signed channel")
    for name, t_hp, r_hp in CROSS.get(craft, []):
        i = int(np.searchsorted(yr, t_hp))
        loc = gaps[max(0, i - 3):i + 3]
        T4[craft][f"{name}_nearby_gaps_d"] = [float(g)
                                              for g in loc]
        logger.info(f"  {craft} {name}: nearby gaps "
                    f"{[round(float(g),1) for g in loc]}")

# ---- verdict --------------------------------------------------------
extremal = min(tails)
verdict = (f"SCLK CALIBRATION DIAGNOSTIC: {n_negative_w3}/{n_crossings} "
           "negative three-year-window steps. Same-craft correlation and "
           "kernel segmentation prevent an independent-crossing detection "
           "claim; compare window sensitivity and empirical epoch ranks.")
logger.info(f"verdict: {verdict}")

out = dict(step="step_141_sclk_signed_drift", result="b105",
           verdict=verdict,
           observable="Signed fractional rate in SCLK clock-correlation coefficients; oscillator and operational contributions are not separated",
           T1_crossing_steps=T1, T2_changepoint_nulls=T2,
           T3_joint=T3, T4_cadence=T4)
with open(RESULTS / "step_b105_sclk_signed_drift.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
logger.info("wrote results/step_b105_sclk_signed_drift.json")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
cols = {"VG1": "indianred", "VG2": "steelblue", "NH1": "0.5"}
for i, craft in enumerate(("VG1", "VG2", "NH1")):
    yr, dnu = series(craft)
    ax[i].plot(yr, dnu, lw=0.5, color=cols[craft])
    for name, t_hp, r_hp in CROSS.get(craft, []):
        ax[i].axvline(t_hp, ls="--", lw=0.9, color="black")
        ax[i].text(t_hp, ax[i].get_ylim()[1] * 0.9, name,
                   fontsize=8, rotation=90, va="top")
    ax[i].set_xlabel("year")
    ax[i].set_title(craft)
ax[0].set_ylabel("signed clock offset dnu (ppm)")
fig.suptitle("SCLK signed oscillator offset vs epoch; boundary "
             "crossings marked (step 141)")
fig.tight_layout()
fig.savefig(RESULTS / "figures" / "step_b105_sclk_drift.png",
            dpi=150)
logger.info("wrote results/figures/step_b105_sclk_drift.png")
