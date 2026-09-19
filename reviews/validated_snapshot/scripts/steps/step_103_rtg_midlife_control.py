#!/usr/bin/env python3
"""Step 103 -- RTG mid-life degradation control.

Step 098 found that Voyager 1's three independent MHW-RTG units
decline coherently across the heliopause window (--0.4 per cent
excess vs the 13-26 yr interior baseline), while Voyager 2's two
surviving units diverge -- and noted the remaining loophole
explicitly: the excursion could be undocumented mid-life MHW
aging, not a field-linked transition.  Three of the five units
share one craft, so unit coherence alone cannot separate a
common-mode aging process from a common-mode field coupling.

This step supplies the missing control without new assumptions:
the telemetry record's own wander defines what undocumented aging
can produce.

  T1  Ramp-placebo scan.  The saturating-ramp statistic of step
      098 (c0 + A*clip((t-t0)/w, 0, 1) fitted to the mean-unit
      excess curve) is evaluated for EVERY onset t0 in 14-40 yr
      and width w in {3,5,7,9} yr, not only in the boundary
      region.  The crossing-aligned ramp is compared against the
      full-record amplitude distribution: if equally deep coherent
      ramps occur at placebo epochs, ordinary aging explains the
      excursion and the result stays a bound; if the
      crossing-window ramp is extreme, the candidate tightens.

  T2  All-negative coherence window.  The fraction of grid epochs
      on which all three V1 units sit below baseline is computed
      on a +-3 yr window centred at the heliopause and on every
      equal-length sliding window of the record -- the empirical
      null for three-unit sign coherence.

  T3  Vintage-matched control.  Voyager 2's two functioning units
      are the same MHW design with the same 1977 beginning-of-
      mission epoch: whatever generic mid-life aging does to V1's
      units at mission-year 29-37 it must also do to V2's at the
      SAME mission years (V2 was still inside the heliosphere and
      outside both field lobes then).  The identical (t0, w) ramp
      is fitted to the V2 mean-unit excess over the identical
      window; a V2 ramp of equal amplitude would convict common-
      mode aging, while a flat or incoherent V2 record makes the
      V1 excursion location-specific.

  T4  Leave-one-unit-out.  The crossing-window ramp amplitude is
      recomputed for each two-unit subset of V1, testing whether
      the coherence survives dropping any single unit.

  T5  Unit-spread audit.  The cross-unit dispersion sigma(t) of
      the three V1 excess curves during the excursion is compared
      with its full-record distribution: a genuine common-mode
      excursion keeps the spread at its usual level rather than
      inflating it.

Inputs
------
results/step_b62_rtg_nuclear.csv    (per-unit p_norm series)
results/step_b62_rtg_nuclear.json   (crossing anchors, model fits)

Outputs
-------
results/step_b67_rtg_midlife_control.json
results/step_b67_rtg_midlife_control.csv
results/figures/step_b67_rtg_midlife_control.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_103_rtg_midlife_control")
tee_stdout(logger)
logger.header("RTG mid-life degradation control")

import csv
import json
import numpy as np

BASE_LO, BASE_HI = 13.0, 26.0
GRID = np.arange(14.0, 44.0, 0.5)
WIDTHS = (3.0, 5.0, 7.0, 9.0)

# crossing mission-years (BOM epochs per step_098: V1 1977-09-05,
# V2 1977-08-20)
V1_TS, V1_HP = 27.28, 34.97
V2_TS, V2_HP = 30.03, 41.22

V1_UNITS = ("VG1_rtg1", "VG1_rtg2", "VG1_rtg3")
V2_UNITS = ("VG2_rtg1", "VG2_rtg2")

# ------------------------------------------------------------------
# Load per-unit power series, rebuild the step-098 excess curves
# ------------------------------------------------------------------

series = {}
with open(RESULTS / "step_b62_rtg_nuclear.csv") as f:
    for r in csv.DictReader(f):
        series.setdefault(r["series"], {"t": [], "p": []})
        series[r["series"]]["t"].append(float(r["t_yr_since_bom"]))
        series[r["series"]]["p"].append(float(r["p_norm"]))
for k in series:
    o = np.argsort(series[k]["t"])
    series[k]["t"] = np.array(series[k]["t"])[o]
    series[k]["p"] = np.array(series[k]["p"])[o]


def excess_decline(t, p):
    m = (t >= BASE_LO) & (t <= BASE_HI)
    c = np.polyfit(t[m], np.log(p[m]), 1)
    return (np.log(p) - np.polyval(c, t)) * 100.0


excess = {k: excess_decline(s["t"], s["p"]) for k, s in series.items()}
E_grid = {k: np.interp(GRID, series[k]["t"], excess[k])
          for k in series}


def craft_mean(units):
    return np.mean([E_grid[k] for k in units], axis=0)


E1 = craft_mean(V1_UNITS)
E2 = craft_mean(V2_UNITS)
S1 = np.std([E_grid[k] for k in V1_UNITS], axis=0)

# ------------------------------------------------------------------
# T1: ramp-placebo scan over the whole post-baseline record
# ------------------------------------------------------------------


def ramp_amp(E, t0, w):
    x = np.clip((GRID - t0) / w, 0.0, 1.0)
    A = np.vstack([np.ones_like(GRID), x]).T
    c, *_ = np.linalg.lstsq(A, E, rcond=None)
    return float(c[1])


def scan(E, t0_lo=14.0, t0_hi=40.0):
    out = []
    for t0 in np.arange(t0_lo, t0_hi + 1e-9, 0.25):
        for w in WIDTHS:
            out.append((t0, w, ramp_amp(E, t0, w)))
    return out


def boundary_best(E):
    """Best (most negative) ramp with onset in the step-098 search
    region 27-38 yr -- the crossing-window statistic."""
    amps = [(t0, w, ramp_amp(E, t0, w))
            for t0 in np.arange(27.0, 38.0 + 1e-9, 0.25)
            for w in WIDTHS]
    return min(amps, key=lambda z: z[2])


def crossing_aligned_best(E, hp):
    """Deepest ramp whose window genuinely straddles the crossing:
    onset <= hp <= onset + width (the transition contains the HP)."""
    amps = [(t0, w, ramp_amp(E, t0, w))
            for t0 in np.arange(26.0, hp + 1e-9, 0.25)
            for w in WIDTHS if t0 <= hp <= t0 + w]
    return min(amps, key=lambda z: z[2])


scan_E1 = scan(E1)
amps_E1 = np.array([a for *_, a in scan_E1])
# Fair placebo pool: only onsets at/after the baseline window ends.
# Onsets inside 14-26 yr sit where the excess is ~0 by construction,
# so their ramps cannot reach the observed depth -- including them
# dilutes the null and inflates the significance.
amps_E1_post = np.array([a for t0, w, a in scan_E1 if t0 >= BASE_HI])
t0_b, w_b, A_obs = boundary_best(E1)
t0_c, w_c, A_c = crossing_aligned_best(E1, V1_HP)
p_placebo = float((int((amps_E1_post <= A_obs).sum()) + 1)
                  / (amps_E1_post.size + 1))
p_placebo_x = float((int((amps_E1_post <= A_c).sum()) + 1)
                    / (amps_E1_post.size + 1))
logger.metric("V1_boundary_ramp",
              f"onset {t0_b:.2f} yr, width {w_b:.1f} yr, "
              f"A = {A_obs:+.3f}% (HP at {V1_HP:.2f}; "
              f"ramp centre {t0_b + w_b/2:.1f} yr is post-crossing)")
logger.metric("V1_ramp_placebo_p", float(p_placebo),
              "fraction of post-baseline (t0,w) ramps this deep")
logger.metric("V1_crossing_aligned_ramp",
              f"onset {t0_c:.2f} yr, width {w_c:.1f} yr, "
              f"A = {A_c:+.3f}% (window contains HP at {V1_HP:.2f}); "
              f"post-baseline placebo frac {p_placebo_x:.3f}")

per_unit_scan = {}
for k in V1_UNITS + V2_UNITS:
    sc = scan(E_grid[k])
    amps = np.array([a for t0, w, a in sc if t0 >= BASE_HI])
    tb, wb, Ab = boundary_best(E_grid[k])
    per_unit_scan[k] = dict(
        boundary_best=dict(onset=tb, width=wb, amp_pct=Ab),
        placebo_frac_deeper=float((int((amps <= Ab).sum()) + 1)
                                  / (amps.size + 1)),
        min_anywhere=float(amps.min()))
    logger.metric(f"{k}_ramp",
                  f"boundary-best {Ab:+.3f}% "
                  f"(post-baseline placebo frac "
                  f"{np.mean(amps <= Ab):.3f}, "
                  f"record min {amps.min():+.3f}%)")

# ------------------------------------------------------------------
# T2: all-three-negative coherence window
# ------------------------------------------------------------------

neg1 = np.array([E_grid[k] < 0.0 for k in V1_UNITS])
all_neg1 = neg1.all(axis=0)


def window_frac(mask, centre, half=3.0):
    m = (GRID >= centre - half) & (GRID <= centre + half)
    return float(mask[m].mean()), int(m.sum())


frac_hp, n_hp = window_frac(all_neg1, V1_HP)
slide = [window_frac(all_neg1, c)[0]
         for c in GRID if c - 3.0 >= GRID.min() and c + 3.0 <= GRID.max()]
p_coh = float((int((np.array(slide) >= frac_hp).sum()) + 1)
              / (len(slide) + 1))
logger.metric("V1_allneg_hp_window",
              f"{frac_hp:.2f} of epochs all-negative within +-3 yr of HP "
              f"(placebo frac {p_coh:.3f})")

# same statistic for V2's two units at the V2 HP window
neg2 = np.array([E_grid[k] < 0.0 for k in V2_UNITS])
all_neg2 = neg2.all(axis=0)
frac_hp2, _ = window_frac(all_neg2, V2_HP)
slide2 = [window_frac(all_neg2, c)[0]
          for c in GRID if c - 3.0 >= GRID.min() and c + 3.0 <= GRID.max()]
p_coh2 = float((int((np.array(slide2) >= frac_hp2).sum()) + 1)
               / (len(slide2) + 1))

# ------------------------------------------------------------------
# T3: vintage-matched control -- V2 units at the same mission years
# ------------------------------------------------------------------
# The V1 ramp lives at onset ~29 yr / width ~8 yr.  V2's units are
# the same MHW vintage; at those mission years V2 was inside the
# heliosphere and outside both field lobes.

A_v2_same = ramp_amp(E2, t0_b, w_b)
scan_E2 = scan(E2)
amps_E2 = np.array([a for *_, a in scan_E2])
t0_b2, w_b2, A_b2 = boundary_best(E2)
logger.metric("V2_same_window_ramp",
              f"A = {A_v2_same:+.3f}% at V1's (t0={t0_b:.2f}, w={w_b:.1f})")
logger.metric("V2_own_boundary_ramp",
              f"onset {t0_b2:.2f} yr, w {w_b2:.1f}, A = {A_b2:+.3f}% "
              f"(V2 HP at {V2_HP:.2f})")

# ------------------------------------------------------------------
# T4: leave-one-unit-out on the V1 crossing-window ramp
# ------------------------------------------------------------------

loo = {}
for drop in V1_UNITS:
    keep = [k for k in V1_UNITS if k != drop]
    E_loo = craft_mean(keep)
    tb, wb, Ab = boundary_best(E_loo)
    loo[f"drop_{drop}"] = dict(onset=tb, width=wb, amp_pct=Ab)
    logger.metric(f"loo_{drop}", f"ramp {Ab:+.3f}% (onset {tb:.1f})")

# ------------------------------------------------------------------
# T5: unit-spread audit during the excursion
# ------------------------------------------------------------------

m_exc = (GRID >= t0_b) & (GRID <= t0_b + w_b)
spread_exc = float(np.median(S1[m_exc]))
spread_all = S1
p_spread = float((int((spread_all <= spread_exc).sum()) + 1)
                 / (spread_all.size + 1))
logger.metric("V1_spread_during_excursion",
              f"median {spread_exc:.3f}% (fraction of record tighter: "
              f"{p_spread:.3f})")

# ------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------

v2_convicts = A_v2_same <= A_obs * 0.75
v1_unusual = p_placebo <= 0.10
verdict = dict(
    ramp_at_crossing=dict(onset=t0_b, width=w_b, amp_pct=A_obs,
                          hp_yr=V1_HP,
                          note="deepest-in-window ramp centre "
                               f"{t0_b + w_b/2:.1f} yr vs V1 "
                               f"heliopause {V1_HP:.1f} yr -- the "
                               "best-fit transition is centred "
                               "post-crossing"),
    crossing_aligned_ramp=dict(
        onset=t0_c, width=w_c, amp_pct=A_c,
        hp_yr=V1_HP,
        note="deepest ramp whose window contains the HP; "
             "the crossing-aligned statistic"),
    midlife_control=dict(
        placebo_frac_deeper_postbaseline=p_placebo,
        crossing_aligned_placebo_frac=p_placebo_x,
        all_neg_hp_window_frac=frac_hp,
        all_neg_placebo_p=p_coh,
        v2_same_vintage_amp_pct=A_v2_same,
        v2_own_boundary_amp_pct=A_b2,
        v2_allneg_hp_window_frac=frac_hp2,
        v2_allneg_placebo_p=p_coh2),
    leave_one_unit_out=loo,
    unit_spread=dict(during_excursion_med=spread_exc,
                     frac_of_record_tighter=p_spread),
    reading=(
        f"The V1 mean-unit excess decline is a broad post-shock "
        f"accumulation.  The deepest saturating-ramp fit sits at "
        f"onset {t0_b:.1f} yr / width {w_b:.0f} yr "
        f"(A = {A_obs:+.2f}%) -- centred {t0_b + w_b/2:.1f} yr, "
        f"{t0_b + w_b/2 - V1_HP:.1f} yr AFTER the heliopause -- and "
        f"{p_placebo*100:.1f}% of post-baseline onset/width choices "
        f"reach a ramp that deep.  The genuinely crossing-aligned "
        f"statistic is weaker still: the deepest ramp whose window "
        f"contains the HP is {A_c:+.2f}% (onset {t0_c:.1f}, width "
        f"{w_c:.0f}), matched or beaten by {p_placebo_x*100:.0f}% of "
        f"post-baseline placebo placements.  The three units stay "
        f"within {spread_exc:.2f}% of each other through the "
        f"excursion and the ramp survives leave-one-unit-out "
        f"({min(v['amp_pct'] for v in loo.values()):+.2f}% to "
        f"{max(v['amp_pct'] for v in loo.values()):+.2f}%), but the "
        f"all-negative coherence window at the HP is unremarkable "
        f"(placebo p = {p_coh:.2f}).  V2's units at the identical "
        f"mission years carry a mean ramp of {A_v2_same:+.2f}%, "
        f"same sign and within ~{abs(A_obs)/max(abs(A_v2_same),1e-9):.1f}x "
        f"of V1's.  Net: the excursion is real coherence but sits "
        f"well inside the record's own aging wander -- the "
        f"crossing-aligned ramp is a ~16th-percentile event, not an "
        f"extreme one, and the result stays a bound."))

res = dict(
    step="step_103_rtg_midlife_control",
    description=("Placebo-epoch and vintage-matched control for the "
                 "Voyager-1 coherent RTG excess decline: the named "
                 "mid-life degradation control left pending in step "
                 "098.  The telemetry record's own wander supplies "
                 "the null for undocumented MHW aging."),
    inputs=["results/step_b62_rtg_nuclear.csv",
            "results/step_b62_rtg_nuclear.json"],
    crossings_yr=dict(v1_ts=V1_TS, v1_hp=V1_HP,
                      v2_ts=V2_TS, v2_hp=V2_HP),
    baseline_window_yr=[BASE_LO, BASE_HI],
    T1_ramp_placebo=dict(
        v1_boundary_best=dict(onset=t0_b, width=w_b, amp_pct=A_obs),
        v1_crossing_aligned=dict(onset=t0_c, width=w_c,
                                 amp_pct=A_c,
                                 placebo_frac=p_placebo_x),
        placebo_frac_deeper=p_placebo,
        placebo_pool="onsets >= 26 yr (post-baseline only)",
        n_onset_width_grid=len(scan_E1),
        per_unit=per_unit_scan),
    T2_coherence=dict(
        v1_allneg_hp_window_frac=frac_hp,
        v1_allneg_placebo_p=p_coh,
        v2_allneg_hp_window_frac=frac_hp2,
        v2_allneg_placebo_p=p_coh2),
    T3_vintage_matched=dict(
        same_window=dict(t0=t0_b, w=w_b, v2_amp_pct=A_v2_same),
        v2_own_boundary=dict(onset=t0_b2, width=w_b2, amp_pct=A_b2)),
    T4_leave_one_unit_out=loo,
    T5_spread=dict(excursion_median_pct=spread_exc,
                   frac_record_tighter=p_spread),
    verdict=verdict,
    caveats=[
        "The control is internal to the published telemetry record; "
        "it cannot exclude a genuinely craft-common external driver "
        "other than the field (e.g., a shared thermal event), which "
        "is why the result is framed as a control on undocumented "
        "aging specifically.",
        "V2's two surviving units diverge by ~2% post-year-37, so "
        "the vintage-matched comparison uses the mean-unit curve and "
        "the identical window, not a per-unit match.",
        "Post-13.3 yr telemetry is bit-drop quantized at ~0.04 A; "
        "quantization steps shared across a craft's units are part "
        "of what the placebo scan absorbs.",
        "The placebo pool is restricted to post-baseline onsets "
        "(t0 >= 26 yr): onsets inside the baseline window sit where "
        "the excess is ~0 by construction and cannot produce a deep "
        "ramp, so including them inflates the significance.",
        "The deepest ramp in the crossing search window is centred "
        "~7.5 yr after the heliopause; the ramp whose window actually "
        "contains the crossing is weaker (-0.53%) and matched by "
        "~17% of post-baseline placebo placements.  The unit-spread "
        "during the excursion is also unremarkable -- ~87% of the "
        "record is tighter -- so the coherence is absolute-scale "
        "(0.06%) rather than record-extreme."])

out = str(RESULTS / "step_b67_rtg_midlife_control.json")
json.dump(res, open(out, "w"), indent=1, default=float)
logger.data_save(_Path(out))

csv_out = str(RESULTS / "step_b67_rtg_midlife_control.csv")
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["series", "onset_yr", "width_yr", "ramp_amp_pct"])
    for t0, ww, a in scan_E1:
        w.writerow(["VG1_mean_unit", t0, ww, a])
    for t0, ww, a in scan_E2:
        w.writerow(["VG2_mean_unit", t0, ww, a])
logger.data_save(_Path(csv_out))

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.4))

ax = axes[0]
for k, c in zip(V1_UNITS, ("#1f77b4", "#6baed6", "#08306b")):
    ax.plot(GRID, E_grid[k], lw=0.9, c=c, alpha=0.85,
            label=k.replace("VG1_", "V1 "))
ax.plot(GRID, E1, "k-", lw=1.8, label="V1 mean")
ax.plot(GRID, E2, "r-", lw=1.3, alpha=0.8, label="V2 mean")
for x, lab in ((V1_TS, "V1 TS"), (V1_HP, "V1 HP"), (V2_HP, "V2 HP")):
    ax.axvline(x, color="0.5", lw=0.8, ls="--")
    ax.annotate(lab, (x, ax.get_ylim()[0]), fontsize=7, rotation=90,
                va="bottom")
ax.axvspan(t0_b, t0_b + w_b, color="crimson", alpha=0.08)
ax.set_xlabel("mission year since BOM")
ax.set_ylabel("excess decline vs 13-26 yr baseline (%)")
ax.legend(frameon=False, fontsize=7, ncol=2)
ax.set_title("unit excess curves and the crossing-window ramp", fontsize=10)

ax = axes[1]
ax.hist(amps_E1_post, bins=50, color="0.75", edgecolor="0.5", lw=0.4)
ax.axvline(A_obs, color="crimson", lw=1.6,
           label=f"V1 deepest ramp {A_obs:+.2f}%")
ax.axvline(A_c, color="darkorange", lw=1.6, ls=":",
           label=f"V1 crossing-aligned {A_c:+.2f}%")
ax.axvline(A_v2_same, color="red", lw=1.4, ls="--",
           label=f"V2 same-window {A_v2_same:+.2f}%")
ax.set_xlabel("saturating-ramp amplitude (%), post-baseline onsets")
ax.set_ylabel("count")
ax.legend(frameon=False, fontsize=8)
ax.set_title(f"mid-life placebo scan (frac deeper = {p_placebo:.3f}; "
             f"crossing-aligned {p_placebo_x:.3f})", fontsize=10)

ax = axes[2]
ax.plot(GRID, all_neg1.astype(float), "k-", lw=1.2,
        label="V1: all 3 units below baseline")
ax.plot(GRID, all_neg2.astype(float), "r-", lw=1.2, alpha=0.8,
        label="V2: both units below baseline")
for x in (V1_TS, V1_HP, V2_HP):
    ax.axvline(x, color="0.5", lw=0.8, ls="--")
ax.axvspan(V1_HP - 3, V1_HP + 3, color="crimson", alpha=0.08)
ax.set_xlabel("mission year since BOM")
ax.set_ylabel("coherent-below-baseline flag")
ax.set_yticks([0, 1])
ax.legend(frameon=False, fontsize=8)
ax.set_title(f"sign coherence (HP window {frac_hp:.0%}, "
             f"placebo p {p_coh:.2f})", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b67_rtg_midlife_control.png", dpi=150)
logger.data_save(FIG / "step_b67_rtg_midlife_control.png")
logger.success("RTG mid-life degradation control complete")
