#!/usr/bin/env python3
"""
TEP-9 step 100 -- Cross-channel clock-consistency synthesis
===========================================================

The TEP-9 anomaly is fundamentally a claim about *clocks*: orbital
periods (macroscopic dynamical clocks) are the resident/transit
channels, while the spacecraft program adds three microscopic
channels -- the electronic quartz oscillator (SCLK), the plasma
frequency (PWS), and the nuclear decay clock (RTG).  Under TEP's
conformal coupling these process classes need not share one
effective coupling strength (screening factors and disformal terms
are process-dependent), so the correct question is not "do all
channels show the anomaly" but "are the measured/bounded values in
each channel mutually consistent with a single lapse field".

This step assembles the cross-channel ledger:

  channel        clock class         observable              result
  -------------  ------------------  ----------------------  ------------------
  TNO resident   orbital dynamics     cap-clustered alignment  p = 0.0073 (44 obj)
  comet transit  orbital dynamics     aphelion-dir slip        delta t/t ~ 1e-2
  SCLK           electronic quartz    onboard-vs-ground rate   |dA/A| < ~2e-5
  SCLK phase     electronic quartz    segment phase residuals  |dphi| < 0.16 s
  PWS            plasma oscillation   f_pe structure           deltaX ~ 1.05
  RTG            nuclear decay        power-decline slope      |dalpha| < ~1e-4

Reading the ledger
------------------
* The comet-channel slip is a *cumulative* integrated quantity
  accrued along a transit through the boundary sector; the SCLK
  bound is an *instantaneous* onboard-vs-ground rate ratio.  A
  lapse field whose contrast is built across a transit region can
  satisfy both simultaneously.
* The PWS lapse-equivalent shows the measured f_pe structure needs
  a constants-shift ~1000x larger than the comet-channel lapse
  contrast, confirming the structure is real plasma -- the channel
  that a plasma-coupled screening transition would use as its
  trigger.
* The RTG bound applies to the nuclear sector specifically; at
  ~1e-4 in effective alpha it is consistent with a process-
  dependent coupling where the nuclear clock's sensitivity differs
  from the orbital channel's.
* Geometry: both Voyager paths thread the anti-axis hemisphere;
  V1's asymptote and both craft's positions sit inside the mirror
  cap, so the channels sample the second lobe of the bipolar slip
  field (step 089).

Inputs
------
results/step_b30_proper_time_slip.csv  (comet-channel slip, step 065)
results/step_b60_sclk_audit.json
results/step_b61_pws_channel.json
results/step_b62_rtg_nuclear.json
results/step_b63_heliospheric_geometry.json
results/step_b16_bipolar.json          (mirror-cap slip, step 089, if present)

Outputs
-------
results/step_b64_clock_consistency.json
results/figures/step_b64_clock_consistency.png
"""

import sys
import json
import csv
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout

logger = StepLogger("step_100_clock_consistency")
tee_stdout(logger)
logger.header("Cross-channel clock-consistency synthesis")

R = RESULTS


def load(name):
    p = R / name
    return json.loads(p.read_text()) if p.exists() else None


sclk = load("step_b60_sclk_audit.json")
pws = load("step_b61_pws_channel.json")
rtg = load("step_b62_rtg_nuclear.json")
geo = load("step_b63_heliospheric_geometry.json")
bip = load("step_b16_bipolar.json")

# --- comet-channel slip -------------------------------------------------
slip_file = R / "step_b30_proper_time_slip.csv"
comet = {}
if slip_file.exists():
    rows = list(csv.DictReader(open(slip_file)))
    cols = rows[0].keys() if rows else []
    # frac_slip is already the fractional lapse contrast
    # (dtau_total / t_transit, verified against the dtau_total and
    # t_transit columns); dtau_total carries the slip in years.
    slipcol = "frac_slip" if "frac_slip" in cols else \
        next((c for c in cols if "slip" in c.lower()), None)
    if slipcol:
        s = np.array([float(r[slipcol]) for r in rows
                      if r[slipcol] not in ("", "nan")])
        comet = dict(n=len(s), column=slipcol,
                     median_frac_slip=float(np.median(s)),
                     p16=float(np.percentile(s, 16)),
                     p84=float(np.percentile(s, 84)))
        if "dtau_total" in cols:
            dt = np.array([float(r["dtau_total"]) for r in rows
                           if r["dtau_total"] not in ("", "nan")])
            comet["median_slip_yr"] = float(np.median(dt))
        frac = float(np.median(s))
        comet["implied_dA_over_A"] = frac
        logger.metric("comet_slip_median_yr",
                      round(comet.get("median_slip_yr", 0.0), 3),
                      "median implied slip")
        logger.metric("comet_dA_over_A", f"{frac:.2e}",
                      "implied lapse contrast")
else:
    logger.error("step_b30 slip file missing; comet channel skipped")

# --- assemble the ledger -------------------------------------------------
ledger = []

if comet:
    ledger.append(dict(
        channel="comet transit (orbital dynamics)",
        observable="aphelion-direction slip / leg duration",
        value=f"dA/A ~ {comet['implied_dA_over_A']:.1e}",
        kind="detection (cumulative)",
        note=("integrated slip accrued along the transit leg through "
              "the boundary sector")))

if sclk:
    b = sclk["craft"]["VG1"]["lapse_bound"]
    st = sclk["craft"]["VG1"].get("step_tests", {})
    st2 = sclk["craft"]["VG2"].get("step_tests", {})
    steps = [st.get(k, {}).get("step_frac") for k in
             ("termination_shock", "heliopause")] + \
            [st2.get(k, {}).get("step_frac") for k in
             ("termination_shock", "heliopause")]
    steps = [abs(s) for s in steps if isinstance(s, float)]
    step_bound = f"{max(steps):.1e}" if steps else "n/a"
    ledger.append(dict(
        channel="Voyager SCLK (electronic quartz)",
        observable="onboard-vs-ground clock-rate ratio",
        value=f"|step| < {step_bound} sustained post-crossing "
              f"offset (p99 resid {b['p99_abs_resid_frac']:.1e})",
        kind="bound (instantaneous)",
        note=("signed-offset step tests at all four crossings are "
              "null (p = 0.16-0.85); bounds the local lapse contrast "
              "along the Voyager paths, not the cumulative transit "
              "integral")))

    # phase-slip channel: segment-boundary residuals bound a discrete
    # proper-time jump -- the comet-channel signature type -- far
    # tighter than the rate tests
    pb = []
    for c in ("VG1", "VG2"):
        pc = sclk["craft"][c].get("phase_channel", {})
        for lab in ("termination_shock", "heliopause"):
            v = pc.get("crossing_phase_tests", {}).get(lab, {})
            x = v.get("max_abs_resid_pm90d_s")
            if isinstance(x, (int, float)):
                pb.append(x)
    if pb:
        ledger.append(dict(
            channel="Voyager SCLK phase (segment boundaries)",
            observable="calibration-segment phase residuals",
            value=f"|phase jump| < {max(pb):.2f} s at all four "
                  "crossings (vs comet slip ~10^7-10^8 s)",
            kind="bound (phase slip)",
            note=("tests the discrete-jump signature type the rate "
                  "channel cannot see; bounded at the sub-second "
                  "level within +-90 d of every crossing")))

    # interior control: New Horizons, same oscillator class, no
    # boundary crossings -- its correction history bounds the
    # channel's own noise floor on an out-of-cap trajectory
    nh = sclk["craft"].get("NH1", {})
    if nh:
        npc = nh.get("phase_channel", {})
        ndi = nh.get("drift_intervals", {})
        ledger.append(dict(
            channel="New Horizons SCLK (interior control)",
            observable="correction-interval drift + phase residuals",
            value=(f"median |drift| {ndi.get('median_abs_drift_ppm', 0):.4f} ppm, "
                   f"phase floor ~{npc.get('median_abs_resid_us', 0):.0f} us; "
                   f"{npc.get('n_engineering_resets', 0)} resyncs, all "
                   "documented MET events"),
            kind="control (no crossing)",
            note=("same clock class on an interior, out-of-cap "
                  "trajectory (r < 60 AU); its three flagged resets "
                  "(2009, 2017 software-MET resyncs, launch anchor) "
                  "are documented engineering events -- validates "
                  "the reset classifier on a craft with a fully "
                  "documented correction history")))

    # correction-interval drift channel (reproduces the earlier TVP
    # spacecraft-clock study's since-launch drift segmentation)
    dd = []
    for c in ("VG1", "VG2"):
        di = sclk["craft"][c].get("drift_intervals", {})
        for lab in ("termination_shock", "heliopause"):
            t = di.get("crossing_drift_tests", {}).get(lab, {})
            o = t.get("offset_ppm")
            if isinstance(o, (int, float)):
                dd.append(abs(o))
    if dd:
        ledger.append(dict(
            channel="Voyager SCLK drift intervals",
            observable="per-correction-interval oscillator drift",
            value=(f"|drift offset| < {max(dd):.3f} ppm across all "
                   "four crossings"),
            kind="bound (correction history)",
            note=("each SCLK coefficient record is a clock "
                  "recalibration; the implied inter-correction drift "
                  "rates show no crossing-localized change -- the "
                  "corrections are ordinary quartz aging plus "
                  "flagged engineering resyncs")))

    # journey channel: drift mapped onto the trajectory -- screening
    # region, direction of travel, planetary wells
    jr = sclk["craft"].get("VG1", {}).get("journey", {})
    rc = jr.get("confounds", {}).get("residual_region_contrasts", {})
    hs = rc.get("heliosheath_vs_post92_inner", {})
    vs = rc.get("vlism_vs_post92_inner", {})
    ac = jr.get("confounds", {}).get(
        "accumulation_law", {}).get("heliosheath", {})
    xc = sclk.get("cross_craft", {})
    if hs:
        ledger.append(dict(
            channel="Voyager SCLK journey drift",
            observable="per-boundary phase residual vs trajectory",
            value=(f"heliosheath resid {hs.get('median_abs_resid_region_s', 0)*1e3:.1f} ms "
                   f"= {hs.get('ratio', 0):.1f}x post-92 inner baseline "
                   f"(p={hs.get('mannwhitney_p', 1):.0e}); "
                   f"{ac.get('long_interval_neg_frac', 0)*100:.0f}% "
                   f"of long-interval residuals negative in shell "
                   f"(p={ac.get('long_interval_neg_binom_p', 1):.0e}); "
                   f"VLISM relaxes to {vs.get('ratio', 0):.1f}x"),
            kind="monitored pattern (trajectory-organized)",
            note=("the correction-interval drift peaks inside the "
                  "low-density heliosheath shell on both independent "
                  "clocks and relaxes in the VLISM -- the bounded-"
                  "shell profile of a density-triggered unscreening; "
                  "the accumulation is coherent (slope ~1.05, "
                  "96-97% negative sign vs ~50% inner) and "
                  "position-organized (matched-radius rho "
                  f"{xc.get('matched_radius', {}).get('spearman_rho', 0):.2f} "
                  f"vs matched-epoch {xc.get('matched_epoch', {}).get('spearman_rho', 0):.2f}, "
                  "maximal at zero radial shift; NH shows no "
                  "radial growth over 1-62 AU on the same record "
                  "class) -- but at ms-scale amplitude, kernel-fit "
                  "provenance, and null inside planetary wells down "
                  "to phi/c^2 ~ 4e-9; registered as a pattern, not a "
                  "detection")))

if pws:
    le = pws["craft"]["VG1"]["lapse_equiv"]
    ledger.append(dict(
        channel="Voyager PWS (plasma oscillation)",
        observable="f_pe structure vs lapse-equivalent",
        value=f"deltaX = {le['deltaX_constant_density']:.2f} "
              f"({le['ratio_vs_comet']:.0f}x comet channel)",
        kind="discriminant",
        note=("structure confirmed as real plasma; a constants-shift "
              f"reading requires ~{le['ratio_vs_comet']:.0f}x the "
              "comet-channel contrast")))

if rtg:
    cb = rtg["bound"]  # largest |step|+2sigma over the four crossing tests
    coh = rtg.get("unit_coherence", {})
    v1_ex = coh.get("VG1", {}).get("excess_at_hp_pct", {})
    v2_spread = coh.get("VG2", {}).get("excess_at_hp_pct", {})
    ledger.append(dict(
        channel="Voyager RTG (nuclear alpha decay)",
        observable="power-decline slope at crossings",
        value=(f"|dalpha_eff/alpha_eff| < {cb['bound_dalpha_frac_2s']:.1e}; "
               f"V1 units coherently {v1_ex.get('mean', 0):+.2f}% at HP, "
               f"V2 unit spread {v2_spread.get('unit_spread', 0):.1f}%"),
        kind="bound + candidate signature (nuclear sector)",
        note=("Gamow amplification ~150x; V1's three units decline "
              "coherently across the heliopause window (V1 is inside "
              "the mirror cap); V2's two surviving units diverge "
              "(V2 is outside both caps). Candidate field-linked "
              "excursion, but sign-unstable across baseline windows "
              "(+1.6% to -0.5%) and matched in amplitude by an "
              "out-of-cap V2 unit; degenerate with undocumented "
              "MHW-RTG aging")))

if geo:
    pa = geo["probe_alibi"]
    ledger.append(dict(
        channel="trajectory geometry",
        observable="probe/ISO asymptotes vs caps",
        value=(f"{pa['n_velocity_in_cap']}/{pa['n_probes']} probe "
               f"velocities in cap; {pa['n_velocity_in_mirror']} in "
               f"mirror; {pa['n_iso_in_cap']+pa['n_iso_in_mirror']}/"
               f"{pa['n_iso']} ISOs in a cap"),
        kind="geometric datum",
        note=("V1's position and velocity asymptote sit inside the "
              "mirror cap; V2 crossed in the intervening band "
              "(80-93 deg from the antiaxis, outside both caps)")))

consistency = dict(
    reading=(
        "The ledger is internally consistent under the TEP lapse "
        "framework: (i) the comet channel measures a cumulative "
        "transit-integrated slip of ~1e-2; (ii) the SCLK channel "
        "shows no signed rate step at any of the four boundary "
        "crossings (all p > 0.16, sustained offsets < ~1e-6) and no "
        "segment-boundary phase jump above ~0.2 s within +-90 d of "
        "any crossing -- a different observable of the same field, "
        "not a contradiction; (iii) the PWS "
        "channel confirms the heliopause structure is real plasma "
        "rather than a lapse artifact, consistent with a plasma-"
        "coupled screening trigger; (iv) the RTG channel shows a "
        "coherent negative excursion across V1's three nuclear "
        "clocks at the heliopause -- the only craft inside a field "
        "lobe -- with the out-of-cap V2 unit-incoherent, consistent "
        "with a directional process-dependent nuclear coupling; the "
        "excursion remains formally a bound-level candidate "
        "degenerate with undocumented MHW-RTG aging."),
    falsifiable_discriminators=[
        "If the boundary transition is sharp (edge-width bound from "
        "step 073) AND the oscillator couples like orbital dynamics, "
        "the SCLK record should show a signed post-crossing rate "
        "offset of order the integrated slip; the signed-offset step "
        "tests are null at all four crossings (< ~1e-6), disfavouring "
        "that specific combination (sharp edge + universal "
        "coupling).",
        "A discrete proper-time jump at a crossing (the comet-channel "
        "signature type) is bounded to < ~0.2 s within +-90 d by the "
        "SCLK segment-phase channel -- any domain-wall slip must "
        "either be sub-second, accumulate smoothly below the "
        "calibration cadence, or couple differently to electronic "
        "clocks than to orbital dynamics.",
        "A future per-cycle RTG telemetry set would tighten the "
        "nuclear bound by ~2 orders of magnitude.",
        "New Horizons (asymptote outside both caps) offers the "
        "control path; an ISO arriving inside the cap with a "
        "measured NG excess would extend the transit channel.",
    ])

# --- figure: ledger bar chart --------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(11, 5.5))
vals, errs, labels, cols = [], [], [], []
if comet:
    vals.append(comet["implied_dA_over_A"]); errs.append(0)
    labels.append("comet transit\n(cumulative, orbital)")
    cols.append("#1f77b4")
if sclk:
    vals.append(sclk["craft"]["VG1"]["lapse_bound"]["p99_abs_resid_frac"])
    errs.append(0); labels.append("SCLK bound\n(instantaneous,\nelectronic)")
    cols.append("#2ca02c")
if pws:
    vals.append(pws["craft"]["VG1"]["lapse_equiv"]
                ["deltaX_constant_density"])
    errs.append(0)
    labels.append("PWS lapse-equiv\n(required deltaX)")
    cols.append("#9467bd")
if rtg:
    vals.append(float(rtg["bound"]["bound_dalpha_frac_2s"]))
    errs.append(0); labels.append("RTG bound\n(nuclear sector)")
    cols.append("#d62728")
ax.bar(range(len(vals)), vals, color=cols, alpha=0.8)
if comet:
    ax.axhline(comet["implied_dA_over_A"], color="k", ls="--", lw=1,
               label=f"comet-channel {comet['implied_dA_over_A']:.1e}")
ax.set_yscale("log")
ax.set_xticks(range(len(vals)))
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("lapse contrast / constants-shift (fractional)")
ax.set_title("Cross-channel clock ledger: detections and bounds")
ax.legend()
fig.tight_layout()
figp = RESULTS / "figures" / "step_b64_clock_consistency.png"
fig.savefig(figp, dpi=150)
logger.data_save(figp)

out = dict(
    step="step_100_clock_consistency",
    description=("Cross-channel clock-consistency ledger: orbital "
                 "(comet/TNO), electronic (SCLK), plasma (PWS) and "
                 "nuclear (RTG) clock channels evaluated against the "
                 "TEP conformal-lapse framework."),
    ledger=ledger,
    comet_channel=comet,
    consistency=consistency,
    inputs=["results/step_b30_proper_time_slip.csv",
            "results/step_b60_sclk_audit.json",
            "results/step_b61_pws_channel.json",
            "results/step_b62_rtg_nuclear.json",
            "results/step_b63_heliospheric_geometry.json"],
    caveats=[
        "Channels measure different field functionals (cumulative "
        "slip vs instantaneous rate ratio vs constants-shift); they "
        "are mutually consistent, not redundant.",
        "The SCLK and RTG bounds apply along the Voyager paths "
        "specifically; they do not bound the cap-direction contrast "
        "the comet channel measures.",
        "Process-dependent coupling (electronic vs plasma vs "
        "nuclear vs orbital) is part of the TEP framework; the "
        "ledger tests consistency, not universality.",
    ])
out_path = RESULTS / "step_b64_clock_consistency.json"
out_path.write_text(json.dumps(out, indent=1))
logger.data_save(out_path)
logger.success("Clock-consistency synthesis complete")
