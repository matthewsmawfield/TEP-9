"""Step 138: One-way downlink channel audit (step_b102).

Registers the one-way-transport channel the two-way cancellation
theorem requires.  In a two-way radio link traversing a
proper-time (lapse) gradient, the round trip traverses the same
gradient twice with opposite sign in the observable transport and
the lapse term cancels; a one-way downlink referenced to an
onboard oscillator does not cancel.  This step audits which
one-way deep-space datasets exist, what distances they cover, and
whether a TEP-sensitive residual product can be built from them.

T1  Transport-class mapping: the onboard-oscillator-referenced
    one-way carrier is the required observable class.  The GNSS
    satellite-to-receiver downlink is the same class -- atomic-
    clock referenced, uplink-free -- so the corpus already
    realizes the one-way test through the GNSS-II clock
    correlation channel (step 135 / b99: free axis 0.6 deg from
    the comet meridian plane).

T2  Deep-space archive inventory:
    - Voyager delta-VLBI demonstration, 1979 Jupiter approach
      (22 of 25 scheduled passes usable; goal 50 nrad, achieved
      ~0.5 urad scatter; AIAA 82-1471; IPN PR 42-60) and the
      Saturn-era delta-DOR measurements via the 360 kHz telemetry
      subcarrier.  Downlink-only mode, 5-10 AU.  Observables not
      in a public machine-readable archive.
    - Cassini VLBA astrometry (published position solutions;
      raw delay observables not publicly archived).
    - PRIDE Venus Express / Mars Express open-loop one-way
      tracking (EVN, AuScope, NICT; X-band USO-referenced;
      mHz-level Doppler products; raw correlator-level data at
      the EVN archive; processed products in the published
      record).  Inner-system paths (<<2 AU).

T3  Boundary-distance coverage: the candidate lapse boundary sits
    at heliospheric distances and beyond (>=100 AU).  None of the
    public one-way archives reaches beyond ~10 AU; at 150+ AU the
    Voyager carrier is receivable only by 70-m-class apertures and
    no public VLBI record exists there.  Registered: no public
    one-way dataset covers the boundary-distance regime.

T4  Sensitivity ordering: the 1979 delta-VLBI pass scatter
    (~0.5 urad, ~few ns differential delay) against the
    GNSS-class one-way clock sensitivity already exploited in the
    corpus -- ordering the channels by realized reach.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_138_oneway_downlink_channel")
tee_stdout(logger)
logger.header("One-way downlink channel audit")

import json

# ---- T1 transport-class mapping ------------------------------------
logger.info("T1 transport-class mapping")
T1 = dict(
    requirement=("onboard-oscillator-referenced one-way carrier "
                 "transport across the lapse gradient; the lapse "
                 "term cannot cancel as in two-way links"),
    realized_in_corpus=dict(
        channel="GNSS satellite->receiver downlink clocks",
        reference="onboard atomic clocks (Rb/Cs/H-maser)",
        result="step 135 / b99: free 2701-direction clock axis "
               "0.6 deg from the comet meridian plane, 22.0 deg "
               "from the CMB apex (out-of-sample p=0.010)"),
    note=("GNSS downlink is the same photon-transport class as "
          "VLBI spacecraft downlink: one-way, uplink-free, "
          "onboard-oscillator referenced.  The non-cancelling "
          "channel the cancellation theorem demands is therefore "
          "already realized in the corpus and returns the "
          "programme's strongest cross-channel directional hit."))
logger.info("  GNSS downlink = same transport class; the realized "
            "one-way test is the GNSS-II axis (b99)")

# ---- T2 archive inventory ------------------------------------------
logger.info("T2 deep-space archive inventory")
T2 = dict(
    voyager_delta_vlbi=dict(
        epochs="1979 Jupiter approach; Saturn-era delta-DOR",
        mode="downlink-only delta-VLBI / delta-DOR (360 kHz "
             "telemetry subcarrier harmonics)",
        passes="22 of 25 scheduled usable (1979 demo)",
        achieved_accuracy="~0.5 urad scatter vs 50 nrad goal",
        references=["AIAA 82-1471 (1982)",
                    "IPN Progress Report 42-60",
                    "Border 2009 ISSFD review"],
        distance_au="5-10",
        public_observables=False,
        note="observables not in a public machine-readable "
             "archive; pre-interstellar distances only"),
    cassini_vlba=dict(
        epochs="published astrometric campaigns",
        mode="one-way downlink, phase-referenced to ICRF quasars",
        distance_au="8-10",
        public_observables=False,
        note="position solutions published; raw delay "
             "observables not publicly archived"),
    pride_vex_mex=dict(
        epochs="VEX 2009-2014; MEX campaigns",
        mode="open-loop one-way X-band, USO-referenced, EVN/"
             "AuScope/NICT",
        sensitivity="mHz-level open-loop Doppler (Bocanegra "
                    "Bahamon et al. 2019; Duev et al. 2012)",
        distance_au="<1.5",
        public_observables="raw correlator-level data at the EVN "
                           "archive; processed Doppler products "
                           "in the published record",
        note="inner-system downlink paths only"))
for k, v in T2.items():
    logger.info(f"  {k}: {v['distance_au']} AU, "
                f"public={v.get('public_observables')}")

# ---- T3 boundary-distance coverage ----------------------------------
logger.info("T3 boundary-distance coverage")
T3 = dict(
    boundary_distance_au=">=100 (heliopause class distances)",
    max_public_oneway_au=10.0,
    voyager_at_150au=("carrier receivable only by 70-m-class "
                      "apertures; no public VLBI record exists "
                      "at boundary distances"),
    conclusion=("no public one-way dataset covers the "
                "boundary-distance regime; the deep-space "
                "one-way channel is identified and scoped but "
                "not populated at the relevant distance"))
logger.info("  no public one-way archive reaches >10 AU")

# ---- T4 sensitivity ordering ---------------------------------------
logger.info("T4 sensitivity ordering")
T4 = dict(
    delta_vlbi_1979=dict(scatter_urad=0.5,
                         note="~few-ns differential-delay class "
                              "per pass"),
    gnss_clock_channel=dict(
        note=("sub-ns carrier-phase clock comparisons over a "
              "global network, 25 yr baseline -- the realized "
              "highest-sensitivity one-way channel in the "
              "corpus")),
    ordering=("GNSS clock network >> PRIDE open-loop Doppler > "
              "1979 delta-VLBI demonstration, in realized "
              "one-way sensitivity"))
logger.info("  GNSS clock channel is the realized deepest "
            "one-way channel")

verdict = ("CHANNEL SCOPED: the non-cancelling one-way transport "
           "requirement is already realized in the corpus through "
           "the GNSS downlink clock channel (step 135/b99 axis "
           "hit); deep-space public one-way archives (Voyager "
           "delta-VLBI 1979-81, Cassini VLBA, PRIDE VEX/MEX) are "
           "inner-system or not publicly reduced at observable "
           "level; no public one-way dataset covers the "
           "boundary-distance regime -- registered as an "
           "identified, open channel with named holdings, not a "
           "tested one")
logger.info(f"verdict: {verdict}")

out = dict(step="step_138_oneway_downlink_channel", result="b102",
           verdict=verdict,
           T1_transport_class=T1, T2_archives=T2,
           T3_boundary_coverage=T3, T4_sensitivity=T4)
with open(RESULTS / "step_b102_oneway_channel.json", "w") as f:
    json.dump(out, f, indent=1)
logger.data_save(RESULTS / "step_b102_oneway_channel.json")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(9, 4.6))
rows = [("GNSS clock downlink\n(corpus, realized)", 1.0,
         "steelblue"),
        ("PRIDE VEX/MEX\nopen-loop", 1.5, "gold"),
        ("Voyager delta-VLBI\n1979-81", 10.0, "indianred"),
        ("Cassini VLBA\nastrometry", 10.0, "0.5"),
        ("boundary-distance\nrequirement", 100.0, "white")]
for i, (lbl, d, c) in enumerate(rows):
    edge = "black" if c == "white" else c
    ax.barh(i, d, color=c, edgecolor=edge,
            hatch="//" if c == "white" else None)
    ax.text(d * 1.15, i, f"{d:g} AU", va="center", fontsize=9)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows], fontsize=8)
ax.set_xscale("log")
ax.set_xlabel("heliocentric distance of one-way observable (AU)")
ax.set_xlim(0.5, 400)
ax.axvline(100, color="black", ls=":", lw=0.8)
ax.set_title("One-way downlink channels vs boundary distance "
             "(step 138)")
fig.tight_layout()
fig.savefig(RESULTS / "figures" / "supplementary" / "step_b102_oneway_channel.png",
            dpi=300)
logger.data_save(RESULTS / "figures/supplementary/step_b102_oneway_channel.png")
