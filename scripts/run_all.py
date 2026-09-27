#!/usr/bin/env python3
"""
TEP-9 Pipeline - Run All Steps
==============================

Runs every script listed in ``CORE_STEPS`` in order. The order is
explicitly defined in this file.

Phases
------
1. Data acquisition (001-003): JPL SBDB, CDS/OSSOS+Warsaw, CODE.
2. Resident signature (010-028): detached-TNO clustering, bias
   gauntlet, footprint nulls, structure diagnostics.
3. Transit signature (030-041): Warsaw/CODE comet clock channel,
   three-leg decomposition, planetary baseline, alternatives.
4. Synthesis (050-051): one-map convergence, evidence summary.
5. Deepened controls (058-063): leverage audits, axis consistency,
   aphelion dipole, Planet Nine insertion, bidirectional validation.
6. Replication & TEP observables (064-094): Warsaw N-body
   replication, implied proper-time slip, detached-cluster secular
   stability with/without Planet Nine, prospective detection power,
   perturber-zoo closure, MPCORB replication, joint axis tests.
7. Spacecraft clock & nuclear channels (095-100): Voyager SCLK
   clock-correction audit, PWS plasma-frequency channel, RTG
   nuclear-decay channel, heliospheric geometry and probe/ISO
   asymptotes, cross-channel clock-consistency synthesis.
8. Predictive validation & census channels (101-116): ISO
   non-gravitational channel, slip-map cross-validation and
   spatial holdout, RTG mid-life control, third-catalogue
   bidirectional replication, harmonic axis transfer, named-direction
   audit, CometEls aphelion-dipole census and bias audits,
   directional-vs-axisymmetric decomposition, detachment-coherence
   map, OSSOS ensemble audit, resident ledger, claims-trace audit,
   mixture posterior and empirical discovery-coupling audits.
9. Prospective lineage holdout and independent-refit synthesis
   (117-154): post-2017 SBDB LPC prospective cohort through the
   identical bidirectional instrument, displaced-dipole axis audit,
   dipole-mechanism identification, pre-2018 era-by-lineage factorial,
   masked-signal discrimination, anomaly localization, leg-coherence
   decomposition, the reconstruction-free lineage-independent
   aphelion dipole, displaced-structure temporal characterization,
   the ten-channel global cross-survey synthesis, the non-Warsaw
   two-leg refits of the transit record on raw MPC astrometry
   (pre-2018, post-2017, and the LPC 1902-1950 signature cohort),
   frame-anchored bipolar unification, the epoch-resolved axis
   trajectory, bipolar cone-shell/meridian/harmonic coherence,
   cross-channel axis convergence against the sibling clock and
   lunar pipelines, the OMNI solar-wind screening-epoch test, the
   RTG cross-generator regression, one-way downlink and SCLK
   signed-drift/radial-gradient channels, the frozen-geometry
   predictive battery, astrometric-catalogue stratification, the
   Gaia covariance-sink audit, the comet non-gravitational channel,
   the giant-planet secular apsidal channel, the slip-injection
   transfer function, the arrival-latitude conditioned nulls, the
   zonal-distortion null, and the SCLK link-plasma regression.
10. Forward-model and generative closure (155-164): the 1I/'Oumuamua
   slip-injection to non-gravitational-term transfer test, the
   boundary-radius forward-model scan, the conventional generative
   population model (isotropic Oort-spike source with secular
   Galactic-tide injection and empirical discovery selection), the
   metric-derived winding-holonomy slip field tested against the
   measured equivalent-time offsets, the perihelion-time leg-slip
   channel, the independent two-leg refit of the interstellar
   objects, the APDB planetary astrometric O-C channel, the NIMA
   detached-resident O-C pilot, the forward-modelled comet discovery
   footprint, and the TEP-side generative model for the arrival
   anisotropy.
11. Publication validation (step 114, run last): the exact
   claims-trace audit binding every cited step, result, figure and
   headline number to a leaf in ``results/``.

Each step writes ``logs/<step_name>.log``; ``logs/pipeline.log``
records the orchestrator run.  Machine-readable outputs land in
``results/``; figures in ``results/figures/``.
"""

import sys
import subprocess
import time
import json
import fcntl
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Tuple, Dict

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STEPS_DIR = PROJECT_ROOT / 'scripts' / 'steps'
LOGS_DIR = PROJECT_ROOT / 'logs'

LOGS_DIR.mkdir(parents=True, exist_ok=True)

LOCK_FILE = LOGS_DIR / "pipeline.lock"


def acquire_pipeline_lock():
    """Exclusive run lock: a second concurrent run_all.py would
    interleave pipeline.log and overwrite results/provenance
    mid-flight, corrupting the audit trail.  The lock is held for the
    life of this process (flock releases on exit, even on a crash).
    Returns the open fd, or None when another run holds the lock."""
    fd = open(LOCK_FILE, "w")
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        fd.close()
        return None
    fd.write(f"pid={os.getpid()} started_utc="
             f"{datetime.now(timezone.utc).isoformat()}\n")
    fd.flush()
    return fd

CORE_STEPS: List[Tuple[str, str]] = [
    # Phase 1: Data acquisition & provenance
    ('step_001_download_sbdb.py', 'Step 001: JPL SBDB download'),
    ('step_002_download_cds.py', 'Step 002: CDS download (OSSOS + Warsaw)'),
    ('step_003_download_code.py', 'Step 003: CODE catalogue acquisition'),
    ('step_004_download_mpc.py', 'Step 004: MPC / JPL element files'),

    # Phase 2: Resident signature -- extreme-TNO spatial alignment
    ('step_010_tno_clustering.py', 'Step 010: Detached-TNO orientation clustering'),
    ('step_011_discriminators.py', 'Step 011: Field-wall vs point-mass discriminators'),
    ('step_012_radial_profile.py', 'Step 012: Radial turn-on profile'),
    ('step_013_empirical_null.py', 'Step 013: Empirical discovery-filter null'),
    ('step_014_discovery_footprint.py', 'Step 014: Per-object discovery footprint'),
    ('step_015_proxy_validation.py', 'Step 015: Footprint proxy validation'),
    ('step_016_footprint_conditional.py', 'Step 016: Distance-aware conditional null'),
    ('step_017_decorrelation_radius.py', 'Step 017: Axis-footprint decorrelation radius'),
    ('step_018_plane_warp.py', 'Step 018: Coherent mean-plane warp'),
    ('step_019_pole_geometry.py', 'Step 019: Orbit-pole eigenstructure'),
    ('step_020_axis_coupling.py', 'Step 020: Axis-element coupling'),
    ('step_021_environment_checks.py', 'Step 021: Galactic-tide and omega decomposition'),
    ('step_022_edges.py', 'Step 022: Radial edge structure'),
    ('step_023_robustness.py', 'Step 023: Jackknife/era/quality robustness'),
    ('step_024_population_selectivity.py', 'Step 024: Resident-vs-transit selectivity'),
    ('step_025_patch_scale.py', 'Step 025: Intrinsic patch scale'),
    ('step_026_antipode_deficit.py', 'Step 026: Anti-axis deficit'),
    ('step_027_antipodal_populations.py', 'Step 027: JFC bipolar population'),
    ('step_028_ossos_bias.py', 'Step 028: OSSOS pointing coupling'),

    # Phase 3: Transit signature -- cometary clock channel
    ('step_030_warsaw_direction.py', 'Step 030: Warsaw comet directional test'),
    ('step_031_clock_channel.py', 'Step 031: Orbit-reconstruction clock channel'),
    ('step_032_nested_field.py', 'Step 032: Nested-field multi-axis test'),
    ('step_033_three_leg.py', 'Step 033: Three-leg orbit decomposition'),
    ('step_034_code_confirm.py', 'Step 034: CODE independent confirmation'),
    ('step_035_consolidation.py', 'Step 035: Confirmation consolidation'),
    ('step_036_axis_convergence.py', 'Step 036: Independent axis recovery'),
    ('step_037_blind_validation.py', 'Step 037: Blind validation battery'),
    ('step_038_model_comparison.py', 'Step 038: Model comparison and controls'),
    ('step_039_signed_channels.py', 'Step 039: Signed channels and matched control'),
    ('step_040_planetary_baseline.py', 'Step 040: Planetary baseline residual'),
    ('step_041_alternatives.py', 'Step 041: Conventional alternatives dismantled'),
    ('step_042_nbody_planetary_baseline.py', 'Step 042: N-body planetary baseline (REBOUND/DE440s)'),

    # Phase 4: Synthesis -- convergence on a single domain boundary
    ('step_050_one_map.py', 'Step 050: One-map axis coincidence'),
    ('step_051_evidence_summary.py', 'Step 051: Multiple-testing evidence summary'),
    ('step_052_injection_chain.py', 'Step 052: Injected-population depth chain'),
    ('step_053_catalog_concordance.py', 'Step 053: MPC catalogue concordance'),
    ('step_054_prospective_predictions.py', 'Step 054: Dated falsifiable predictions'),
    ('step_055_morphology.py', 'Step 055: Cluster threshold morphology'),
    ('step_056_morphology_conditional.py', 'Step 056: Morphology vs footprint null'),
    ('step_057_incap_audit.py', 'Step 057: In-cap property audit and synthesis map'),

    # Phase 5: Deepened controls -- observational leverage, axis
    # consistency, uncertainty-normalized and third-lineage channels
    ('step_058_comet_obs_audit.py', 'Step 058: Comet observational-leverage audit'),
    ('step_059_axis_consistency.py', 'Step 059: Axis-consistency test'),
    ('step_060_uncertainty_channel.py', 'Step 060: Uncertainty channel and MW08 third lineage'),
    ('step_061_aphelion_dipole.py', 'Step 061: Aphelion-sky dipole and tide calibration'),
    ('step_062_planet_nine_insertion.py', 'Step 062: Planet Nine insertion test (REBOUND/DE440s)'),
    ('step_063_bidirectional_rotation.py', 'Step 063: Bidirectional rotation residual (REBOUND/DE440s)'),

    # Phase 6: Independent-catalogue replication, TEP-native observables,
    # dynamical stability, and prospective power
    ('step_064_warsaw_bidirectional.py', 'Step 064: Warsaw bidirectional validation (REBOUND/DE440s)'),
    ('step_065_proper_time_slip.py', 'Step 065: Implied proper-time slip'),
    ('step_066_tno_secular_stability.py', 'Step 066: Detached-cluster secular stability +/- P9 (REBOUND/DE440s)'),
    ('step_067_detection_power.py', 'Step 067: Prospective detection power'),
    ('step_068_wall_vs_field.py', 'Step 068: Wall-vs-field slip scaling discriminator'),
    ('step_069_joint_axis_test.py', 'Step 069: Joint axis-permutation global test'),
    ('step_070_galactic_tide.py', 'Step 070: Galactic-tide insertion control (REBOUND/DE440s)'),
    ('step_071_radial_slip_profile.py', 'Step 071: Radial slip profile, multi-shell legs (REBOUND/DE440s)'),
    ('step_072_epoch_stability.py', 'Step 072: Epoch stability of the implied slip'),
    ('step_073_edge_width.py', 'Step 073: Boundary edge width, logistic fit'),
    ('step_074_loco_joint.py', 'Step 074: Leave-one-channel-out joint axis test'),
    ('step_075_signed_slip.py', 'Step 075: Signed slip direction, lag vs advance'),
    ('step_076_inner_slip_profile.py', 'Step 076: Inner-region slip profile 8-100 AU (REBOUND/DE440s)'),
    ('step_077_dipole_era.py', 'Step 077: Era stability of the aphelion dipole'),
    ('step_078_tno_era.py', 'Step 078: Discovery-era stability of the TNO cluster'),
    ('step_079_cross_fitter.py', 'Step 079: Cross-fitter element audit, MPC vs CODE/Warsaw'),
    ('step_080_signed_residual.py', 'Step 080: Signed residual coherence, planets subtracted (REBOUND/DE440s)'),
    ('step_081_geocentric_null.py', 'Step 081: Geocentric observation-direction null'),
    ('step_082_sky_matched.py', 'Step 082: Same-patch geocentric control'),
    ('step_083_mass_ladder.py', 'Step 083: Perturber mass ladder, required-mass bound (REBOUND/DE440s)'),
    ('step_084_newest_cohort.py', 'Step 084: Newest detached-cohort audit and Ammonite jackknife'),
    ('step_085_des_resident.py', 'Step 085: DES independent resident-cohort audit (Bernardinelli+ catalogue)'),
    ('step_086_slip_physics.py', 'Step 086: Slip-physics discriminator, angular vs boundary vs arc-epoch'),
    ('step_087_distance_ladder.py', 'Step 087: Perturber distance ladder, required-mass map (REBOUND/DE440s)'),
    ('step_088_tide_resident.py', 'Step 088: Resident-side galactic-tide insertion, maintenance + regeneration (REBOUND/DE440s)'),
    ('step_089_bipolar_morphology.py', 'Step 089: Bipolar morphology audit, mirror-cap second-lobe test'),
    ('step_090_mean_plane.py', 'Step 090: Bias-free mean-plane warp (SCT25 estimator, MPCORB lineage)'),
    ('step_091_inner_perturber.py', 'Step 091: Inner-perturber (Planet-Y) ladder vs ephemeris bound (REBOUND/DE440s)'),
    ('step_092_mpcorb_clustering.py', 'Step 092: MPCORB-lineage replication of detached-TNO clustering'),
    ('step_093_perturber_zoo.py', 'Step 093: Perturber-zoo closure, candidates vs required-mass map (REBOUND/DE440s)'),
    ('step_094_joint_axis_mp.py', 'Step 094: Six-channel joint axis-permutation test (mean-plane channel)'),

    # Phase 7: Spacecraft clock and nuclear-process channels
    ('step_095_download_spacecraft.py', 'Step 095: Spacecraft clock/plasma/power data acquisition (NAIF SCLK+SPK+LSK, PDS PPI VLISM)'),
    ('step_096_sclk_clock_audit.py', 'Step 096: Voyager SCLK clock-correction audit (resync/rate history vs crossings and radius)'),
    ('step_097_pws_plasma_channel.py', 'Step 097: Voyager PWS plasma-frequency channel (density structure vs lapse-equivalent)'),
    ('step_098_rtg_nuclear_channel.py', 'Step 098: RTG nuclear-decay channel (Pu-238 clock, Gamow-amplified constants bound)'),
    ('step_099_heliospheric_geometry.py', 'Step 099: Heliospheric geometry (probe/ISO asymptotes vs axis and caps)'),
    ('step_100_clock_consistency.py', 'Step 100: Cross-channel clock-consistency synthesis (orbital/electronic/plasma/nuclear)'),

    # Phase 8: Predictive validation & census channels
    ('step_101_iso_ng_channel.py', 'Step 101: Interstellar-object non-gravitational channel (SBDB NG solutions vs cap-transit geometry)'),
    ('step_102_slip_map_validation.py', 'Step 102: Out-of-sample validation of the bipolar slip map (5-fold CV + CODE<->Warsaw transfer)'),
    ('step_103_rtg_midlife_control.py', 'Step 103: RTG mid-life degradation control (placebo-epoch ramp scan + vintage-matched V2 reference)'),
    ('step_104_spatial_holdout.py', 'Step 104: Spatial holdout and axis localization of the slip map (azimuthal halves, band extrapolation, displacement curve)'),
    ('step_105_lpc_bidirectional.py', 'Step 105: Third-catalogue bidirectional replication (one-apparition J/A+A/571/A63) and out-of-sample slip-map transfer'),
    ('step_106_harmonic_axis_transfer.py', 'Step 106: Harmonic-order out-of-sample discrimination, third-cohort axis recovery, amplitude transfer, and cross-solution concordance'),
    ('step_107_named_directions.py', 'Step 107: Named-direction audit of the recovered axes (galactic/solar-apex/CMB/ecliptic/P9 separations + coincidence null)'),
    ('step_108_cometels_dipole.py', 'Step 108: Full-population aphelion dipole on the MPC CometEls near-parabolic census (tide-aware null, membership split, free dipole recovery)'),
    ('step_109_census_bias_audit.py', 'Step 109: Discovery-bias audit of the CometEls census dipole (ecliptic-footprint null, magnitude/era/q stability, heterogeneity test)'),
    ('step_110_aphelion_direction_vs_axis.py', 'Step 110: Directional vs axisymmetric decomposition of the aphelion asymmetry (cos, |cos|, cos2 moments under the tide-aware null, all comet samples)'),
    ('step_111_detachment_coherence.py', 'Step 111: Detachment-cut coherence map of the resident axis (varpi significance + recovered direction across the a_min x q_min plane, conditioned null)'),
    ('step_112_ossos_ensemble.py', 'Step 112: OSSOS characterized-ensemble audit (class-aware footprint baselines, selectivity profile, inter-fitter floor, plutino positive control)'),
    ('step_113_resident_ledger.py', 'Step 113: Independent-survey resident ledger (SBDB+DES+OSSOS vs own baselines, Fisher combination, deduplicated union)'),
    ('step_115_mixture_posterior.py', 'Step 115: Mixture posterior-predictive audit (each resident cohort scored vs footprint-only and TEP-mixture posteriors on its own discovery longitude)'),
    ('step_116_discovery_coupling.py', 'Step 116: Empirical discovery-coupling kernel and closed-channel audit (LOO measured varpi-lam kernel replaces assumed WN50; decisive-window split; cross-fitter instability; forward predictions)'),

    # Phase 9: Prospective lineage holdout & independent-refit synthesis
    ('step_117_prospective_lpc.py', 'Step 117: Prospective post-2017 SBDB LPC cohort through the identical bidirectional instrument (out-of-time, out-of-lineage replication)'),
    ('step_118_prospective_axis_audit.py', 'Step 118: Displaced-dipole audit of the prospective cohort (free axis recovery, global nulls, longitude harmonics, cross-cohort coverage)'),
    ('step_119_dipole_mechanism.py', 'Step 119: Mechanism identification for the post-2017 displaced dipole (observing-geometry covariates, survey lineage, quality stratification)'),
    ('step_120_pre2018_sbdb_factorial.py', 'Step 120: Era x lineage factorial -- pre-2018 SBDB LPC cohort (declared vs displaced axis recovery, fit-vintage split, CODE overlap)'),
    ('step_121_masked_signal.py', 'Step 121: Masked-signal discrimination -- joint declared-cap + private-dipole decomposition of the post-2017 cohort (Freedman-Lane added-term null, dipole-placement scan, clean strata, CODE/pre-2018 controls)'),
    ('step_122_pre2018_localization.py', 'Step 122: Anomaly localization on the pre-2018 record -- CODE-membership lock vs broad-population readings, era decomposition, fat-tail audit'),
    ('step_123_leg_coherence.py', 'Step 123: Leg-coherence decomposition -- is the overlap anomaly a single-leg reconstruction floor or an azimuthally-twisted inter-leg mismatch; noise-floor audit, cross-fitter concordance, class-constraint gradient'),
    ('step_124_lineage_aphelion.py', 'Step 124: Lineage-independent aphelion dipole -- reconstruction-free spatial channel on the post-2017 and pre-2018 SBDB cohorts (tide-aware + ecliptic nulls, quality strata, displaced-axis control)'),
    ('step_125_displaced_temporal.py', 'Step 125: Displaced-structure temporal characterization -- raw-field free-scan extremeness, era x cap factorial, transit-epoch gradient with quality/energy controls, geometry-matched pairs, named-direction audit incl. CMB rest frame'),
    ('step_126_global_cross_survey_synthesis.py', 'Step 126: Ten-channel global cross-survey and multi-lineage synthesis -- S_10 omnibus stat, 20k permutation null, pairwise angular concordance matrix, axis-localization bootstrap'),
    ('step_127_nonwarsaw_refit.py', 'Step 127: Non-Warsaw two-leg refit -- independent LM differential-correction fits on raw MPC astrometry propagated through the same boundary machinery; carrier-degeneracy close (drot concordance, d_in carrier replication, cross-arc audit)'),
    ('step_128_post2017_refit.py', 'Step 128: Post-2017 two-leg refit masking discriminator -- independent LM refit of the prospective cohort on raw MPC astrometry; leg-fit vs single-fit vs catalogue carrier constructions at declared and displaced axes'),
    ('step_129_bipolar_unification.py', 'Step 129: Frame-anchored bipolar unification test -- per-era CMB-frame dipoles, joint opposite-polarity statistic, direction shuffle and bipolar-cap audit on the independent two-leg refit record'),
    ('step_130_epoch_axis_trajectory.py', 'Step 130: Epoch-resolved axis trajectory on the pooled independent dual-leg record -- per-bin contrasts, arc-trajectory scan, drift-vs-switch model comparison, in-cap epoch interactions'),
    ('step_131_cone_shell_coherence.py', 'Step 131: Bipolar cone-shell coherence -- per-era cone-angle profiles, common-cone joint test on the measured axes, ring-vs-pole morphology, and band-resident contrasts on the independent refit record'),
    ('step_132_bipolar_meridian.py', 'Step 132: Bipolar meridian and polarity-mirror coherence -- equator-mirror pair of the era axes, azimuthal coplanarity of the measured axes about the CMB bipolar axis, and meridian-resident residual cells on the independent refit record'),
    ('step_133_harmonic_bipolar.py', 'Step 133: Bipolar harmonic decomposition and cross-population geometry -- per-era Legendre power spectrum, free-dipole mirror-pair test, ISO meridian ledger, resident radial-gradient control'),
    ('step_134_lpc_signature_refit.py', 'Step 134: Signature-cohort independent refit -- non-Warsaw two-leg LM refit of the LPC 1902-1950 signature-bearing sample (step_127 machinery); carrier-degeneracy close on the objects that carry the anomaly'),
    ('step_135_cross_channel.py', 'Step 135: Cross-channel axis convergence -- joint bipolar-frame ledger of independently recovered directions across TEP channels (GNSS-II clock axis, MGEX projection, LLR Planck rank), out-of-sample meridian coplanarity, Fisher-combined frame evidence, and clock-epoch coincidence at the comet switch boundary'),
    ('step_136_solar_wind_screening.py', 'Step 136: Solar-wind screening-epoch test -- OMNI daily solar-wind dynamic pressure matched to per-comet perihelion epochs; pooled and within-era residual-pressure coupling, era pressure contrast, and regime-split dipole polarity'),
    ('step_137_rtg_ensemble_regression.py', 'Step 137: RTG cross-generator degradation regression -- all unit-level P/P0 series across six RTG families; pooled 9-yr wander null, Voyager crossing-window residuals, sibling-coherence thermal control, Gamow bound translation'),
    ('step_138_oneway_downlink_channel.py', 'Step 138: One-way downlink channel audit -- non-cancelling transport-class mapping to the GNSS clock channel, deep-space archive inventory (Voyager delta-VLBI, Cassini VLBA, PRIDE VEX/MEX), boundary-distance coverage, sensitivity ordering'),
    ('step_139_frozen_geometry.py', 'Step 139: Frozen-geometry predictive battery -- post-hoc cone/meridian characterization converted to scored predictions: parameter-free mirror prediction, frozen meridian plane and cone band evaluated on held-out directions (displaced, resident, ISO legs, GNSS-II, MGEX, LLR) with provenance partition'),
    ('step_140_astcat_stratification.py', 'Step 140: Astrometric-catalogue stratification -- per-observation astcat leg-composition profiles on the dual-leg refit record; uniform-Gaia and uniform-legacy stratum controls on the declared/displaced carriers, composition-vs-signal association, matched-composition audit of the mixed stratum'),
    ('step_141_sclk_signed_drift.py', 'Step 141: SCLK signed-drift crossing channel -- signed oscillator fractional-offset level steps at the four boundary crossings, raw and smooth-ageing detrended, changepoint extremeness, joint sign-consistency and Fisher tails, cadence/segmentation audit, NH control note'),
    ('step_142_gaia_sink.py', 'Step 142: Gaia covariance-sink audit -- cross-arc leg misfit, leg-disagreement, NG-solution and element-sigma channels on the uniform-Gaia vs sub-threshold strata; seed-lock and time-domain-carrier tests'),
    ('step_143_comet_ng_channel.py', 'Step 143: Comet non-gravitational acceleration channel -- Warsaw tableb4 A1/A2/A3/tau vectors and tablea1 NG incidence read as a directional observable against the declared cap; transverse-term leakage test with radial/normal controls and n_obs confound audit'),
    ('step_144_giant_planet_apsidal.py', 'Step 144: Giant-planet secular apsidal channel -- 5 Myr REBOUND/DE440s integration of the four giants; forced apsidal azimuths, joint axis-proximity test, dwell fractions, present-epoch 3-D geometry; secular context of the registered Neptune mean-element coincidence'),
    ('step_145_transfer_function.py', 'Step 145: Impulse-versus-holonomy transfer-function discriminator -- crossing-speed scaling and energy-channel morphology of the slip on the independent dual-leg refit records; selects the microscopic realization of the equivalent-time observable'),
    ('step_146_sclk_radial_gradient.py', 'Step 146: SCLK radial lapse-gradient and bipolar-sign channel -- signed dnu vs radius per craft, outer-region drift sign split vs trajectory projection on the CMB and comet bipolar axes, radius-vs-ageing discrimination, NH control'),
    ('step_147_ephemeris_audit.py', 'Step 147: Planetary-ephemeris invariance audit -- bidirectional boundary integration of the class-1 CODE sample under DE440s and DE430; per-leg boundary periapsis displacements and the kick-and-geometry-regressed in-cap gap compared across ephemeris generations'),
    ('step_148_pws_anchored_sclk.py', 'Step 148: PWS-anchored SCLK epoch-lock test -- strongest negative signed step within each crossing window vs canonical epoch, random-target null, joint Fisher, and V1 PWS plasma-step coincidence'),
    ('step_149_newtonian_directional_null.py', 'Step 149: Newtonian directional null -- the complete directional battery (cap contrast, free sky scan, look-elsewhere, shell morphology, permutation and axis-coincidence nulls) applied to the pure N-body rotation field, the catalogued field, the observed-minus-model residual, the energy kick and rotation-per-kick'),
    ('step_150_channel_dissociation.py', 'Step 150: Channel-dissociation sensitivity audit -- per-comet Monte-Carlo element perturbations from SBDB sigmas propagated through the boundary machinery; rotation-channel vs aphelion-channel noise sensitivity, element-group attribution and quality-gradient comparison to the observed stratification'),
    ('step_151_injection_transfer.py', 'Step 151: Slip-injection transfer function -- synthesized astrometry on the real observing chain (epochs, stations, noise) with injected slip realizations (integrable control, position/velocity translations, transverse kick, within-apparition holonomy) refit by the standard LM/DE440s machinery; measures which realizations survive standard orbit fitting into the record'),
    ('step_152_arrival_latitude_null.py', 'Step 152: Arrival-anisotropy conditioned nulls -- ecliptic- and Galactic-latitude-conditioned longitude shuffles and axis-longitude specificity scan for the in-cap arrival excess, replacing the uniform-sky binomial null of step 149 T6'),
    ('step_153_zonal_distortion.py', 'Step 153: Zonal-distortion null -- legacy star-catalogue warp classes (zonal bands, regional tiles, rigid frame rotation) injected into real MPC astrometry, refit by the identical LM/DE440s machinery and propagated to the boundary sphere; tests whether a catalogue error field can reproduce the anomaly simultaneously in element selectivity, energy flatness, axis organization and amplitude'),
    ('step_154_sclk_plasma_geometry.py', 'Step 154: SCLK link-plasma and solar-activity regression -- per-boundary phase residuals and drift vs interval-matched OMNI P_dyn/|B|/Vsw/Tp/SSN/F10.7/ap and Sun-Earth-probe elongation; partial radius-vs-plasma correlations and heliosheath activity stratification'),

    # Phase 10: Forward-model closure -- generative source models,
    # field-level derivations, and injected-realization transfer tests
    ('step_155_iso_injection.py', 'Step 155: 1I/Oumuamua slip-injection -> NG-term transfer test -- synthetic astrometry on the real MPC observing chain with boundary-slip and sustained lapse-rate realizations, refit by gravity-only and Marsden-NG-augmented LM/DE440s machinery; measures the fitted A1 transfer function and validates against the SBDB 7c solution'),
    ('step_156_boundary_radius_scan.py', 'Step 156: Boundary-radius forward-model scan -- position-slip holonomy injected at controlled inbound crossing radii {8,25,50,100,150,250} AU through the real MPC observing chain, refit by the standard LM/DE440s machinery, with shell-resolved fitted-vs-true reconstruction testing whether a distant crossing produces a fitted offset already developed at 8 AU'),
    ('step_157_generative_population.py', 'Step 157: Conventional generative population model -- isotropic Oort-spike source, secular Galactic-tide injection (orbit-integrated torque), stellar-impulse remixing and empirical ecliptic-latitude discovery selection, tested against the observed cap excess and arrival-resultant alignment'),
    ('step_158_metric_winding.py', 'Step 158: Metric-derived TEP slip -- winding-holonomy field phi = m*theta_A*sigma(r) about the boundary axis; path-integrated proper-time offsets computed on the true dual-leg trajectories and tested against the measured equivalent-time offsets, cap partition, velocity flatness, control axes and boundary-radius sensitivity'),
    ('step_159_tp_slip.py', 'Step 159: Perihelion-time leg-slip channel -- osculating perihelion-epoch disagreement between independently fitted inbound/outbound legs on the step_b91 dual-leg record; validity concordance vs catalogue dtau/ddirf, cap contrast, bipolar polarity, noise dissociation'),
    ('step_160_iso_legrefit.py', 'Step 160: Independent two-leg refit of the interstellar objects -- raw-MPC leg-split LM/DE440s fits of 1I/2I/3I with per-leg boundary asymptotes, cap classification, inter-leg non-closure, perihelion-epoch slip, SBDB NG-ledger confound context and post-2017 cohort percentiles'),
    ('step_161_apdb_planet.py', 'Step 161: APDB planetary astrometric O-C channel -- O\'Handley-format Pluto/Uranus/Neptune optical record (1914-1998) reduced against DE440s with FK4/FK5/TETE equinox handling, light-time and topocentric correction; per-source sanity floor, direction organisation, within-source robustness, normal-point ledger'),
    ('step_162_nima_pilot.py', 'Step 162: NIMA detached-resident O-C pilot -- Lucky Star/NIMA per-observation residuals for the resident-ledger intersection; coverage, per-object residual level and secular drift vs cap membership; small-N pilot, not a detection channel'),
    ('step_163_footprint_forward.py', 'Step 163: Forward-modelled comet discovery footprint -- isotropic reorientation of each real orbit propagated through the bright window with elongation/declination/geocentric cuts and brightness weighting; tests whether physical discovery selection can produce the observed aphelion dipole direction, with selection-only null, axis+footprint decomposition, scenario robustness and longitude-marginal diagnostic'),
    ('step_164_tep_generative.py', 'Step 164: TEP-side generative model for the arrival anisotropy -- symmetric counterpart of the conventional generative test; quantifies the slip-redirection bound (measured fitted-aphelion displacement vs coherent requirement), solves the required source-level modulation amplitude on the isotropic and tidal baselines, and discriminates lobe vs bipolar source morphology via the anti-cap fraction across all five cohorts'),
    ('step_165_disformal_slip.py', 'Step 165: Canonical disformal transport -- admissible-branch benchmark for the comet boundary slip; integrates the synchronisation connection dsigma = -(B/A^2)(u.grad phi) P grad phi (Jakarta A3.2) along comet trajectories through a boundary-layer field excursion, evaluates the B0 > 0 causal branch explicitly, separates the unsigned equivalent-time amplitude from signed open-path transport, solves the excursion depth u_b required at the benchmark scale, and tests antipodal double-crossing bipolarity, per-comet correlation, tracking baseline and channel decomposition'),
    # Phase 11: Publication validation -- claims-trace audit.
    # Depends on every upstream result.
    ('step_114_claims_trace.py', 'Step 114: Manuscript claims-trace audit (every cited step/result/figure/number traced to results/)'),

]


class PipelineLogger:
    """Handles consistent logging for pipeline execution."""

    def __init__(self, append: bool = False):
        self.start_time = datetime.now(timezone.utc)
        self.log_file = LOGS_DIR / "pipeline.log"
        self.step_results: List[Dict] = []
        # resume runs append so the aggregate log keeps the full-session
        # audit trail rather than truncating the earlier positions
        self._append = append

    def _write(self, message: str, level: str = "INFO"):
        timestamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        log_line = f"[{timestamp}] [{level:8}] {message}"
        print(message)
        mode = 'a' if (hasattr(self, '_log_initialized') or self._append) else 'w'
        with open(self.log_file, mode, encoding='utf-8') as f:
            f.write(log_line + '\n')
        self._log_initialized = True

    def header(self, title: str, width: int = 70):
        self._write("=" * width)
        self._write(title)
        self._write("=" * width)

    def subheader(self, title: str, width: int = 70):
        self._write("-" * width)
        self._write(title)
        self._write("-" * width)

    def info(self, message: str):
        self._write(message, "INFO")

    def success(self, message: str):
        self._write(f"[OK] {message}", "SUCCESS")

    def error(self, message: str):
        self._write(f"[FAIL] {message}", "ERROR")

    def warning(self, message: str):
        self._write(f"[WARN] {message}", "WARNING")

    def step_start(self, step_num: int, total_steps: int, description: str):
        progress = f"[{step_num}/{total_steps}]"
        self._write("")
        self.header(f"{progress} {description}")
        return time.time()

    def step_complete(self, description: str, start_time: float, status: str = "SUCCESS"):
        duration = time.time() - start_time
        duration_str = f"{duration:.2f}s" if duration < 60 else f"{duration/60:.2f}m"
        if status == "SUCCESS":
            self.success(f"{description} completed in {duration_str}")
        else:
            self.error(f"{description} failed after {duration_str}")
        return duration

    def add_step_result(self, step_name: str, description: str,
                        status: str, duration: float, exit_code: int = 0):
        self.step_results.append({
            'step': step_name,
            'description': description,
            'status': status,
            'duration_seconds': duration,
            'exit_code': exit_code,
            'timestamp': datetime.now(timezone.utc).isoformat()
        })

    def generate_research_audit(self):
        """Generate a research audit for reproducibility."""
        import platform
        import hashlib

        def get_file_hash(path):
            if not path.exists():
                return "MISSING"
            h = hashlib.sha256()
            with open(path, "rb") as f:
                for block in iter(lambda: f.read(4096), b""):
                    h.update(block)
            return h.hexdigest()

        telemetry = {
            "os": f"{platform.system()} {platform.release()}",
            "python_version": platform.python_version(),
            "processor": platform.processor(),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "project_root": str(PROJECT_ROOT)
        }

        integrity = {fn: get_file_hash(STEPS_DIR / fn)
                     for fn, _ in CORE_STEPS}

        audit_data = {
            "audit_type": "TEP-9 Research-Grade Audit",
            "project": "TEP-9 (outer solar system temporal boundary)",
            "telemetry": telemetry,
            "pipeline_results": self.step_results,
            "script_integrity": integrity,
            "shared_source_integrity": {str(p.relative_to(PROJECT_ROOT)): get_file_hash(p)
                for folder in ("scripts/utils", "core", "config")
                for p in (PROJECT_ROOT / folder).rglob("*")
                if p.is_file() and p.suffix in {".py", ".json", ".yaml"}}
        }

        audit_dir = PROJECT_ROOT / "results" / "audits"
        audit_dir.mkdir(parents=True, exist_ok=True)
        audit_path = audit_dir / f"RESEARCH_AUDIT_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(audit_path, "w", encoding="utf-8") as f:
            json.dump(audit_data, f, indent=4)
        self.info(f"Research audit generated: {audit_path.name}")
        return audit_path

    def final_summary(self, total_steps: int) -> bool:
        total_duration = (datetime.now(timezone.utc) - self.start_time).total_seconds()
        success_count = sum(1 for r in self.step_results if r['status'] == 'SUCCESS')
        fail_count = sum(1 for r in self.step_results if r['status'] == 'FAILED')

        self._write("")
        self.header("PIPELINE EXECUTION SUMMARY")
        self.subheader("Execution Timeline")
        self.info(f"Pipeline started: {self.start_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
        self.info(f"Pipeline ended:   {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}")
        self.info(f"Total duration:   {total_duration/60:.2f} minutes")

        self.subheader("Step-by-Step Results")
        self.info(f"{'Step':<10} {'Description':<58} {'Status':<10} {'Duration':<12}")
        self.info("-" * 90)
        for result in self.step_results:
            symbol = "OK" if result['status'] == 'SUCCESS' else "X"
            desc = result['description'][:56]
            self.info(f"{result['step']:<10} {desc:<58} {symbol} {result['status']:<8} {result['duration_seconds']:.2f}s")

        self.subheader("Statistics")
        attempted = len(self.step_results)
        self.info(f"Total steps:     {total_steps}")
        self.info(f"Attempted:       {attempted}")
        self.info(f"Successful:      {success_count}")
        self.info(f"Failed:          {fail_count}")
        if attempted:
            self.info(f"Success rate:    {100*success_count/attempted:.1f}%")

        self.subheader("Research Accountability & Reproducibility")
        self.generate_research_audit()

        if fail_count == 0 and attempted == total_steps:
            self.subheader("Output Locations")
            self.info(f"Log file:         {self.log_file}")
            self.info(f"Results:          {PROJECT_ROOT / 'results'}")
            self.info(f"Figures:          {PROJECT_ROOT / 'results' / 'figures'}")
            self._write("")
            self.success("PIPELINE COMPLETED SUCCESSFULLY")
            return True
        self._write("")
        if fail_count:
            self.error(f"PIPELINE FAILED - {fail_count} step(s) did not complete")
        else:
            self.success(
                f"PARTIAL RUN COMPLETED - {success_count}/{total_steps} "
                f"positions ran, none failed")
        return fail_count == 0


def step_outputs(step_path: Path) -> List[Path]:
    """Result JSONs a step declares it writes (RESULTS / "name.json"
    or "results/name.json" literal paths)."""
    import re
    src = step_path.read_text(errors='ignore')
    names = set(re.findall(
        r'["\'](?:results[/\\])?(step_[A-Za-z0-9_]+\.json)["\']', src))
    return [PROJECT_ROOT / 'results' / n for n in sorted(names)]


def step_is_fresh(filename: str, logger: PipelineLogger) -> bool:
    """True when every declared output exists, is newer than its declared
    inputs, and is newer than the producing script itself."""
    step_path = STEPS_DIR / filename
    outs = [p for p in step_outputs(step_path) if p.exists()]
    if not outs:
        return False
    script_m = step_path.stat().st_mtime
    for out in outs:
        try:
            rec = json.loads(out.read_text())
        except Exception:
            return False
        out_m = out.stat().st_mtime
        if out_m < script_m:
            return False
        for rel in (rec.get('inputs') or []):
            if not isinstance(rel, str) or '*' in rel or '(' in rel:
                continue  # annotated/globbed input description, not a path
            inp = PROJECT_ROOT / rel
            if not inp.exists() or inp.stat().st_mtime > out_m:
                return False
    logger.info(f"skip {filename}: outputs fresh")
    return True


def run_step(filename: str, description: str, step_num: int, total_steps: int,
             logger: PipelineLogger) -> bool:
    step_path = STEPS_DIR / filename
    if not step_path.exists():
        logger.error(f"File not found: {step_path}")
        logger.add_step_result(f"Step {step_num:03d}", description, "FAILED", 0.0, -1)
        return False

    start_time = logger.step_start(step_num, total_steps, description)
    try:
        result = subprocess.run(
            [sys.executable, str(step_path)],
            cwd=str(PROJECT_ROOT),
            capture_output=False,
            text=True
        )
        duration = logger.step_complete(
            description, start_time,
            "SUCCESS" if result.returncode == 0 else "FAILED")
        logger.add_step_result(
            f"Step {step_num:03d}", description,
            "SUCCESS" if result.returncode == 0 else "FAILED",
            duration, result.returncode)
        return result.returncode == 0
    except Exception as e:
        duration = time.time() - start_time
        logger.error(f"Exception running {description}: {e}")
        logger.add_step_result(f"Step {step_num:03d}", description, "FAILED", duration, -1)
        return False


def main():
    # Optional positional slicing: ``--from N`` (1-based first position)
    # and ``--to M`` (inclusive last position) run a contiguous subset,
    # e.g. resuming after an acquisition failure without re-downloading.
    first, last = 1, len(CORE_STEPS)
    argv = sys.argv[1:]
    if "--from" in argv:
        i = argv.index("--from")
        first = max(1, int(argv[i + 1]))
    if "--to" in argv:
        i = argv.index("--to")
        last = min(len(CORE_STEPS), int(argv[i + 1]))
    skip_fresh = "--skip-fresh" in argv

    lock_fd = acquire_pipeline_lock()
    if lock_fd is None:
        holder = LOCK_FILE.read_text(errors="ignore").strip()
        print(f"[FAIL] another pipeline run holds {LOCK_FILE} "
              f"({holder}); refusing to start a concurrent run")
        return 1

    logger = PipelineLogger(append=(first > 1))
    title = ("TEP-9 PIPELINE - FULL EXECUTION" if (first, last) == (1, len(CORE_STEPS))
             else f"TEP-9 PIPELINE - RESUME (positions {first}-{last})")
    logger.header(title, 80)
    logger.info(f"Project root: {PROJECT_ROOT}")
    logger.info(f"Steps dir:     {STEPS_DIR}")
    logger.info(f"Log file:      {logger.log_file}")

    total_steps = len(CORE_STEPS)
    logger.subheader("EXECUTING PIPELINE STEPS", 80)

    for i, (filename, description) in enumerate(CORE_STEPS, 1):
        if not (first <= i <= last):
            continue
        if skip_fresh and step_is_fresh(filename, logger):
            continue
        if not run_step(filename, description, i, total_steps, logger):
            logger.error("Stopping pipeline due to failure")
            break

    success = logger.final_summary(total_steps)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
