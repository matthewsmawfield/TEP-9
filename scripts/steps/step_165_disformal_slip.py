"""step_165: canonical disformal transport -- amplitude closure for
the comet boundary slip (b131).

Step 158 (b124) established that a single-valued scalar winding
realisation cannot produce the slip: d(phi) is exact, the physical
endpoint term is leg-symmetric, and the tested open-path statistic
was cap-blind.  The canonical theory paper (0-TEP, Jakarta revision,
Sec. 3.3 / 6 / A3.2) supplies the structure that step was missing:
the invariant observable is the synchronisation connection

    dsigma_mu = -(B/A^2)(u . grad phi) P_mu^nu grad_nu phi ,

whose scalar (conformal) part is exact and whose leading NON-exact
part is the projected disformal term.  A comet is a transported
clock: along its worldline u.grad phi = du/dt, so the open-path
desynchronisation relative to the barycentric network is

    dt_slip = (1/c) int dsigma_i dx^i
            = -(1/c^2) int b(u) R_H^2 / (A^2 N^2) (du/dt)^2 dt

with the corpus field-space envelope (excursion convention: u is the
field measured from the contemporary ambient; the absolute-u convention
is macroscopically excluded in the Paper-0 triangle calculation)

    b(u) = B0 * u^2/(1+u^2) * exp(-u^4 / 2 sigma_B^4),
    B0 = +3.2e-3   (admissible-branch benchmark).  The magnitude is
    mirrored from the Paper-0 volume-balance reconstruction only to
    set an instrumental scale; that reconstruction found B0 < 0 and
    is excluded by the Paper-0 null-cone condition.  It is therefore
    not a calibration of the admissible normalization.

The primary step-065 datum is an equivalent-time amplitude obtained
from an angular separation, delta t_eq = |delta theta| r_b^2/h.  It is
unsigned.  Consequently the fit below must use |delta t_slip|.  The
sign of the open-path synchronization transport is retained separately
and is never inferred from that amplitude statistic.

Properties that follow from the form alone (all measured
independently in this pipeline):
  * sign-definite in B: one sign of slip, no tuned sign;
  * quadratic in the traversal rate -> the slip concentrates where
    |grad u| is concentrated: a boundary transition gives a discrete
    crossing offset (step 071 shell flatness, step 068 transit-time
    independence);
  * the comet crosses the boundary twice, at antipodal points, so an
    axis-anchored field G(th) contributes G(th)^4 + G(180-th)^4 --
    every odd harmonic cancels identically and the open-path slip is
    automatically bipolar (cos 2th dominated) for ANY axis-anchored
    field -- the morphology the residual slip field actually carries
    (step 086: cos2th rho=+0.17, p=0.010; cos th flat at p=0.30),
    while the resident-cloud source response is linear and stays
    dipolar (step 164);
  * the theta-uniform part of the lag is absorbed into the fitted
    elements (parameter absorption, step 151) -- the observable is
    the theta-dependent contrast, which is what the harmonic map
    measures.

The open question tested quantitatively: at a positive-B admissible
benchmark of the same magnitude as the excluded reconstruction scale,
what boundary field excursion u_b does the measured amplitude contrast
require?  The answer is a scale relation, not a normalization
prediction, and is compared against the corpus anchors
(Newtonian-tracking 4e-11, Galactic ambient ~7e-7, solar surface ~2e-6,
cosmological u ~ 0.75).

Channels:
  M (motion/quadratic): u.grad phi = du/dt along the comet
      worldline -> per-crossing lag ~ (|B0| R_H^2/c^2) v_r *
      int env_b(u) (du/dr)^2 dr; bipolar response via the
      antipodal sum (exact for every axis-anchored G).
  R (network roll): u.grad phi = -Pi W (partitioned roll of the
      static clock congruence, Pi = H0) -> net transport
      ~ (|B0| R_H^2/c^2) Pi * [U(G_in) - U(G_out)], signed by the
      antipodal field difference; dipolar response, subdominant.

Inputs
  results/step_b30_proper_time_slip.csv   (per-comet residual slip,
      dtau_unexplained, and axis separation theta; 229 pooled)
  results/step_b124_metric_winding.csv    (per-comet boundary
      crossing speed v_cross)

Outputs
  results/step_b131_disformal_slip.json
  results/step_b131_disformal_slip.csv
  results/figures/supplementary/step_b131_disformal_slip.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout
logger = StepLogger("step_165_disformal_slip")
tee_stdout(logger)
logger.header("Canonical disformal transport -- comet-slip "
              "amplitude closure")

import csv
import json
import math
import numpy as np
from scipy.optimize import brentq
from scipy.stats import spearmanr

SEED = 20261103
rng = np.random.default_rng(SEED)

# ------------------------------------------------------------------
# constants -- admissible branch and benchmark scale
# ------------------------------------------------------------------
C = 2.998e8                       # m/s
H0 = 2.27e-18                     # s^-1 (partitioned roll Pi_bar)
R_H = C / H0                      # ~1.32e26 m
B0 = +3.2e-3                      # admissible B >= 0 benchmark
B0_EXCLUDED = -3.2e-3             # Paper-0 volume-balance diagnostic
SIG_B = 1.0
AU_M = 1.496e11
YR_S = 3.156e7
R_B = 250.0 * AU_M
CAP = 60.0

U_TRACK = 1.32712440018e20 / (C**2 * R_B)   # GM_sun/(c^2 r_b) ~4.4e-11
U_GAL = 7.0e-7                              # Galactic ambient
U_SUN_SURF = 2.1e-6                         # solar surface
U_COSMO = 0.75                              # cosmological ambient


def env_b(u):
    """b(u)/B0 = u^2/(1+u^2) exp(-u^4/2 s_B^4) -- array-safe."""
    u = np.asarray(u, dtype=float)
    return u * u / (1.0 + u * u) * np.exp(-u**4 / (2 * SIG_B**4))


# ------------------------------------------------------------------
# cohort: axis separation + crossing speed + measured residual slip
# ------------------------------------------------------------------
slip_rows = {r["desig"]: r for r in csv.DictReader(
    open(RESULTS / "step_b30_proper_time_slip.csv"))}
wind_rows = {r["desig"]: r for r in csv.DictReader(
    open(RESULTS / "step_b124_metric_winding.csv"))}

comets = []
for d, r in slip_rows.items():
    if r.get("dtau_unexplained") in ("", None) or r.get("theta") in (
            "", None):
        continue
    wr = wind_rows.get(d, {})
    vc = float(wr["v_cross"]) if wr.get("v_cross") not in (
        "", None) else 0.88
    comets.append(dict(desig=d, theta=float(r["theta"]),
                       dtau=float(r["dtau_unexplained"]),
                       daa=abs(float(r["daa_sim"]))
                       if r.get("daa_sim") not in ("", None) else np.nan,
                       v_cross=vc, cohort=r["cohort"]))
logger.info(f"cohort: {len(comets)} comets with residual slip, "
            f"theta, crossing speed")

TH = np.array([c["theta"] for c in comets])
DT_OBS = np.array([c["dtau"] for c in comets])          # yr
VR = np.array([c["v_cross"] for c in comets]) * AU_M / YR_S  # m/s

in_pole = (TH < CAP) | (TH > 180.0 - CAP)
midband = ~in_pole
meas_contrast = float(np.median(DT_OBS[in_pole])
                      - np.median(DT_OBS[midband]))
meas_pole = float(np.median(DT_OBS[in_pole]))
meas_mid = float(np.median(DT_OBS[midband]))
logger.info(f"measured residual: pole {meas_pole:.2f} yr "
            f"(n={in_pole.sum()}), midband {meas_mid:.2f} yr "
            f"(n={midband.sum()}), contrast {meas_contrast:.2f} yr")

# ------------------------------------------------------------------
# transport integrals -- vectorised over the per-comet G values
# ------------------------------------------------------------------
# Excursion convention: u(r,th) = u_b G(th) sigma(r) is the field
# measured from ambient; sigma logistic, 1 inside, 0 outside.
# Inbound crossing at theta_in = theta; outbound at the antipode
# theta_out = 180 - theta.  Crossing layer s = r - r_b.
N_S = 8001


def _layer(ell_m):
    s = np.linspace(-6.0 * ell_m, 6.0 * ell_m, N_S)
    sig = 1.0 / (1.0 + np.exp(s / ell_m))
    dsig = sig * (1.0 - sig) / ell_m
    return s, sig, dsig


def K_matrix(Gvals, u_b, ell_m):
    """int env_b(u) (du/dr)^2 dr [m^-1] for each G in Gvals."""
    s, sig, dsig = _layer(ell_m)
    u = u_b * np.outer(sig, Gvals)                  # (N_S, n_G)
    integrand = env_b(u) * (u_b * Gvals[None, :]
                            * dsig[:, None]) ** 2
    return np.trapz(integrand, s, axis=0)


def U_matrix(Gvals, u_b, n=2001):
    """int_0^{u_b G} env_b(u) du  (roll channel)."""
    uu = np.linspace(0.0, 1.0, n)[None, :] * (u_b * Gvals)[:, None]
    return np.trapz(env_b(uu), uu, axis=1)


def slip_profiles(Gfun, u_b, ell_m):
    """Signed per-comet transport [yr], channels M and R.

    The leading connection carries an explicit minus sign, so the
    transport orientation reverses with B0.  Observable amplitude
    comparisons must be made outside this function.
    """
    Gin = Gfun(TH)
    Gout = Gfun(180.0 - TH)
    Kin = K_matrix(Gin, u_b, ell_m)
    Kout = K_matrix(Gout, u_b, ell_m)
    dtM = (-B0 * R_H**2 / C**2) * VR * (Kin + Kout)
    dtR = (-B0 * R_H**2 / C**2) * H0 * (
        U_matrix(Gin, u_b) - U_matrix(Gout, u_b))
    return dtM / YR_S, dtR / YR_S


def slip_amplitude(Gfun, u_b, ell_m):
    """Unsigned observable corresponding to the step-065 conversion."""
    dtM, dtR = slip_profiles(Gfun, u_b, ell_m)
    return np.abs(dtM), np.abs(dtR)


def harmonics_uniform(y, th_deg):
    """cos th / cos 2th content on a uniform theta grid."""
    c1 = np.cos(np.radians(th_deg))
    c2 = np.cos(2 * np.radians(th_deg))
    return dict(cos1=float(np.dot(y, c1) / np.dot(c1, c1)),
                cos2=float(np.dot(y, c2) / np.dot(c2, c2)))


def P2(x):
    return 0.5 * (3.0 * x * x - 1.0)


MORPHS = {
    "iso":     lambda t: np.ones_like(np.asarray(t, dtype=float)),
    "dip_a1":  lambda t: 1.0 + np.cos(np.radians(t)),
    "dip_a2":  lambda t: 1.0 + 2.0 * np.cos(np.radians(t)),
    "quad_P2": lambda t: 1.0 + P2(np.cos(np.radians(t))),
}

TH_GRID = np.linspace(0.0, 180.0, 721)   # uniform grid for
                                         # morphology harmonics

ELL_LIST = [4.0, 8.0, 16.0, 32.0, 64.0]
U_B_GRID = np.logspace(-8, -1, 60)

results = {}
csv_rows = []
for mname, Gf in MORPHS.items():
    for ell_au in ELL_LIST:
        ell_m = ell_au * AU_M
        tag = f"{mname}|ell{ell_au:g}"
        grid_c = np.empty(len(U_B_GRID))
        for j, ub in enumerate(U_B_GRID):
            pM, _ = slip_amplitude(Gf, ub, ell_m)
            grid_c[j] = (np.median(pM[in_pole])
                         - np.median(pM[midband]))
        ub_req = np.nan
        if meas_contrast <= grid_c[-1]:
            ub_req = float(np.interp(meas_contrast, grid_c,
                                     U_B_GRID))
        ub_req_out = None if np.isnan(ub_req) else ub_req
        # morphology harmonics on the uniform grid
        Gin = Gf(TH_GRID)
        Gout = Gf(180.0 - TH_GRID)
        prof = K_matrix(Gin, 1e-4, ell_m) + K_matrix(Gout, 1e-4,
                                                   ell_m)
        h = harmonics_uniform(prof, TH_GRID)
        rec = dict(tag=tag, morph=mname, ell_au=ell_au,
                   cos1_uniform=h["cos1"], cos2_uniform=h["cos2"],
                   ub_for_contrast=ub_req_out,
                   contrast_max_yr=float(grid_c[-1]))
        results[tag] = rec
        csv_rows.append(rec)
        logger.info(f"  {tag}: uniform cos1={h['cos1']:.2e} "
                    f"cos2={h['cos2']:.2e} | u_b(contrast) = "
                    f"{ub_req:.2e}")

# ------------------------------------------------------------------
# headline: dipolar field, ell = 8 AU (the step_158 envelope width)
# ------------------------------------------------------------------
Gf = MORPHS["dip_a1"]
ell_m = 8.0 * AU_M
grid = np.logspace(-7, -1, 120)
gc = np.empty(len(grid))
for j, ub in enumerate(grid):
    pM, _ = slip_amplitude(Gf, ub, ell_m)
    gc[j] = np.median(pM[in_pole]) - np.median(pM[midband])


def _headline_contrast(log_ub):
    pM, _ = slip_amplitude(Gf, math.exp(log_ub), ell_m)
    return (float(np.median(pM[in_pole]) - np.median(pM[midband]))
            - meas_contrast)


ub_head = float(math.exp(brentq(
    _headline_contrast, math.log(grid[0]), math.log(grid[-1]),
    xtol=1e-12, rtol=1e-12)))
prof_M_signed, prof_R_signed = slip_profiles(Gf, ub_head, ell_m)
prof_M, prof_R = np.abs(prof_M_signed), np.abs(prof_R_signed)
rho_pred = spearmanr(prof_M, DT_OBS)
rho_in = spearmanr(prof_M[TH < CAP], DT_OBS[TH < CAP])
null_rho = np.array([spearmanr(prof_M[rng.permutation(len(prof_M))],
                               DT_OBS).statistic
                     for _ in range(2000)])
p_rho = float((np.abs(null_rho) >= abs(rho_pred.statistic)).mean())
# uniform-grid morphology for the headline
prof_grid = K_matrix(Gf(TH_GRID), ub_head, ell_m) + K_matrix(
    Gf(180.0 - TH_GRID), ub_head, ell_m)
h_head = harmonics_uniform(prof_grid, TH_GRID)

headline = dict(
    morph="dip_a1: field 1 + cos th; antipodal double crossing",
    ell_au=8.0, u_b_required=ub_head,
    harmonics_uniform=dict(cos1=h_head["cos1"], cos2=h_head["cos2"]),
    rho_pred_vs_resid=float(rho_pred.statistic),
    p_pred_vs_resid=float(p_rho),
    rho_in_cap=float(rho_in.statistic),
    p_in_cap=float(rho_in.pvalue),
    pred_pole_med=float(np.median(prof_M[in_pole])),
    pred_mid_med=float(np.median(prof_M[midband])),
    pred_mean_absorbed=float(prof_M.mean()),
    signed_transport_pole_med=float(np.median(prof_M_signed[in_pole])),
    signed_transport_mid_med=float(np.median(prof_M_signed[midband])),
    signed_transport_contrast=float(
        np.median(prof_M_signed[in_pole])
        - np.median(prof_M_signed[midband])),
    roll_pole_med=float(np.median(prof_R[in_pole])),
    roll_mid_med=float(np.median(prof_R[midband])),
    roll_cos1_uniform=float(harmonics_uniform(
        U_matrix(Gf(TH_GRID), ub_head)
        - U_matrix(Gf(180.0 - TH_GRID), ub_head), TH_GRID)["cos1"]),
    measured_contrast_yr=meas_contrast)
logger.info("headline: " + json.dumps(
    {k: (round(v, 4) if isinstance(v, float) else v)
     for k, v in headline.items()}))

# Newtonian-tracking baseline
dtM_track, _ = slip_amplitude(Gf, U_TRACK, ell_m)
track_pred = dict(
    u_b=U_TRACK,
    pred_contrast_s=float((np.median(dtM_track[in_pole])
                           - np.median(dtM_track[midband])) * YR_S),
    statement="potential-tracking amplitude under-produces the "
              "measured slip by ~24 orders -- a genuine boundary "
              "field structure is required, not local tracking")

anchors = {"Newtonian_tracking_250AU": U_TRACK,
           "Galactic_ambient": U_GAL,
           "solar_surface": U_SUN_SURF,
           "required_u_b_dip_ell8": ub_head,
           "cosmological_ambient": U_COSMO}

verdict = {
    "required_u_b": ub_head,
    "statement": (
        "On the admissible B0 = +3.2e-3 benchmark, whose magnitude "
        "is mirrored from the sign-excluded Paper-0 reconstruction "
        "only as an instrumental scale, the canonical "
        "synchronisation connection "
        "integrated along comet trajectories reproduces the "
        "measured bipolar slip contrast (+%.1f yr pole-versus-"
        "midband) for a boundary field excursion u_b ~ %.1e.  The "
        "morphology is not tuned: an axis-anchored field crossed "
        "at antipodal points gives a bipolar (cos 2th) open-path "
        "slip identically -- every odd harmonic cancels in the "
        "double crossing -- matching the measured harmonic "
        "structure (cos 2th dominated, cos th flat), while the "
        "resident-cloud source response stays linear/dipolar "
        "(step 164).  The required excursion sits ~10^3 above "
        "Galactic ambient tracking and ~10^3 below the "
        "cosmological value: the 250 AU boundary would have to be "
        "a genuine temporal-structure boundary -- the edge of the "
        "Solar System's screening domain -- rather than a "
        "kinematic marker.  Potential-tracking amplitudes "
        "(u ~ 4e-11) under-produce the slip by ~24 orders."
        % (meas_contrast, ub_head)),
    "sign_check": (
        "The fitted step-065 equivalent time is constructed from an "
        "unsigned angular separation and therefore cannot select the "
        "sign of B.  The admissible B0 > 0 and excluded B0 < 0 "
        "benchmarks give identical amplitude profiles and opposite "
        "signed open-path transport.  A signed transfer function is "
        "required before the independent signed-rotation statistic "
        "can be mapped to the synchronization-connection orientation."),
    "morphology_check": (
        "any axis-anchored field -> bipolar slip via the antipodal "
        "double crossing (odd harmonics cancel identically); the "
        "network-roll channel is dipolar and subdominant"),
    "absorption_note": (
        "the theta-uniform part of the lag is absorbed into the "
        "fitted elements (parameter absorption, step 151); the "
        "observable is the theta-dependent contrast"),
    "remaining_gap": (
        "u_b is measured, not derived: closure requires solving "
        "the screening-domain field profile at the 250 AU "
        "boundary and checking whether the solved depth reaches "
        "~1e-4"),
    "admissibility": (
        "see admissibility_ledger: the required excursion is "
        "unsourceable by the local field hierarchy, carries a "
        "mandatory ~120 ppm conformal clock step and a ~10^10x "
        "over-bound shear impulse, and the admissible B >= 0 "
        "branch predicts the wrong sign; the result is a "
        "falsification condition on the solved field profile, "
        "not mechanism closure")}

# ------------------------------------------------------------------
# admissibility ledger -- does the corpus's own architecture admit
# the required excursion?  Computed on real anchors, not asserted.
# ------------------------------------------------------------------
BETA_A = -1.0                          # universal conformal coupling
A_AMB = 3.9e-10                        # ambient Galactic shear floor
                                       # m/s^2 (Paper 0)
SCLK_SUST_PPM = 4.0                    # spacecraft sustained-step
                                       # reach (candidates -0.9..-4.4 ppm)
SCLK_INST = 2.0e-5                     # instantaneous onboard-vs-ground
                                       # contrast bound (step 096/145)
DALPHA_HOL = 6.04e-13                  # step_145 equivalent conformal
                                       # step, pre-2018 record
G_SUN_RB = 1.32712440018e20 / R_B**2   # solar g at r_b ~ 9.5e-8 m/s^2
V_CROSS_MED = float(np.median(VR))     # m/s

ledger = {}
ledger["required_u_b"] = ub_head
ledger["sourceable_budget"] = {
    "galactic_ambient_u": U_GAL,
    "ratio_vs_galactic": ub_head / U_GAL,
    "solar_tracking_250AU_u": U_TRACK,
    "ratio_vs_solar": ub_head / U_TRACK,
    "note": ("an excursion intermediate between the Galactic ambient "
             "and the cosmological field is landscape-level structure, "
             "not a solar screening transition -- a screening edge joins "
             "the solar solution to the ambient, it cannot exceed it")}

# conformal correlate: any real phi excursion carries a clock-rate step
dlnA_ppm = ub_head * 1.0e6             # |d ln A| = |beta_A| u_b
ledger["conformal_correlate"] = {
    "delta_lnA_ppm": dlnA_ppm,
    "sclk_sustained_bound_ppm": SCLK_SUST_PPM,
    "ratio_vs_sustained_bound": dlnA_ppm / SCLK_SUST_PPM,
    "sclk_instantaneous_bound": SCLK_INST,
    "ratio_vs_instantaneous_bound": ub_head / SCLK_INST,
    "voyager_note": ("at r_b <= 165 AU the step is already excluded by "
                     "the Voyager SCLK record; at the step-156 best fit "
                     "(~250 AU) no spacecraft has crossed, so the "
                     "~120 ppm crossing step is a dated prediction "
                     "the text must state, not omit")}

# shear correlate: same gradient sources c^2 grad ln A; the integrated
# impulse across the wall is width-independent, c^2 u_b / v_cross
shear = {}
for w_au in ELL_LIST + [250.0]:
    a_shear = C * C * ub_head / (w_au * AU_M)
    shear[f"W_{w_au:g}AU"] = {
        "a_m_s2": a_shear,
        "ratio_vs_ambient_floor": a_shear / A_AMB,
        "ratio_vs_solar_g_at_rb": a_shear / G_SUN_RB}
ledger["shear_correlate"] = shear
# energy-channel bound on the wall-normal impulse, computed from the
# measured |D(1/a)| channel rather than hardcoded: the along-track
# impulse satisfies |dv_par| <= GM |D(1/a)| / (2 v); against the
# crossing impulse c^2 u_b / v the speed cancels, leaving
# u_b <= GM |daa| / (2 c^2) per comet (daa in 1e-6 AU^-1).
GM_SI = 1.32712440018e20          # m^3 s^-2
DAA = np.array([c["daa"] for c in comets])
DAA_MED = float(np.nanmedian(DAA))          # 1e-6 AU^-1
DV_ENERGY = (GM_SI * DAA_MED * 1e-6 / AU_M
             / (2.0 * V_CROSS_MED))          # m/s
U_B_ENERGY = GM_SI * DAA_MED * 1e-6 / AU_M / (2.0 * C**2)
ledger["shear_impulse"] = {
    "delta_v_m_s": C * C * ub_head / V_CROSS_MED,
    "energy_channel_bound_m_s": DV_ENERGY,
    "energy_channel_bound_source": (
        "median measured |D(1/a)| = %.1f x1e-6 AU^-1 over the slip "
        "cohort; total channel (encounter kicks included), hence "
        "conservative.  The earlier hardcoded 0.3 m/s traced to a "
        "step-145 expression that substituted the planetary-approach "
        "distance (AU) for the energy residual (1e-6 AU^-1) -- a "
        "wrong-variable defect corrected in step_145 T2."
        % DAA_MED),
    "ratio": C * C * ub_head / V_CROSS_MED / DV_ENERGY,
    "u_b_bound_from_energy_channel": U_B_ENERGY,
    "note": ("the wall-normal impulse c^2 du/dr integrated across the "
             "layer is c^2 u_b / v -- independent of W; exceeds the "
             "along-track energy-channel bound by ~10^8")}

# internal realization consistency
ledger["realization_consistency"] = {
    "holonomy_equiv_conformal_step": DALPHA_HOL,
    "disformal_required_u_b": ub_head,
    "ratio": ub_head / DALPHA_HOL,
    "note": ("the two quoted realizations differ by ~8 orders; they "
             "cannot describe the same field")}


def _contrast_at(ub):
    pM, _ = slip_amplitude(Gf, ub, ell_m)
    return float(np.median(pM[in_pole]) - np.median(pM[midband]))


ub_sust = 4.0e-6
slip_sust = _contrast_at(ub_sust)
slip_inst = _contrast_at(SCLK_INST)
ledger["admissible_sector"] = {
    "admissible_branch": "B >= 0 (matter cone inside gravitational "
                         "cone, Paper 0 Section 4)",
    "sign_under_admissible_branch": (
        "the admissible B >= 0 envelope gives the opposite signed "
        "open-path transport to the excluded B < 0 reconstruction. "
        "The step-065 amplitude used here is unsigned, so it neither "
        "selects nor falsifies either orientation; the earlier "
        "wrong-sign conclusion was a sign-observability error"),
    "u_b_at_sustained_clock_bound": ub_sust,
    "max_slip_contrast_sustained_yr": slip_sust,
    "u_b_at_instantaneous_clock_bound": SCLK_INST,
    "max_slip_contrast_instantaneous_yr": slip_inst,
    "deficit_vs_measured_sustained": meas_contrast / slip_sust,
    "deficit_vs_measured_instantaneous": meas_contrast / slip_inst}

# joint (B0, u_b) requirement plane -- the falsification condition is
# two-dimensional, not a property of the single mirrored benchmark.
# The slip is exactly linear in B0, so the B0 that would close the
# measured contrast at any candidate excursion is
#   B0_req(u_b) = B0 * meas_contrast / contrast(u_b).
# Each channel bound on u_b therefore maps to a required admissible
# normalization.
u_b_energy = U_B_ENERGY   # energy-channel crossing bound (data-derived)
b0_plane = {}
for name, ub in [
        ("galactic_ambient", U_GAL),
        ("sclk_sustained_bound", 4.0e-6),
        ("sclk_instantaneous_bound", SCLK_INST),
        ("energy_channel_bound", u_b_energy),
        ("solar_tracking_250AU", U_TRACK)]:
    c_ub = _contrast_at(ub)
    b0_plane[name] = {
        "u_b": ub,
        "slip_contrast_at_benchmark_yr": c_ub,
        "B0_required_for_measured": B0 * meas_contrast / c_ub}
ledger["b0_requirement_plane"] = {
    "scaling": "slip proportional to B0 exactly; contrast ~ B0 u_b^4 "
               "in the small-u regime (env_b ~ u^2)",
    "entries": b0_plane,
    "note": ("the energy-channel bound (u_b <= %.1e, from the "
             "wall-normal impulse c^2 u_b / v <= %.1f m/s, the "
             "measured |D(1/a)| floor) is the binding constraint: "
             "closure at that excursion requires an admissible "
             "B0 ~ %.1e.  The verdict 'unsourceable excursion' is "
             "equivalently 'underived normalization': the "
             "falsification condition is the (B0, u_b) plane, and it "
             "stays closed unless either a sourceable ~1e-4 wall or "
             "an admissible B0 many orders above the mirrored "
             "benchmark is derived."
             % (u_b_energy, DV_ENERGY,
                b0_plane["energy_channel_bound"]
                ["B0_required_for_measured"]))}

ledger["verdict"] = (
    "falsification condition, not mechanism closure: u_b ~ %.1e is a "
    "required input the solved admissible field profile must supply. "
    "It exceeds the entire ambient field hierarchy by ~10^2-10^6, "
    "carries a mandatory %.0f ppm clock step (~%.0fx the sustained "
    "spacecraft bound) and a %.1e x over-floor shear impulse "
    "that violates the flat energy channel. The branch-sign audit "
    "does not add an exclusion: the fitted equivalent-time datum is "
    "unsigned, and the admissible B >= 0 branch has now been evaluated "
    "explicitly. The channel is thereby converted to a bound: the "
    "admissible boundary-localised "
    "disformal slip is <= %.1e yr (sustained clock bound) to "
    "%.1e yr (instantaneous bound), versus the measured %.1f yr."
    % (ub_head, dlnA_ppm, dlnA_ppm / SCLK_SUST_PPM,
       ledger["shear_impulse"]["ratio"],
       slip_sust, slip_inst, meas_contrast))

logger.info("admissibility ledger: " + json.dumps(
    {k: (round(v, 4) if isinstance(v, float) else v)
     for k, v in ledger["sourceable_budget"].items()
     if isinstance(v, float)}))
logger.info("  conformal step %.0f ppm = %.0fx sustained bound; "
            "shear impulse %.2e m/s = %.1e x bound; "
            "admissible slip <= %.1e-%.1e yr vs measured %.1f yr"
            % (dlnA_ppm, dlnA_ppm / SCLK_SUST_PPM,
               ledger["shear_impulse"]["delta_v_m_s"],
               ledger["shear_impulse"]["ratio"],
               slip_sust, slip_inst, meas_contrast))

branch_comparison = {
    "admissible_B0": B0,
    "excluded_reference_B0": B0_EXCLUDED,
    "normalization_provenance": (
        "|B0| = 3.2e-3 is inherited from the excluded Paper-0 "
        "volume-balance reconstruction only as an amplitude benchmark; "
        "the admissible normalization is not derived"),
    "admissible_signed_transport_contrast_yr": float(
        headline["signed_transport_contrast"]),
    "excluded_signed_transport_contrast_yr": float(
        -headline["signed_transport_contrast"]),
    "observable_amplitude_contrast_yr": float(
        np.median(prof_M[in_pole]) - np.median(prof_M[midband])),
    "sign_observability": (
        "step_065 uses angular separations and measures |delta t_eq|; "
        "it cannot adjudicate the transport sign.  The sign retained "
        "by step_080 has not yet been propagated through a signed "
        "plus/minus disformal injection and is not used to select B.")}

# ------------------------------------------------------------------
# outputs
# ------------------------------------------------------------------
out_json = RESULTS / "step_b131_disformal_slip.json"
out_csv = RESULTS / "step_b131_disformal_slip.csv"
out_fig = RESULTS / "figures/supplementary/step_b131_disformal_slip.png"

json.dump(dict(step="step_165_disformal_slip", seed=SEED,
               model="dt_slip = (1/c) int dsigma.dx; dsigma_mu = "
                     "-(B/A^2)(u.grad phi) P grad phi; admissible "
                     "positive-B benchmark; excursion convention for u",
               constants=dict(B0=B0, B0_excluded_reference=B0_EXCLUDED,
                              R_H_m=R_H, H0=H0,
                              r_b_AU=250.0, cap_deg=CAP),
               anchors=anchors,
               headline=headline,
               tracking_baseline=track_pred,
               scan=results,
               verdict=verdict,
               branch_comparison=branch_comparison,
               admissibility_ledger=ledger),
          open(out_json, "w"), indent=1)

with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
    w.writeheader()
    for r in csv_rows:
        w.writerow(r)

# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

ax = axes[0]
ub_scan = np.logspace(-6, -2, 60)
contr = np.empty(len(ub_scan))
for j, ub in enumerate(ub_scan):
    pM, _ = slip_amplitude(Gf, ub, ell_m)
    contr[j] = np.median(pM[in_pole]) - np.median(pM[midband])
ax.loglog(ub_scan, np.abs(contr), lw=2)
ax.axhline(meas_contrast, color="r", ls="--",
           label=f"measured {meas_contrast:.1f} yr")
ax.axvline(ub_head, color="k", ls=":", label=f"u_b = {ub_head:.1e}")
for nm, uu, c in (("tracking", U_TRACK, "gray"),
                  ("Galactic", U_GAL, "tab:blue"),
                  ("solar", U_SUN_SURF, "tab:green")):
    ax.axvline(uu, color=c, alpha=0.35)
ax.set_xlabel("boundary field excursion $u_b$")
ax.set_ylabel("pole-midband slip contrast (yr)")
ax.legend(fontsize=8)
ax.set_title("Amplitude closure")

ax = axes[1]
for mname, Gfm in MORPHS.items():
    prof = K_matrix(Gfm(TH_GRID), ub_head, ell_m) + K_matrix(
        Gfm(180.0 - TH_GRID), ub_head, ell_m)
    prof = (abs(B0) * R_H**2 / C**2) * np.median(VR) * prof / YR_S
    ax.plot(TH_GRID, prof - prof.mean(), label=mname)
ax.scatter(TH, DT_OBS - np.median(DT_OBS), s=6, alpha=0.25,
           color="gray", label="measured (centred)")
ax.set_xlabel(r"axis separation $\theta_A$ (deg)")
ax.set_ylabel("slip pattern (yr)")
ax.legend(fontsize=8)
ax.set_title("Predicted theta-structure")

ax = axes[2]
ax.scatter(prof_M, DT_OBS, s=8, alpha=0.4)
ax.set_xlabel("predicted slip (yr)")
ax.set_ylabel("measured residual slip (yr)")
ax.set_title(f"per-comet: rho = {rho_pred.statistic:.3f}, "
             f"p = {p_rho:.3f}")
fig.tight_layout()
fig.savefig(out_fig, dpi=140)
logger.info(f"outputs: {out_json}, {out_csv}, {out_fig}")
