"""step_165: canonical disformal transport benchmark and geometry controls.

Fits the step-065 unsigned equivalent-time amplitude using the positive-B
benchmark envelope b(u)=B0*u^2/(1+u^2)*exp(-u^4/(2*sigma_B^4)). The leading
connection term is integrated in a radial boundary layer. This is a scale
requirement, not a derived admissible normalization or field solution.

Per-comet radial speeds are computed from the residual cohort's inverse
semimajor axes and angular momenta at 250 AU. Spatial entry and exit directions
come from the published original/future conic planes. Near-parabolic orbits
cross near the same aphelion-side region; they are not antipodal. The previous
antipodal construction imposed generic odd-mode cancellation artificially.

The exploratory signed-dipole excursion G(theta)=cos(theta) has an even
quartic amplitude in the coincident-aphelion limit because the envelope is
even and the traversal gradient is squared. Offset dipoles and quadrupolar
profiles are retained as comparison cases. Shape correlations are conditional
on this model exploration, not calibrated discovery probabilities.

The signed transport and unsigned fitted statistic are stored separately.
Clock and one-crossing impulse comparisons are conditional diagnostics; the
full metric, trajectory and observation response remains to be validated.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, DATA_RAW, tee_stdout, parse_code, perih_dir, lv, sep
from scripts.utils.statistics import monte_carlo_tail
from scripts.utils.boundary_transport import (radial_crossing_speed, angular_harmonics,
                                                conic_boundary_direction)
logger = StepLogger("step_165_disformal_slip")
tee_stdout(logger)
logger.header("Disformal transport: amplitude and spatial crossing geometry")

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
logger.data_load(RESULTS / "step_b30_proper_time_slip.csv",
                 "Residual slip, inverse semimajor axis and angular momentum")

# Use the published orbital planes to sample the field at the spatial
# entry/exit locations, rather than imposing antipodal positions.
code_in = parse_code(str(DATA_RAW / "code/code_original.html"))
code_out = parse_code(str(DATA_RAW / "code/code_future.html"))

def _warsaw(path, preference):
    rows = {}
    for line in path.read_text().splitlines():
        if len(line) < 115:
            continue
        try:
            name, model = line[5:17].strip(), line[3].strip()
            record = dict(q=float(line[42:56]), e=float(line[56:70]),
                          w=float(line[70:82]), Om=float(line[82:94]),
                          i=float(line[94:106]), model=model)
        except ValueError:
            continue
        if name not in rows or preference.get(model, 9) < preference.get(rows[name]["model"], 9):
            rows[name] = record
    return rows

warsaw_in = _warsaw(DATA_RAW / "warsaw/warsaw_tablec.dat", {"a":0,"h":0,"e":1,"b":2})
warsaw_out = _warsaw(DATA_RAW / "warsaw/warsaw_tabled.dat", {"i":0,"l":0,"j":2,"k":2})
FIELD_AXIS = lv(34.0, -13.0)

def _crossing(record, outbound):
    w, node, inc = map(math.radians, (record["w"], record["Om"], record["i"]))
    peri = perih_dir(w, node, inc)
    normal = np.array([math.sin(inc)*math.sin(node),
                       -math.sin(inc)*math.cos(node), math.cos(inc)])
    return conic_boundary_direction(record["q"], record["e"], peri, normal,
                                    R_B / AU_M, outbound)

comets = []
for d, r in slip_rows.items():
    if r.get("dtau_unexplained") in ("", None) or r.get("theta") in (
            "", None):
        continue
    vc = float(radial_crossing_speed(
        float(r["aa_back"]), float(r["h"]), R_B / AU_M))
    incoming, outgoing = ((code_in, code_out) if r["cohort"] == "code"
                          else (warsaw_in, warsaw_out))
    entry, exit_ = _crossing(incoming[d], False), _crossing(outgoing[d], True)
    comets.append(dict(desig=d, theta=float(r["theta"]),
                       theta_entry=sep(entry, FIELD_AXIS), theta_exit=sep(exit_, FIELD_AXIS),
                       crossing_separation=sep(entry, exit_),
                       dtau=float(r["dtau_unexplained"]),
                       daa=abs(float(r["daa_sim"]))
                       if r.get("daa_sim") not in ("", None) else np.nan,
                       v_cross=vc, cohort=r["cohort"]))
logger.info(f"cohort: {len(comets)} comets with residual slip, "
            f"theta, crossing speed")

TH = np.array([c["theta"] for c in comets])
TH_ENTRY = np.array([c["theta_entry"] for c in comets])
TH_EXIT = np.array([c["theta_exit"] for c in comets])
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
# Entry/exit use published conic spatial directions. The radial-layer
# approximation freezes each direction across the layer.
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
    return np.trapezoid(integrand, s, axis=0)


def U_matrix(Gvals, u_b, n=2001):
    """int_0^{u_b G} env_b(u) du  (roll channel)."""
    uu = np.linspace(0.0, 1.0, n)[None, :] * (u_b * Gvals)[:, None]
    return np.trapezoid(env_b(uu), uu, axis=1)


def slip_profiles(Gfun, u_b, ell_m):
    """Signed per-comet transport [yr], channels M and R.

    The leading connection carries an explicit minus sign, so the
    transport orientation reverses with B0.  Observable amplitude
    comparisons must be made outside this function.
    """
    Gin = Gfun(TH_ENTRY)
    Gout = Gfun(TH_EXIT)
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
    return angular_harmonics(y, th_deg)


def P2(x):
    return 0.5 * (3.0 * x * x - 1.0)


MORPHS = {
    "signed_dip": lambda t: np.cos(np.radians(t)),
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
        Gout = Gf(TH_GRID)
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
# headline: exploratory signed-dipole excursion, ell = 8 AU
# ------------------------------------------------------------------
Gf = MORPHS["signed_dip"]
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
p_rho = monte_carlo_tail(null_rho, rho_pred.statistic, "two-sided")
# uniform-grid morphology for the headline
prof_grid = K_matrix(Gf(TH_GRID), ub_head, ell_m) + K_matrix(
    Gf(TH_GRID), ub_head, ell_m)
h_head = harmonics_uniform(prof_grid, TH_GRID)

headline = dict(
    morph="signed_dip: excursion cos(theta); published conic spatial crossings",
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
        - U_matrix(Gf(TH_GRID), ub_head), TH_GRID)["cos1"]),
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
           "required_u_b_signed_dip_ell8": ub_head,
           "cosmological_ambient": U_COSMO}

verdict = {
    "required_u_b": ub_head,
    "statement": (
        "At the positive-B benchmark, the exploratory signed-dipole excursion "
        "cos(theta) evaluated at published conic spatial entry/exit directions "
        "fits the %.3g-year amplitude contrast with u_b = %.3g. The observed "
        "contrast is fitted, not predicted. A signed dipole gives an even "
        "quartic amplitude in the coincident-aphelion limit because env_b is "
        "even and the traversal term is squared. Offset dipoles do not share "
        "that symmetry. Nearly parabolic entry and exit are nearby, rather "
        "than antipodal, so the previous generic odd-mode cancellation claim "
        "does not apply to these comet trajectories. Source-profile, metric "
        "and observation-transfer closure remain required."
        % (meas_contrast, ub_head)),
    "sign_check": (
        "The fitted step-065 equivalent time is constructed from an "
        "unsigned angular separation and therefore cannot select the "
        "sign of B.  The admissible B0 > 0 and excluded B0 < 0 "
        "benchmarks give identical amplitude profiles and opposite "
        "signed open-path transport.  A signed transfer function is "
        "required before the independent signed-rotation statistic "
        "can be mapped to the synchronization-connection orientation."),
    "morphology_check": "Even response is specific to the tested signed-dipole excursion in the symmetric limit; actual crossing locations are retained per comet. Profile choice is exploratory.",
    "absorption_note": (
        "the theta-uniform part of the lag is absorbed into the "
        "fitted elements (parameter absorption, step 151); the "
        "observable is the theta-dependent contrast"),
    "remaining_gap": (
        "u_b is measured, not derived: closure requires solving "
        "the screening-domain field profile at the 250 AU "
        "boundary and checking whether the solved depth reaches "
        "~%.3g" % ub_head),
    "admissibility": "See the computed admissibility_ledger below."}

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
V_CROSS_MED = float(np.median(VR))     # radial m/s
V_TOTAL_MED = float(np.median([
    math.hypot(c["v_cross"], float(slip_rows[c["desig"]]["h"]) / 250.0)
    for c in comets])) * AU_M / YR_S

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
             / (2.0 * V_TOTAL_MED))          # m/s, along-track
U_B_ENERGY = GM_SI * DAA_MED * 1e-6 / AU_M / (2.0 * C**2)
ledger["shear_impulse"] = {
    "delta_v_m_s": C * C * ub_head / V_CROSS_MED,
    "along_track_delta_v_m_s": C * C * ub_head / V_TOTAL_MED,
    "energy_channel_bound_m_s": DV_ENERGY,
    "energy_channel_bound_source": (
        "median measured |D(1/a)| = %.1f x1e-6 AU^-1 over the slip "
        "cohort; total channel (encounter kicks included), hence "
        "conservative.  The earlier hardcoded 0.3 m/s traced to a "
        "step-145 expression that substituted the planetary-approach "
        "distance (AU) for the energy residual (1e-6 AU^-1) -- a "
        "wrong-variable defect corrected in step_145 T2."
        % DAA_MED),
    "ratio": C * C * ub_head / V_TOTAL_MED / DV_ENERGY,
    "u_b_bound_from_energy_channel": U_B_ENERGY,
    "note": ("One-crossing, weak-impulse diagnostic only. The radial kick is "
             "c^2 u_b/v_r; its along-track projection is c^2 u_b/v. "
             "A static conservative field can cancel the energy change on "
             "the return crossing, so an inter-leg net-energy floor is not "
             "a direct bound on each crossing without a transfer model. "
             "The fitted benchmark violates the small-impulse approximation."),
    "bound_status": "conditional on an uncancelled along-track impulse"}

# internal realization consistency
ledger["realization_consistency"] = {
    "holonomy_equiv_conformal_step": DALPHA_HOL,
    "disformal_required_u_b": ub_head,
    "ratio": ub_head / DALPHA_HOL,
    "note": ("The equivalent impulse scale and fitted disformal excursion "
             "use different transfer models. Their ratio is diagnostic; "
             "equality is not required without deriving their relation.")}


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
    "note": ("Conditional scale map, not a calibrated exclusion. The energy "
             "row assumes the one-crossing impulse survives in the measured "
             "net-energy channel; clock rows apply only where a trajectory "
             "samples the assumed boundary. None supplies an independently "
             "derived admissible B0 or validates the leading-order transport.")}

# Track where the approximations used to obtain this benchmark fail.
_s, _sig, _dsig = _layer(ell_m)
_u = ub_head * _sig
_grad = ub_head * _dsig
ledger["approximation_diagnostics"] = {
    "max_disformal_spatial_term": float(np.max(B0 * env_b(_u) * R_H**2 * _grad**2)),
    "formal_radial_kick_over_crossing_speed": C**2 * ub_head / V_CROSS_MED**2,
    "status": "leading-order benchmark requires full metric/trajectory validation",
    "note": "The leading connection uses A^2 N^2 approximately one and the impulse uses fixed unperturbed speed. Large diagnostics invalidate treating the benchmark as a closed physical prediction."}


ledger["verdict"] = (
    "The fitted leading-order benchmark requires u_b = %.3g and implies "
    "a %.1f ppm conformal step. This is an amplitude requirement, not a "
    "derived field profile. The formal one-crossing along-track impulse "
    "is %.3g times the measured energy floor, but the small-impulse "
    "approximation fails and the net two-crossing energy response is not "
    "computed here. Spacecraft clock limits do not directly exclude a "
    "250-AU wall that no spacecraft has crossed. The unsigned equivalent "
    "time cannot exclude either transport sign. A sourceable profile, "
    "admissible normalization and full trajectory/observation response "
    "remain required."
    % (ub_head, dlnA_ppm, ledger["shear_impulse"]["ratio"]))

# One source of truth for the human-readable summary.
verdict["admissibility"] = ledger["verdict"]
logger.info("admissibility ledger: " + json.dumps(
    {k: (round(v, 4) if isinstance(v, float) else v)
     for k, v in ledger["sourceable_budget"].items()
     if isinstance(v, float)}))
logger.info("  conformal step %.0f ppm = %.0fx sustained bound; "
            "formal shear impulse %.2e m/s = %.1e x conditional floor; "
            "clock-threshold benchmark %.1e-%.1e yr vs measured %.1f yr"
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
               inputs=["results/step_b30_proper_time_slip.csv",
                       "data/raw/code/code_original.html", "data/raw/code/code_future.html",
                       "data/raw/warsaw/warsaw_tablec.dat", "data/raw/warsaw/warsaw_tabled.dat"],
               crossing_geometry=dict(
                   approximation="Published original/future two-body conics at 250 AU; not full metric trajectories",
                   median_separation_deg=float(np.median([c["crossing_separation"] for c in comets])),
                   min_separation_deg=float(min(c["crossing_separation"] for c in comets)),
                   max_separation_deg=float(max(c["crossing_separation"] for c in comets)),
                   field_axis_ecliptic_deg=[34.0,-13.0],
                   profile_selection="Exploratory signed-dipole excursion; not an independently calibrated profile",
                   note="Spatial entry and exit are not generally antipodal. Uniform-grid harmonics use the coincident-aphelion limit; cos(theta)^4 is even because the chosen excursion is odd, not because all odd modes cancel on real crossings."),
               velocity_provenance=dict(
                   source="step_b30 aa_back and h; two-body radial vis-viva",
                   radius_AU=R_B / AU_M, n_comets=len(comets),
                   n_imputed=0,
                   median_AU_yr=float(np.median(VR) * YR_S / AU_M),
                   note="Replaces failed exact-name join to a 100-AU total-speed table; no fallback speeds."),
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
ax.set_title("Required amplitude (benchmark)")

ax = axes[1]
for mname, Gfm in MORPHS.items():
    prof = K_matrix(Gfm(TH_GRID), ub_head, ell_m) + K_matrix(
        Gfm(TH_GRID), ub_head, ell_m)
    prof = (abs(B0) * R_H**2 / C**2) * np.median(VR) * prof / YR_S
    ax.plot(TH_GRID, prof - prof.mean(), label=mname)
ax.scatter(TH, DT_OBS - np.median(DT_OBS), s=6, alpha=0.25,
           color="gray", label="measured (centred)")
ax.set_xlabel(r"axis separation $\theta_A$ (deg)")
ax.set_ylabel("slip pattern (yr)")
ax.legend(fontsize=8)
ax.set_title("Profiles at shared fitted amplitude")

ax = axes[2]
ax.scatter(prof_M, DT_OBS, s=8, alpha=0.4)
ax.set_xlabel("predicted slip (yr)")
ax.set_ylabel("measured residual slip (yr)")
ax.set_title(f"per-comet: rho = {rho_pred.statistic:.3f}, "
             f"p = {p_rho:.3f}")
fig.tight_layout()
fig.savefig(out_fig, dpi=140)
for output in (out_json, out_csv, out_fig):
    logger.data_save(output)
logger.info(f"outputs: {out_json}, {out_csv}, {out_fig}")
