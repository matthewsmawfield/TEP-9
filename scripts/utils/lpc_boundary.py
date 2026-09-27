"""Bidirectional boundary integrator for near-parabolic comets.

Shared by step_117 (post-2017 prospective cohort), step_120 (pre-2018
SBDB era x lineage factorial) and any later LPC cohort step.  The
instrument is identical to step_063/064/065: each comet's published
osculating elements are converted to a Cartesian state at the
osculation epoch (anchored at perihelion so the e ~ 1 mean-motion
limit never enters), then integrated both directions to the 250 AU
barycentric sphere under REBOUND/IAS15 with the Sun plus planetary-
system barycentres on DE440s.

All functions are pure given the catalogue inputs; the SPICE kernel
must be furnished in every process that calls them (see worker_init).
"""

import math

import numpy as np
import rebound
import spiceypy as sp

from scripts.utils.tep9_common import perih_dir, sep, lv

AU_KM = 149597870.7
DAY_YR = 365.25
EPS = math.radians(23.4392911)
RX = np.array([[1, 0, 0],
               [0, math.cos(EPS), math.sin(EPS)],
               [0, -math.sin(EPS), math.cos(EPS)]])

GM_SUN = 1.32712440018e11
GM = {"1": 2.2031868551e4, "2": 3.2485859200e5, "3": 4.0350323562e5,
      "4": 4.2828375814e4, "5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6, "9": 1.0868657e3}
PLANET_IDS = list(GM.keys())
MU = 4 * math.pi ** 2
R_STOP = 250.0
T_MAX = 20000.0
DT_OUT = 1.0
V_BND_MAX = 10.0  # |v| ceiling at the 250 AU sphere (escape ~0.6 AU/yr)
R_MIN_ENC = 0.001  # exit threshold: comet driven inside 0.001 AU of a
                   # massive body is inside its Roche zone and is
                   # destroyed -- bounds IAS15 against adaptive-step
                   # collapse without ever firing on a real approach
DT_MIN_YR = 1e-8  # ~0.3 s: below this the integrator is numerically
                  # degenerate even absent a registered encounter;
                  # real deep-encounter transients stay above it

TNO = lv(34.0, -13.0)          # pre-declared transit axis (steps 030-086)


def worker_init(spk_path):
    """Pool initializer: fork children inherit the parent's kernel-
    table entries but not their DAF file records, so the inherited
    handles are dead; reset the pool and load a private kernel."""
    sp.kclear()
    sp.furnsh(str(spk_path))


def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return (RX @ np.array(st[:3]) / AU_KM,
            RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR)


def comet_state(rec):
    """Heliocentric ecliptic-J2000 state at the SBDB osculation epoch.

    sp.conics converts the published osculating elements to a Cartesian
    state -- the elements are ecliptic J2000, so the state lands in the
    same frame the DE440s barycentres are rotated into.  Anchoring the
    conic at perihelion (M0 = 0 at t0 = tp) is exact regardless of the
    e ~ 1 mean-motion limit.
    """
    et0 = (rec["epoch"] - 2451545.0) * 86400.0
    tp_et = (rec["tp"] - 2451545.0) * 86400.0
    elts = [rec["q"] * AU_KM, rec["e"], math.radians(rec["i"]),
            math.radians(rec["om"]), math.radians(rec["w"]),
            0.0, tp_et, GM_SUN]
    st = np.array(sp.conics(elts, et0))
    return (st[:3] / AU_KM, st[3:] / AU_KM * 86400.0 * DAY_YR, et0)


def init_sim(et):
    sim = rebound.Simulation()
    sim.G = MU
    ps, vs = body_state("10", et)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2],
            m=1.0)
    for b in PLANET_IDS:
        pp, vv = body_state(b, et)
        sim.add(x=pp[0], y=pp[1], z=pp[2], vx=vv[0], vy=vv[1], vz=vv[2],
                m=GM[b] / GM_SUN)
    return sim


def boundary_orbit(r_rel, v_rel, mtot):
    mu = MU * mtot
    r = np.linalg.norm(r_rel)
    h = np.cross(r_rel, v_rel)
    evec = (np.cross(v_rel, h) / mu) - r_rel / r
    en = np.linalg.norm(evec)
    phat = evec / en if en > 1e-12 else r_rel / r
    E = v_rel.dot(v_rel) / 2 - mu / r
    return phat, -2 * E / mu * 1e6


def integrate_leg(r0, v0, et0, direction, t_max=None,
                  v_max=V_BND_MAX):
    """One leg to the +/-250 AU barycentric sphere. direction=-1 gives
    the inbound (original-analogue) asymptote, +1 the outbound.
    t_max overrides the global T_MAX iteration budget; callers that
    perturb orbits can pass a smaller bound so pathological draws
    (deeply bound solutions that never reach the sphere) cannot stall
    a worker for the full 20 kyr march.  v_max overrides the boundary
    speed sanity ceiling; interstellar objects carry real hyperbolic
    excess (~5-13 AU/yr), so ISO refits pass a wider bound while the
    comet default stays at the escape-speed scale."""
    if t_max is None:
        t_max = T_MAX
    sim = init_sim(et0)
    p = sim.particles
    ps = np.array([p[0].x, p[0].y, p[0].z])
    vs = np.array([p[0].vx, p[0].vy, p[0].vz])
    sim.add(x=r0[0] + ps[0], y=r0[1] + ps[1], z=r0[2] + ps[2],
            vx=v0[0] + vs[0], vy=v0[1] + vs[1], vz=v0[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = R_MIN_ENC
    denc = np.full(sim.N - 1, np.inf)
    t = direction * DT_OUT
    while abs(t) < t_max:
        sim.integrate(t, exact_finish_time=0)
        if abs(sim.dt) < DT_MIN_YR:
            raise RuntimeError("IAS15 timestep collapse "
                               f"(dt={sim.dt:.2e} yr)")
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y,
                          p[nc].z - p[0].z])
        for j in range(1, sim.N - 1):
            d = math.sqrt((p[nc].x - p[j].x) ** 2 +
                          (p[nc].y - p[j].y) ** 2 +
                          (p[nc].z - p[j].z) ** 2)
            if d < denc[j - 1]:
                denc[j - 1] = d
        if np.linalg.norm(r_rel) >= R_STOP:
            break
        t += direction * DT_OUT
    else:
        return None
    mtot = sum(pp.m for pp in sim.particles)
    rb = np.zeros(3); vb = np.zeros(3)
    for pp in sim.particles:
        rb += pp.m * np.array([pp.x, pp.y, pp.z])
        vb += pp.m * np.array([pp.vx, pp.vy, pp.vz])
    rb /= mtot; vb /= mtot
    p = sim.particles
    r_rel = np.array([p[nc].x, p[nc].y, p[nc].z]) - rb
    v_rel = np.array([p[nc].vx, p[nc].vy, p[nc].vz]) - vb
    if np.linalg.norm(v_rel) > v_max:
        raise RuntimeError(
            f"unphysical boundary speed |v|={np.linalg.norm(v_rel):.1f} AU/yr")
    phat, aa = boundary_orbit(r_rel, v_rel, mtot)
    return dict(phat=phat, aa=aa, denc=float(denc.min()),
                t_years=float(t))


def process_comet(rec, axis=TNO):
    """Integrate one comet through both boundary legs -> row dict.

    Returns (row, None) on success or (None, reason) on failure.
    rec carries name, yr, q, e, i, om, w, epoch, tp, cc, arc, nobs, ng.
    """
    try:
        r0, v0, et0 = comet_state(rec)
        rb = integrate_leg(r0, v0, et0, -1)
        rf = integrate_leg(r0, v0, et0, +1)
    except Exception as exc:
        return None, str(exc)
    if rb is None or rf is None:
        return None, "boundary not reached"
    p_osc = perih_dir(math.radians(rec["w"]), math.radians(rec["om"]),
                      math.radians(rec["i"]))
    drot = sep(rb["phat"], rf["phat"])
    d_in = sep(rb["phat"], p_osc)
    d_out = sep(rf["phat"], p_osc)
    aph = -rb["phat"]
    theta = sep(aph, axis)
    a_loc = 1e6 / rb["aa"] if rb["aa"] != 0 else float("inf")
    e_loc = 1.0 - rec["q"] / a_loc
    h = math.sqrt(MU * rec["q"] * (1.0 + e_loc))
    om_b = h / R_STOP ** 2
    row = dict(
        desig=rec["name"], yr=rec["yr"], q=rec["q"], i=rec["i"],
        e=rec["e"], cc=rec["cc"], arc=rec["arc"], nobs=rec["nobs"],
        ng=rec["ng"],
        theta=theta, drot=drot, d_in=d_in, d_out=d_out,
        daa=rf["aa"] - rb["aa"],
        denc=min(rb["denc"], rf["denc"]),
        aa_back=rb["aa"], aa_fwd=rf["aa"],
        t_back=rb["t_years"], t_fwd=rf["t_years"],
        aph_lb=[float(x) for x in np.round(aph, 6)],
        dtau=math.radians(drot) / om_b,
        dtau_in=math.radians(d_in) / om_b,
        dtau_out=math.radians(d_out) / om_b)
    return row, None
