#!/usr/bin/env python3
"""step_156: boundary-radius forward-model scan (result b122).

The radial-shell record (step_076, b41) finds the in-cap rotation
excess already fully developed at the 8 AU shell and flat outward,
while the putative physical boundary sits at ~50-150 AU.  The
open question is whether that is a contradiction: a crossing slip
injected far outside the observed arc is absorbed by the orbit fit
and redistributed, so a fitted angular offset could in principle be
present at 8 AU even when the physical event sat at 100+ AU.  This
step measures that directly.

For each dual-leg cohort comet (the step_151 sample) a position-slip
holonomy (the realization step_151/step_145 select) is injected into
the truth trajectory at a controlled inbound crossing radius
R_INJ in {8, 25, 50, 100, 150, 250} AU.  Synthetic astrometry is
generated on the REAL observing chain (same MPC epochs, stations and
noise floor) and refit with the identical LM/DE440s machinery.  The
fitted solutions are then propagated through the shell set and the
reconstruction error is resolved per shell:

  d_fit_vs_true_in[s]   all-fit backward leg vs true PRE-slip inbound
                        leg, osc-periapsis separation at shell s
  d_fit_vs_true_out[s]  all-fit forward leg vs true POST-slip outbound
                        leg at shell s
  d_legfits[s]          in-leg fit (backward) vs out-leg fit
                        (forward) -- the synthetic analogue of the
                        CODE orig-vs-fut catalogue discrepancy

If a slip injected at R_INJ = 100-250 AU produces fitted offsets that
are already fully developed at s = 8 AU, the observed flat-at-8
profile is consistent with a distant boundary; if the offset only
develops for s >= R_INJ, the observed profile would require the
imprint interior to 8 AU.

Also recorded per comet x radius: fitted element shifts, boundary
products (drot, d_in, d_out, daa) and leg disagreement (ddirf), the
same transfer vector as step_151 -- so the radius dependence of the
fitted-observable amplitudes is measured on one construction.

Outputs:
  results/step_b122_boundary_radius.json
  results/step_b122_boundary_radius.csv
  results/figures/supplementary/step_b122_boundary_radius.png
"""

import csv
import json
import math
import sys as _sys
from pathlib import Path as _Path

import numpy as np
import rebound
import spiceypy as sp

_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, sep, tee_stdout
from scripts.utils.parallel import cli_workers as _cli_workers

logger = StepLogger("step_156_boundary_radius_scan")
tee_stdout(logger)
logger.header("Boundary-radius forward-model scan: "
              "slip injected at R_INJ in {8,25,50,100,150,250} AU, "
              "refit through the real observing chain, shell-resolved")

from scripts.utils import mpc_refit as R
from scripts.utils.lpc_boundary import (
    AU_KM, DAY_YR, GM, GM_SUN, PLANET_IDS, R_STOP, MU,
    body_state, init_sim, boundary_orbit, integrate_leg)

SPK = DATA_RAW / "spice" / "de440s.bsp"
LSK = DATA_RAW / "naif" / "naif0012.tls"
MPC_DIR = DATA_RAW / "mpc"
OBS_DIR = MPC_DIR / "obs"
SBDB_DIR = MPC_DIR / "sbdb_fp"
OBSC_PATH = MPC_DIR / "obscodes.json"
PROV_PATH = MPC_DIR / "provenance.json"
CKPT = RESULTS / "step_b122_radius_scan.jsonl"

SEED = 20260920
N_IN, N_OUT = 10, 10
MIN_OBS_LEG = R.MIN_OBS_LEG
MIN_LEG_SPAN_D = R.MIN_LEG_SPAN_D
AS_RAD = R.AS_RAD
RX = R.RX
C_AU_DAY = R.C_AU_DAY

R_INJ_LIST = (8.0, 25.0, 50.0, 100.0, 150.0, 250.0)
SHELLS = (8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0,
          60.0, 80.0, 100.0, 150.0, 200.0, 250.0)

WORKERS = _cli_workers(_sys.argv)
LIMIT = None
if "--limit" in _sys.argv:
    _i = _sys.argv.index("--limit")
    LIMIT = int(_sys.argv[_i + 1])
REDO = "--redo" in _sys.argv

R.configure(OBS_DIR, SBDB_DIR, OBSC_PATH, PROV_PATH, SPK, LSK,
            "step_156_boundary_radius_scan",
            extra_spks=(DATA_RAW / "spice" / "de430.bsp",))

SEC_YR = 86400.0 * DAY_YR
T_MAX = 20000.0
DT_OUT = 1.0
R_MIN_ENC = 0.001
DT_MIN_YR = 1e-8


# ------------------------------------------------------------------
# shell-resolved leg integrator (step_076 construction on a Cartesian
# state): march from (r,v,et) outward (direction=-1 past, +1 future),
# recording the osculating periapsis direction at each shell crossing
# ------------------------------------------------------------------
def leg_shells(r0, v0, et0, direction, shells=SHELLS, t_max=T_MAX):
    sim = init_sim(et0)
    p = sim.particles
    ps = np.array([p[0].x, p[0].y, p[0].z])
    vs = np.array([p[0].vx, p[0].vy, p[0].vz])
    sim.add(x=r0[0] + ps[0], y=r0[1] + ps[1], z=r0[2] + ps[2],
            vx=v0[0] + vs[0], vy=v0[1] + vs[1], vz=v0[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = R_MIN_ENC
    hits = {}
    nxt = 0
    t = direction * DT_OUT
    while abs(t) < t_max and nxt < len(shells):
        sim.integrate(t, exact_finish_time=0)
        if abs(sim.dt) < DT_MIN_YR:
            raise RuntimeError("IAS15 timestep collapse")
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y,
                          p[nc].z - p[0].z])
        rr = np.linalg.norm(r_rel)
        if rr >= shells[nxt]:
            mtot = sum(pp.m for pp in sim.particles)
            rb = np.zeros(3); vb = np.zeros(3)
            for pp in sim.particles:
                rb += pp.m * np.array([pp.x, pp.y, pp.z])
                vb += pp.m * np.array([pp.vx, pp.vy, pp.vz])
            rb /= mtot; vb /= mtot
            r_b = np.array([p[nc].x, p[nc].y, p[nc].z]) - rb
            v_b = np.array([p[nc].vx, p[nc].vy, p[nc].vz]) - vb
            phat, aa = boundary_orbit(r_b, v_b, mtot)
            while nxt < len(shells) and rr >= shells[nxt]:
                hits[shells[nxt]] = (phat, float(t))
                nxt += 1
        t += direction * DT_OUT
    if len(hits) < len(shells):
        return None
    return hits


def crossing_state(r0, v0, et0, r_target):
    """March backward from et0 until |r_rel| first reaches r_target.
    Returns (et_x, r_x, v_x) or None."""
    sim = init_sim(et0)
    p = sim.particles
    ps = np.array([p[0].x, p[0].y, p[0].z])
    vs = np.array([p[0].vx, p[0].vy, p[0].vz])
    sim.add(x=r0[0] + ps[0], y=r0[1] + ps[1], z=r0[2] + ps[2],
            vx=v0[0] + vs[0], vy=v0[1] + vs[1], vz=v0[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = R_MIN_ENC
    t = -DT_OUT
    while abs(t) < T_MAX:
        sim.integrate(t, exact_finish_time=0)
        if abs(sim.dt) < DT_MIN_YR:
            return None
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y,
                          p[nc].z - p[0].z])
        if np.linalg.norm(r_rel) >= r_target:
            return (et0 + t * SEC_YR, r_rel,
                    np.array([p[nc].vx - p[0].vx, p[nc].vy - p[0].vy,
                              p[nc].vz - p[0].vz]))
        t -= DT_OUT
    return None


def prop2b_years(r, v, dt_yr):
    pv = np.r_[r * AU_KM, v * AU_KM / (DAY_YR * 86400.0)]
    pv2 = sp.prop2b(GM_SUN, pv, dt_yr * 86400.0 * DAY_YR)
    return pv2[:3] / AU_KM, pv2[3:] / AU_KM * 86400 * DAY_YR


def slip_propagate(r_t, v_t, et0, ets, et_slip=None, slip_fn=None):
    """Identical construction to step_151: propagate truth to each
    epoch, applying slip_fn when the march passes et_slip."""
    sim = R.init_sim(et0)
    p = sim.particles
    ps = np.array([p[0].x, p[0].y, p[0].z])
    vs = np.array([p[0].vx, p[0].vy, p[0].vz])
    sim.add(x=r_t[0] + ps[0], y=r_t[1] + ps[1], z=r_t[2] + ps[2],
            vx=v_t[0] + vs[0], vy=v_t[1] + vs[1], vz=v_t[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = R_MIN_ENC
    out = np.empty((len(ets), 6))
    order = np.argsort(ets)
    slipped = slip_fn is None or et_slip is None
    for ei in order:
        et = ets[ei]
        if not slipped and et >= et_slip:
            sim.integrate((et_slip - et0) / SEC_YR,
                          exact_finish_time=1)
            rx = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y,
                           p[nc].z - p[0].z])
            vx = np.array([p[nc].vx - p[0].vx, p[nc].vy - p[0].vy,
                           p[nc].vz - p[0].vz])
            rx, vx = slip_fn(rx, vx)
            p[nc].x, p[nc].y, p[nc].z = rx + np.array(
                [p[0].x, p[0].y, p[0].z])
            p[nc].vx, p[nc].vy, p[nc].vz = vx + np.array(
                [p[0].vx, p[0].vy, p[0].vz])
            slipped = True
        sim.integrate((et - et0) / SEC_YR, exact_finish_time=1)
        out[ei] = [p[nc].x - p[0].x, p[nc].y - p[0].y,
                   p[nc].z - p[0].z, p[nc].vx - p[0].vx,
                   p[nc].vy - p[0].vy, p[nc].vz - p[0].vz]
    return out


def observer_ecl(obs):
    out = np.empty((len(obs), 3))
    ok = np.ones(len(obs), bool)
    for i, o in enumerate(obs):
        topo = R.topo_ecl(o["stn"], o["et"])
        if topo is None:
            ok[i] = False
            continue
        out[i] = (R._body_state("399", o["et"])[0]
                  - R._body_state("10", o["et"])[0] + topo)
    return out, ok


def synth_obs(r_t, v_t, et0, obs, sigma_as, rng,
              et_slip=None, slip_fn=None):
    ets = np.array([o["et"] for o in obs])
    robs, ok = observer_ecl(obs)
    st = slip_propagate(r_t, v_t, et0, ets, et_slip, slip_fn)
    rc0 = st[:, :3]
    tau = np.linalg.norm(rc0 - robs, axis=1) / C_AU_DAY
    st2 = slip_propagate(r_t, v_t, et0, ets - tau * 86400.0,
                         et_slip, slip_fn)
    rc = st2[:, :3]
    u = rc - robs
    un = np.linalg.norm(u, axis=1)
    un[un == 0] = 1.0
    u /= un[:, None]
    with np.errstate(all="ignore"):
        u_eq = (RX.T @ u.T).T
    ra_p = np.arctan2(u_eq[:, 1], u_eq[:, 0])
    dec_p = np.arcsin(np.clip(u_eq[:, 2], -1, 1))
    sig = sigma_as / AS_RAD
    syn = []
    for i, o in enumerate(obs):
        if not ok[i]:
            continue
        syn.append(dict(et=o["et"], stn=o["stn"],
                        ra=float(ra_p[i] + rng.normal(0, sig)),
                        dec=float(dec_p[i] + rng.normal(0, sig))))
    return syn


def leg_ok(leg):
    return (len(leg) >= MIN_OBS_LEG and
            (leg[-1]["et"] - leg[0]["et"]) / 86400.0 >= MIN_LEG_SPAN_D)


def refit_synthetic(r_t, v_t, et0, syn, tp_et, et_slip, slip_fn):
    inb = [o for o in syn if o["et"] < tp_et]
    ob = [o for o in syn if o["et"] >= tp_et]
    legs = {"in": inb if leg_ok(inb) else None,
            "out": ob if leg_ok(ob) else None,
            "all": syn if leg_ok(syn) else None}
    fits = {}
    for lab in ("in", "out", "all"):
        leg = legs[lab]
        if leg is None:
            continue
        mid = 0.5 * (leg[0]["et"] + leg[-1]["et"])
        a = slip_propagate(r_t, v_t, et0, np.array([mid]),
                           et_slip, slip_fn)[0]
        rl, vl, rmsl, kl, _ = R.fit(a[:3], a[3:], mid, leg)
        fits[lab] = dict(r=rl, v=vl, et=mid, rms=rmsl,
                         n_keep=int(kl.sum()), n=len(leg))
    return fits


def truth_state(des):
    sb = R.get_sbdb(des)
    el = {e["name"]: float(e["value"])
          for e in sb["orbit"]["elements"] if e.get("value")}
    el["epoch"] = float(sb["orbit"]["epoch"])
    et0 = (el["epoch"] - 2451545.0) * 86400.0
    tp_et = (el["tp"] - 2451545.0) * 86400.0
    st = np.array(sp.conics(
        [el["q"] * AU_KM, el["e"], math.radians(el["i"]),
         math.radians(el["om"]), math.radians(el["w"]),
         0.0, tp_et, GM_SUN], et0))
    return (st[:3] / AU_KM, st[3:] / AU_KM * 86400 * DAY_YR,
            et0, tp_et, el)


def _pos_slip(dtau):
    def _f(r, v):
        return r + v * dtau, v
    return _f


def process_comet(rec):
    des = rec["des"]
    import zlib
    rng = np.random.default_rng(
        SEED + zlib.crc32(des.encode()) % 100000)
    obs = R.parse_obs(R.get_obs(des))
    if len(obs) < 2 * MIN_OBS_LEG:
        return None, f"only {len(obs)} obs"
    r_t, v_t, et0, tp_et, el = truth_state(des)
    sigma = float(np.clip(rec.get("fit_all_rms") or 1.0, 0.3, 5.0))
    dtau = rec.get("our_dtau") or rec.get("cat_dtau") or 20.0

    # shell-resolved TRUE legs: pre-slip inbound (march backward) and
    # reference post-slip forward legs are per-radius below
    try:
        hits_true_in = leg_shells(r_t, v_t, et0, -1)
    except Exception as exc:
        return None, f"truth leg failed: {exc}"
    if hits_true_in is None:
        return None, "truth leg did not reach 250 AU"
    # crossing epochs per radius
    cross = {}
    for r_inj in R_INJ_LIST:
        if r_inj <= el["q"]:
            cross[r_inj] = None
            continue
        cross[r_inj] = crossing_state(r_t, v_t, et0, r_inj)
    row = dict(desig=rec["desig"], des=des, yr=rec["yr"],
               theta=rec["cat_theta"], dtau=dtau, sigma=sigma,
               n_obs=len(obs), radii={})

    for r_inj in R_INJ_LIST:
        cx = cross[r_inj]
        if cx is None:
            row["radii"][str(r_inj)] = {"failed": "no crossing"}
            continue
        et_x, r_x, v_x = cx
        et_slip = et_x
        fn = _pos_slip(dtau)
        try:
            syn = synth_obs(r_t, v_t, et0, obs, sigma, rng,
                            et_slip, fn)
            fits = refit_synthetic(r_t, v_t, et0, syn, tp_et,
                                   et_slip, fn)
        except Exception as exc:
            row["radii"][str(r_inj)] = {"failed": str(exc)[:150]}
            continue
        ent = {"et_slip": et_x, "n_syn": len(syn)}
        # post-slip truth leg: slipped state at et_x marched forward
        try:
            r_s, v_s = fn(np.array(r_x), np.array(v_x))
            hits_true_out = leg_shells(r_s, v_s, et_x, +1)
            if not hits_true_out:
                ent["shell_err"] = "post-slip truth leg incomplete"
        except Exception as exc:
            hits_true_out = None
            ent["shell_err"] = f"post-slip truth leg: {str(exc)[:120]}"
        # shell-resolved fitted legs
        if "all" in fits:
            fa = fits["all"]
            try:
                hb = leg_shells(fa["r"], fa["v"], fa["et"], -1)
                hf = leg_shells(fa["r"], fa["v"], fa["et"], +1)
                if not (hb and hf):
                    ent["shell_err"] = (ent.get("shell_err", "") +
                                        "|fitted leg incomplete")
            except Exception as exc:
                hb = hf = None
                ent["shell_err"] = (ent.get("shell_err", "") +
                                    f"|fitted leg: {str(exc)[:120]}")
            if hb and hf and hits_true_out:
                ent["d_fit_vs_true_in"] = {
                    str(s): round(sep(hb[s][0], hits_true_in[s][0]), 5)
                    for s in SHELLS}
                ent["d_fit_vs_true_out"] = {
                    str(s): round(sep(hf[s][0], hits_true_out[s][0]),
                                  5)
                    for s in SHELLS}
                ent["d_fit_legs"] = {
                    str(s): round(sep(hb[s][0], hf[s][0]), 5)
                    for s in SHELLS}
            ent["all_rms"] = fa["rms"]
        if "in" in fits and "out" in fits:
            fi, fo = fits["in"], fits["out"]
            try:
                hbi = leg_shells(fi["r"], fi["v"], fi["et"], -1)
                hfo = leg_shells(fo["r"], fo["v"], fo["et"], +1)
            except Exception:
                hbi = hfo = None
            if hbi and hfo:
                ent["d_legfits"] = {
                    str(s): round(sep(hbi[s][0], hfo[s][0]), 5)
                    for s in SHELLS}
            # boundary products of the leg fits (ddirf analogue)
            try:
                rbi = integrate_leg(fi["r"], fi["v"], fi["et"], -1)
                rfo = integrate_leg(fo["r"], fo["v"], fo["et"], +1)
                if rbi and rfo:
                    ent["ddirf"] = sep(rbi["phat"], rfo["phat"])
            except Exception:
                pass
        # joint-fit boundary products
        if "all" in fits:
            try:
                fa = fits["all"]
                rb = integrate_leg(fa["r"], fa["v"], fa["et"], -1)
                rf = integrate_leg(fa["r"], fa["v"], fa["et"], +1)
                if rb and rf:
                    ent["drot"] = sep(rb["phat"], rf["phat"])
                    ent["daa"] = rf["aa"] - rb["aa"]
            except Exception:
                pass
        row["radii"][str(r_inj)] = ent
    return row, None


# ------------------------------------------------------------------
# cohort: same dual-leg sample as step_151
# ------------------------------------------------------------------
cohort = []
with open(RESULTS / "step_b91_refit.jsonl") as _fh:
    for line in _fh:
        r = json.loads(line)
        if r.get("failed") or "our_ddirf" not in r:
            continue
        cohort.append(r)
inc = sorted([r for r in cohort if r["cat_theta"] < 60],
             key=lambda r: r["cat_theta"])[:N_IN]
out = sorted([r for r in cohort if r["cat_theta"] >= 60],
             key=lambda r: -r["cat_theta"])[:N_OUT]
todo0 = inc + out
if LIMIT:
    todo0 = todo0[:LIMIT]
logger.info(f"radius-scan cohort: {len(inc)} in-cap + "
            f"{len(out)} out-of-cap dual-leg comets")

done = {}
if CKPT.exists() and not REDO:
    for line in CKPT.read_text().splitlines():
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not r.get("failed"):
            done[r["desig"]] = r
elif CKPT.exists():
    CKPT.unlink()

todo = [r for r in todo0 if r["desig"] not in done]
rows = [done[r["desig"]] for r in todo0 if r["desig"] in done]
logger.info(f"{len(rows)} checkpointed, {len(todo)} to scan")


def _worker_init():
    R.refit_worker_init()


def _process(rec):
    try:
        return process_comet(rec)
    except Exception as exc:
        return None, str(exc)[:200]


def _account(rec, rowres, err):
    def _f(o):
        return float(o) if isinstance(o, np.floating) else str(o)
    with open(CKPT, "a") as ck:
        if rowres is None:
            ck.write(json.dumps({"desig": rec["desig"],
                                 "failed": True,
                                 "err": str(err)[:200]}) + "\n")
        else:
            ck.write(json.dumps(rowres, default=_f) + "\n")
            rows.append(rowres)


if todo:
    if WORKERS > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(WORKERS, initializer=_worker_init) as pool:
            for i, (rec, (rr, err)) in enumerate(zip(
                    todo, pool.imap(_process, todo))):
                _account(rec, rr, err)
                logger.info(f"[{i+1}/{len(todo)}] {rec['desig']} "
                            f"{'ok' if rr else err}")
    else:
        for i, rec in enumerate(todo):
            rr, err = _process(rec)
            _account(rec, rr, err)
            logger.info(f"[{i+1}/{len(todo)}] {rec['desig']} "
                        f"{'ok' if rr else err}")

# ------------------------------------------------------------------
# aggregate: median fitted-error profile per shell per radius
# ------------------------------------------------------------------
agg = {}
for r_inj in R_INJ_LIST:
    key = str(r_inj)
    prof = {"d_fit_vs_true_in": {str(s): [] for s in SHELLS},
            "d_fit_vs_true_out": {str(s): [] for s in SHELLS},
            "d_fit_legs": {str(s): [] for s in SHELLS},
            "d_legfits": {str(s): [] for s in SHELLS},
            "drot": [], "ddirf": [], "n": 0}
    n_disrupted = 0
    for row in rows:
        ent = row["radii"].get(key)
        if not ent or "failed" in ent:
            continue
        prof["n"] += 1
        if ent.get("shell_err"):
            n_disrupted += 1
        for chan in ("d_fit_vs_true_in", "d_fit_vs_true_out",
                     "d_fit_legs", "d_legfits"):
            if chan in ent:
                for s in SHELLS:
                    prof[chan][str(s)].append(ent[chan][str(s)])
        if "drot" in ent:
            prof["drot"].append(ent["drot"])
        if "ddirf" in ent:
            prof["ddirf"].append(ent["ddirf"])
    out_p = {"n": prof["n"], "n_disrupted": n_disrupted}
    for chan in ("d_fit_vs_true_in", "d_fit_vs_true_out",
                 "d_fit_legs", "d_legfits"):
        out_p[chan] = {
            s: {"med": float(np.median(v)),
                "p16": float(np.percentile(v, 16)),
                "p84": float(np.percentile(v, 84)), "n": len(v)}
            for s, v in prof[chan].items() if v}
    if prof["drot"]:
        out_p["drot_med"] = float(np.median(prof["drot"]))
    if prof["ddirf"]:
        out_p["ddirf_med"] = float(np.median(prof["ddirf"]))
    agg[key] = out_p

# the decisive comparison: fitted-inbound error at the innermost shell
# vs at the injection radius -- does the offset develop below R_INJ?
radial_transfer = {}
for r_inj in R_INJ_LIST:
    key = str(r_inj)
    p = agg[key].get("d_fit_vs_true_in", {})
    if str(8.0) in p and key in p:
        s_in = "8.0"
        # shell nearest the injection radius
        s_at = str(min(SHELLS, key=lambda s: abs(s - r_inj)))
        if s_at in p:
            radial_transfer[key] = {
                "med_at_8AU": p[s_in]["med"],
                "med_at_rinj_shell": p[s_at]["med"],
                "ratio_8AU_to_rinj": (p[s_in]["med"] /
                                      p[s_at]["med"]
                                      if p[s_at]["med"] > 0 else
                                      None)}

result = {
    "step": "step_156_boundary_radius_scan",
    "seed": SEED,
    "r_inj_list": list(R_INJ_LIST),
    "shells": list(SHELLS),
    "n_comets": len(rows),
    "aggregate": agg,
    "radial_transfer": radial_transfer,
    "note": ("position-slip holonomy injected at controlled inbound "
             "crossing radii; synthetic astrometry on the real MPC "
             "chain; identical LM/DE440s refit; fitted error resolved "
             "per boundary shell.  d_fit_vs_true_in[s] measures where "
             "the fitted orbit's implied inbound periapsis direction "
             "departs from the true pre-slip leg"),
}
RESULTS.mkdir(exist_ok=True)
(RESULTS / "step_b122_boundary_radius.json").write_text(
    json.dumps(result, indent=1))
logger.data_save(RESULTS / "step_b122_boundary_radius.json")

# flat csv
with open(RESULTS / "step_b122_boundary_radius.csv", "w",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["desig", "theta", "r_inj", "shell",
                "d_fit_vs_true_in", "d_fit_vs_true_out",
                "d_legfits", "drot", "ddirf"])
    for row in rows:
        for rk, ent in row["radii"].items():
            if "failed" in ent:
                continue
            for s in SHELLS:
                w.writerow([row["desig"], row["theta"], rk, s,
                            ent.get("d_fit_vs_true_in", {}).get(
                                str(s), ""),
                            ent.get("d_fit_vs_true_out", {}).get(
                                str(s), ""),
                            ent.get("d_legfits", {}).get(str(s), ""),
                            ent.get("drot", ""),
                            ent.get("ddirf", "")])
logger.data_save(RESULTS / "step_b122_boundary_radius.csv")

logger.info("aggregate median d_fit_vs_true_in (deg) per shell:")
hdr = "  r_inj | " + " ".join(f"{s:6.0f}" for s in SHELLS[:8])
logger.info(hdr)
for r_inj in R_INJ_LIST:
    p = agg[str(r_inj)].get("d_fit_vs_true_in", {})
    line = "  " + f"{r_inj:5.0f} | " + " ".join(
        f"{p[str(s)]['med']:6.3f}" if str(s) in p else "     -"
        for s in SHELLS[:8])
    logger.info(line)
logger.info(f"wrote results/step_b122_boundary_radius.json "
            f"({len(rows)} comets x {len(R_INJ_LIST)} radii)")

# ---------------- figure ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.utils.tep9_style import apply_style
apply_style()

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axes[0]
cmap = plt.cm.viridis(np.linspace(0, 0.9, len(R_INJ_LIST)))
for c, r_inj in zip(cmap, R_INJ_LIST):
    p = agg[str(r_inj)].get("d_fit_vs_true_in", {})
    xs = [float(s) for s in p]
    ys = [p[s]["med"] for s in p]
    ax.plot(xs, ys, "o-", color=c, ms=3, lw=1,
            label=f"$R_{{inj}}$ = {r_inj:.0f} AU")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("shell radius s (AU)")
ax.set_ylabel("median fitted-vs-true inbound offset (deg)")
ax.set_title("Where the absorbed slip enters the fitted orbit")
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

ax = axes[1]
for c, r_inj in zip(cmap, R_INJ_LIST):
    p = agg[str(r_inj)].get("d_legfits", {})
    xs = [float(s) for s in p]
    ys = [p[s]["med"] for s in p]
    ax.plot(xs, ys, "o-", color=c, ms=3, lw=1,
            label=f"{r_inj:.0f} AU")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("shell radius s (AU)")
ax.set_ylabel("median leg-fit disagreement (deg)")
ax.set_title("orig/fut analogue: leg-fit discrepancy per shell")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(RESULTS / "figures/supplementary/"
            "step_b122_boundary_radius.png", dpi=200)
logger.data_save(RESULTS / "figures/supplementary/"
                 "step_b122_boundary_radius.png")
logger.info("figure written")
