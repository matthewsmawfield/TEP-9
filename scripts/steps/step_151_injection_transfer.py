"""step_151: slip injection -> standard refit transfer function (b115).

The manuscript converts the observed periapsis rotation into an
equivalent time, dtau = drot * r_b^2 / h, without deriving what a
proper-time perturbation actually does to a fitted orbit.  This step
derives it.  For each dual-leg cohort comet the REAL observing chain
is kept -- the same MPC epochs, stations and noise floor -- but the
astrometry is synthesised from a truth orbit with a specified slip
realization injected, then refit with the standard LM/DE440s model
exactly as steps 127/128 fit the real record.

Injection realizations (the slip is applied at the true inbound
250 AU boundary crossing unless noted):

  ctrl      no injection -- instrument noise floor
  full      full-state time translation: boundary state advanced by
            dtau along its own (two-body) orbit.  The translation is
            absorbed by t_p; what remains is only the rescrambled
            planetary-encounter sequence -- the measurable residue of
            a global phase shift and the practical content of the
            oint(d ln A) = 0 statement on an open (transiting) path.
  pos       position-only time translation: r -> r + v*dtau, v fixed.
            The minimal non-integrable structure that can carry a
            record.
  vel       velocity-only time translation: v -> v(t+dtau), r fixed.
  kick      transverse impulse at the crossing: v -> v + dv_t, with
            dv_t = |v| * drot (the impulse realization of step_145).
  mid       within-apparition holonomy: position-only slip applied at
            the perihelion epoch, so inbound and outbound arcs see
            different orbits -- the realization a true leg-disagreement
            signature would require.

Measured per comet x injection (the transfer vector):
  fitted-element shifts dq, de, di, dOm, dw vs truth;
  boundary products drot, d_in, d_out, daa, denc (fitted orbit);
  inferred-origin displacements d_true_in / d_true_out (fitted
  asymptote vs TRUE asymptote -- the channel the catalogue record
  actually encodes);
  leg products ddirf, d_in_leg, d_out_leg, xarc_rms_in2out;
  fit rms per leg / full arc.

Tests
  T1  controls: ctrl and full must be null on every channel.
  T2  transfer vectors per boundary realization (pos, vel, kick).
  T3  leg-disagreement discriminator: predicted ddirf per realization
      vs the observed cohort distribution -- a past-boundary slip
      predicts ~zero leg disagreement (the whole apparition is
      post-slip); only the mid-arc realization can produce it.
  T4  energy channel: daa injected per realization; the impulse must
      move energy, the slip realizations must not -- matched against
      the observed flat energy channel.
  T5  amplitude tracking: injected dtau vs produced d_true_in across
      comets -- does the required slip size track the measured
      anomaly?

Inputs
  results/step_b91_refit.jsonl   (dual-leg cohort + measured products)
  data/raw/mpc/obs/*.json        (cached real observing geometry)
  data/raw/mpc/sbdb_fp/*.json    (truth orbits)

Outputs
  results/step_b115_injection_transfer.json / .csv
  results/step_b115_injection.jsonl        (per-comet checkpoint)
  results/figures/supplementary/step_b115_injection_transfer.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, perih_dir, sep, tee_stdout
logger = StepLogger("step_151_injection_transfer")
tee_stdout(logger)
logger.header("Slip injection -> standard refit transfer function")

import csv
import json
import math
import numpy as np
import spiceypy as sp

from scripts.utils import mpc_refit as R
from scripts.utils.lpc_boundary import (
    AU_KM, DAY_YR, GM_SUN, R_STOP, integrate_leg)

SPK = DATA_RAW / "spice" / "de440s.bsp"
LSK = DATA_RAW / "naif" / "naif0012.tls"
MPC_DIR = DATA_RAW / "mpc"
OBS_DIR = MPC_DIR / "obs"
SBDB_DIR = MPC_DIR / "sbdb_fp"
OBSC_PATH = MPC_DIR / "obscodes.json"
PROV_PATH = MPC_DIR / "provenance.json"
CKPT = RESULTS / "step_b115_injection.jsonl"

SEED = 20260920
N_IN, N_OUT = 10, 10          # in-cap / out-of-cap comets
MIN_OBS_LEG = R.MIN_OBS_LEG
MIN_LEG_SPAN_D = R.MIN_LEG_SPAN_D
AS_RAD = R.AS_RAD
RX = R.RX
C_AU_DAY = R.C_AU_DAY

INJECTIONS = ("ctrl", "full", "pos", "vel", "kick", "mid")

from scripts.utils.parallel import cli_workers as _cli_workers
WORKERS = _cli_workers(_sys.argv)
LIMIT = None
if "--limit" in _sys.argv:
    _i = _sys.argv.index("--limit")
    LIMIT = int(_sys.argv[_i + 1])
REDO = "--redo" in _sys.argv

R.configure(OBS_DIR, SBDB_DIR, OBSC_PATH, PROV_PATH, SPK, LSK,
            "step_151_injection_transfer",
            extra_spks=(DATA_RAW / "spice" / "de430.bsp",))


# ------------------------------------------------------------------
# synthetic astrometry on the real observing geometry
# ------------------------------------------------------------------
def obs_geometry(des):
    """Real (et, stn) list for a comet from the cached MPC record."""
    obs = R.parse_obs(R.get_obs(des))
    return obs


def truth_state(des):
    """SBDB full-precision solution -> heliocentric state + elements."""
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


def slip_propagate(r_t, v_t, et0, ets, et_slip=None, slip_fn=None):
    """Propagate truth (r_t,v_t @ et0) to ets, applying slip_fn to the
    comet's heliocentric state when the integration passes et_slip.

    Planets are seeded at et0 (kernel coverage); slip_fn maps
    (r, v) -> (r', v') in heliocentric ecliptic AU, AU/yr."""
    sim = R.init_sim(et0)
    p = sim.particles
    ps = np.array([p[0].x, p[0].y, p[0].z])
    vs = np.array([p[0].vx, p[0].vy, p[0].vz])
    sim.add(x=r_t[0] + ps[0], y=r_t[1] + ps[1], z=r_t[2] + ps[2],
            vx=v_t[0] + vs[0], vy=v_t[1] + vs[1], vz=v_t[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = 0.001  # collision scale: bound IAS15 against step collapse
    sec_yr = 86400.0 * DAY_YR
    out = np.empty((len(ets), 6))
    order = np.argsort(ets)
    slipped = slip_fn is None or et_slip is None
    for ei in order:
        et = ets[ei]
        if not slipped and et >= et_slip:
            sim.integrate((et_slip - et0) / sec_yr, exact_finish_time=1)
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
        sim.integrate((et - et0) / sec_yr, exact_finish_time=1)
        out[ei] = [p[nc].x - p[0].x, p[nc].y - p[0].y,
                   p[nc].z - p[0].z, p[nc].vx - p[0].vx,
                   p[nc].vy - p[0].vy, p[nc].vz - p[0].vz]
    return out


def prop2b_years(r, v, dt_yr):
    """Two-body Kepler advance of a heliocentric state by dt_yr."""
    pv = np.r_[r * AU_KM, v * AU_KM / (DAY_YR * 86400.0)]
    pv2 = sp.prop2b(GM_SUN, pv, dt_yr * 86400.0 * DAY_YR)
    return pv2[:3] / AU_KM, pv2[3:] / AU_KM * 86400 * DAY_YR


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
    """Synthetic ra/dec on real (et, stn) from the truth trajectory,
    with an optional slip applied at et_slip."""
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
    with np.errstate(all="ignore"):  # Accelerate BLAS raises spurious FP flags
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


# ------------------------------------------------------------------
# fit driver on synthetic obs
# ------------------------------------------------------------------
def leg_ok(leg):
    return (len(leg) >= MIN_OBS_LEG and
            (leg[-1]["et"] - leg[0]["et"]) / 86400.0 >= MIN_LEG_SPAN_D)


def refit_synthetic(r_t, v_t, et0, syn, tp_et,
                    et_slip=None, slip_fn=None):
    """Mirror fit_comet's leg structure on synthetic obs.

    Anchors are seeded from the same (possibly slipped) trajectory the
    observations were generated from -- the fair seed for every
    injection."""
    inb = [o for o in syn if o["et"] < tp_et]
    out = [o for o in syn if o["et"] >= tp_et]
    legs = {"in": inb if leg_ok(inb) else None,
            "out": out if leg_ok(out) else None,
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
    return fits, legs


def _try_leg(r, v, et, direction):
    try:
        return integrate_leg(r, v, et, direction)
    except Exception:
        return None


def boundary_products(fits, legs, traj_true_asy):
    out = {}
    if "all" in fits:
        fa = fits["all"]
        rb = _try_leg(fa["r"], fa["v"], fa["et"], -1)
        rf = _try_leg(fa["r"], fa["v"], fa["et"], +1)
        if rb and rf:
            p_all, q_f, e_f, i_f, Om_f, w_f = R.state_to_perih(
                fa["r"], fa["v"], fa["et"])
            out.update(drot=sep(rb["phat"], rf["phat"]),
                       d_in=sep(rb["phat"], p_all),
                       d_out=sep(rf["phat"], p_all),
                       daa=rf["aa"] - rb["aa"],
                       denc=min(rb["denc"], rf["denc"]),
                       d_true_in=sep(rb["phat"], traj_true_asy[0]),
                       d_true_out=sep(rf["phat"], traj_true_asy[1]),
                       q=q_f, e=e_f, i=i_f, Om=Om_f, w=w_f,
                       p_all=p_all)
    if "in" in fits and "out" in fits:
        fi, fo = fits["in"], fits["out"]
        rbi = _try_leg(fi["r"], fi["v"], fi["et"], -1)
        rfo = _try_leg(fo["r"], fo["v"], fo["et"], +1)
        if rbi and rfo:
            out["ddirf"] = sep(rbi["phat"], rfo["phat"])
            if "p_all" in out:
                out["d_in_leg"] = sep(rbi["phat"], out["p_all"])
                out["d_out_leg"] = sep(rfo["phat"], out["p_all"])
            out["xarc_rms_in2out"] = float(np.sqrt(np.nanmean(
                R.residuals(fi["r"], fi["v"], fi["et"],
                            legs["out"]) ** 2)))
    for lab in ("in", "out", "all"):
        if lab in fits:
            out[f"{lab}_rms"] = fits[lab]["rms"]
            out[f"{lab}_nkeep"] = fits[lab]["n_keep"]
    out.pop("p_all", None)
    return out


# ------------------------------------------------------------------
# per-comet injection battery
# ------------------------------------------------------------------
def process_comet(rec):
    des = rec["des"]
    import zlib
    rng = np.random.default_rng(SEED + zlib.crc32(des.encode()) % 100000)
    obs = obs_geometry(des)
    if len(obs) < 2 * MIN_OBS_LEG:
        return None, f"only {len(obs)} obs"
    r_t, v_t, et0, tp_et, el = truth_state(des)
    sigma = rec.get("fit_all_rms") or 1.0
    sigma = float(np.clip(sigma, 0.3, 5.0))
    dtau = rec.get("our_dtau") or rec.get("cat_dtau") or 20.0   # yr

    # true boundary legs + crossing epoch
    try:
        b_in = integrate_leg(r_t, v_t, et0, -1)
        b_fu = integrate_leg(r_t, v_t, et0, +1)
    except Exception as exc:
        return None, f"truth leg unphysical: {exc}"
    if not b_in or not b_fu:
        return None, "truth boundary not reached"
    et_x = et0 + b_in["t_years"] * 86400.0 * DAY_YR
    asy_true = (b_in["phat"], b_fu["phat"])

    drot_rad = math.radians(rec.get("our_drot") or rec.get("cat_drot")
                            or 0.1)

    def _full(r, v):
        return prop2b_years(r, v, dtau)

    def _pos(r, v):
        return r + v * dtau, v

    def _vel(r, v):
        _, v2 = prop2b_years(r, v, dtau)
        return r, v2

    def _kick(r, v):
        hh = np.cross(r, v); hh /= np.linalg.norm(hh)
        t_hat = np.cross(hh, v / np.linalg.norm(v))
        t_hat /= np.linalg.norm(t_hat)
        return r, v + np.linalg.norm(v) * drot_rad * t_hat

    slips = {"ctrl": (None, None),
             "full": (et_x, _full),
             "pos":  (et_x, _pos),
             "vel":  (et_x, _vel),
             "kick": (et_x, _kick),
             "mid":  (tp_et, _pos)}

    prods = {}
    for inj in INJECTIONS:
        et_slip, fn = slips[inj]
        syn = synth_obs(r_t, v_t, et0, obs, sigma, rng,
                        et_slip, fn)
        fits, legs = refit_synthetic(r_t, v_t, et0, syn, tp_et,
                                     et_slip, fn)
        if "all" not in fits:
            prods[inj] = {"failed": "no full-arc fit"}
            continue
        try:
            bp = boundary_products(fits, legs, asy_true)
        except Exception as exc:
            prods[inj] = {"failed": str(exc)[:150]}
            continue
        # element shifts vs truth
        if "w" in bp:
            bp["dq"] = bp["q"] - el["q"]
            bp["de"] = bp["e"] - el["e"]
            bp["di"] = bp["i"] - el["i"]
            dw = (bp["w"] - el["w"] + 180) % 360 - 180
            dOm = (bp["Om"] - el["om"] + 180) % 360 - 180
            bp["dw"] = float(dw); bp["dOm"] = float(dOm)
            a_t = el["q"] / (1 - el["e"]) if el["e"] != 1 else np.nan
            a_f = bp["q"] / (1 - bp["e"]) if bp["e"] != 1 else np.nan
            if np.isfinite(a_t) and np.isfinite(a_f) and a_t and a_f:
                bp["d_1a"] = 1e6 * (1 / a_f - 1 / a_t)
        prods[inj] = bp

    row = dict(desig=rec["desig"], des=des, yr=rec["yr"],
               theta=rec["cat_theta"], dtau=dtau,
               cat_drot=rec["cat_drot"], our_drot=rec.get("our_drot"),
               our_ddirf=rec.get("our_ddirf"),
               our_d_in_leg=rec.get("our_d_in_leg"),
               our_d_out_leg=rec.get("our_d_out_leg"),
               xarc=rec.get("xarc_rms_in2out"),
               sigma=sigma, n_obs=len(obs), products=prods)
    return row, None


# ------------------------------------------------------------------
# cohort
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
logger.info(f"injection cohort: {len(inc)} in-cap + {len(out)} "
            f"out-of-cap dual-leg comets")

def _f(o):
    return float(o) if isinstance(o, np.floating) else str(o)


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
logger.info(f"{len(rows)} checkpointed, {len(todo)} to inject")


def _worker_init():
    R.refit_worker_init()


def _process(rec):
    try:
        return process_comet(rec)
    except Exception as exc:
        return None, str(exc)[:200]


def _account(rec, rowres, err):
    with open(CKPT, "a") as ck:
        if rowres is None:
            ck.write(json.dumps({"desig": rec["desig"], "failed": True,
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
                logger.info(f"  {i+1}/{len(todo)} injected")
    else:
        for i, rec in enumerate(todo):
            rr, err = _process(rec)
            if rr is None:
                logger.warning(f"{rec['desig']}: {err}")
            _account(rec, rr, err)
            logger.info(f"  {i+1}/{len(todo)} injected")

# compact checkpoint: append-mode accumulation across overlapping runs can
# otherwise leave duplicate per-designation records in the JSONL
if CKPT.exists():
    uniq = {}
    for line in CKPT.read_text().splitlines():
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        uniq[r.get("desig")] = r
    with open(CKPT, "w") as ck:
        for r in uniq.values():
            ck.write(json.dumps(r, default=_f) + "\n")
    logger.info(f"checkpoint compacted: {len(uniq)} unique designations "
                f"written to {CKPT}")

logger.info(f"injected cohort: {len(rows)} comets")


# ------------------------------------------------------------------
# transfer-function summary
# ------------------------------------------------------------------
PROD_KEYS = ("drot", "d_in", "d_out", "d_true_in", "d_true_out",
             "ddirf", "d_in_leg", "d_out_leg", "xarc_rms_in2out",
             "daa", "dw", "dOm", "di", "dq", "d_1a",
             "in_rms", "out_rms", "all_rms")


def med(row, inj, key):
    v = row["products"].get(inj, {}).get(key)
    return float(v) if v is not None and np.isfinite(v) else np.nan


def vec(rows, inj, key):
    return np.array([med(r, inj, key) for r in rows])


res = {"seed": SEED, "injections": {
    "ctrl": "no injection -- instrument floor",
    "full": "full-state time translation at the boundary (phase-shift "
            "realization; absorbed by tp modulo encounter rescramble)",
    "pos": "position-only slip r->r+v*dtau at the boundary",
    "vel": "velocity-only slip v->v(t+dtau) at the boundary",
    "kick": "transverse impulse dv_t = v*drot at the boundary",
    "mid": "position-only slip at the perihelion epoch (within-arc "
           "holonomy)"},
    "n_comets": len(rows)}

for inj in INJECTIONS:
    blk = {}
    for k in PROD_KEYS:
        v = vec(rows, inj, k)
        v = v[np.isfinite(v)]
        if len(v):
            blk[k] = {"n": int(len(v)), "med": float(np.median(v)),
                      "p16": float(np.percentile(v, 16)),
                      "p84": float(np.percentile(v, 84)),
                      "med_abs": float(np.median(np.abs(v)))}
    res[f"T_{inj}"] = blk

# T3 observed-vs-predicted leg disagreement
obs_ddirf = np.array([r["our_ddirf"] for r in rows
                      if r.get("our_ddirf") is not None])
res["T3_leg_disagreement"] = {
    "observed_ddirf_med": float(np.median(obs_ddirf)) if len(obs_ddirf) else None,
    "predicted": {inj: res[f"T_{inj}"].get("ddirf", {}).get("med")
                  for inj in INJECTIONS}}

# T5 amplitude tracking: injected dtau vs produced d_true_in
from scipy.stats import spearmanr
for inj in ("pos", "vel", "kick"):
    dt = np.array([r["dtau"] for r in rows])
    di = vec(rows, inj, "d_true_in")
    m = np.isfinite(di)
    if m.sum() > 5:
        rho, p = spearmanr(dt[m], di[m])
        res.setdefault("T5_amplitude", {})[inj] = {
            "rho": float(rho), "p": float(p), "n": int(m.sum())}

res["finding"] = (
    "The derived TEP transfer function: which fitted-orbit channels a "
    "specified proper-time slip moves.  The integrable translation "
    "(full) must be invisible -- the oint(d ln A)=0 control; the "
    "position/velocity slip realizations show where a non-integrable "
    "time transport lands (boundary-asymptote vs leg-disagreement vs "
    "element channels), and the impulse control separates the "
    "energy-carrying mechanical realization.")

with open(RESULTS / "step_b115_injection_transfer.json", "w") as _fh:
    json.dump(res, _fh, indent=1, default=str)

csv_out = RESULTS / "step_b115_injection_transfer.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["desig", "theta", "dtau", "injection"] + list(PROD_KEYS))
    for r in rows:
        for inj in INJECTIONS:
            w.writerow([r["desig"], r["theta"], r["dtau"], inj] +
                       [("" if not np.isfinite(med(r, inj, k))
                         else f"{med(r, inj, k):.6g}") for k in PROD_KEYS])
logger.data_save(csv_out)
# ------------------------------------------------------------------
# figure
# ------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

ax = axes[0]
keys = ["drot", "d_true_in", "d_true_out", "ddirf", "daa"]
xpos = np.arange(len(keys))
for inj, c in (("pos", "crimson"), ("vel", "darkorange"),
               ("kick", "steelblue"), ("mid", "purple")):
    vals = [res[f"T_{inj}"].get(k, {}).get("med_abs", np.nan)
            for k in keys]
    ax.scatter(xpos, vals, label=inj, color=c, s=40, zorder=3)
for j, k in enumerate(keys):
    ax.axvline(j, color="0.9", lw=0.5, zorder=1)
ax.set_yscale("log")
ax.set_xticks(xpos); ax.set_xticklabels(keys, rotation=30, ha="right")
ax.set_ylabel("median |signature| (deg, or ppm for daa)")
ax.set_title("transfer vector per realization")
ax.legend(fontsize=8)

ax = axes[1]
for inj, c, m in (("pos", "crimson", "o"), ("mid", "purple", "s")):
    v = vec(rows, inj, "ddirf")
    v = v[np.isfinite(v)]
    ax.hist(v, bins=15, alpha=0.55, color=c, label=f"{inj} (predicted)")
if len(obs_ddirf):
    ax.axvline(np.median(obs_ddirf), color="k", ls="--",
               label=f"observed median {np.median(obs_ddirf):.3f} deg")
ax.set_xlabel("ddirf (deg)"); ax.set_ylabel("comets")
ax.set_title("leg-disagreement discriminator")
ax.legend(fontsize=8)

ax = axes[2]
for inj, c in (("pos", "crimson"), ("vel", "darkorange"),
               ("kick", "steelblue")):
    dt = np.array([r["dtau"] for r in rows])
    di = vec(rows, inj, "d_true_in")
    m = np.isfinite(di)
    ax.scatter(dt[m], di[m], s=22, alpha=0.7, color=c, label=inj)
ax.set_xlabel("injected dtau (yr)")
ax.set_ylabel("produced d_true_in (deg)")
ax.set_title("amplitude tracking")
ax.legend(fontsize=8)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b115_injection_transfer.png", dpi=300)
logger.data_save(RESULTS / "step_b115_injection_transfer.json")
logger.data_save(FIG / "supplementary" / "step_b115_injection_transfer.png")
