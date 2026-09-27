"""Shared machinery for independent two-leg orbit refits on raw MPC
astrometry (steps 127+).

The module provides everything needed to rebuild a comet boundary
record without any catalogue fit products:

  * raw-observation retrieval from the MPC get-obs service (ADES_DF)
    with on-disk caching and provenance;
  * full-precision JPL SBDB orbit retrieval used purely as an
    integration seed, likewise cached;
  * the MPC observatory-code table (cached) for topocentric parallax;
  * a Levenberg-Marquardt differential-correction fitter for the six
    Cartesian state parameters at an anchor epoch, under the same
    REBOUND/IAS15 + DE440s force model the boundary integrator uses,
    with one-way light-time correction and a MAD outlier clip;
  * propagation of fitted states to the +/-250 AU barycentric sphere
    through the shared integrate_leg machinery.

A step calls ``configure(paths...)`` once at start-up (furnishes the
SPICE kernels, creates the cache directories, and loads the
observatory table), then calls ``fit_comet(desig)`` per comet.  The
returned dict carries the fitted states, the SBDB seed elements, and
every boundary product both constructions support:

  single-fit construction (any full-arc fit):
      drot, d_in, d_out, daa, denc, aph, dtau
  leg-fit construction (requires both legs fittable):
      ddirf, xarc_rms_in2out, d_in_leg, d_out_leg

``d_in`` is the pipeline's registered carrier convention: the
separation between the propagated inbound boundary asymptote and the
fitted osculating periapsis direction.  ``d_in_leg`` is the
fit-vs-fit analogue available only on dual-leg objects: the
independently fitted inbound leg's boundary asymptote versus the
joint (full-arc) osculating periapsis direction -- the construction a
three-leg catalogue record realizes, and the one a single-solution
fit can never express.
"""

import datetime as _dt
import hashlib
import json
import math
import re
import time
import urllib.parse

import numpy as np
import rebound
import requests
import spiceypy as sp

from scripts.utils.tep9_common import perih_dir, sep
from scripts.utils.lpc_boundary import (
    AU_KM, DAY_YR, GM, GM_SUN, MU, PLANET_IDS, R_STOP, TNO,
    integrate_leg, worker_init as _lpc_worker_init)

# ------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------

MPC_OBS_URL = "https://data.minorplanetcenter.net/api/get-obs"
SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"
OBSC_URL = "https://www.minorplanetcenter.net/iau/lists/ObsCodesF.html"

MIN_OBS_LEG = 12
MIN_LEG_SPAN_D = 5.0
CLIP = 4.0
ITERS = 6

C_AU_DAY = 173.144632674
R_EARTH_AU = 6378.137 / AU_KM
AS_RAD = 206264.806247

EPS = math.radians(23.4392911)
RX = np.array([[1, 0, 0],
               [0, math.cos(EPS), math.sin(EPS)],
               [0, -math.sin(EPS), math.cos(EPS)]])

OBS_DIR = None
SBDB_DIR = None
OBSC_PATH = None
PROV_PATH = None
SPK_PATH = None
LSK_PATH = None
EXTRA_SPKS = ()
OBSC = {}


def configure(obs_dir, sbdb_dir, obsc_path, prov_path, spk, lsk,
              prov_step, extra_spks=()):
    """Point the caches + kernels at the given paths, furnish SPICE,
    and load the observatory table.  Must run before any fit call
    (and before spawning workers under fork).  ``extra_spks`` are
    furnished before the primary SPK so they only serve epochs the
    primary kernel does not cover."""
    global OBS_DIR, SBDB_DIR, OBSC_PATH, PROV_PATH
    global SPK_PATH, LSK_PATH, EXTRA_SPKS, OBSC, _PROV_STEP
    OBS_DIR = obs_dir
    SBDB_DIR = sbdb_dir
    OBSC_PATH = obsc_path
    PROV_PATH = prov_path
    SPK_PATH = spk
    LSK_PATH = lsk
    EXTRA_SPKS = tuple(extra_spks)
    _PROV_STEP = prov_step
    for d in (OBS_DIR, SBDB_DIR):
        d.mkdir(parents=True, exist_ok=True)
    sp.furnsh(str(LSK_PATH))
    for x in EXTRA_SPKS:
        sp.furnsh(str(x))
    sp.furnsh(str(SPK_PATH))
    OBSC = load_obscodes()


def refit_worker_init():
    """Pool initializer: furnish kernels inside a worker process.
    Fork children inherit kernel-table entries but not DAF records,
    so the pool is reset and every kernel is refurnished; extras go
    before the primary SPK so it keeps priority on overlap epochs."""
    sp.kclear()
    for x in EXTRA_SPKS:
        sp.furnsh(str(x))
    sp.furnsh(str(SPK_PATH))
    sp.furnsh(str(LSK_PATH))


# ------------------------------------------------------------------
# provenance + downloads
# ------------------------------------------------------------------

_PROV_STEP = "mpc_refit"


def _prov():
    if PROV_PATH.exists():
        try:
            return json.loads(PROV_PATH.read_text())
        except json.JSONDecodeError:
            # interleaved write from a racing process; recover the
            # first complete document rather than losing the ledger
            try:
                d, _ = json.JSONDecoder().raw_decode(
                    PROV_PATH.read_text().lstrip())
                return d
            except json.JSONDecodeError:
                pass
    return {"step": _PROV_STEP,
            "generated_utc": _dt.datetime.now(
                _dt.timezone.utc).isoformat(),
            "files": {}}


def _prov_write(p):
    tmp = PROV_PATH.with_suffix(".tmp")
    tmp.write_text(json.dumps(p, indent=1))
    tmp.replace(PROV_PATH)


def _prov_add(key, url, data: bytes, note):
    """Append a fetch record under an exclusive file lock so
    concurrent workers cannot interleave the ledger."""
    import fcntl
    lock_path = PROV_PATH.with_suffix(".lock")
    with open(lock_path, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        p = _prov()
        p["files"][key] = {
            "url": url,
            "retrieved_utc": _dt.datetime.now(
                _dt.timezone.utc).isoformat(),
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "note": note}
        _prov_write(p)


def _get(url, binary=False, tries=4, json_payload=None):
    err = None
    for k in range(tries):
        try:
            if json_payload is not None:
                r = requests.get(url, json=json_payload, timeout=120,
                                 headers={"User-Agent": "tep9/1.0"})
            else:
                r = requests.get(url, timeout=120,
                                 headers={"User-Agent": "tep9/1.0"})
            if r.status_code == 200:
                return r.content if binary else r.text
            err = f"HTTP {r.status_code}"
        except Exception as exc:
            err = str(exc)
        time.sleep(1.5 * (k + 1))
    raise RuntimeError(err or "fetch failed")


def load_obscodes():
    if OBSC_PATH.exists():
        return json.loads(OBSC_PATH.read_text())
    html = _get(OBSC_URL)
    m = re.search(r"<pre>(.*?)</pre>", html, re.S)
    txt = m.group(1) if m else html
    obsc = {}
    for line in txt.splitlines():
        mm = re.match(r"\s*(\S{3})\s+(-?\d+\.?\d*)\s+"
                      r"(-?\d+\.?\d*)\s+([+-]?\d+\.?\d*)\s+(.*)", line)
        if not mm:
            continue
        code = mm.group(1)
        if not re.fullmatch(r"[0-9A-Z]{3}", code):
            continue
        name = re.sub(r"<[^>]+>", "", mm.group(5)).strip()
        obsc[code] = {"lon": float(mm.group(2)),
                      "cospar": float(mm.group(3)),
                      "sinpar": float(mm.group(4)),
                      "name": name}
    OBSC_PATH.write_text(json.dumps(obsc))
    _prov_add("obscodes.json", OBSC_URL,
              html.encode(), "MPC observatory code table "
              "(longitude deg, rho*cos/sin phi' parallax constants)")
    return obsc


def _safe(des):
    return re.sub(r"[^A-Za-z0-9]+", "_", des).strip("_")


def get_obs(des):
    path = OBS_DIR / f"{_safe(des)}.json"
    if path.exists():
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            path.unlink()  # corrupt cache (interrupted write): refetch
    err = None
    d = None
    body = None
    for k in range(4):
        try:
            r = requests.get(
                MPC_OBS_URL,
                json={"desigs": [des], "output_format": ["ADES_DF"]},
                timeout=180, headers={"User-Agent": "tep9/1.0"})
            if r.status_code == 200:
                # the service intermittently appends a second payload
                # to the body under load; take the first complete
                # JSON document and ignore the trailer
                d, _end = json.JSONDecoder().raw_decode(
                    r.text.lstrip())
                body = r.content
                break
            err = f"get-obs HTTP {r.status_code}"
        except Exception as exc:
            err = str(exc)
        time.sleep(1.5 * (k + 1))
    if d is None:
        raise RuntimeError(err or "get-obs failed")
    rows = d[0].get("ADES_DF") if isinstance(d, list) and d else None
    rows = rows or []
    path.write_text(json.dumps(rows))
    _prov_add(f"obs/{_safe(des)}.json", MPC_OBS_URL +
              f"?desigs={urllib.parse.quote(des)}",
              body or b"", f"MPC get-obs ADES_DF for {des}")
    return rows


def get_sbdb(des):
    path = SBDB_DIR / f"{_safe(des)}.json"
    if path.exists():
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            path.unlink()
    url = (SBDB_URL + "?" +
           urllib.parse.urlencode({"des": des, "full-prec": "true"}))
    txt = _get(url)
    stripped = txt.lstrip()
    d, _end = json.JSONDecoder().raw_decode(stripped)
    if "orbit" not in d:
        raise RuntimeError("no orbit block")
    path.write_text(stripped[:_end])
    _prov_add(f"sbdb_fp/{_safe(des)}.json", url, txt.encode(),
              f"JPL SBDB full-precision solution for {des}")
    return d


# ------------------------------------------------------------------
# astrometry + dynamics
# ------------------------------------------------------------------

def gmst_deg(et):
    jd_ut = (et - sp.deltet(et, "ET")) / 86400.0 + 2451545.0
    return (280.46061837 + 360.98564736629 * (jd_ut - 2451545.0)) % 360.0


def topo_ecl(stn, et):
    o = OBSC.get(stn)
    if o is None:
        return None
    lst = math.radians(gmst_deg(et) + o["lon"])
    v = np.array([o["cospar"] * math.cos(lst),
                  o["cospar"] * math.sin(lst),
                  o["sinpar"]]) * R_EARTH_AU
    return RX @ v


def parse_obs(rows):
    obs = []
    for o in rows:
        if o.get("Obstype") != "optical":
            continue
        ra, dec = o.get("ra"), o.get("dec")
        stn = o.get("stn")
        if ra is None or dec is None or stn is None or stn not in OBSC:
            continue
        try:
            et = sp.utc2et(str(o["obstime"]).rstrip("Z"))
        except Exception:
            continue
        obs.append(dict(et=et, ra=math.radians(float(ra)),
                        dec=math.radians(float(dec)), stn=stn))
    obs.sort(key=lambda o: o["et"])
    return obs


def _body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return (RX @ np.array(st[:3]) / AU_KM,
            RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR)


def init_sim(et):
    sim = rebound.Simulation()
    sim.G = MU
    ps, vs = _body_state("10", et)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2],
            m=1.0)
    for b in PLANET_IDS:
        pp, vv = _body_state(b, et)
        sim.add(x=pp[0], y=pp[1], z=pp[2], vx=vv[0], vy=vv[1], vz=vv[2],
                m=GM[b] / GM_SUN)
    return sim


def propagate_states(r0, v0, et0, ets):
    """Comet heliocentric ecliptic states (pos+vel) at each et."""
    sim = init_sim(et0)
    p = sim.particles
    ps = np.array([p[0].x, p[0].y, p[0].z])
    vs = np.array([p[0].vx, p[0].vy, p[0].vz])
    sim.add(x=r0[0] + ps[0], y=r0[1] + ps[1], z=r0[2] + ps[2],
            vx=v0[0] + vs[0], vy=v0[1] + vs[1], vz=v0[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = 0.001  # collision scale: bound IAS15 against step collapse
    out = np.empty((len(ets), 6))
    order = np.argsort(ets)
    for ei in order:
        t_yr = (ets[ei] - et0) / 86400.0 / DAY_YR
        sim.integrate(t_yr, exact_finish_time=1)
        out[ei] = [p[nc].x - p[0].x, p[nc].y - p[0].y,
                      p[nc].z - p[0].z, p[nc].vx - p[0].vx,
                      p[nc].vy - p[0].vy, p[nc].vz - p[0].vz]
    return out


def observer_geom(obs):
    """Orbit-independent observation geometry: epochs, observatory
    ecliptic positions, invalid row mask.  Computed once per obs list;
    every LM residual evaluation reuses it instead of repeating the
    per-observation SPICE lookups."""
    ets = np.array([o["et"] for o in obs])
    robs = np.empty((len(obs), 3))
    bad = np.zeros(len(obs), bool)
    for i, o in enumerate(obs):
        topo = topo_ecl(o["stn"], o["et"])
        if topo is None:
            bad[i] = True
            continue
        robs[i] = (_body_state("399", o["et"])[0]
                   - _body_state("10", o["et"])[0] + topo)
    ra = np.array([o["ra"] for o in obs])
    cosdec = np.array([math.cos(o["dec"]) for o in obs])
    dec = np.array([o["dec"] for o in obs])
    return {"ets": ets, "robs": robs, "bad": bad,
            "ra": ra, "dec": dec, "cosdec": cosdec}


def residuals(r0, v0, et0, obs, geom=None):
    """Computed-minus-observed residuals, arcsec (dRa*cosDec, dDec).
    geom: precomputed observer_geom for this obs list (same order); if
    omitted it is built here, which is the slow path for repeated
    evaluation."""
    if geom is None:
        geom = observer_geom(obs)
    ets, robs, bad = geom["ets"], geom["robs"], geom["bad"]
    st = propagate_states(r0, v0, et0, ets)
    rc0 = st[:, :3]
    tau = np.linalg.norm(rc0 - robs, axis=1) / C_AU_DAY
    rc = propagate_states(r0, v0, et0, ets - tau * 86400.0)[:, :3]
    u = rc - robs
    # 'over' included: rejected LM proposals can propagate to absurd
    # radii; the resulting inf/nan residuals are filtered below
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        u /= np.linalg.norm(u, axis=1)[:, None]
        u_eq = (RX.T @ u.T).T
        ra_p = np.arctan2(u_eq[:, 1], u_eq[:, 0])
        dec_p = np.arcsin(np.clip(u_eq[:, 2], -1, 1))
    res = np.empty((len(obs), 2))
    dra = (ra_p - geom["ra"] + math.pi) % (2 * math.pi) - math.pi
    res[:, 0] = dra * geom["cosdec"] * AS_RAD
    res[:, 1] = (dec_p - geom["dec"]) * AS_RAD
    res[bad] = np.nan
    return res


def _geom_sub(geom, mask):
    return {k: v[mask] for k, v in geom.items()}


def fit(r0, v0, et0, obs, iters=ITERS, clip=CLIP):
    """LM differential correction of (r,v) at et0 against obs."""
    r0 = np.array(r0, float)
    v0 = np.array(v0, float)
    keep = np.ones(len(obs), bool)
    geom_all = observer_geom(obs)
    hs = [1e-5, 1e-5, 1e-5, 1e-7, 1e-7, 1e-7]
    for _it in range(iters):
        ob = [o for k, o in enumerate(obs) if keep[k]]
        gb = _geom_sub(geom_all, keep)
        res = residuals(r0, v0, et0, ob, gb)
        finite = np.isfinite(res).all(axis=1)
        ob = [o for f, o in zip(finite, ob) if f]
        res = res[finite]
        gf = _geom_sub(gb, finite)
        if len(res) < 4:
            break
        J = np.empty((len(res) * 2, 6))
        for j in range(6):
            dr, dv = r0.copy(), v0.copy()
            if j < 3:
                dr[j] += hs[j]
            else:
                dv[j - 3] += hs[j]
            rp = residuals(dr, dv, et0, ob, gf)
            J[:, j] = ((rp - res).ravel() / hs[j])
        bad = ~np.isfinite(J).all(axis=1)
        J[bad] = 0.0
        y = res.ravel()
        A = J.T @ J
        b = J.T @ y
        lam = 1e-6 * np.trace(A) / 6.0
        try:
            d = np.linalg.solve(A + lam * np.eye(6), -b)
        except np.linalg.LinAlgError:
            break
        # accept step only if it reduces rms; else damp harder
        r_try, v_try = r0 + d[:3], v0 + d[3:]
        rms_old = float(np.sqrt(np.nanmean(res ** 2)))
        rms_new = float(np.sqrt(np.nanmean(
            residuals(r_try, v_try, et0, ob, gf) ** 2)))
        tries = 0
        while not np.isfinite(rms_new) or rms_new > rms_old:
            d *= 0.25
            r_try, v_try = r0 + d[:3], v0 + d[3:]
            rms_new = float(np.sqrt(np.nanmean(
                residuals(r_try, v_try, et0, ob, gf) ** 2)))
            tries += 1
            if tries > 8 or np.linalg.norm(d[:3]) < 1e-12:
                break
        # Exhausting damping is a rejected proposal, not permission to
        # accept a worse or nonfinite orbit. Keep the last accepted state.
        if not np.isfinite(rms_new) or rms_new > rms_old:
            break
        r0, v0 = r_try, v_try
        res2 = residuals(r0, v0, et0, obs, geom_all)
        r2 = np.linalg.norm(res2, axis=1)
        r2[~np.isfinite(r2)] = np.nan
        med = np.nanmedian(r2)
        mad = np.nanmedian(np.abs(r2 - med)) * 1.4826 + 1e-9
        newkeep = (np.abs(r2 - med) < clip * mad) & np.isfinite(r2)
        if (newkeep == keep).all() and np.linalg.norm(d[:3]) < 1e-10:
            keep = newkeep
            break
        keep = newkeep
    res = residuals(r0, v0, et0, obs, geom_all)
    rms_all = float(np.sqrt(np.nanmean(res[keep] ** 2))) \
        if keep.any() else float("nan")
    return r0, v0, rms_all, keep, res


def state_to_perih(r, v, et):
    el = sp.oscelt(np.r_[r * AU_KM, v * AU_KM / (DAY_YR * 86400.0)],
                   et, GM_SUN)
    q, e, i, Om, w = (el[0] / AU_KM, el[1], math.degrees(el[2]),
                      math.degrees(el[3]), math.degrees(el[4]))
    p = perih_dir(math.radians(w), math.radians(Om), math.radians(i))
    return p, q, e, i, Om, w


# ------------------------------------------------------------------
# per-comet refit driver
# ------------------------------------------------------------------

def fit_comet(des, min_obs_leg=MIN_OBS_LEG,
              min_leg_span_d=MIN_LEG_SPAN_D, sbdb_des=None):
    """Fetch, fit, and propagate one comet.

    Returns (result_dict, None) or (None, error_string).  The result
    carries the SBDB seed elements, the fitted states per leg label
    ("in", "out", "all"), the single-fit boundary products, and --
    when both legs are fittable -- the leg-fit products (ddirf,
    cross-arc rms, d_in_leg / d_out_leg).

    ``sbdb_des`` overrides the designator used for the SBDB seed
    fetch when the MPC observation designator differs (e.g. ISOs:
    MPC astrometry under "C/2025 N1", SBDB orbit under "3I")."""
    try:
        raw = get_obs(des)
    except Exception as exc:
        return None, f"mpc fetch: {exc}"
    obs = parse_obs(raw)
    if len(obs) < 2 * min_obs_leg:
        return None, f"only {len(obs)} usable obs"
    try:
        sb = get_sbdb(sbdb_des or des)
        el = {e["name"]: float(e["value"])
              for e in sb["orbit"]["elements"] if e.get("value")}
        el["epoch"] = float(sb["orbit"]["epoch"])
    except Exception as exc:
        return None, f"sbdb: {exc}"
    if not all(k in el for k in ("q", "e", "i", "om", "w", "tp")):
        return None, "sbdb elements incomplete"
    et0 = (el["epoch"] - 2451545.0) * 86400.0
    tp_et = (el["tp"] - 2451545.0) * 86400.0
    st = np.array(sp.conics(
        [el["q"] * AU_KM, el["e"], math.radians(el["i"]),
         math.radians(el["om"]), math.radians(el["w"]),
         0.0, tp_et, GM_SUN], et0))
    r_seed = st[:3] / AU_KM
    v_seed = st[3:] / AU_KM * 86400 * DAY_YR

    inb = [o for o in obs if o["et"] < tp_et]
    out = [o for o in obs if o["et"] >= tp_et]

    def leg_ok(leg):
        return (len(leg) >= min_obs_leg and
                (leg[-1]["et"] - leg[0]["et"]) / 86400.0
                >= min_leg_span_d)

    legs = {"in": inb if leg_ok(inb) else None,
            "out": out if leg_ok(out) else None,
            "all": obs if leg_ok(obs) else None}
    if legs["in"] is None and legs["out"] is None:
        return None, ("no leg with "
                      f">={min_obs_leg} obs/{min_leg_span_d}d "
                      f"(in={len(inb)}, out={len(out)})")

    def anchor(r_s, v_s, et_s, et_t):
        st6 = propagate_states(r_s, v_s, et_s, np.array([et_t]))[0]
        return st6[:3], st6[3:]

    fits = {}
    try:
        # fit the leg nearer the seed epoch first; its solution seeds
        # the other leg (sequential-leg initialisation)
        if et0 >= tp_et:
            todo = ["out", "in"]
        else:
            todo = ["in", "out"]
        r_s, v_s, et_s = r_seed, v_seed, et0
        for lab in todo:
            leg = legs[lab]
            if leg is None:
                continue
            mid = 0.5 * (leg[0]["et"] + leg[-1]["et"])
            ra, va = anchor(r_s, v_s, et_s, mid)
            rl, vl, rmsl, kl, _ = fit(ra, va, mid, leg)
            fits[lab] = dict(r=rl, v=vl, et=mid, rms=rmsl,
                             n_keep=int(kl.sum()), n=len(leg))
            r_s, v_s, et_s = rl, vl, mid
        if legs["all"] is not None:
            leg = legs["all"]
            mid = 0.5 * (leg[0]["et"] + leg[-1]["et"])
            ra, va = anchor(r_seed, v_seed, et0, mid)
            rl, vl, rmsl, kl, res_all = fit(ra, va, mid, leg)
            fits["all"] = dict(r=rl, v=vl, et=mid, rms=rmsl,
                               n_keep=int(kl.sum()), n=len(leg))
    except Exception as exc:
        return None, f"fit: {exc}"

    res = dict(des=des,
               n_obs=len(obs), n_in=len(inb), n_out=len(out),
               tp_jd=el["tp"], elements=el, fits=fits)

    # element agreement of the full-arc fit vs its SBDB seed
    p_all = None
    if "all" in fits:
        fa = fits["all"]
        p_all, q_f, e_f, i_f, Om_f, w_f = state_to_perih(
            fa["r"], fa["v"], fa["et"])
        p_seed = perih_dir(math.radians(el["w"]),
                           math.radians(el["om"]),
                           math.radians(el["i"]))
        res.update(fit_all_rms=fa["rms"], fit_all_n=fa["n_keep"],
                   q_fit=q_f, e_fit=e_f, i_fit=i_f,
                   om_fit=Om_f, w_fit=w_f, p_osc=p_all,
                   del_q=q_f - el["q"], del_e=e_f - el["e"],
                   del_i=i_f - el["i"], del_pdeg=sep(p_all, p_seed))
        for lab in ("in", "out"):
            if lab in fits:
                res[f"{lab}_rms"] = fits[lab]["rms"]
                res[f"{lab}_nkeep"] = fits[lab]["n_keep"]

    # boundary asymptotes; hyperbolic seeds (ISOs) carry real escape
    # speed at the 250 AU sphere, so the speed ceiling is relaxed
    v_bnd = 40.0 if el["e"] > 1.0 else None
    _kw = {} if v_bnd is None else {"v_max": v_bnd}
    try:
        if "all" in fits:
            fa = fits["all"]
            rb = integrate_leg(fa["r"], fa["v"], fa["et"], -1, **_kw)
            rf = integrate_leg(fa["r"], fa["v"], fa["et"], +1, **_kw)
            if rb and rf:
                drot = sep(rb["phat"], rf["phat"])
                aph = -rb["phat"]
                a_loc = 1e6 / rb["aa"] if rb["aa"] != 0 else np.inf
                e_loc = 1.0 - q_f / a_loc
                h_loc = math.sqrt(MU * q_f * (1.0 + e_loc))
                om_b = h_loc / R_STOP ** 2
                res.update(drot=drot,
                           d_in=sep(rb["phat"], p_all),
                           d_out=sep(rf["phat"], p_all),
                           daa=rf["aa"] - rb["aa"],
                           denc=min(rb["denc"], rf["denc"]),
                           theta=sep(aph, TNO),
                           aph=[float(x) for x in np.round(aph, 6)],
                           dtau=math.radians(drot) / om_b)
        if "in" in fits and "out" in fits:
            fi, fo = fits["in"], fits["out"]
            rbi = integrate_leg(fi["r"], fi["v"], fi["et"], -1, **_kw)
            rfo = integrate_leg(fo["r"], fo["v"], fo["et"], +1, **_kw)
            if rbi and rfo:
                res["ddirf"] = sep(rbi["phat"], rfo["phat"])
                if p_all is not None:
                    res["d_in_leg"] = sep(rbi["phat"], p_all)
                    res["d_out_leg"] = sep(rfo["phat"], p_all)
                # cross-arc prediction: inbound fit -> withheld
                # outbound astrometry (never fitted)
                rms_x = float(np.sqrt(np.nanmean(
                    residuals(fi["r"], fi["v"], fi["et"], out) ** 2)))
                res["xarc_rms_in2out"] = rms_x
    except Exception as exc:
        res["boundary_err"] = str(exc)[:150]

    if "drot" not in res and "ddirf" not in res:
        return None, "boundary propagation failed"
    return res, None
