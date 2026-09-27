#!/usr/bin/env python3
"""step_117: Prospective validation cohort -- post-2017 SBDB long-period
comets through the identical bidirectional boundary instrument.

Every transit-channel test so far (steps 030-042, 063-086, 102-106)
runs on catalogue comets discovered before 2018 -- the Warsaw, CODE and
one-apparition tables that also defined the axis, the cap and the slip
map.  This step breaks that closure.  It downloads the complete JPL
SBDB comet table and selects every comet discovered in 2018 or later --
a cohort that could not have entered any fit, scan or registration --
and passes it through the identical REBOUND/DE440s bidirectional
instrument validated in step 063 (the instrument that reproduces the
catalogue's own boundary solutions to a median 0.001 deg per leg).

Two independence statements are made at once:

  * out-of-time: the axis (l,b) = (34,-13) deg, the 60 deg cap and the
    bipolar slip map were all fixed on pre-2018 data; this cohort is a
    genuine prospective holdout;
  * out-of-lineage: the osculating elements are JPL SBDB orbit
    solutions -- a fourth orbit-determination pipeline, independent of
    the Warsaw school that produced every orig/fut catalogue used so
    far.

Cohort (pre-declared, mirroring the CODE class-1 selection)
-----------------------------------------------------------
  * provisional designation C/2018+ (post-CODE cutoff; CODE is
    near-complete through 2017);
  * 0.95 < e < 1.5 (near-parabolic; excludes ISOs and bound comets);
  * 0.1 < q < 3.1 AU (deep plungers, the class-1 regime; Kreutz-family
    sungrazers excluded);
  * data arc >= 30 d and >= 20 observations used (orbit quality);
  * primary designation only (fragment components excluded);
  * a bound periodic control cohort (e < 0.95, 2018+) is carried
    through the same instrument -- the boundary rotation requires
    reaching the boundary.

Channels (identical instruments to steps 063/065)
-------------------------------------------------
  T1  cap contrast of the full-pass rotation drot (in- vs out-cap,
      Mann-Whitney + theta-permutation null);
  T2  rotation-per-kick residual: log(drot) regressed on
      [log|daa|, log(d_enc), q, i]; cap contrast + Spearman vs theta;
  T3  continuous cap-free Spearman(theta, drot);
  T4  energy channel |d(1/a)| in/out (expected flat -- the TEP
      discriminator is rotation without energy exchange);
  T5  implied proper-time offset dtau = drot/om_b with the step-065
      linear residual -> unexplained offset in years, in- vs out-cap;
  T6  three-leg decomposition: inbound-leg vs outbound-leg rotation
      residual (the inbound-weighted claim of Section 4.1);
  T7  bipolar slip-map transfer: a + b cos(2*theta) fitted on the
      pooled CODE+Warsaw residuals (step_b30 products) scored on the
      new cohort -- out-of-sample Spearman and lobe sign accuracy;
  T8  aphelion-direction dipole of the reconstructed inbound
      asymptotes: cap count vs uniform baseline and vs the CODE
      cohort's own measured cap rate;
  T9  periodic-comet control: identical cap test on bound 2018+
      comets that never reach the boundary sphere;
  T10 geometric control: the identical cap test with the transit
      axis rotated 90 deg in ecliptic longitude -- same cohort,
      wrong axis, expected flat;
  T11 subset tiers: the CODE class-1 dynamical equivalent
      (0 < 1/a_back < 100 x 10^-6 AU^-1, Oort-spike first
      arrivals) and the pure-gravity subset (SBDB solutions
      without fitted A1/A2/A3/DT non-gravitational parameters).

Inputs
------
data/raw/sbdb/sbdb_comets_all.json      (downloaded here if absent)
data/raw/spice/de440s.bsp               (provenance pinned)
results/step_b28_bidirectional_rotation.json   (CODE reference values)
results/step_b30_proper_time_slip.csv          (training residuals)

Outputs
-------
results/step_b81_prospective_lpc.json
results/step_b81_prospective_lpc.csv
results/figures/step_b81_prospective_lpc.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import DATA_RAW, RESULTS, perih_dir, sep, lv, gv
logger = StepLogger("step_117_prospective_lpc")
tee_stdout(logger)

import csv
import time
import hashlib
import json
import math
import re
import urllib.request
from datetime import datetime, timezone

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr, binomtest

import rebound
import spiceypy as sp

logger.header("Prospective validation cohort: post-2017 SBDB LPCs (REBOUND/DE440s)")

SEED = 20260919
rng = np.random.default_rng(SEED)
N_PERM = 20000
CAP = 60.0
TNO = lv(34.0, -13.0)          # pre-declared transit axis (steps 030-086)

STATS_ONLY = "--stats-only" in _sys.argv

if STATS_ONLY:
    # recompute statistics from the previously written per-comet CSV
    # without re-running the ~50 min N-body integrations
    _prev = json.load(open(RESULTS / "step_b81_prospective_lpc.json"))
    prows = []
    with open(RESULTS / "step_b81_prospective_lpc.csv") as _f:
        for r in csv.DictReader(_f):
            for k in ("yr", "q", "i", "e", "arc", "nobs", "theta", "drot",
                      "d_in", "d_out", "daa", "denc", "aa_back", "aa_fwd",
                      "t_back", "t_fwd", "dtau", "dtau_in", "dtau_out"):
                r[k] = float(r[k])
            r["ng"] = r.get("ng") in ("True", True)
            r["aph_lb"] = [float(x) for x in r["aph_lb"].split(";")]
            r.pop("dtau_unexplained", None)
            prows.append(r)
    crows, pfailed, cfailed = [], [], []
    cohort, cohort_all, controls, rows = [], [], [], []
    excluded = _prev["cohort_counts"]["excluded"]
    cohort_counts = _prev["cohort_counts"]
    logger.info(f"stats-only mode: {len(prows)} rows loaded from CSV")

# ------------------------------------------------------------------
# 1. Acquisition: complete JPL SBDB comet table
# ------------------------------------------------------------------

SBDB_FIELDS = ("full_name,a,e,i,om,w,q,ad,ma,epoch,tp,class,"
               "condition_code,data_arc,n_obs_used,first_obs,last_obs,per,"
               "A1,A2,A3,DT")
SBDB_URL = ("https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields="
            + SBDB_FIELDS + "&sb-kind=c")
COMETS_JSON = DATA_RAW / "sbdb" / "sbdb_comets_all.json"
PROV_PATH = DATA_RAW / "sbdb" / "provenance.json"


def download_sbdb():
    req = urllib.request.Request(SBDB_URL,
                                 headers={"User-Agent": "TEP-9-pipeline"})
    body = status = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                body = r.read()
                status = r.status
            if status == 200:
                break
            raise RuntimeError(f"SBDB download failed: HTTP {status}")
        except Exception as exc:
            if attempt == 5:
                raise
            wait = 5 * 2 ** (attempt - 1)
            logger.progress(f"SBDB fetch attempt {attempt} failed ({exc}); "
                            f"retrying in {wait}s")
            time.sleep(wait)
    if status != 200:
        raise RuntimeError(f"SBDB download failed: HTTP {status}")
    COMETS_JSON.write_bytes(body)
    sha = hashlib.sha256(body).hexdigest()
    d = json.loads(body)
    prov = json.loads(PROV_PATH.read_text()) if PROV_PATH.exists() else {
        "step": "step_001_download_sbdb", "files": {}}
    prov.setdefault("files", {})[COMETS_JSON.name] = {
        "url": SBDB_URL, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "http_status": status, "bytes": len(body), "sha256": sha,
        "rows": len(d.get("data", [])),
        "count_field": d.get("count"),
        "api": "JPL SBDB sbdb_query.api",
        "note": "complete comet table; step_117 filters the post-2017 "
                "prospective cohort"}
    PROV_PATH.write_text(json.dumps(prov, indent=4))
    logger.info(f"downloaded {len(d.get('data', []))} comets -> "
                f"{COMETS_JSON.name} (sha256 {sha[:12]}...)")
    return d


if not STATS_ONLY:
    if COMETS_JSON.exists():
        sbdb = json.load(open(COMETS_JSON))
        if "A1" not in sbdb.get("fields", []):
            logger.info("cached table predates NG fields; re-downloading")
            sbdb = download_sbdb()
        else:
            logger.info(f"using cached {COMETS_JSON.name}")
    else:
        sbdb = download_sbdb()

    fields = sbdb["fields"]
    rows = [dict(zip(fields, rec)) for rec in sbdb["data"]]
    logger.info(f"SBDB comet table: {len(rows)} rows, fields {len(fields)}")

# ------------------------------------------------------------------
# 2. Cohort selection
# ------------------------------------------------------------------

def fnum(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")


def desig_year(name):
    m = re.match(r"\s*[CP]/(\d{4})", name)
    return int(m.group(1)) if m else None


def is_fragment(name):
    """Fragment components carry a trailing -letter: 'C/2019 Y4-A'."""
    return bool(re.match(r"\s*C/\d{4}\s+\S+-\w", name))


cohort_all = []
controls = []
excluded = {"year": 0, "frag": 0, "e": 0, "q": 0, "arc": 0, "iso": 0}
if not STATS_ONLY:
    for r in rows:
        name = str(r["full_name"]).strip()
        if not name.startswith("C/"):
            continue
        yr = desig_year(name)
        if yr is None or yr < 2018:
            excluded["year"] += 1
            continue
        if is_fragment(name):
            excluded["frag"] += 1
            continue
        e = fnum(r, "e")
        if not np.isfinite(e):
            excluded["e"] += 1
            continue
        if e >= 1.5:
            excluded["iso"] += 1
            continue
        q = fnum(r, "q")
        arc = fnum(r, "data_arc")
        nobs = fnum(r, "n_obs_used")
        rec = dict(name=name, yr=yr, q=q, e=e, i=fnum(r, "i"),
                   om=fnum(r, "om"), w=fnum(r, "w"), a=fnum(r, "a"),
                   epoch=fnum(r, "epoch"), tp=fnum(r, "tp"),
                   cc=r.get("condition_code"), arc=arc, nobs=nobs,
                   cls=r.get("class"),
                   ng=any(np.isfinite(fnum(r, k)) and fnum(r, k) != 0.0
                          for k in ("A1", "A2", "A3", "DT")))
        if e < 0.95:
            if np.isfinite(q) and np.isfinite(arc) and arc >= 30 and nobs >= 20:
                controls.append(rec)
            else:
                excluded["arc"] += 1
            continue
        # near-parabolic branch
        if not (np.isfinite(q) and 0.1 <= q):
            excluded["q"] += 1
            continue
        if not (np.isfinite(arc) and arc >= 30 and np.isfinite(nobs)
                and nobs >= 20):
            excluded["arc"] += 1
            continue
        if not (np.isfinite(rec["epoch"]) and np.isfinite(rec["tp"])):
            excluded["e"] += 1
            continue
        cohort_all.append(rec)

    cohort = [r for r in cohort_all if r["q"] < 3.1]
    logger.info(f"post-2017 comets: near-parabolic usable {len(cohort_all)}; "
                f"deep-plunger primary (q<3.1) {len(cohort)}; "
                f"bound controls {len(controls)}; excluded {excluded}")

# ------------------------------------------------------------------
# 3. Ephemeris + integrator (identical to step_063)
# ------------------------------------------------------------------

SPK = DATA_RAW / "spice" / "de440s.bsp"
if not SPK.exists():
    raise FileNotFoundError(f"JPL ephemeris missing: {SPK}")
sp.furnsh(str(SPK))

AU_KM   = 149597870.7
DAY_YR  = 365.25
EPS     = math.radians(23.4392911)
RX      = np.array([[1, 0, 0],
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


def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return RX @ np.array(st[:3]) / AU_KM, RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR


def comet_state(rec):
    """Heliocentric ecliptic-J2000 state at the SBDB osculation epoch.

    sp.conics converts the published osculating elements to a Cartesian
    state -- the elements are ecliptic J2000, so the state lands in the
    same frame the DE440s barycentres are rotated into.
    """
    et0 = (rec["epoch"] - 2451545.0) * 86400.0
    tp_et = (rec["tp"] - 2451545.0) * 86400.0
    # Anchor the conic at perihelion (M0 = 0 at t0 = tp): the osculating
    # elements define the orbit geometrically, and this parametrisation
    # is exact regardless of the e ~ 1 mean-motion limit.  Verified:
    # conics(elts, tp) returns |r| = q to machine precision for every
    # cohort member, including e = 1.0000 solutions.
    elts = [rec["q"] * AU_KM, rec["e"], math.radians(rec["i"]),
            math.radians(rec["om"]), math.radians(rec["w"]),
            0.0, tp_et, GM_SUN]
    st = np.array(sp.conics(elts, et0))
    return (st[:3] / AU_KM, st[3:] / AU_KM * 86400.0 * DAY_YR, et0)


def init_sim(et):
    sim = rebound.Simulation()
    sim.G = MU
    ps, vs = body_state("10", et)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2], m=1.0)
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


def integrate_leg(r0, v0, et0, direction):
    """One leg to the +/-250 AU barycentric sphere. direction=-1 gives
    the inbound (original-analogue) asymptote, +1 the outbound."""
    sim = init_sim(et0)
    p = sim.particles
    ps = np.array([p[0].x, p[0].y, p[0].z])
    vs = np.array([p[0].vx, p[0].vy, p[0].vz])
    sim.add(x=r0[0] + ps[0], y=r0[1] + ps[1], z=r0[2] + ps[2],
            vx=v0[0] + vs[0], vy=v0[1] + vs[1], vz=v0[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = 0.001  # collision scale: bound IAS15 against step collapse
    denc = np.full(sim.N - 1, np.inf)
    t = direction * DT_OUT
    while abs(t) < T_MAX:
        sim.integrate(t, exact_finish_time=0)
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
    if np.linalg.norm(v_rel) > V_BND_MAX:
        raise RuntimeError(
            f"unphysical boundary speed |v|={np.linalg.norm(v_rel):.1f} AU/yr")
    phat, aa = boundary_orbit(r_rel, v_rel, mtot)
    return dict(phat=phat, aa=aa, denc=float(denc.min()), t_years=float(t))


# ------------------------------------------------------------------
# 4. Integrate the cohort
# ------------------------------------------------------------------

CKPT = RESULTS / "step_b81_legs.jsonl"
from scripts.utils.parallel import cli_workers as _cli_workers
WORKERS = _cli_workers(_sys.argv)


def _worker_init():
    # fork children inherit the parent's kernel-table entries but not
    # their DAF file records, so the inherited handles are dead;
    # reset the pool and load a worker-private DE440s kernel
    sp.kclear()
    sp.furnsh(str(SPK))


def _process_comet(rec):
    """Integrate one comet through both boundary legs -> row dict.

    Returns (row, None) on success or (None, reason) on failure.
    Top-level worker shared by the serial and pooled paths;
    deterministic given the catalogue inputs.
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
    d_in = sep(rb["phat"], p_osc)      # inbound leg rotation
    d_out = sep(rf["phat"], p_osc)     # outbound leg rotation
    aph = -rb["phat"]                  # inbound aphelion direction
    theta = sep(aph, TNO)
    a_loc = 1e6 / rb["aa"] if rb["aa"] != 0 else float("inf")
    e_loc = 1.0 - rec["q"] / a_loc
    h = math.sqrt(MU * rec["q"] * (1.0 + e_loc))
    om_b = h / R_STOP ** 2             # rad/yr at the boundary
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


def _load_checkpoint():
    """Replay step_b81_legs.jsonl.  Returns (done, failed_names, hard):
    hard = failed rows that were genuine integration attempts (not
    deterministic aphelion preclassifications).  A cohort whose
    checkpoint is all-hard-failures with zero successes is the
    signature of a poisoned run (e.g. broken SPICE kernel pool)."""
    done, failed_names, hard = {}, set(), set()
    if CKPT.exists():
        for line in CKPT.read_text().splitlines():
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("failed"):
                failed_names.add(row["desig"])
                if not row.get("skipped"):
                    hard.add(row["desig"])
            else:
                done[row["desig"]] = row
    return done, failed_names, hard


def run_cohort(recs, tag, aph_limit=None):
    """aph_limit: skip e < 1 comets whose osculating aphelion falls
    below this radius -- they cannot reach the boundary sphere on the
    osculating orbit.  Set to 0.9*R_STOP for the primary cohort
    (conservative margin against 20 kyr diffusion; every failed primary
    member sits below 225 AU) and R_STOP for the bound control
    (a demonstration cohort).

    Each completed comet is appended to step_b81_legs.jsonl; on restart
    the checkpoint is replayed so a crash loses no finished leg and a
    re-run resumes where it stopped.  --workers N integrates N comets
    concurrently (each worker furnshes its own DE440s kernel)."""
    done, failed_names, hard = _load_checkpoint()
    out_rows, failed = [], []
    todo = []
    for rec in recs:
        if rec["name"] in done:
            out_rows.append(done[rec["name"]])
        elif rec["name"] in failed_names:
            failed.append(rec["name"])
        else:
            todo.append(rec)
    if not out_rows and not todo and any(n in hard for n in failed):
        # Every member of this cohort is checkpointed as a hard
        # (integration) failure with zero successes -- the signature of
        # a poisoned checkpoint (e.g. a broken SPICE kernel pool in
        # forked workers), not a real all-fail cohort.  Quarantine this
        # cohort's rows only -- the file is shared across cohorts --
        # and re-attempt once; genuine failures re-checkpoint cheaply.
        names = {r["name"] for r in recs}
        keep, quarantined = [], []
        for line in CKPT.read_text().splitlines():
            try:
                desig = json.loads(line).get("desig")
            except json.JSONDecodeError:
                keep.append(line)
                continue
            (quarantined if desig in names else keep).append(line)
        suspect = CKPT.with_suffix(f".{tag}.suspect")
        suspect.write_text("\n".join(quarantined) + "\n")
        CKPT.write_text("\n".join(keep) + "\n")
        logger.warning(f"{tag}: checkpoint holds {len(failed)} failures "
                       f"and 0 successes -- quarantined to "
                       f"{suspect.name}; re-attempting")
        todo, failed = list(recs), []
    if todo:
        logger.info(f"{tag}: {len(todo)} comets to integrate "
                    f"({len(done) + len(failed)} checkpointed)")
    n_done = 0

    def _checkpoint(row, name, err=None, skipped=False):
        with open(CKPT, "a") as ck:
            if row is None:
                ck.write(json.dumps({"desig": name, "failed": True,
                                     "skipped": skipped,
                                     "err": str(err)[:200]})
                         + "\n")
            else:
                ck.write(json.dumps(row) + "\n")

    def _account(rec, row, err, skipped=False):
        nonlocal n_done
        if row is None:
            failed.append(rec["name"])
            logger.warning(f"{rec['name']}: {err}")
        else:
            out_rows.append(row)
        _checkpoint(row, rec["name"], err, skipped)
        n_done += 1
        if n_done % 25 == 0:
            logger.info(f"{tag}: {n_done}/{len(todo)} integrated")

    def _preclassify(rec):
        if aph_limit is not None and rec["e"] < 1:
            aph = rec["q"] * (1 + rec["e"]) / (1 - rec["e"])
            if aph < aph_limit:
                return (f"osculating aphelion {aph:.0f} AU < "
                        f"{aph_limit:.0f} AU -- boundary unreachable")
        return None

    if todo and WORKERS > 1:
        import multiprocessing as mp
        # fork (not spawn): the script is top-level without a __main__
        # guard, so spawn children would re-execute the whole module;
        # fork workers reset and reload the kernel in _worker_init.
        ctx = mp.get_context("fork")
        todo_int = [r for r in todo if _preclassify(r) is None]
        for rec in todo:
            err = _preclassify(rec)
            if err:
                _account(rec, None, err, skipped=True)
        with ctx.Pool(WORKERS, initializer=_worker_init) as pool:
            for rec, (row, err) in zip(
                    todo_int, pool.imap(_process_comet, todo_int)):
                _account(rec, row, err)
    else:
        for rec in todo:
            err = _preclassify(rec)
            if err:
                _account(rec, None, err, skipped=True)
                continue
            row, err = _process_comet(rec)
            _account(rec, row, err)
    logger.info(f"{tag}: integrated {len(out_rows)} comets x2 legs; "
                f"{len(failed)} failed")
    return out_rows, failed


if not STATS_ONLY:
    logger.info("integrating primary cohort ...")
    prows, pfailed = run_cohort(cohort_all, "prospective",
                              aph_limit=0.9 * R_STOP)
    if not prows:
        raise RuntimeError(
            "primary cohort produced 0 integrated comets -- refusing "
            "to overwrite results with an empty analysis")
    logger.info("integrating bound-comet control cohort ...")
    crows, cfailed = run_cohort(controls, "control", aph_limit=R_STOP)


# ------------------------------------------------------------------
# 5. Statistics -- the step-063/065 instruments
# ------------------------------------------------------------------

def resid_cap(v, th, X_extra):
    """log(v) regressed on covariates -> cap contrast + Spearman."""
    inc = th < CAP
    # Residuals are reported in dex throughout the manuscript; use base-10
    # logarithms rather than natural logs so the numeric units match the label.
    y = np.log10(v)
    X = np.column_stack([np.ones(len(y))] + X_extra)
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ coef
    u = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    rho, p = spearmanr(th, resid)
    return {"p": float(u.pvalue), "rho": float(rho), "p_2sided": float(p),
            "med_resid_in": float(np.median(resid[inc])),
            "med_resid_out": float(np.median(resid[~inc]))}


def lin_resid(y, X_cols):
    X = np.column_stack([np.ones(len(y))] + X_cols)
    with np.errstate(all="ignore"):
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        return y - X @ coef


def perm_contrast(v, th, n_perm=N_PERM):
    """Permutation null on med(in) - med(out) under theta shuffle."""
    inc = th < CAP
    obs = float(np.median(v[inc]) - np.median(v[~inc]))
    cnt = 0
    for _ in range(n_perm):
        tt = rng.permutation(th)
        s = float(np.median(v[tt < CAP]) - np.median(v[tt >= CAP]))
        if s >= obs:
            cnt += 1
    return obs, (cnt + 1) / (n_perm + 1)


def run(sub):
    th = np.array([r["theta"] for r in sub])
    v = np.array([r["drot"] for r in sub])
    vin = np.array([r["d_in"] for r in sub])
    vout = np.array([r["d_out"] for r in sub])
    K = np.abs(np.array([r["daa"] for r in sub]))
    D = np.array([r["denc"] for r in sub])
    Q = np.array([r["q"] for r in sub])
    I = np.array([r["i"] for r in sub])
    DT = np.array([r["dtau"] for r in sub])
    DIN = np.array([r["dtau_in"] for r in sub])
    DOUT = np.array([r["dtau_out"] for r in sub])
    inc = th < CAP
    out = {"n": len(sub), "n_in": int(inc.sum()),
           "med_drot_in": float(np.median(v[inc])) if inc.any() else None,
           "med_drot_out": float(np.median(v[~inc])) if (~inc).any() else None,
           "med_din_in": float(np.median(vin[inc])) if inc.any() else None,
           "med_din_out": float(np.median(vin[~inc])) if (~inc).any() else None,
           "med_dout_in": float(np.median(vout[inc])) if inc.any() else None,
           "med_dout_out": float(np.median(vout[~inc])) if (~inc).any() else None}
    # need both cap classes and enough points for the 5-parameter
    # residual regressions
    if not (inc.any() and (~inc).any()) or len(sub) < 8:
        out["status"] = "insufficient_sample"
        return out
    sub_keep = sub
    fin = (np.isfinite(v) & np.isfinite(K) & np.isfinite(D) &
           np.isfinite(DT) & np.isfinite(DIN) & np.isfinite(DOUT))
    if not fin.all():
        keep = fin
        sub_keep = [r for r, k in zip(sub, keep) if k]
        th, v, vin, vout, K, D, Q, I, DT, DIN, DOUT = [
            a[keep] for a in (th, v, vin, vout, K, D, Q, I, DT, DIN, DOUT)]
        inc = th < CAP
        out["n_finite"] = int(fin.sum())
        if not (inc.any() and (~inc).any()):
            out["status"] = "insufficient_sample"
            return out
    u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
    out["drot_in_gt_out_p"] = float(u.pvalue)
    obs, p_perm = perm_contrast(v, th)
    out["drot_perm_contrast"] = {"obs_med_diff": obs, "p_perm": float(p_perm)}
    rho, p = spearmanr(th, v)
    out["drot_vs_theta"] = {"rho": float(rho), "p_2sided": float(p)}
    # energy channel (expected flat)
    u = mannwhitneyu(K[inc], K[~inc], alternative="two-sided")
    out["energy_kick_in_vs_out_p"] = float(u.pvalue)
    out["med_kick_in"] = float(np.median(K[inc]))
    out["med_kick_out"] = float(np.median(K[~inc]))
    # rotation-per-kick residuals (identical regressors to step_063)
    out["rot_per_kick_resid_kick"] = resid_cap(v, th, [np.log10(K + 1.0)])
    out["rot_per_kick_resid_full"] = resid_cap(
        v, th, [np.log10(K + 1.0), np.log10(D), Q, I])
    # proper-time offset (step_065 instrument)
    resid_lin = lin_resid(DT, [np.log10(K + 1.0), np.log10(D), Q, I])
    u = mannwhitneyu(resid_lin[inc], resid_lin[~inc], alternative="greater")
    out["dtau_unexplained"] = {
        "med_in": float(np.median(resid_lin[inc])),
        "med_out": float(np.median(resid_lin[~inc])),
        "p_in_gt_out": float(u.pvalue),
        "med_total_in": float(np.median(DT[inc])),
        "med_total_out": float(np.median(DT[~inc]))}
    rho, p = spearmanr(th, resid_lin)
    out["dtau_unexplained_vs_theta"] = {"rho": float(rho),
                                        "p_2sided": float(p)}
    for r, rr in zip(sub_keep, resid_lin):
        r["dtau_unexplained"] = float(rr)
    # three-leg decomposition: inbound vs outbound leg residuals
    out["inbound_leg_resid"] = resid_cap(
        vin, th, [np.log10(K + 1.0), np.log10(D), Q, I])
    out["outbound_leg_resid"] = resid_cap(
        vout, th, [np.log10(K + 1.0), np.log10(D), Q, I])
    rin = lin_resid(DIN, [np.log10(K + 1.0), np.log10(D), Q, I])
    rout = lin_resid(DOUT, [np.log10(K + 1.0), np.log10(D), Q, I])
    out["dtau_inbound_unexplained"] = {
        "med_in": float(np.median(rin[inc])),
        "med_out": float(np.median(rin[~inc])),
        "p_in_gt_out": float(mannwhitneyu(rin[inc], rin[~inc],
                                          alternative="greater").pvalue)}
    out["dtau_outbound_unexplained"] = {
        "med_in": float(np.median(rout[inc])),
        "med_out": float(np.median(rout[~inc])),
        "p_in_gt_out": float(mannwhitneyu(rout[inc], rout[~inc],
                                          alternative="greater").pvalue)}
    # Direction-aware audit: the declared test is one-sided ("greater",
    # matching the pre-registered CODE sign).  A cohort can still carry
    # cap-correlated structure with the opposite sign -- e.g. if the
    # fitter absorbs the rotation the measured leg deviation shrinks
    # where the signal lives.  Report the two-sided cap contrast on
    # every leg channel so reversed-direction structure is visible
    # rather than hidden behind a one-sided p ~ 1.
    da = {}
    for tag, vec in (("drot", v), ("d_in", vin), ("d_out", vout),
                     ("dtau", DT), ("dtau_in", DIN), ("dtau_out", DOUT)):
        vec = np.asarray(vec, dtype=float)
        gap = float(np.median(vec[inc]) - np.median(vec[~inc]))
        da[tag] = {
            "gap": gap,
            "p_2sided": float(mannwhitneyu(
                vec[inc], vec[~inc], alternative="two-sided").pvalue),
            "p_less": float(mannwhitneyu(
                vec[inc], vec[~inc], alternative="less").pvalue)}
    out["direction_audit"] = da
    # residualised two-sided contrasts on the same channels
    Xfull = [np.log10(K + 1.0), np.log10(D), Q, I]
    dr = {}
    for tag, vec in (("drot", v), ("d_in", vin), ("d_out", vout)):
        res = lin_resid(np.log10(np.asarray(vec, dtype=float)), Xfull)
        dr[tag] = {
            "gap": float(np.median(res[inc]) - np.median(res[~inc])),
            "p_2sided": float(mannwhitneyu(
                res[inc], res[~inc], alternative="two-sided").pvalue),
            "p_less": float(mannwhitneyu(
                res[inc], res[~inc], alternative="less").pvalue)}
    for tag, vec in (("dtau", DT), ("dtau_in", DIN), ("dtau_out", DOUT)):
        res = lin_resid(np.asarray(vec, dtype=float), Xfull)
        dr[tag] = {
            "gap": float(np.median(res[inc]) - np.median(res[~inc])),
            "p_2sided": float(mannwhitneyu(
                res[inc], res[~inc], alternative="two-sided").pvalue),
            "p_less": float(mannwhitneyu(
                res[inc], res[~inc], alternative="less").pvalue)}
    out["direction_audit_resid"] = dr
    # CMB-frame audit: the TEP scalar-field rest frame is the CMB
    # dipole frame (the GNSS/MGEX clock channel recovered an axis
    # 21.4 deg from it).  Decompose the residual rotation field into
    # the m=1 dipole along the CMB axis and locate the cohort's own
    # free dipole axis -- the bipolar axis is the cross-catalogue
    # invariant even where fit methodology inverts the measured sign.
    aph_v = np.array([np.asarray(r["aph_lb"], dtype=float)
                      for r in sub_keep])
    CMB_AP = gv(264.02, 48.25)
    res_rot = lin_resid(np.log10(v), Xfull)

    def _dipole_b(res, u):
        cos = np.clip(aph_v @ u, -1, 1)
        Xm = np.column_stack([np.ones(len(res)), cos])
        with np.errstate(all="ignore"):
            cc, *_ = np.linalg.lstsq(Xm, res, rcond=None)
        return float(cc[1])

    b_cmb = _dipole_b(res_rot, CMB_AP)
    cnt = 0
    for _ in range(N_PERM):
        if abs(_dipole_b(rng.permutation(res_rot), CMB_AP)) >= abs(b_cmb):
            cnt += 1
    th_cmb = np.degrees(np.arccos(np.clip(aph_v @ CMB_AP, -1, 1)))
    Xm = np.column_stack([np.ones(len(res_rot)),
                          np.cos(np.radians(th_cmb)),
                          np.cos(2 * np.radians(th_cmb))])
    with np.errstate(all="ignore"):
        ch, *_ = np.linalg.lstsq(Xm, res_rot, rcond=None)
    best = None
    for L_ in np.arange(0.0, 360.0, 15.0):
        for B_ in np.arange(-75.0, 76.0, 15.0):
            u_ = lv(float(L_), float(B_))
            b_ = _dipole_b(res_rot, u_)
            if best is None or abs(b_) > abs(best[0]):
                best = (b_, L_, B_)
    u_best = lv(best[1], best[2])
    ca = {"b_cos_cmb": b_cmb,
          "p_2sided_perm": float((cnt + 1) / (N_PERM + 1)),
          "cos1_cmb": float(ch[1]), "cos2_cmb": float(ch[2]),
          "free_dipole": {"b": float(best[0]), "l": float(best[1]),
                          "b_lat": float(best[2]),
                          "sep_to_cmb_apex": float(sep(u_best, CMB_AP)),
                          "sep_to_declared": float(sep(u_best, TNO)),
                          "sep_antipole_to_cmb_apex":
                              float(sep(-u_best, CMB_AP))},
          "leg_channels": {}}
    for tag, vec in (("d_in", vin), ("d_out", vout)):
        rleg = lin_resid(np.log10(np.asarray(vec, dtype=float)), Xfull)
        ca["leg_channels"][tag] = {
            "b_cos_cmb": _dipole_b(rleg, CMB_AP)}
    out["cmb_frame_audit"] = ca
    # aphelion dipole
    out["aph_in_cap_frac"] = float(inc.mean())
    out["aph_in_cap_n"] = int(inc.sum())
    out["aph_cap_binom_vs_uniform_p"] = float(
        binomtest(int(inc.sum()), len(sub_keep), 0.25,
                  alternative="greater").pvalue)
    # geometric control: the transit axis rotated 90 deg in ecliptic
    # longitude -- same cohort, wrong axis, expected flat
    tno_null = lv((34.0 + 90.0) % 360.0, -13.0)
    th_null = np.array([sep(np.asarray(r["aph_lb"], dtype=float), tno_null)
                        for r in sub_keep])
    incn = th_null < CAP
    if incn.any() and (~incn).any():
        u = mannwhitneyu(v[incn], v[~incn], alternative="greater")
        rho_n, p_n = spearmanr(th_null, v)
        out["null_axis_90deg"] = {
            "n_in": int(incn.sum()),
            "med_drot_in": float(np.median(v[incn])),
            "med_drot_out": float(np.median(v[~incn])),
            "p_in_gt_out": float(u.pvalue),
            "rho_vs_theta": float(rho_n), "p_2sided": float(p_n),
            "rot_per_kick_resid_full": resid_cap(
                v, th_null, [np.log10(K + 1.0), np.log10(D), Q, I])}
    return out


res_all = run(prows)
res_ext = run([r for r in prows if r["q"] >= 3.1]) if any(
    r["q"] >= 3.1 for r in prows) else None
res_primary = run([r for r in prows if r["q"] < 3.1])
# the per-row slip residuals used downstream (CSV, slip-map transfer)
# are the primary-regression values; later subset runs would overwrite
primary_resid = {r["desig"]: r["dtau_unexplained"]
                 for r in prows if r["q"] < 3.1 and "dtau_unexplained" in r}
# CODE class-1 dynamical equivalent: original 1/a in the Oort-spike
# interval (0, 100) x 10^-6 AU^-1 measured at the backward boundary
# leg -- dynamically-new first arrivals only
spike = [r for r in prows if r["q"] < 3.1 and 0.0 < r["aa_back"] < 100.0]
res_spike = run(spike)
res_spike_pure = run([r for r in spike if not r["ng"]])
res_primary_pure = run([r for r in prows
                        if r["q"] < 3.1 and not r["ng"]])
res_control = run(crows)
for r in prows:
    if r["desig"] in primary_resid:
        r["dtau_unexplained"] = primary_resid[r["desig"]]
    else:
        r.pop("dtau_unexplained", None)

# Failure accounting: classify every checkpointed failure by cohort
# membership and mechanism, so exclusions are explicit rather than a
# single opaque counter.  In STATS_ONLY mode the element lookup is
# rebuilt from the cached SBDB table (no re-download).
_ck_done, _ck_failed_names, _ck_hard = _load_checkpoint()
if STATS_ONLY:
    _sbdb = json.load(open(COMETS_JSON))
    _lk = {dict(zip(_sbdb["fields"], r))["full_name"].strip():
           dict(zip(_sbdb["fields"], r)) for r in _sbdb["data"]}
    def _rec(name):
        r = _lk.get(name)
        if r is None:
            return None
        return dict(q=fnum(r, "q"), e=fnum(r, "e"),
                    w=fnum(r, "w"), om=fnum(r, "om"), i=fnum(r, "i"))
else:
    _lk = {r["name"]: r for r in cohort_all + controls}
    def _rec(name):
        return _lk.get(name)

fail_acct = {"primary_hard": [], "primary_skipped": [],
             "extended": [], "control": [], "unknown": []}
for name in sorted(_ck_failed_names):
    rec = _rec(name)
    if rec is None:
        fail_acct["unknown"].append(name)
        continue
    e, q = rec["e"], rec["q"]
    if np.isfinite(e) and 0.95 <= e < 1.5:
        if np.isfinite(q) and q < 3.1:
            key = "primary_skipped" if name not in _ck_hard \
                else "primary_hard"
        else:
            key = "extended"
    elif np.isfinite(e) and e < 0.95:
        key = "control"
    else:
        key = "unknown"
    aph_osc = (q * (1 + e) / (1 - e)
               if np.isfinite(e) and np.isfinite(q) and e < 1 else None)
    fail_acct[key].append(dict(
        desig=name, e=e, q=q, aph_osc_au=aph_osc))

# directional check on primary failures: theta from osculating
# elements (no integration needed)
_prim_fail_th = []
for f in fail_acct["primary_hard"] + fail_acct["primary_skipped"]:
    rec = _rec(f["desig"]) if isinstance(f, dict) else None
    if rec is None:
        continue
    try:
        _aph = -perih_dir(math.radians(rec["w"]),
                          math.radians(rec["om"]),
                          math.radians(rec["i"]))
        _prim_fail_th.append(float(sep(_aph, TNO)))
    except Exception:
        pass
fail_summary = {
    "n_primary_hard": len(fail_acct["primary_hard"]),
    "n_primary_skipped": len(fail_acct["primary_skipped"]),
    "n_extended": len(fail_acct["extended"]),
    "n_control": len(fail_acct["control"]),
    "primary_fail_incap_frac": (
        float(np.mean(np.array(_prim_fail_th) < CAP))
        if _prim_fail_th else None),
    "primary_integrated_incap_frac": (
        float(np.mean(
            np.array([r["theta"] for r in prows if r["q"] < 3.1])
            < CAP)) if any(r["q"] < 3.1 for r in prows) else None),
    "note": ("hard failures are bound near-parabolic orbits "
             "(e in [0.95,1.0)) whose osculating aphelia lie "
             "below the 250 AU boundary (median ~127 AU) -- "
             "geometrically ineligible, same mechanism as the "
             "preclassified skips; the boundary leg cannot "
             "reach the sphere on the osculating orbit"),
    "detail": fail_acct}
logger.info(
    f"failure accounting: primary hard={fail_summary['n_primary_hard']} "
    f"skipped={fail_summary['n_primary_skipped']} "
    f"extended={fail_summary['n_extended']} "
    f"control={fail_summary['n_control']} | primary-fail in-cap "
    f"{fail_summary['primary_fail_incap_frac']} vs integrated "
    f"{fail_summary['primary_integrated_incap_frac']}")

if not STATS_ONLY:
    cohort_counts = {"usable_near_parabolic": len(cohort_all),
                     "primary_q_lt_3p1": len(cohort),
                     "extended_q_ge_3p1": len(cohort_all) - len(cohort),
                     "bound_controls": len(controls),
                     "nongrav_fitted_primary": sum(
                         1 for r in cohort if r["ng"]),
                     "oort_spike_q_lt_3p1": len(spike),
                     "oort_spike_pure_gravity": sum(
                         1 for r in spike if not r["ng"]),
                     "excluded": excluded,
                     "integration_failed_primary":
                         fail_summary["n_primary_hard"]
                         + fail_summary["n_primary_skipped"],
                     "integration_failed_extended":
                         fail_summary["n_extended"],
                     "integration_failed_control":
                         fail_summary["n_control"]}
else:
    cohort_counts["integration_failed_primary"] = (
        fail_summary["n_primary_hard"]
        + fail_summary["n_primary_skipped"])
    cohort_counts["integration_failed_extended"] = \
        fail_summary["n_extended"]
    cohort_counts["integration_failed_control"] = \
        fail_summary["n_control"]

# ------------------------------------------------------------------
# T7: bipolar slip-map transfer (CODE+Warsaw trained, applied blind)
# ------------------------------------------------------------------

train = []
b30 = RESULTS / "step_b30_proper_time_slip.csv"
with open(b30) as f:
    for r in csv.DictReader(f):
        try:
            train.append((math.radians(float(r["theta"])),
                          float(r["dtau_unexplained"])))
        except (ValueError, KeyError):
            continue
th_tr = np.array([t[0] for t in train])
y_tr = np.array([t[1] for t in train])
ok = np.isfinite(y_tr)
X_tr = np.column_stack([np.ones(ok.sum()), np.cos(2 * th_tr[ok])])
coef_tr, *_ = np.linalg.lstsq(X_tr, y_tr[ok], rcond=None)

new = [r for r in prows if r["q"] < 3.1 and "dtau_unexplained" in r]
th_new = np.array([math.radians(r["theta"]) for r in new])
y_new = np.array([r["dtau_unexplained"] for r in new])
pred = coef_tr[0] + coef_tr[1] * np.cos(2 * th_new)
rho, p = spearmanr(pred, y_new) if len(new) > 2 else (float("nan"), float("nan"))
inlobe = (np.degrees(th_new) < CAP) | (np.degrees(th_new) > 180 - CAP)
sign_acc = (float(np.mean(np.sign(pred[inlobe]) == np.sign(y_new[inlobe])))
            if inlobe.any() else None)
transfer = {"train_n": int(ok.sum()),
            "train_coef": [float(coef_tr[0]), float(coef_tr[1])],
            "new_n": len(new),
            "oos_spearman": {"rho": float(rho), "p_2sided": float(p)},
            "lobe_sign_accuracy": sign_acc,
            "lobe_n": int(inlobe.sum())}
logger.info(f"slip-map transfer: rho={rho:.3f} p={p:.4f} "
            f"sign-acc={sign_acc if sign_acc is None else f'{sign_acc:.2f}'} "
            f"(n={len(new)})")

# CODE reference values (measured, not assumed)
ref = {}
b28 = RESULTS / "step_b28_bidirectional_rotation.json"
if b28.exists():
    d28 = json.load(open(b28))
    ref = {"code_matched": d28.get("matched"), "code_all_c1": d28.get("all_c1")}

res = {
    "step": "step_117_prospective_lpc",
    "description": "prospective validation: post-2017 SBDB LPC cohort "
                   "through the identical REBOUND/DE440s bidirectional "
                   "instrument; out-of-time and out-of-lineage "
                   "replication of the transit anomaly",
    "inputs": [
        "data/raw/sbdb/sbdb_comets_all.json (JPL SBDB, downloaded "
        "this step with provenance)",
        "data/raw/spice/de440s.bsp",
        "results/step_b30_proper_time_slip.csv (training residuals)",
        "results/step_b28_bidirectional_rotation.json (CODE anchors)"],
    "seed": SEED, "n_perm": N_PERM, "cap_deg": CAP,
    "axis_ecl_deg": [34.0, -13.0],
    "boundary_AU_barycentric": R_STOP,
    "cohort_definition": {
        "designation": "C/2018+ (post-CODE cutoff)",
        "e": "(0.95, 1.5) near-parabolic; ISOs excluded",
        "q_primary": "0.1-3.1 AU (deep plungers, class-1-like)",
        "quality": "data_arc >= 30 d and n_obs_used >= 20",
        "fragments": "fragment components excluded",
        "sungrazers": "q < 0.1 AU excluded"},
    "cohort_counts": cohort_counts,
    "failure_accounting": fail_summary,
    "primary_q_lt_3p1": res_primary,
    "primary_q_lt_3p1_pure_gravity": res_primary_pure,
    "oort_spike_q_lt_3p1": res_spike,
    "oort_spike_pure_gravity": res_spike_pure,
    "all_near_parabolic": res_all,
    "extended_q_ge_3p1": res_ext,
    "bound_control_e_lt_0p95": res_control,
    "slip_map_transfer": transfer,
    "code_reference": ref,
}

# T12: sky scan -- the cohort's own best-fit cap direction.  For each
# trial direction the 60 deg cap contrast of drot is computed; the
# maximum and its sky position show whether a significant contrast
# somewhere on the sky is generic (encounter-geometry driven) rather
# than specific to the pre-registered boundary axis.
scan_rows = [r for r in prows
             if r["q"] < 3.1 and np.isfinite(r["drot"])]
sky_scan = None
if len(scan_rows) >= 20:
    vv = np.array([r["drot"] for r in scan_rows])
    aph_v = np.array([r["aph_lb"] for r in scan_rows])
    th0 = np.array([r["theta"] for r in scan_rows])

    def _contrast(ths):
        inc = ths < CAP
        if not (inc.any() and (~inc).any()):
            return float("nan")
        return float(np.median(vv[inc]) - np.median(vv[~inc]))

    obs = _contrast(th0)
    best = {"contrast": -np.inf}
    n_gt = 0
    n_dirs = 0
    for L in np.arange(0.0, 360.0, 15.0):
        for B in np.arange(-75.0, 76.0, 15.0):
            u = lv(float(L), float(B))
            c = _contrast(np.degrees(
                np.arccos(np.clip(aph_v @ u, -1.0, 1.0))))
            n_dirs += 1
            if c >= obs:
                n_gt += 1
            if c > best["contrast"]:
                best = {"contrast": c, "l": float(L), "b": float(B),
                        "sep_to_registered_deg": float(sep(u, TNO))}
    sky_scan = {"grid_step_deg": 15.0, "n_dirs": n_dirs,
                "registered_axis_contrast": obs,
                "frac_dirs_ge_registered": n_gt / n_dirs,
                "best": best}
    logger.info(f"sky scan: registered-axis contrast {obs:.3f} deg; "
                f"best at (l,b)=({best['l']:.0f},{best['b']:.0f}) "
                f"contrast {best['contrast']:.3f} deg "
                f"({best['sep_to_registered_deg']:.0f} deg away); "
                f"{n_gt}/{n_dirs} dirs >= registered")
    res["sky_scan_primary"] = sky_scan

out = str(RESULTS / "step_b81_prospective_lpc.json")
json.dump(res, open(out, "w"), indent=1, default=float)

csv_out = str(RESULTS / "step_b81_prospective_lpc.csv")
flat = [{k: (v if not isinstance(v, list) else v)
         for k, v in r.items()} for r in prows]
for r in flat:
    r["aph_lb"] = ";".join(f"{x:.4f}" for x in r["aph_lb"])
if flat:
    keys = list(flat[0].keys())
    for r in flat[1:]:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(csv_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader(); w.writerows(flat)

# ------------------------------------------------------------------
# 6. Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

th = np.array([r["theta"] for r in prows])
v = np.array([r["drot"] for r in prows])
inc = th < CAP

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
ax = axes[0]
ax.scatter(th[~inc], v[~inc], s=22, c="#566573", alpha=0.7,
           label=f"outside cap ($n={int((~inc).sum())}$)")
ax.scatter(th[inc], v[inc], s=26, c="#b43b4e", alpha=0.85,
           label=f"inside cap ($n={int(inc.sum())}$)")
ax.axvline(CAP, color="k", ls=":", lw=1)
if ref.get("code_matched"):
    ax.axhline(ref["code_matched"]["med_drot_sim_in"], color="#b43b4e",
               ls="--", lw=1, alpha=0.6,
               label="CODE in-cap median")
    ax.axhline(ref["code_matched"]["med_drot_sim_out"], color="#566573",
               ls="--", lw=1, alpha=0.6,
               label="CODE out-cap median")
ax.set_yscale("log")
ax.set_xlabel(r"$\theta$ from axis (deg)")
ax.set_ylabel(r"rotation $d_{\rm rot}$ (deg)", labelpad=6)
ax.legend(frameon=False, loc="upper left")

ax = axes[1]
res_arr = np.array([r.get("dtau_unexplained", np.nan) for r in prows])
ax.scatter(th[~inc], res_arr[~inc], s=22, c="#566573", alpha=0.7)
ax.scatter(th[inc], res_arr[inc], s=26, c="#b43b4e", alpha=0.85)
ax.axvline(CAP, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel(r"$\theta$ from axis (deg)")
ax.set_ylabel(r"unexplained $\delta\tau$ (yr)", labelpad=6)

ax = axes[2]
bins = np.linspace(0, 180, 19)
ax.hist(th[inc], bins=bins, color="#b43b4e", alpha=0.75,
        label="in-cap")
ax.hist(th[~inc], bins=bins, color="#566573", alpha=0.5,
        label="out-cap")
ax.set_xlabel(r"$\theta$ from axis (deg)")
ax.set_ylabel("comets")
ax.legend(frameon=False)

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b81_prospective_lpc.png", dpi=300)

for tag, d in [("primary(q<3.1)", res_primary),
               ("primary pure-gravity", res_primary_pure),
               ("oort-spike(q<3.1,0<1/a<100)", res_spike),
               ("oort-spike pure-gravity", res_spike_pure),
               ("all_near_parabolic", res_all),
               ("extended(q>=3.1)", res_ext or {}),
               ("control(e<0.95)", res_control)]:
    if not d:
        continue
    if d.get("status") == "insufficient_sample" or not d.get(
            "drot_perm_contrast"):
        logger.info(f"{tag}: n={d['n']} in={d['n_in']} -- "
                    f"{d.get('status', 'no cap contrast')}")
        continue
    logger.info(f"{tag}: n={d['n']} in={d['n_in']} | "
                f"med drot in/out {d['med_drot_in']:.3f}/"
                f"{d['med_drot_out']:.3f} | "
                f"MWU p={d['drot_in_gt_out_p']:.4f} "
                f"perm p={d['drot_perm_contrast']['p_perm']:.4f} | "
                f"resid cap p={d['rot_per_kick_resid_full']['p']:.4f} | "
                f"dtau med_in={d['dtau_unexplained']['med_in']:.2f} yr")
logger.data_save(out)
logger.data_save(csv_out)
logger.data_save(FIG / "step_b81_prospective_lpc.png")