#!/usr/bin/env python3
"""Step 093 -- perturber-zoo closure: every proposed body against the
required-mass map and the observational bounds
====================================================================

Steps 083/087/091 measured, with IAS15 + DE440s integrations of the
CODE class-1 comet sample, the point mass each candidate perturber
would need to reproduce the observed boundary rotation.  This step
closes the zoo: it places EVERY published candidate on one map --
required mass vs exclusion bound -- and adds the one real candidate
body reported to date.

Candidates / bounds
-------------------
* BB21 Planet Nine median geometry (a9 = 461 AU, step 083/087):
  required ~10^3 M_earth -- deep inside the IRAS/WISE all-sky
  exclusion (>~1 M_Jup).
* Distance ladder (step 087, a9 = 250-1000 AU): required-mass scaling
  m_req ~ a9^gamma measured on the same integrations.
* SCT25 Planet Y (a9 = 100-200 AU, step 091): required mass vs the
  Gomes+23 5-sigma ephemeris line 1.1 M_earth (d/400)^3 and vs the
  SCT25 warp-producing box 0.06-1 M_earth.
* IRAS-AKARI candidate (Phan et al. 2025, arXiv:2504.17288): a real
  source pair with candidate separation (RA, Dec) ~ (35.5, -48.9)
  deg -> ecliptic (lam, beta) ~ (0, -57) deg, inferred distance
  500-700 AU and mass 7-17 M_earth IF it is a planet.  Its sky
  position is measured against (i) the detached-TNO axis
  (lam, beta) = (49, -17) deg, (ii) the required anti-aligned
  shepherding direction (229, 17) deg, and (iii) the BB21 perihelion
  direction.  Its capacity is then measured directly: the candidate
  is inserted into the same REBOUND machinery at 600 AU toward its
  catalogued position at 7 and 17 M_earth and the injected rotation
  is compared with the observed boundary rotation.

Outputs
-------
results/step_b58_perturber_zoo.json
results/step_b58_perturber_zoo.csv
results/figures/supplementary/step_b58_perturber_zoo.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_093_perturber_zoo")
tee_stdout(logger)
logger.header("Perturber-zoo closure: candidates vs required-mass map")

import csv
import json
import math
import re
import numpy as np
from html.parser import HTMLParser

import rebound
import spiceypy as sp

rng = np.random.default_rng(20260918)

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l),
                     math.sin(b)])

def sep(a, b):
    return math.degrees(math.acos(np.clip(float(np.dot(a, b)), -1, 1)))

TNO_AXIS = lv(49.0, -17.0)
ANTIPODE = -TNO_AXIS
BB21_PERIH = lv(241.0, -15.0)     # BB21 median perihelion direction
INVAR_POLE = lv(107.0, 88.4)
TNO_TRANSIT = lv(34.0, -13.0)     # transit-axis convention (Phase-3)

# ------------------------------------------------------------------
# A. IRAS-AKARI candidate geometry (Phan et al. 2025)
# ------------------------------------------------------------------

# RA=35.5, Dec=-48.9 -> ecliptic (J2000, eps=23.4393 deg)
ra, dec = math.radians(35.5), math.radians(-48.9)
eq = np.array([math.cos(dec)*math.cos(ra), math.cos(dec)*math.sin(ra),
               math.sin(dec)])
EPS = math.radians(23.4392911)
RX = np.array([[1, 0, 0],
               [0, math.cos(EPS), math.sin(EPS)],
               [0, -math.sin(EPS), math.cos(EPS)]])
ec = RX @ eq
cand_lam = math.degrees(math.atan2(ec[1], ec[0])) % 360
cand_bet = math.degrees(math.asin(ec[2]))
CAND = lv(cand_lam, cand_bet)

geom = dict(
    candidate_eq=dict(ra_deg=35.5, dec_deg=-48.9),
    candidate_ecl=dict(lam_deg=round(cand_lam, 1),
                       bet_deg=round(cand_bet, 1)),
    dist_range_AU=[500, 700], mass_range_earth=[7, 17],
    sep_vs_tno_axis_deg=round(sep(CAND, TNO_AXIS), 1),
    sep_vs_antipode_deg=round(sep(CAND, ANTIPODE), 1),
    sep_vs_bb21_perih_deg=round(sep(CAND, BB21_PERIH), 1),
    sep_vs_transit_axis_deg=round(sep(CAND, TNO_TRANSIT), 1),
    note="IRAS/AKARI slow-mover pair (Phan+25); not confirmed -- "
         "the AKARI-epoch counterpart was absent at the predicted "
         "position in follow-up inspection")
logger.info("IRAS-AKARI candidate at ecliptic "
            f"(lam,beta)=({cand_lam:.1f},{cand_bet:.1f}); "
            f"{geom['sep_vs_tno_axis_deg']} deg from the TNO axis, "
            f"{geom['sep_vs_antipode_deg']} deg from the required "
            "anti-aligned direction")

# ------------------------------------------------------------------
# B. Capacity of the real candidate: REBOUND insertion at 600 AU
# ------------------------------------------------------------------

class TP(HTMLParser):
    def __init__(self):
        super().__init__(); self.rows = []; self.cur = []; self.buf = ""; self.in_td = False
    def handle_starttag(self, t, a):
        if t == "tr": self.cur = []
        elif t == "td": self.in_td = True; self.buf = ""
    def handle_endtag(self, t):
        if t == "td": self.in_td = False; self.cur.append(self.buf.strip())
        elif t == "tr" and self.cur: self.rows.append(self.cur)
    def handle_data(self, d):
        if self.in_td: self.buf += d

def parse_code(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            out[r[0].strip()] = dict(desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                T=r[7].strip(),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

SPK = DATA_RAW / "spice" / "de440s.bsp"
if not SPK.exists():
    raise FileNotFoundError(f"JPL ephemeris missing: {SPK}")
sp.furnsh(str(SPK))

AU_KM   = 149597870.7
DAY_YR  = 365.25
GM_SUN = 1.32712440018e11
GM = {"1": 2.2031868551e4, "2": 3.2485859200e5, "3": 4.0350323562e5,
      "4": 4.2828375814e4, "5": 1.2671276480e8, "6": 3.7940626000e7,
      "7": 5.7945490100e6, "8": 6.8365271006e6, "9": 1.0868657e3}
PLANET_IDS = list(GM.keys())
MU = 4 * math.pi ** 2
M_EARTH_SUN = 3.0035e-6
R_STOP, T_MAX, DT_OUT = 255.0, -20000.0, 1.0

def jd_tt(y, m, d):
    if m <= 2: y -= 1; m += 12
    A = y // 100; B = 2 - A + A // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5

def perihelion_et(T):
    m = re.match(r"(\d{4})\s+(\d{1,2})\s+(\d+\.?\d*)", T)
    if not m:
        return None
    return (jd_tt(int(m.group(1)), int(m.group(2)), float(m.group(3))) - 2451545.0) * 86400.0

def body_state(body, et):
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "0")
    return RX @ np.array(st[:3]) / AU_KM, RX @ np.array(st[3:]) / AU_KM * 86400 * DAY_YR

def init_sim(et):
    sim = rebound.Simulation(); sim.G = MU
    ps, vs = body_state("10", et)
    sim.add(x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2], m=1.0)
    for b in PLANET_IDS:
        pp, vv = body_state(b, et)
        sim.add(x=pp[0], y=pp[1], z=pp[2], vx=vv[0], vy=vv[1], vz=vv[2],
                m=GM[b] / GM_SUN)
    return sim, ps, vs

def state_at_periapsis(ro):
    w_, O_, i_ = map(math.radians, (ro["w"], ro["Om"], ro["i"]))
    ph = perih_dir(w_, O_, i_)
    hh = np.array([math.sin(i_) * math.sin(O_), -math.sin(i_) * math.cos(O_), math.cos(i_)])
    th = np.cross(hh, ph); th /= np.linalg.norm(th)
    return ro["q"] * ph, math.sqrt(MU * (1 + ro["e"]) / ro["q"]) * th

def boundary_orbit(r_rel, v_rel, mtot):
    mu = MU * mtot
    r = np.linalg.norm(r_rel)
    h = np.cross(r_rel, v_rel)
    evec = (np.cross(v_rel, h) / mu) - r_rel / r
    en = np.linalg.norm(evec)
    phat = evec / en if en > 1e-12 else r_rel / r
    return phat, -2 * (v_rel.dot(v_rel) / 2 - mu / r) / mu * 1e6

def integrate(ro, et, p9=None):
    sim, ps, vs = init_sim(et)
    if p9 is not None:
        sim.add(x=p9[0][0], y=p9[0][1], z=p9[0][2],
                vx=p9[1][0], vy=p9[1][1], vz=p9[1][2], m=p9[2])
    rvec, vvec = state_at_periapsis(ro)
    sim.add(x=rvec[0] + ps[0], y=rvec[1] + ps[1], z=rvec[2] + ps[2],
            vx=vvec[0] + vs[0], vy=vvec[1] + vs[1], vz=vvec[2] + vs[2])
    nc = sim.N - 1
    sim.integrator = "ias15"
    sim.exit_min_distance = 0.001  # collision scale: bound IAS15 against step collapse
    t = -DT_OUT
    # t is negative for the backward leg; compare elapsed time, not signed t.
    while abs(t) < abs(T_MAX):
        sim.integrate(t, exact_finish_time=0)
        p = sim.particles
        r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y,
                          p[nc].z - p[0].z])
        if np.linalg.norm(r_rel) >= R_STOP:
            break
        t -= DT_OUT
    else:
        return None
    mtot = sum(pp.m for pp in sim.particles)
    p = sim.particles
    r_rel = np.array([p[nc].x - p[0].x, p[nc].y - p[0].y, p[nc].z - p[0].z])
    v_rel = np.array([p[nc].vx - p[0].vx, p[nc].vy - p[0].vy,
                      p[nc].vz - p[0].vz])
    phat, aa_out = boundary_orbit(r_rel, v_rel, mtot)
    return dict(phat=phat, aa_out=aa_out)

# Candidate state: 600 AU toward the catalogued position, on an
# orbit whose apoapsis sits there (e=0.3 -> a=462 AU, residence
# weighted to apocentre -- the generic configuration for a distant
# body observed near apoapsis).
D_CAND = 600.0
E_CAND = 0.30
A_CAND = D_CAND / (1.0 + E_CAND)

def cand_state(m_earth):
    """Candidate at apoapsis toward CAND: r = D*CAND, v tangential."""
    r = D_CAND * CAND
    # tangential velocity at apoapsis in an arbitrary but fixed plane:
    # take the plane containing CAND and the ecliptic pole tilted to
    # i=30 deg toward the node at lam=270 deg.
    n = lv(270.0, 60.0)                    # orbit pole (i=30, Om=0)
    t = np.cross(n, CAND); t /= np.linalg.norm(t)
    v_apo = math.sqrt(MU * (1.0 - E_CAND) / (A_CAND * (1.0 + E_CAND)))
    return r, v_apo * t

osc  = parse_code(str(DATA_RAW / "code" / "code_osculating.html"))
orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse_code(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip() for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")) if len(l) > 115}

sample = []
for k, ro in orig.items():
    if k not in fut: continue
    if not (0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b")) or k in warsaw:
        continue
    if k not in osc: continue
    et = perihelion_et(osc[k]["T"])
    if et is None: continue
    sample.append((k, ro, osc[k], et))

logger.info(f"sample: {len(sample)} class-1 CODE comets; candidate "
            f"inserted at {D_CAND:.0f} AU, e={E_CAND}, "
            "m=7 and 17 M_earth")

def _run_comet(k, ro, oo, et):
    base = integrate(oo, et, None)
    if base is None:
        return None
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                   math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]),
                   math.radians(fut[k]["i"]))
    rec = dict(desig=k, theta=sep(-po, TNO_TRANSIT),
               d_of=sep(-po, -pf))
    for m_e in (7.0, 17.0):
        rc, vc = cand_state(m_e)
        rp = integrate(oo, et, (rc, vc, m_e * M_EARTH_SUN))
        key = f"m{int(m_e)}"
        rec[f"{key}_drot"] = sep(base["phat"], rp["phat"]) if rp else float("nan")
        rec[f"{key}_mreq"] = (m_e * rec["d_of"] / rec[f"{key}_drot"]
                              if np.isfinite(rec[f"{key}_drot"])
                              and rec[f"{key}_drot"] > 0
                              else float("nan"))
    return rec


def _work_comet(job):
    k, ro, oo, et = job
    try:
        return _run_comet(k, ro, oo, et)
    except Exception as e:
        return ("__error__", k, str(e))


from scripts.utils.parallel import default_workers as _default_workers

if "--workers" in _sys.argv:
    N_WORK = max(1, int(_sys.argv[_sys.argv.index("--workers") + 1]))
else:
    N_WORK = _default_workers()

if N_WORK > 1 and len(sample) > 1:
    import multiprocessing as mp
    ctx = mp.get_context("fork") if _sys.platform != "win32" \
        else mp.get_context("spawn")

    def _init():
        # forked children inherit the parent's BSP fd; concurrent spkezr
        # reads through a shared descriptor corrupt each other -- reopen
        # the kernel so each worker holds its own file handle.
        sp.kclear()
        sp.furnsh(str(SPK))

    with ctx.Pool(min(N_WORK, len(sample)), initializer=_init) as pool:
        recs = pool.map(_work_comet, sample)
else:
    recs = [_work_comet(s) for s in sample]

rows = []
failed = []
for s, rec in zip(sample, recs):
    if rec is None:
        failed.append(s[0])
    elif isinstance(rec, tuple) and rec[0] == "__error__":
        failed.append(rec[1]); logger.warning(f"{rec[1]}: {rec[2]}")
    else:
        rows.append(rec)

logger.info(f"integrated {len(rows)} comets x 3 runs; "
            f"{len(failed)} failed")

inc = np.array([r["theta"] < 60 for r in rows])
cand_run = {}
for m_e in (7, 17):
    dr = np.array([r[f"m{m_e}_drot"] for r in rows])
    mr = np.array([r[f"m{m_e}_mreq"] for r in rows])
    ok = np.isfinite(mr)
    cand_run[f"m{m_e}"] = dict(
        med_drot=float(np.nanmedian(dr)),
        m_req_median_in=float(np.nanmedian(mr[ok & inc])),
        m_req_p05_in=float(np.nanpercentile(mr[ok & inc], 5)),
        frac_in_req_gt_318=float(np.mean(mr[ok & inc] > 318)))
    logger.info(f"candidate at {m_e} M_earth: med drot "
                f"{cand_run[f'm{m_e}']['med_drot']:.5f} deg | required "
                f"mass in-cap med {cand_run[f'm{m_e}']['m_req_median_in']:.0f} "
                "M_earth")

# ------------------------------------------------------------------
# C. Closure table: hypothesis x required mass x exclusion bound
# ------------------------------------------------------------------

def loadj(name):
    p = RESULTS / name
    return json.load(open(p)) if p.exists() else None

b52 = loadj("step_b52_distance_ladder.json")
b48 = loadj("step_b48_mass_ladder.json")
b56 = loadj("step_b56_inner_perturber.json")

closure = []
if b48:
    closure.append(dict(
        hypothesis="BB21 Planet Nine (a=461 AU, nominal 6.9 M_earth)",
        source="step_b48_mass_ladder",
        required_mass_earth=b48["required_mass"]["m_req_median_in"],
        bound="IRAS/WISE all-sky ~1 M_Jup = 318 M_earth",
        excluded=None))
if b52:
    for key, r in b52["rungs"].items():
        closure.append(dict(
            hypothesis=f"anti-aligned perturber a={r['a9_AU']:.0f} AU",
            source="step_b52_distance_ladder",
            required_mass_earth=r["m_req_median_in"],
            bound="IRAS/WISE all-sky ~1 M_Jup = 318 M_earth",
            excluded=bool(r["m_req_p05_in"] > 318)))
if b56:
    for key, r in b56["rungs"].items():
        closure.append(dict(
            hypothesis=f"{r['geometry']} inner perturber "
                       f"a={r['a9_AU']:.0f} AU",
            source="step_b56_inner_perturber",
            required_mass_earth=r["m_req_median_in"],
            bound=f"ephemeris 5-sigma {r['eph_bound_5sig_earth']:.3f} "
                  "M_earth (Gomes+23); SCT25 box 0.06-1 M_earth",
            excluded=bool(r["frac_in_below_eph_bound"] < 0.05)))
for m_e in (7, 17):
    cr = cand_run[f"m{m_e}"]
    closure.append(dict(
        hypothesis=f"IRAS-AKARI candidate ({m_e} M_earth, 600 AU, "
                   f"lam~{cand_lam:.0f} bet~{cand_bet:.0f})",
        source="this step",
        required_mass_earth=cr["m_req_median_in"],
        bound="candidate's own mass estimate 7-17 M_earth",
        excluded=bool(cr["m_req_p05_in"] > 17)))

out = dict(meta=dict(
    purpose="single-map closure of all proposed point-mass "
            "explanations against the required-mass ladder",
    axes=dict(tno_axis="(lam,beta)=(49,-17)", antipode="(229,17)",
              bb21_perih="(241,-15)", transit_axis="(34,-13)"),
    references=["Phan et al. 2025, arXiv:2504.17288 (IRAS-AKARI)",
                "Gomes et al. 2023, PSJ 4:66 (ephemeris bound)",
                "Siraj, Chyba & Tremaine 2025, MNRAS Lett slaf091",
                "Brown & Batygin 2021, AJ 161:27"]),
    iras_akari_geometry=geom,
    iras_akari_capacity=cand_run,
    closure=closure)

json.dump(out, open(RESULTS / "step_b58_perturber_zoo.json", "w"),
          indent=1, default=float)
with open(RESULTS / "step_b58_perturber_zoo.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
logger.data_save(RESULTS / "step_b58_perturber_zoo.json")
logger.data_save(RESULTS / "step_b58_perturber_zoo.csv")
# ------------------------------------------------------------------
# Figure: required-mass map
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))

ax = axes[0]
# sky positions of the named directions
dirs = [("TNO axis (49,-17)", 49, -17, "crimson", "*"),
        ("antipode (229,17)", 229, 17, "0.5", "*"),
        ("IRAS-AKARI cand.", cand_lam, cand_bet, "navy", "D"),
        ("BB21 perih (241,-15)", 241, -15, "darkgreen", "s"),
        ("transit axis (34,-13)", 34, -13, "purple", "^")]
for lb, l_, b_, c_, m_ in dirs:
    ax.plot(l_, b_, m_, ms=11, color=c_, label=lb)
ax.set(xlabel="ecliptic longitude (deg)", ylabel="ecliptic latitude (deg)",
       xlim=(0, 360), ylim=(-90, 90),
       title="candidate sky positions vs the boundary axes")
ax.grid(alpha=0.25); ax.legend(fontsize=7, loc="upper right")

ax = axes[1]
d = np.logspace(np.log10(80), np.log10(1200), 200)
ax.plot(d, 1.1 * (d / 400.0) ** 3, "b-", lw=1.5,
        label="ephemeris 5$\\sigma$ line (Gomes+23)")
ax.fill_between(d, 318, 1e6, color="steelblue", alpha=0.08)
ax.axhline(318, color="steelblue", ls=":", lw=1.2,
           label="1 $M_{Jup}$ (IRAS/WISE all-sky)")
ax.axhspan(0.06, 1.0, color="orange", alpha=0.12,
           label="SCT25 Planet-Y box")
if b52:
    a9 = [r["a9_AU"] for r in b52["rungs"].values()]
    mr = [r["m_req_median_in"] for r in b52["rungs"].values()]
    ax.plot(a9, mr, "o-", color="crimson", ms=7,
            label="required mass (anti-aligned, in-cap median)")
if b56:
    for gname, col in (("anti-aligned", "crimson"),
                       ("warp-plane", "darkgreen")):
        aa = [r["a9_AU"] for r in b56["rungs"].values()
              if r["geometry"] == gname]
        mm = [r["m_req_median_in"] for r in b56["rungs"].values()
              if r["geometry"] == gname]
        ax.plot(aa, mm, "s--" if gname == "warp-plane" else "o-",
                color=col, ms=6, alpha=0.8,
                label=f"required mass ({gname}, inner)")
for m_e, mk in ((7, "D"), (17, "^")):
    ax.plot(D_CAND, cand_run[f"m{m_e}"]["m_req_median_in"], mk,
            color="navy", ms=9,
            label=f"IRAS cand. required ({m_e} $M_\\oplus$ inserted)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set(xlabel="perturber distance / semimajor axis (AU)",
       ylabel="mass ($M_\\oplus$)",
       title="required mass vs exclusion bounds")
ax.legend(fontsize=6.5, loc="upper left")

fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b58_perturber_zoo.png", dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b58_perturber_zoo.png')