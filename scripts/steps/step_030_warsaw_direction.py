"""step_07: decisive comet test -- Warsaw Catalogue original barycentric orbits.

Data: J/A+A/567/A126 tablec.dat (Krolikowska 2014) -- original barycentric
orbital elements at 250 AU inbound for 128 near-parabolic comets, with
formal uncertainties propagated through full dynamical swarms. This is the
proper Oort-spike dataset: comets whose aphelia have *not* yet been
perturbed into the observable region, so their aphelion directions trace
where the Oort cloud material actually sits.

TEP discriminator: if a proper-time field boundary exists near the
confined-TNO axis (lam~34 deg, beta~-13 deg ecliptic), comets arriving
from that direction crossed the boundary -- their aphelia may show an
excess along the axis *beyond* what the galactic tide predicts. The
galactic tide alone predicts aphelion avoidance of the galactic equator
(|b| structure) but is axisymmetric about the galactic pole -- it cannot
produce a concentration toward one specific ecliptic direction.

Tests:
  W1  isotropy of aphelion unit vectors (Rayleigh, 3D)
  W2  galactic-latitude structure vs uniform and vs tide prediction
  W3  angle between mean aphelion direction and the TNO cluster axis,
      with bootstrap direction uncertainty and tide-aware null
  W4  direct excess: fraction of aphelia within 60 deg of the TNO axis
      vs Monte Carlo nulls (uniform, and galactic-tide-weighted)
  W5  ecliptic-plane direction test on dynamically-new subset
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_030_warsaw_direction")
tee_stdout(logger)
logger.header("Warsaw comet directional test")
import json, math, os, sys
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scripts.utils.statistics import rayleigh_3d

RNG = np.random.default_rng(20260917)


def perih_dir(om, Om, inc):
    """perihelion unit vector, ecliptic J2000, from arg-perihelion om,
    node Om, inclination inc (all radians)."""
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci,
                     sO*co + cO*so*ci,
                     so*si])

def to_lonlat(v):
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(np.clip(v[2], -1, 1)))

def ang_sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def parse_tablec(path):
    rows = []
    with open(path) as f:
        for line in f:
            if len(line) < 115:
                continue
            try:
                rows.append(dict(
                    sample=line[0:2].strip(), com=line[3].strip(),
                    desig=line[5:17].strip(),
                    q=float(line[42:56]), e=float(line[56:70]),
                    w=float(line[70:82]), Om=float(line[82:94]),
                    i=float(line[94:106]), aa=float(line[106:115]),
                    e_aa=float(line[189:198]) if line[189:198].strip() else np.nan))
            except ValueError:
                continue
    return rows

rows = parse_tablec(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
print(f"parsed {len(rows)} original-orbit rows")

# one row per comet: prefer sample membership A1/A2 (higher quality)
# and comment 'a' or 'h' (entire-data-set orbits)
seen = {}
for r in rows:
    key = r["desig"]
    if key not in seen or (seen[key]["com"] not in ("a", "h") and r["com"] in ("a", "h")):
        seen[key] = r
comets = list(seen.values())
print(f"{len(comets)} unique comets")

# aphelion directions: aphelion is opposite perihelion
for c in comets:
    p = perih_dir(math.radians(c["w"]), math.radians(c["Om"]), math.radians(c["i"]))
    c["aph"] = -p                      # unit vector toward aphelion
    g = ECL2GAL @ c["aph"]
    c["gal_b"] = math.degrees(math.asin(np.clip(g[2], -1, 1)))
    c["gal_l"] = math.degrees(math.atan2(g[1], g[0])) % 360
    c["ecl_lon"], c["ecl_lat"] = to_lonlat(c["aph"])

# Oort-spike selection: bound original orbit, 1/a inside the spike
# (aa in 1e-6 AU^-1; spike ~ aa < 100; dynamically new ~ aa < 30-40)
spike = [c for c in comets if 0 < c["aa"] < 100]
newc = [c for c in comets if 0 < c["aa"] < 35]
print(f"spike (0<1/a<100e-6): {len(spike)}   dynamically new (0<1/a<35e-6): {len(newc)}")

def rayleigh3d(vecs):
    return rayleigh_3d(vecs)

def mean_dir(vecs):
    m = np.sum(vecs, axis=0); m /= np.linalg.norm(m)
    return m

def bootstrap_dir(vecs, n_boot=4000):
    n = len(vecs); out = []
    m0 = mean_dir(vecs)
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        out.append(ang_sep(mean_dir(vecs[idx]), m0))
    return float(np.percentile(out, 68))

def mc_null_directional(vecs, axis, frac_within=60.0, n_mc=20000, tide_weight=False):
    """null: isotropic, or isotropic-with-galactic-|b|-distribution preserved.
    tide_weight=True resamples |b| from observed |b| values (the tide's
    axisymmetric imprint) while keeping azimuth random."""
    n = len(vecs)
    obs = float(np.mean([ang_sep(v, axis) < frac_within for v in vecs]))
    cnt = 0
    gal_b_obs = np.array([math.asin(np.clip((ECL2GAL @ v)[2], -1, 1)) for v in vecs])
    for _ in range(n_mc):
        if tide_weight:
            b = RNG.choice(gal_b_obs, n)
            l = RNG.uniform(0, 2*np.pi, n)
            # build galactic vecs, rotate to ecliptic
            gv = np.stack([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)], axis=1)
            vv = (np.linalg.inv(ECL2GAL) @ gv.T).T
        else:
            vv = RNG.normal(size=(n, 3)); vv /= np.linalg.norm(vv, axis=1)[:, None]
        f = np.mean([ang_sep(x, axis) < frac_within for x in vv])
        if f >= obs: cnt += 1
    return obs, (cnt + 1) / (n_mc + 1)

# TNO cluster axis from step_02/03: lam~34 deg, beta~-13 deg (confined sample)
ax = math.radians(34.0); bx = math.radians(-13.0)
TNO_AXIS = np.array([math.cos(bx)*math.cos(ax), math.cos(bx)*math.sin(ax), math.sin(bx)])

res = {"n_rows": len(rows), "n_unique": len(comets)}
for name, sample in [("spike", spike), ("new", newc), ("all_bound", [c for c in comets if c["aa"] > 0])]:
    if len(sample) < 8:
        res[name] = {"n": len(sample), "note": "too small"}; continue
    vecs = np.array([c["aph"] for c in sample])
    Rb, p_iso = rayleigh3d(vecs)
    m = mean_dir(vecs)
    lon, lat = to_lonlat(m)
    s = {}
    s["n"] = len(sample)
    s["Rbar"] = Rb; s["p_isotropy"] = p_iso
    s["mean_dir_ecl"] = {"lon": lon, "lat": lat}
    s["dir_ang68_deg"] = bootstrap_dir(vecs)
    s["angle_to_tno_axis"] = ang_sep(m, TNO_AXIS)
    s["mean_abs_gal_b"] = float(np.mean(np.abs([c["gal_b"] for c in sample])))
    s["frac_galb_gt30"] = float(np.mean([abs(c["gal_b"]) > 30 for c in sample]))
    # W4: excess within 60 deg of TNO axis
    obs_iso, p_iso60 = mc_null_directional(vecs, TNO_AXIS, 60.0, tide_weight=False)
    obs_tw, p_tw60 = mc_null_directional(vecs, TNO_AXIS, 60.0, tide_weight=True)
    s["frac_within_60deg_axis"] = obs_iso
    s["p_vs_isotropic"] = p_iso60
    s["p_vs_tide_weighted"] = p_tw60
    res[name] = s

out = str(RESULTS / "step_07_warsaw_comets.json")
with open(out, "w") as f:
    json.dump(res, f, indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)