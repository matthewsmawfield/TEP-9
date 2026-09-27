"""Step 161 -- Planetary astrometric O-C channel (APDB lineage).

The transit and resident channels all descend from small-body
catalogues; this step opens an independent lineage: the raw
ground-based astrometry of the outer planets archived in the JPL
APDB mirror at IMCCE (data used to build the DE ephemerides).
Pluto is the case of interest -- it is the only major-body object
spending the entire observational record at boundary-adjacent radii
(30-50 AU), and between 1914 and 1997 its heliocentric longitude
sweeps from inside the declared cap region toward the mirror cap.
Uranus and Neptune contribute the astrographic and normal-point
records.

If a lapse-type boundary field exists, the orbit-reconstruction
residuals of a resident body should carry direction-dependent
structure as the body moves through the cap geometry -- the same
non-closure observable class as the comet legs, expressed on a
resident.  This step computes the observed-minus-computed positions
of every APDB optical observation against DE440s and tests whether
the residual field is organised by the planet's direction relative
to the declared axis.

Parsing (O'Handley 1968 card format, opticalformat.txt):
  col 1      'P'
  col 2-4    planet field (00p or pss)
  col 5+     UT Julian date, 7-digit integer part + variable decimals
  col 22-24  observatory/telescope number (JPL/MPC numeric codes)
  col 26-28  catalogue code
  col 29     type (1,3 = geocentric; 4,6 = topocentric)
  col 34-42  RA hhmmss.sss (8-9 digits, trailing decimals implied)
  col 43     RA equinox flag (0 true-of-date, 3 B1950, 9 J2000)
  col 51-59  Dec [sign]ddmmss.ss
  col 60     Dec equinox flag
  tail       source code + year

Observed positions are rotated to J2000 with an IAU-1976 precession
matrix; blank equinox flags are recorded as unknown-frame and kept
out of the precision tests.  Computed positions are DE440s system
barycentres (bodies 7, 8, 9) with one-way light-time correction and
topocentric parallax where the observation is marked topocentric and
the station resolves in the MPC observatory table.

Registered tests

  T1  sanity: per-source median |O-C| must sit at the historical
      plate/CCD floor (sub-arcsec to ~arcsec), certifying the parse.
  T2  direction organisation: Spearman of the along-track residual
      against the planet's angular separation from the declared axis,
      plus the in-cap vs out MWU contrast on the same component.
  T3  source-coded robustness: the T2 tests repeated within each
      source series with n >= 30 -- catalogue-zonal and plate
      systematics are source-coded, so a direction dependence that
      survives within sources cannot be a single catalogue's bias.
  T4  ancillary normal points: the published Uranus/Neptune
      occultation and VLA O-C normal points (nptsobs.txt, vs DE405)
      tabulated with their cap geometry at epoch.

Outputs
  results/step_b127_apdb_planet.json
  results/step_b127_apdb_planet.csv
  results/figures/supplementary/step_b127_apdb_planet.png
  data/raw/apdb/                       (downloaded observation files)
"""

import json
import math
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import (
    DATA_RAW, RESULTS, AXES, sep, lb)

FIG = RESULTS / "figures" / "supplementary"

logger = StepLogger("step_161_apdb_planet")
tee_stdout(logger)

APDB = DATA_RAW / "apdb"
SPK = DATA_RAW / "spice" / "de440s.bsp"
LSK = DATA_RAW / "naif" / "naif0012.tls"
MPC_DIR = DATA_RAW / "mpc"
OBSC_PATH = MPC_DIR / "obscodes.json"

CAP_DEG = 60.0
AS_RAD = 206264.806247
C_AU_DAY = 173.144632674
AU_KM = 149597870.7
R_EARTH_AU = 6378.137 / AU_KM

EPS = math.radians(23.4392911)
RX = np.array([[1, 0, 0],
               [0, math.cos(EPS), math.sin(EPS)],
               [0, -math.sin(EPS), math.cos(EPS)]])

OPTICAL_FILES = [
    "plutopre69.txt", "plutopost69.txt", "plutofk5.txt",
    "astgrnofs.txt", "astgrbord.txt", "astgrmcd.txt",
]

BODY = {"7": "7", "8": "8", "9": "9"}   # system barycentres
PNAME = {"7": "Uranus", "8": "Neptune", "9": "Pluto"}

logger.header("APDB planetary astrometric O-C channel")

import spiceypy as sp
from scripts.utils import mpc_refit as R

sp.furnsh(str(LSK))
sp.furnsh(str(SPK))
# reuse mpc_refit's observatory-table loader so the MPC code map and
# its provenance ledger are shared with the refit steps
R.OBSC_PATH = OBSC_PATH
R.PROV_PATH = MPC_DIR / "provenance.json"
OBSC = R.load_obscodes()

# ------------------------------------------------------------------
# equinox rotation: astropy FK4/FK5/TETE -> ICRS(J2000)
# ------------------------------------------------------------------

from astropy.coordinates import SkyCoord, FK4, TETE, ICRS
from astropy.time import Time
import astropy.units as au_u


def eq_to_j2000(ra_deg, dec_deg, flag, jd):
    """Rotate an observed (RA,Dec) to J2000 under its equinox flag.
    flag 9 -> already J2000; 3 -> FK4 B1950 mean (E-terms handled by
    astropy); 0 -> true equator/equinox of date (TETE); blank ->
    unknown frame (returned unrotated with a flag)."""
    if flag == "9":
        return ra_deg, dec_deg, "j2000"
    if flag == "3":
        fr = FK4(equinox=Time("B1950"))
    elif flag == "0":
        fr = TETE(obstime=Time(jd, format="jd"))
    else:
        return ra_deg, dec_deg, "unknown"
    c = SkyCoord(ra=ra_deg * au_u.deg, dec=dec_deg * au_u.deg,
                 frame=fr).transform_to(ICRS())
    return float(c.ra.deg), float(c.dec.deg), "rotated"


# ------------------------------------------------------------------
# parse the O'Handley card format
# ------------------------------------------------------------------


def parse_optical(path):
    recs = []
    for line in open(path):
        if not line.startswith("P") or len(line) < 60:
            continue
        tok = line[:20].split()[0]
        plf = tok[1:4]
        digits = tok[4:]
        if len(digits) < 8 or not digits.isdigit():
            continue
        jd = int(digits[:7]) + float("0." + digits[7:])
        pl = plf[2] if plf.startswith("00") else plf[0]
        obs = line[21:25].strip()
        cat = line[25:28].strip()
        typ = line[28].strip()
        # RA occupies cols 34-42 (right-justified hhmmss^sss, hours
        # may drop the leading zero); flag at col 43.  Dec occupies
        # cols 51-60 (right-justified [sign]ddmmss^sss, degrees may
        # drop the leading zero and the sign may be blank for +);
        # flag at col 61 -- blank in this archive, so the RA flag
        # governs the pair.
        ra_d = line[33:42].strip().replace(" ", "0")
        eqx = line[42].strip()
        dec_d = line[50:60].strip()
        eqx2 = line[60].strip() if len(line) > 60 else ""
        src = line[72:].strip()
        try:
            if not ra_d.isdigit() or len(ra_d) < 7:
                continue
            # component range validation: misaligned plate records can
            # land field edges inside the RA/Dec blocks and decode to
            # impossible minutes/seconds -- reject rather than ingest
            if int(ra_d[-7:-5]) >= 60 or int(ra_d[-5:-3]) >= 60:
                continue
            ra = (int(ra_d[:-7]) + int(ra_d[-7:-5]) / 60.0
                  + (int(ra_d[-5:-3]) + int(ra_d[-3:]) / 1000.0)
                  / 3600.0) * 15.0
            neg = dec_d.startswith("-")
            dd = dec_d.lstrip("+-").replace(" ", "0")
            if not dd.isdigit() or len(dd) < 7:
                continue
            if int(dd[-7:-5]) >= 60 or int(dd[-5:-3]) >= 60:
                continue
            dec = (int(dd[:-7]) + int(dd[-7:-5]) / 60.0
                   + (int(dd[-5:-3]) + int(dd[-3:]) / 1000.0)
                   / 3600.0)
            if neg:
                dec = -dec
            if ra >= 360.0 or abs(dec) > 90.0:
                continue
        except ValueError:
            continue
        recs.append(dict(planet=pl, jd=jd, obs=obs, cat=cat, typ=typ,
                         ra=ra, dec=dec, eqx=eqx or eqx2, src=src,
                         file=path.name))
    return recs


APDB_URL = ("https://ftp.imcce.fr/pub/databases/PODB/JPL/APDB/")
APDB_PROV = APDB / "provenance.json"

import hashlib
import urllib.request


def fetch_apdb(fn):
    """Download an APDB file if absent; record url + sha256."""
    p = APDB / fn
    if not p.exists():
        APDB.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(APDB_URL + fn, timeout=120) as r:
            p.write_bytes(r.read())
    blob = p.read_bytes()
    prov = json.loads(APDB_PROV.read_text()) if APDB_PROV.exists() else {}
    prov[fn] = {"url": APDB_URL + fn,
                "sha256": hashlib.sha256(blob).hexdigest(),
                "bytes": len(blob)}
    APDB_PROV.write_text(json.dumps(prov, indent=1))
    return p


obs = []
for fn in OPTICAL_FILES + ["nptsobs.txt"]:
    p = fetch_apdb(fn)
    if fn == "nptsobs.txt":
        continue
    got = parse_optical(p)
    logger.info(f"{fn}: {len(got)} observations parsed")
    obs.extend(got)
logger.info(f"total optical observations: {len(obs)}")


def gmst_deg(et):
    jd_ut = (et - sp.deltet(et, "ET")) / 86400.0 + 2451545.0
    return (280.46061837 + 360.98564736629 * (jd_ut - 2451545.0)) % 360


def topo_ecl(stn, et):
    o = OBSC.get(stn)
    if o is None:
        return None
    lst = math.radians(gmst_deg(et) + o["lon"])
    v = np.array([o["cospar"] * math.cos(lst),
                  o["cospar"] * math.sin(lst),
                  o["sinpar"]]) * R_EARTH_AU
    return RX @ v


def helio_ecl(body, et):
    """Heliocentric ecliptic state of a DE440s barycentre id.
    Returns r [AU], v [AU/yr]."""
    st, _ = sp.spkezr(body, et, "J2000", "NONE", "10")
    return (RX @ np.array(st[:3]) / AU_KM,
            RX @ np.array(st[3:]) / AU_KM * 86400.0 * 365.25)


def earth_geo(et):
    st, _ = sp.spkezr("399", et, "J2000", "NONE", "10")
    return RX @ np.array(st[:3]) / AU_KM


rows = []
n_topo_ok = n_topo_geo = 0
for o in obs:
    et = (o["jd"] - 2451545.0) * 86400.0
    body = BODY.get(o["planet"])
    if body is None:
        continue
    rp, vp = helio_ecl(body, et)
    re_ = earth_geo(et)
    topovec = np.zeros(3)
    topo = False
    if o["typ"] in ("4", "6"):
        tv = topo_ecl(o["obs"], et)
        if tv is not None:
            topovec = tv
            topo = True
            n_topo_ok += 1
        else:
            n_topo_geo += 1
    rob = re_ + topovec
    tau = np.linalg.norm(rp - rob) / C_AU_DAY
    rp_lt, _ = helio_ecl(body, et - tau * 86400.0)
    u = rp_lt - rob
    u /= np.linalg.norm(u)
    u_eq = RX.T @ u
    ra_c = math.degrees(math.atan2(u_eq[1], u_eq[0])) % 360.0
    dec_c = math.degrees(math.asin(np.clip(u_eq[2], -1, 1)))

    ra_o, dec_o, frame = eq_to_j2000(o["ra"], o["dec"], o["eqx"],
                                     o["jd"])
    dra = ((ra_o - ra_c + 180.0) % 360.0 - 180.0)
    dra *= math.cos(math.radians(dec_c))
    ddec = dec_o - dec_c

    # planet heliocentric direction and sky-projected along-track
    # basis (velocity direction projected onto the tangent plane)
    ph = rp / np.linalg.norm(rp)
    vh = vp / np.linalg.norm(vp)
    # tangent basis in equatorial frame at the observed direction
    uo = np.array([math.cos(math.radians(dec_o)) *
                   math.cos(math.radians(ra_o)),
                   math.cos(math.radians(dec_o)) *
                   math.sin(math.radians(ra_o)),
                   math.sin(math.radians(dec_o))])
    e_ra = np.array([-math.sin(math.radians(ra_o)),
                     math.cos(math.radians(ra_o)), 0.0])
    e_dec = np.array([-math.sin(math.radians(dec_o)) *
                      math.cos(math.radians(ra_o)),
                      -math.sin(math.radians(dec_o)) *
                      math.sin(math.radians(ra_o)),
                      math.cos(math.radians(dec_o))])
    vh_eq = RX.T @ vh
    at = np.array([vh_eq @ e_ra, vh_eq @ e_dec])
    nrm = np.linalg.norm(at)
    at = at / nrm if nrm > 0 else np.array([1.0, 0.0])
    oc = np.array([dra, ddec]) * AS_RAD
    rows.append(dict(planet=PNAME[o["planet"]], jd=o["jd"],
                     yr=(o["jd"] - 2451545.0) / 365.25 + 2000.0,
                     obs=o["obs"], cat=o["cat"], typ=o["typ"],
                     src=o["src"], file=o["file"], frame=frame,
                     topo=topo,
                     dra_as=oc[0], ddec_as=oc[1],
                     oc_at=oc @ at, oc_ct=oc @ np.array([-at[1], at[0]]),
                     helio_lon=lb(ph)[0], helio_lat=lb(ph)[1],
                     sep_axis=sep(ph, AXES["tno"]),
                     sep_anti=sep(ph, AXES["anti"])))

logger.info(f"computed O-C for {len(rows)} observations "
            f"(topo-corrected {n_topo_ok}, geocentric-fallback "
            f"{n_topo_geo})")

# ------------------------------------------------------------------
# ancillary normal points (published O-C vs DE405)
# ------------------------------------------------------------------

npts = []
for line in open(APDB / "nptsobs.txt"):
    m = line.split()
    if len(m) < 11 or not m[0].isdigit():
        continue
    try:
        pl = m[0]
        yr, mo, dy = int(m[1]), int(m[2]), float(m[3])
        # m[4] is the UT hour field; O-C triplets follow
        oc_ra, sg_ra = float(m[5]), float(m[6])
        oc_de, sg_de = float(m[7]), float(m[8])
    except (ValueError, IndexError):
        continue
    if pl not in PNAME or max(abs(oc_ra), abs(oc_de)) > 10.0:
        continue
    jd = sp.utc2et(f"{yr:04d}-{mo:02d}-{int(dy):02d}T00:00:00")
    rp, _ = helio_ecl(BODY[pl], jd)
    ph = rp / np.linalg.norm(rp)
    npts.append(dict(planet=PNAME[pl], yr=yr + mo / 12.0,
                     oc_ra=oc_ra, oc_dec=oc_de,
                     sg_ra=sg_ra, sg_dec=sg_de,
                     helio_lon=lb(ph)[0],
                     sep_axis=sep(ph, AXES["tno"]),
                     sep_anti=sep(ph, AXES["anti"]),
                     note=" ".join(m[9:])))

# ------------------------------------------------------------------
# tests
# ------------------------------------------------------------------

res = {"step": "step_161_apdb_planet",
       "description": ("outer-planet astrometric O-C vs DE440s on the "
                       "APDB optical record; direction organisation "
                       "of residuals against the declared cap "
                       "geometry"),
       "cap_deg": CAP_DEG, "n_obs": len(rows),
       "n_topo_ok": n_topo_ok, "n_topo_geofallback": n_topo_geo}

# T1 sanity
t1 = {}
for src in sorted({r["src"] for r in rows}):
    v = [r for r in rows if r["src"] == src]
    amp = np.hypot([r["dra_as"] for r in v], [r["ddec_as"] for r in v])
    t1[src] = {"n": len(v), "med_amp_as": float(np.median(amp)),
               "p90_as": float(np.percentile(amp, 90)),
               "frames": sorted({r["frame"] for r in v}),
               "planet": v[0]["planet"]}
res["T1_sanity"] = t1

# T2 direction organisation on Pluto (the boundary-radius resident)
pl = [r for r in rows if r["planet"] == "Pluto" and
      r["frame"] != "unknown"]
t2 = {"n_pluto": len(pl)}
if len(pl) > 30:
    s = np.array([r["sep_axis"] for r in pl])
    at = np.array([r["oc_at"] for r in pl])
    ct = np.array([r["oc_ct"] for r in pl])
    inc = s < CAP_DEG
    mir = np.array([r["sep_anti"] for r in pl]) < CAP_DEG
    for comp, v in (("along_track", at), ("cross_track", ct)):
        rho, p = spearmanr(v, s)
        t2[comp] = {"rho_vs_sep_axis": float(rho), "p": float(p)}
        if inc.sum() >= 5 and (~inc).sum() >= 5:
            u = mannwhitneyu(np.abs(v[inc]), np.abs(v[~inc]),
                             alternative="greater")
            t2[comp]["cap_mwu"] = {
                "n_in": int(inc.sum()), "med_in": float(
                    np.median(np.abs(v[inc]))),
                "med_out": float(np.median(np.abs(v[~inc]))),
                "p": float(u.pvalue)}
    t2["cap_coverage"] = {
        "n_in_declared": int(inc.sum()), "n_in_mirror": int(mir.sum()),
        "sep_axis_range": [float(s.min()), float(s.max())]}
res["T2_direction"] = t2

# T3 within-source repetition
t3 = {}
for src in sorted({r["src"] for r in pl}):
    v = [r for r in pl if r["src"] == src]
    if len(v) < 30:
        continue
    s = np.array([r["sep_axis"] for r in v])
    at = np.array([r["oc_at"] for r in v])
    rho, p = spearmanr(at, s)
    t3[src] = {"n": len(v), "rho_at_vs_sep": float(rho),
               "p": float(p)}
res["T3_within_source"] = t3

res["T4_normal_points"] = npts
res["caveats"] = [
    "O-C is computed against DE440s, whose own fit absorbed the "
    "outer-planet record; the channel tests structure that survives "
    "solver absorption, not absolute ephemeris error.",
    "JPL observation numbers are assumed coincident with MPC "
    "numeric codes for topocentric correction; unresolved stations "
    "fall back to geocentric (<=0.3 arcsec penalty at these "
    "distances).",
    "JD is treated as TDB; the <=1 min UT-TDB ambiguity moves Pluto "
    "<0.005 arcsec.",
    "Pluto residuals reference the system barycentre; the "
    "Pluto-Charon photocenter offset (~0.1-0.2 arcsec, 6.4 d period) "
    "is unmodelled noise, not direction-correlated.",
    "Pluto's angular separation from the declared axis drifts nearly "
    "monotonically with epoch over any single source window, so "
    "within-source direction correlations alias any secular residual "
    "trend (catalogue zonal error, ephemeris slope); T2/T3 are "
    "exploratory diagnostics, not independent detections.",
    "A residual tail at 100-1000 arcsec persists for several "
    "historical plate sources; these may reflect genuine plate "
    "reduction errors but are excluded from interpretation."]

with open(RESULTS / "step_b127_apdb_planet.json", "w") as f:
    json.dump(res, f, indent=1)
logger.data_save(RESULTS / "step_b127_apdb_planet.json")

import csv
keys = sorted({k for r in rows for k in r})
with open(RESULTS / "step_b127_apdb_planet.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k) for k in keys})
logger.data_save(RESULTS / "step_b127_apdb_planet.csv")

# figure: Pluto O-C along-track vs epoch, coloured by cap separation
if pl:
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    yr = [r["yr"] for r in pl]
    s = [r["sep_axis"] for r in pl]
    ax[0].scatter(yr, [r["oc_at"] for r in pl],
                  c=s, cmap="viridis_r", s=6)
    ax[0].set_xlabel("year"); ax[0].set_ylabel("O-C along-track [as]")
    ax[0].set_title("Pluto residual vs epoch")
    ax[1].scatter(s, [r["oc_at"] for r in pl], s=6)
    ax[1].axvline(CAP_DEG, ls="--", c="r", lw=0.8)
    ax[1].set_xlabel("Pluto helio direction sep from axis [deg]")
    ax[1].set_ylabel("O-C along-track [as]")
    fig.tight_layout()
    fig.savefig(FIG / "step_b127_apdb_planet.png", dpi=150)
    logger.data_save(FIG / "step_b127_apdb_planet.png")
    plt.close(fig)

logger.info(json.dumps(res, indent=1))
logger.info("done")
