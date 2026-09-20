#!/usr/bin/env python3
"""
TEP-9 step 099 -- Heliospheric geometry and probe/ISO asymptotes
================================================================

Scores the TEP boundary axis (49 deg, -17 deg ecliptic) and its
60-deg caps against the measured geometry of the heliospheric
interface and of every object that has physically crossed -- or is
crossing -- the relevant boundary sector:

Reference directions (literature anchors, cited):
  * ISM inflow (upwind)          (255.8, +5.16)
  * heliotail apex (downwind)    ( 75.8, -5.16)
  * IBEX ribbon centre           (221.0, +39.0)
  * pristine ISMF direction      (227.28, +34.62)  Zirnstein+16

Escaping spacecraft (heliocentric velocity asymptote at epoch
2026-01-01, JPL Horizons vector ephemeris):
  Voyager 1 (-31), Voyager 2 (-32), Pioneer 10 (-23),
  Pioneer 11 (-24), New Horizons (-98)

Interstellar objects (inbound velocity asymptote from SBDB
hyperbolic elements):
  1I/'Oumuamua, 2I/Borisov, 3I/ATLAS

Tests
-----
1. Separations: every direction vs the axis and anti-axis; cap
   membership at the 60-deg radius used throughout the paper.
2. Axis-vs-heliosphere proximity: angular distance from the axis to
   the ISM-inflow antipode (heliotail apex) and to the ISMF
   direction, with a uniform-direction p-value for each.
3. Probe alibi: how many escaping-probe asymptotes fall inside the
   cap vs the mirror cap -- the geometric fact behind "why has no
   spacecraft seen the boundary directly".  Reported descriptively
   (probe directions are not a random draw).
4. ISO arrivals: inbound asymptotes vs caps; a dated datum for the
   ISO-arrival channel.

Inputs
------
data/raw/literature/literature_anchors.json
JPL Horizons vector ephemeris (live query, provenance recorded)
JPL SBDB (live query, provenance recorded)
results/step_b60_sclk_audit.json (Voyager geometry, cross-check)

Outputs
-------
results/step_b63_heliospheric_geometry.json
results/figures/supplementary/step_b63_heliospheric_geometry.png
data/raw/horizons/probe_vectors.json
"""

import sys
import json
import math
import time
import hashlib
import urllib.request
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (DATA_RAW, RESULTS, AXES, sep,
                                       lv, lb, perih_dir, tee_stdout)

logger = StepLogger("step_099_heliospheric_geometry")
tee_stdout(logger)
logger.header("Heliospheric geometry and probe/ISO asymptotes")

LIT = json.loads((DATA_RAW / "literature" /
                  "literature_anchors.json").read_text())
HDIR = LIT["heliosphere_directions"]

CAP = 60.0
EPOCH = "2026-01-01"

PROBES = {"Voyager 1": "-31", "Voyager 2": "-32", "Pioneer 10": "-23",
          "Pioneer 11": "-24", "New Horizons": "-98"}
ISOS = {"1I/'Oumuamua": "1I", "2I/Borisov": "2I", "3I/ATLAS": "3I"}

HOR = "https://ssd.jpl.nasa.gov/api/horizons.api"
SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"
HOR_DIR = DATA_RAW / "horizons"
HOR_DIR.mkdir(parents=True, exist_ok=True)
prov = {}


def fetch(url, label):
    logger.progress(f"GET {url[:110]}")
    req = urllib.request.Request(url, headers={"User-Agent": "tep9/1.0"})
    t0 = datetime.now(timezone.utc)
    body = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            break
        except Exception as exc:
            if attempt == 5:
                raise
            wait = 5 * 2 ** (attempt - 1)
            logger.progress(f"fetch attempt {attempt} failed ({exc}); "
                            f"retrying in {wait}s")
            time.sleep(wait)
    sha = hashlib.sha256(body).hexdigest()
    prov[label] = dict(url=url, retrieved_utc=t0.isoformat(),
                       bytes=len(body), sha256=sha)
    time.sleep(0.4)
    return body.decode()


# ---------------------------------------------------------------------
# 1. probe asymptotes from Horizons vector ephemeris
# ---------------------------------------------------------------------
def probe_asymptote(name, cmd):
    q = dict(format="text", COMMAND=f"'{cmd}'", OBJ_DATA="'YES'",
             MAKE_EPHEM="'YES'", EPHEM_TYPE="'VECTORS'",
             CENTER="'@10'", START_TIME=f"'{EPOCH}'",
             STOP_TIME="'2026-01-02'", STEP_SIZE="'1d'",
             OUT_UNITS="'AU-D'", VEC_TABLE="'2'")
    txt = fetch(HOR + "?" + urllib.parse.urlencode(q), f"horizons_{cmd}")
    m = txt.split("$$SOE")[1].split("$$EOE")[0]
    import re as _re
    vx, vy, vz = [float(x) for x in
                  _re.search(r"VX=\s*([-\d.E+]+)\s*VY=\s*([-\d.E+]+)"
                             r"\s*VZ=\s*([-\d.E+]+)", m).groups()]
    x, y, z = [float(x) for x in
               _re.search(r"X =\s*([-\d.E+]+)\s*Y =\s*([-\d.E+]+)"
                          r"\s*Z =\s*([-\d.E+]+)", m).groups()]
    v = np.array([vx, vy, vz])
    vhat = v / np.linalg.norm(v)
    pos = np.array([x, y, z])
    r = np.linalg.norm(pos)
    return dict(velocity_dir=vhat, position_dir=pos / r, r_au=float(r))


# ---------------------------------------------------------------------
# 2. ISO inbound asymptote from SBDB hyperbolic elements
# ---------------------------------------------------------------------
def iso_asymptote(name, des):
    txt = fetch(SBDB + "?" + urllib.parse.urlencode({"des": des}),
                f"sbdb_{des}")
    d = json.loads(txt)
    el = {e["name"]: float(e["value"]) for e in d["orbit"]["elements"]
          if e.get("value") is not None}
    e, inc = el["e"], math.radians(el["i"])
    om, Om = math.radians(el["w"]), math.radians(el["om"])
    P = perih_dir(om, Om, inc)
    # orbit normal W and in-plane Q = W x P
    W = np.array([math.sin(Om) * math.sin(inc),
                  -math.cos(Om) * math.sin(inc), math.cos(inc)])
    Q = np.cross(W, P)
    nu_max = math.acos(-1.0 / e)
    # inbound velocity at infinity (nu -> -nu_max):
    #   v_inf- ∝ sin(nu_max) P + (e + cos(nu_max)) Q  -- but velocity
    #   at the inbound asymptote points *toward* the Sun, i.e. the
    #   direction of motion along the incoming asymptote.  The
    #   direction on the sky the object comes from is the position
    #   vector at nu = -nu_max:
    u_arr = math.cos(nu_max) * P - math.sin(nu_max) * Q
    u_arr = u_arr / np.linalg.norm(u_arr)
    return dict(arrival_dir=u_arr, e=e, q_au=el["q"], i_deg=el["i"])


def score(dvec, label):
    s_axis = float(sep(dvec, AXES["tno"]))
    s_anti = float(sep(dvec, AXES["anti"]))
    lon, lat = lb(dvec)
    return dict(label=label, ecl_lon=round(lon, 2), ecl_lat=round(lat, 2),
                sep_axis_deg=round(s_axis, 2),
                sep_antiaxis_deg=round(s_anti, 2),
                in_cap=bool(s_axis < CAP), in_mirror_cap=bool(s_anti < CAP))


def main():
    table = []

    # --- reference directions ------------------------------------------
    refs = {
        "ISM inflow (upwind)": lv(HDIR["ism_inflow"]["lon_deg"],
                                HDIR["ism_inflow"]["lat_deg"]),
        "heliotail apex (downwind)": lv(HDIR["heliotail_apex"]["lon_deg"],
                                       HDIR["heliotail_apex"]["lat_deg"]),
        "IBEX ribbon centre": lv(HDIR["ibex_ribbon_centre"]["lon_deg"],
                               HDIR["ibex_ribbon_centre"]["lat_deg"]),
        "pristine ISMF": lv(HDIR["pristine_ismf"]["lon_deg"],
                          HDIR["pristine_ismf"]["lat_deg"]),
    }
    for k, v in refs.items():
        table.append(dict(kind="reference", **score(v, k)))

    # --- axis geometry vs heliosphere -----------------------------------
    axis_geom = {}
    for k, v in refs.items():
        s = float(sep(AXES["tno"], v))
        axis_geom[k] = dict(
            sep_deg=round(s, 2),
            uniform_p=round(float((1 - math.cos(math.radians(s))) / 2), 4))
        logger.metric(f"axis_vs_{k.split()[0]}", round(s, 1), "deg")

    # --- probes ----------------------------------------------------------
    probe_rows = []
    for name, cmd in PROBES.items():
        try:
            r = probe_asymptote(name, cmd)
        except Exception as ex:
            logger.error(f"{name}: Horizons failed: {ex}")
            continue
        row_v = dict(kind="probe_velocity", r_au=r["r_au"],
                     **score(r["velocity_dir"], name))
        row_p = dict(kind="probe_position", r_au=r["r_au"],
                     **score(r["position_dir"], name))
        probe_rows.append(row_v)
        table.extend([row_v, row_p])
        logger.metric(name.replace(" ", "_") + "_vel_dir",
                      f"({row_v['ecl_lon']}, {row_v['ecl_lat']})",
                      f"r={r['r_au']:.1f} AU")

    # --- ISOs --------------------------------------------------------------
    iso_rows = []
    for name, des in ISOS.items():
        try:
            r = iso_asymptote(name, des)
        except Exception as ex:
            logger.error(f"{name}: SBDB failed: {ex}")
            continue
        row = dict(kind="iso_arrival", e=r["e"], q_au=r["q_au"],
                   i_deg=r["i_deg"], **score(r["arrival_dir"], name))
        iso_rows.append(row)
        table.append(row)
        logger.metric(name.replace(" ", "_").replace("/", ""),
                      f"({row['ecl_lon']}, {row['ecl_lat']})",
                      "inbound asymptote")

    # --- probe alibi -------------------------------------------------------
    # fail loudly on dropped records: the cap fractions are only
    # meaningful over the full probe/ISO sets
    if len(probe_rows) != len(PROBES):
        raise RuntimeError(
            f"probe table incomplete: {len(probe_rows)}/{len(PROBES)} "
            "Horizons queries succeeded -- refusing to score cap "
            "fractions on a partial set")
    if len(iso_rows) != len(ISOS):
        raise RuntimeError(
            f"ISO table incomplete: {len(iso_rows)}/{len(ISOS)} "
            "SBDB queries succeeded -- refusing to score cap "
            "fractions on a partial set")
    n_in = sum(x["in_cap"] for x in probe_rows)
    n_mir = sum(x["in_mirror_cap"] for x in probe_rows)
    n_iso_in = sum(x["in_cap"] for x in iso_rows)
    n_iso_mir = sum(x["in_mirror_cap"] for x in iso_rows)
    alibi = dict(
        n_probes=len(probe_rows), n_velocity_in_cap=n_in,
        n_velocity_in_mirror=n_mir,
        uniform_cap_fraction=round(
            float((1 - math.cos(math.radians(CAP))) / 2), 4),
        n_iso=len(iso_rows), n_iso_in_cap=n_iso_in,
        n_iso_in_mirror=n_iso_mir,
        note=("Probe directions are mission-constrained, not a "
              "random draw; reported as a geometric datum, not an "
              "inferential statistic."))
    logger.metric("probes_in_cap", n_in, "of velocity asymptotes")
    logger.metric("probes_in_mirror", n_mir, "of velocity asymptotes")
    logger.metric("iso_in_cap", n_iso_in, "of arrival asymptotes")

    # --- figure: sky map ----------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.set_xlim(360, 0)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("ecliptic longitude (deg)")
    ax.set_ylabel("ecliptic latitude (deg)")

    # caps
    th = np.linspace(0, 2 * np.pi, 200)
    for ctr, col in ((AXES["tno"], "tab:blue"), (AXES["anti"], "tab:red")):
        # parametric cap circle around centre direction
        elon, elat = lb(ctr)
        z = math.sin(math.radians(CAP))
        rr = math.cos(math.radians(CAP))
        # build orthonormal basis around ctr
        u = np.cross(ctr, [0, 0, 1]); u /= np.linalg.norm(u)
        w = np.cross(ctr, u)
        pts = np.array([rr * ctr + z * (math.cos(t) * u + math.sin(t) * w)
                        for t in th])
        lons, lats = zip(*[lb(p / np.linalg.norm(p)) for p in pts])
        lons = np.array(lons)
        # unwrap for plotting
        lats = np.array(lats)
        order = np.argsort(np.cumsum(np.abs(np.diff(
            np.r_[lons, lons[0]]))) * -1)
        ax.plot(np.r_[lons, lons[0]], np.r_[lats, lats[0]],
                color=col, lw=1.0, alpha=0.6)
    ax.plot(*lb(AXES["tno"]), "*", ms=18, color="tab:blue",
            label="TNO axis")
    ax.plot(*lb(AXES["anti"]), "*", ms=18, color="tab:red",
            label="anti-axis")

    mk = {"reference": ("D", "black"), "probe_velocity": ("^", "tab:green"),
          "probe_position": ("v", "tab:olive"), "iso_arrival": ("s", "tab:purple")}
    for row in table:
        m, c = mk[row["kind"]]
        ax.plot(row["ecl_lon"], row["ecl_lat"], marker=m, ms=7,
                color=c, ls="none")
        ax.annotate(row["label"].replace(" (downwind)", "\n(downwind)"),
                    (row["ecl_lon"], row["ecl_lat"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=7)
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker=m, color=c, ls="none", ms=7,
                      label=k.replace("_", " "))
               for k, (m, c) in mk.items()]
    handles += [Line2D([0], [0], marker="*", color="tab:blue", ls="none",
                       ms=14, label="TNO axis (49,-17)"),
                Line2D([0], [0], marker="*", color="tab:red", ls="none",
                       ms=14, label="anti-axis (229,17)")]
    ax.legend(handles=handles, fontsize=7, loc="lower left")
    ax.set_title("TEP axis, 60-deg caps, heliospheric reference "
                 "directions, escaping-probe and ISO asymptotes")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    figp = RESULTS / "figures" / "supplementary" / "step_b63_heliospheric_geometry.png"
    fig.savefig(figp, dpi=300)
    logger.data_save(figp)

    (HOR_DIR / "probe_vectors.json").write_text(json.dumps(
        {"provenance": prov, "table": table}, indent=1))
    logger.data_save(HOR_DIR / "probe_vectors.json")

    out = dict(
        step="step_099_heliospheric_geometry",
        description=("Score the TEP axis and 60-deg caps against "
                     "heliospheric reference directions, escaping-"
                     "probe asymptotes (Horizons vectors), and ISO "
                     "inbound asymptotes (SBDB hyperbolic elements)."),
        epoch=EPOCH,
        axis=dict(ecl_lon=49.0, ecl_lat=-17.0, cap_deg=CAP),
        axis_geometry=axis_geom,
        directions_table=table,
        probe_alibi=alibi,
        provenance=prov,
        caveats=[
            "Probe directions are mission-constrained (Jupiter-"
            "assist and escape geometry), not a random draw; the "
            "probe-alibi counts are a geometric datum, not an "
            "inferential test.",
            "ISO inbound asymptotes are computed from SBDB "
            "osculating hyperbolic elements at their reference "
            "epochs; planetary perturbations shift the asymptote "
            "direction by <~1 deg for these solutions.",
            "Velocity asymptote != position direction: a probe can "
            "sit inside the mirror cap by position while its "
            "velocity asymptote lies outside it (V2 case).",
            "Reference directions are literature values with "
            "published uncertainties; the ISMF direction has "
            "+/-0.7 deg quoted errors.",
        ])
    out_path = RESULTS / "step_b63_heliospheric_geometry.json"
    out_path.write_text(json.dumps(out, indent=1))
    logger.data_save(out_path)
    logger.success("Heliospheric geometry complete")


if __name__ == "__main__":
    main()
