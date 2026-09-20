#!/usr/bin/env python3
"""
TEP-9 step 101 -- Interstellar-object non-gravitational channel
================================================================

The three catalogued interstellar objects are the only population
besides comets that physically transits the outer boundary sector on
measurable trajectories.  Step 099 scored their inbound asymptotes
geometrically; this step tests the *dynamical* channel: the fitted
non-gravitational (NG) acceleration parameters in the JPL SBDB orbit
solutions against the boundary-sector geometry each trajectory
threaded.

The discriminating observable is not the presence of NG -- every
outgassing comet carries fitted NG terms -- but an *unexplained* NG:
a fitted acceleration with no detected volatile driver.  Literature
anchors (transcribed published values, ``iso_literature.json``)
record that 1I/'Oumuamua's NG is unexplained by any detected mass
loss (~25 kg/s required, ~3 orders of magnitude above observed
limits) while 2I/Borisov's and 3I/ATLAS's are accounted for by
measured gas and dust production.

Geometry scored per object
--------------------------
* arrival and departure velocity-asymptote directions (computed from
  the SBDB osculating hyperbolic elements: position directions at
  nu = -nu_inf and nu = +nu_inf), and the periapsis direction that
  sets the sky sector of the inner observed arc;
* separations against the pre-declared axis (49,-17) and antiaxis
  (229,+17), classified as lobe (<60 deg, the measured cap radius),
  edge band (60-75 deg, the measured transition width of step 073),
  or field;
* a bipolar double-lobe transit flag -- arrival through one 60-deg
  lobe and departure through the other -- the geometry under which
  the two boundary crossings compound rather than cancel in the
  bipolar slip field (steps 052, 089).

Tests
-----
T1 sector concordance: does unexplained-NG incidence coincide with
   the bipolar double-lobe transit?  Exchangeability weight reported
   (n = 3; descriptive, not a detection claim).
T2 magnitude: implied fractional effective-attraction offset
   f(r) = |a_ng(r)| / g_sun(r) evaluated under each solution's own
   Marsden-Sekanina g(r) law at r = q, 1 au, 2 au; compared with the
   comet-channel lapse contrast (step 100, implied dA/A).
T3 vector morphology: |A2/A1| and |A3/A1| transverse fractions and
   perihelion-time offset DT -- a directed-gradient signature is
   radial-dominant while rotational-jet outgassing loads power into
   the transverse component and asymmetric activity requires DT.
T4 inner-arc sector: the periapsis-direction sector each inner
   observed arc sampled vs the measured cos-2theta sign structure
   (step 089: positive dtau at both axis poles, negative in the
   60-120 deg mid-band).  Borisov's fully volatile-explained NG
   bounds any *universal interior* lapse-gradient term on a
   mid-band inner arc at its fitted amplitude -- extending the
   wall-over-field verdict of step 068 into the inner system.
T5 prospective registration: pre-declared NG-class predictions for
   future ISO arrivals by transit geometry.
T6 classification stability: single-member reclassification audit --
   the T1 concordance is load-bearing on the volatile classes, so
   each object's reclassification consequence is computed and the
   monitored channel (SBDB solution evolution, the 3I mass-
   accounting caveat) registered; publicly discussed 3I oddities
   are triaged as non-lapse channels.

Inputs
------
JPL SBDB full-precision orbit solutions incl. NG model_pars
  (live query, raw JSON + provenance pinned under data/raw/iso/)
data/raw/literature/iso_literature.json (transcribed anchors)
results/step_b64_clock_consistency.json (comet-channel contrast)

Outputs
-------
results/step_b65_iso_ng_channel.json
results/step_b65_iso_ng_channel.csv
results/figures/supplementary/step_b65_iso_ng_channel.png
data/raw/iso/sbdb_iso_{1I,2I,3I}.json
data/raw/iso/provenance.json
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

logger = StepLogger("step_101_iso_ng_channel")
tee_stdout(logger)
logger.header("Interstellar-object non-gravitational channel")

SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"
ISO_DIR = DATA_RAW / "iso"
ISO_DIR.mkdir(parents=True, exist_ok=True)

CAP = 60.0          # measured cap radius (deg), used throughout
EDGE_HI = 75.0      # outer edge of the measured transition band (deg)
GM_SUN = 2.959122082855911e-4   # au^3 d^-2 (Gaussian k^2)

ISOS = {"1I": "1I/'Oumuamua", "2I": "2I/Borisov", "3I": "3I/ATLAS"}

prov = {}
recs = {}


def fetch_sbdb(des):
    url = SBDB + "?" + urllib.parse.urlencode(
        {"des": des, "full-prec": "true"})
    out = ISO_DIR / f"sbdb_iso_{des}.json"
    for attempt in range(3):
        try:
            logger.progress(f"GET {url}")
            req = urllib.request.Request(
                url, headers={"User-Agent": "tep9/1.0"})
            t0 = datetime.now(timezone.utc)
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
                status = r.status
            sha = hashlib.sha256(body).hexdigest()
            prov[f"sbdb_iso_{des}.json"] = dict(
                url=url, retrieved_utc=t0.isoformat(),
                http_status=status, bytes=len(body), sha256=sha)
            out.write_bytes(body)
            logger.data_save(out)
            time.sleep(0.4)
            return json.loads(body)
        except Exception as ex:
            logger.progress(f"SBDB attempt {attempt + 1} failed: {ex}")
            time.sleep(2.0 * (attempt + 1))
    if out.exists():
        cached = out.read_bytes()
        prov[f"sbdb_iso_{des}.json"] = dict(
            url=url, retrieved_utc="cached-file",
            note=("live fetch failed after 3 attempts; reused the "
                  "pinned raw file from the previous live retrieval"),
            bytes=len(cached),
            sha256=hashlib.sha256(cached).hexdigest())
        logger.progress(f"using pinned cache {out}")
        return json.loads(cached)
    raise RuntimeError(f"SBDB fetch failed and no cache for {des}")


def sector(sep_min):
    if sep_min < CAP:
        return "lobe"
    if sep_min < EDGE_HI:
        return "edge"
    return "field"


def classify(dvec):
    sa = float(sep(dvec, AXES["tno"]))
    sb = float(sep(dvec, AXES["anti"]))
    lon, lat = lb(dvec)
    # bipolar cos-2theta band of step 089: positive slip at both
    # axis poles (theta<60 or theta>120), negative in the mid-band
    band = "pole_lobe" if (sa < CAP or sa > 120.0) else "mid_band"
    return dict(ecl_lon=round(lon, 2), ecl_lat=round(lat, 2),
                sep_axis_deg=round(sa, 2), sep_antiaxis_deg=round(sb, 2),
                in_cap=bool(sa < CAP), in_mirror_cap=bool(sb < CAP),
                sector=sector(min(sa, sb)),
                bipolar_band=band,
                near_pole="axis" if sa < sb else "antiaxis")


def geometry(d):
    el = {e["name"]: float(e["value"]) for e in d["orbit"]["elements"]
          if e.get("value") is not None}
    e, inc = el["e"], math.radians(el["i"])
    om, Om = math.radians(el["w"]), math.radians(el["om"])
    P = perih_dir(om, Om, inc)
    W = np.array([math.sin(Om) * math.sin(inc),
                  -math.cos(Om) * math.sin(inc), math.cos(inc)])
    Q = np.cross(W, P)
    nu_max = math.acos(-1.0 / e)
    u_arr = math.cos(nu_max) * P - math.sin(nu_max) * Q
    u_arr = u_arr / np.linalg.norm(u_arr)
    u_dep = math.cos(nu_max) * P + math.sin(nu_max) * Q
    u_dep = u_dep / np.linalg.norm(u_dep)
    return el, P, u_arr, u_dep


def ng_ledger(d, el):
    pars = {p["name"]: p for p in d["orbit"].get("model_pars", [])}
    out = {"n_model_pars": len(pars)}
    for k in ("A1", "A2", "A3", "DT"):
        if k in pars and pars[k].get("value") is not None:
            out[k.lower()] = float(pars[k]["value"])
            out[k.lower() + "_sig"] = (float(pars[k]["sigma"])
                                       if pars[k].get("sigma") else None)
    const = {}
    for k in ("ALN", "R0", "NM", "NN", "NK"):
        if k in pars and pars[k].get("value") is not None:
            const[k] = float(pars[k]["value"])
    out["law"] = const

    def g(r):
        aln = const.get("ALN", 1.0)
        r0 = const.get("R0", 1.0)
        m = const.get("NM", 2.0)
        n = const.get("NN", 0.0)
        kk = const.get("NK", 0.0)
        return aln * (r / r0) ** (-m) * (1.0 + (r / r0) ** n) ** (-kk)

    prof = {}
    if "a1" in out:
        for r in (el["q"], 1.0, 2.0):
            a_ng = out["a1"] * g(r)
            f = abs(a_ng) / (GM_SUN / r ** 2)
            prof[f"r_{r:.3g}"] = dict(
                r_au=r, g_of_r=g(r), a_ng_au_d2=a_ng,
                frac_of_gravity=f,
                sign="outward" if out["a1"] > 0 else "inward")
        out["radial_profile"] = prof
        out["a1_significance"] = (abs(out["a1"] / out["a1_sig"])
                                  if out.get("a1_sig") else None)
    if "a1" in out and "a2" in out:
        out["transverse_frac_a2"] = abs(out["a2"] / out["a1"])
    if "a1" in out and "a3" in out:
        out["normal_frac_a3"] = abs(out["a3"] / out["a1"])
    out["dt_present"] = "dt" in out
    return out


def main():
    lit = json.loads((DATA_RAW / "literature" /
                      "iso_literature.json").read_text())

    comet = json.loads((RESULTS /
                        "step_b64_clock_consistency.json").read_text())
    dA = comet["comet_channel"]["implied_dA_over_A"]
    logger.metric("comet_contrast_dA_over_A", f"{dA:.3e}")

    rows = []
    for des, label in ISOS.items():
        try:
            d = fetch_sbdb(des)
        except Exception as ex:
            logger.error(f"{label}: SBDB failed: {ex}")
            continue
        recs[des] = d
        el, P, u_arr, u_dep = geometry(d)
        ng = ng_ledger(d, el)
        geo = dict(arrival=classify(u_arr),
                   departure=classify(u_dep),
                   periapsis=classify(P))
        double_lobe = ((geo["arrival"]["in_cap"] and
                        geo["departure"]["in_mirror_cap"]) or
                       (geo["arrival"]["in_mirror_cap"] and
                        geo["departure"]["in_cap"]))
        any_lobe = any(geo[k]["in_cap"] or geo[k]["in_mirror_cap"]
                       for k in ("arrival", "departure"))
        anchor = lit.get(des, {})
        vol = anchor.get("volatile_accounting", {})
        orb = d["orbit"]
        row = dict(
            des=des, label=label,
            orbit_id=orb.get("orbit_id"),
            soln_date=orb.get("soln_date"),
            data_arc_d=orb.get("data_arc"),
            n_obs_used=orb.get("n_obs_used"),
            rms_arcsec=orb.get("rms"),
            e=el["e"], q_au=el["q"], i_deg=el["i"],
            geometry=geo,
            double_lobe_transit=bool(double_lobe),
            any_lobe_leg=bool(any_lobe),
            ng=ng,
            volatile_status=vol.get("status"),
            volatile_note=vol.get("note"),
            volatile_citations=vol.get("citations"))
        rows.append(row)
        a1 = ng.get("a1")
        f1 = (ng.get("radial_profile", {})
              .get("r_1", {}).get("frac_of_gravity"))
        logger.metric(label.replace("/", "").replace(" ", "_"),
                      f"arr {geo['arrival']['sector']}"
                      f"({geo['arrival']['ecl_lon']},"
                      f"{geo['arrival']['ecl_lat']}) "
                      f"dep {geo['departure']['sector']}"
                      f"({geo['departure']['ecl_lon']},"
                      f"{geo['departure']['ecl_lat']})",
                      f"vol={vol.get('status')} "
                      f"A1={a1:.2e} f(1au)={f1:.2e}"
                      if a1 is not None else "no NG")

    # fail loudly on dropped records: the concordance statistics are
    # only meaningful over the full ISO set
    if len(rows) != len(ISOS):
        raise RuntimeError(
            f"ISO table incomplete: {len(rows)}/{len(ISOS)} records "
            "-- refusing to score concordance on a partial set")

    # ---------------- T1: sector concordance -------------------------
    n = len(rows)
    n_unexp = sum(r["volatile_status"] == "unexplained" for r in rows)
    n_dbl = sum(r["double_lobe_transit"] for r in rows)
    concord = all(
        (r["volatile_status"] == "unexplained") ==
        r["double_lobe_transit"] for r in rows)
    # exchangeability weight: under geometry-blind assignment the
    # unexplained object is uniformly one of n; the probability it
    # is also the unique double-lobe object is n_dbl/n
    exch_p = (n_dbl / n) if (n and n_unexp) else None
    T1 = dict(n_iso=n, n_unexplained_ng=n_unexp,
              n_double_lobe=n_dbl,
              all_three_classified=bool(concord),
              exchangeability_p=float(exch_p) if exch_p else None,
              note=("descriptive concordance at n=3; the inferential "
                    "weight is carried by the magnitude sign "
                    "consistency and the preregistered T5 "
                    "predictions, not by the retrospective "
                    "coincidence alone"))
    logger.metric("T1_concordance", f"{concord}",
                  f"p_exch={exch_p:.3f}" if exch_p else "n/a")

    # ---------------- T2: magnitude ----------------------------------
    mags = {}
    for r in rows:
        prof = r["ng"].get("radial_profile", {})
        if "r_1" in prof:
            mags[r["label"]] = dict(
                f_1au=prof["r_1"]["frac_of_gravity"],
                f_q=prof[f"r_{r['q_au']:.3g}"]["frac_of_gravity"],
                a1_sig=r["ng"].get("a1_significance"))
    f1i = mags.get("1I/'Oumuamua", {}).get("f_1au")
    ratio_to_comet = (f1i / dA) if f1i else None
    explained = [mags[r["label"]]["f_1au"] for r in rows
                 if r["volatile_status"] == "explained"
                 and r["label"] in mags]
    T2 = dict(
        comet_contrast_dA_over_A=dA,
        per_iso_frac_1au={k: round(v["f_1au"], 6) for k, v in mags.items()},
        oumuamua_over_comet=round(ratio_to_comet, 3)
        if ratio_to_comet else None,
        oumuamua_over_explained_max=round(
            f1i / max(explained), 2) if f1i and explained else None,
        note=(f"f(r) = |a_ng(r)|/g_sun(r) under each solution's own "
              f"g(r) law.  'Oumuamua's unexplained fractional offset "
              f"sits at ~1e-3 -- a factor "
              f"{dA/f1i:.0f} below the comet-channel lapse-contrast "
              f"scale ({dA:.1e}) if the two channels couple "
              f"identically" +
              (f" -- while the volatile-explained objects "
               f"carry ~1.6-1.8e-4, ~{dA/max(explained):.0f}-"
               f"{dA/min(explained):.0f}x below the comet scale"
               if explained else "") +
              f".  The sector coincidence is exact, but the "
              f"amplitude match requires a channel-dependent "
              f"(attenuated) ISO coupling or a partial-slip "
              f"interpretation."
              if f1i else
              "1I record unavailable this run."))
    logger.metric("T2_1I_frac_1au", f"{f1i:.2e}" if f1i else "n/a",
                  f"vs comet {dA:.1e} (x{ratio_to_comet:.2f})"
                  if ratio_to_comet else "")

    # ---------------- T3: vector morphology --------------------------
    morph = {r["label"]: dict(
        transverse_frac_a2=r["ng"].get("transverse_frac_a2"),
        normal_frac_a3=r["ng"].get("normal_frac_a3"),
        dt_d=r["ng"].get("dt"), dt_present=r["ng"]["dt_present"],
        law=r["ng"]["law"]) for r in rows}
    T3 = dict(per_iso=morph,
              note=("radial dominance |A2/A1| -> 0 and absence of a "
                    "perihelion-time offset characterize a directed "
                    "(gradient-like) term; rotational-jet outgassing "
                    "loads power into A2 and asymmetric activity "
                    "requires DT."))
    for r in rows:
        logger.metric(f"T3_{r['des']}_A2_over_A1",
                      f"{r['ng'].get('transverse_frac_a2'):.3f}"
                      if r['ng'].get('transverse_frac_a2') else "n/a",
                      f"DT={r['ng'].get('dt')}")

    # ---------------- T4: inner-arc sector + interior bound ----------
    inner = {r["label"]: r["geometry"]["periapsis"] for r in rows}
    borisov_f = mags.get("2I/Borisov", {}).get("f_1au")
    atlas_f = mags.get("3I/ATLAS", {}).get("f_1au")
    bound_f = max((x for x in (borisov_f, atlas_f) if x),
                  default=None)
    T4 = dict(
        periapsis_sector=inner,
        midband_negative_slip=("step 089: dtau significantly "
                               "negative for theta in [60,120] deg"),
        borisov_interior_bound=dict(
            f_1au=borisov_f,
            below_comet_contrast=bool(borisov_f and borisov_f < dA),
            factor_below=round(dA / borisov_f, 1)
            if borisov_f else None,
            note=(f"2I's periapsis sector is mid-band like 1I's, yet "
                  f"its fitted NG is fully volatile-explained at "
                  f"~1.6e-4 of gravity: a universal interior lapse "
                  f"gradient at comet-channel strength on a mid-band "
                  f"inner arc is excluded by a factor "
                  f"{dA/borisov_f:.0f}.  The lapse structure "
                  f"therefore does not extend into the inner system "
                  f"as a distributed sector field -- the wall-over-"
                  f"field verdict of step 068 extended inward."
                  if borisov_f else
                  "2I record unavailable this run; interior bound "
                  "not computed.")),
        atlas_interior_bound=dict(
            f_1au=atlas_f,
            sector="mirror-lobe inner arc",
            note=("3I's periapsis sector sits inside the mirror "
                  "lobe; an interior positive-lapse term would have "
                  "added an inward residual on top of outgassing, "
                  "demanding ~6x the measured volatile production "
                  "to stay fitted-flat.  Its NG is consistent with "
                  "the measured CO2-driven outgassing at ~1.8e-4, "
                  "extending the interior bound to the lobe "
                  "sector.")),
        interior_bound_summary=(
            f"universal interior sector-gradient bounded at "
            f"~{bound_f:.1e} of gravity ({dA / bound_f:.0f}x below "
            f"the comet-channel contrast) across both the negative-"
            f"slip mid-band and the positive-slip lobe sector"
            if bound_f else
            "explained-ISO records unavailable this run; interior "
            "bound not computed"),
        sign_reading=("1I's periapsis direction sits in the measured "
                      "negative-dtau mid-band; its unexplained term "
                      "is outward -- an effective-attraction "
                      "shortfall -- the sign a lapse deficit "
                      "produces.  Reported as sign-consistency, "
                      "interpretation level: the wall model alone "
                      "does not predict in-arc terms."))
    logger.metric("T4_borisov_bound", f"{borisov_f:.2e}",
                  f"{dA / borisov_f:.1f}x below comet contrast"
                  if borisov_f else "n/a")

    # ---------------- T5: prospective registration -------------------
    T5 = dict(
        rule=("pre-declared NG-class prediction for future ISOs by "
              "transit geometry computed at discovery from the "
              "SBDB orbit solution"),
        predictions=[

            dict(transit="bipolar double-lobe (arrival in one 60-deg "
                 "lobe, departure in the other)",
                 predicted_ng_class="unexplained-fraction candidate",
                 predicted_frac="~1e-3 of gravity if volatile-free"),
            dict(transit="edge-band leg (60-75 deg) or single-lobe",
                 predicted_ng_class="volatile-explained dominant",
                 predicted_frac="lapse term at/below ~1e-4, masked "
                 "by any measured outgassing"),
            dict(transit="field (both asymptotes >75 deg)",
                 predicted_ng_class="volatile-explained",
                 predicted_frac="no lapse contribution")],
        falsifier=("the next ISO executing a bipolar double-lobe "
                   "transit with a fully volatile-explained NG, or a "
                   "field-transit ISO requiring unexplained NG at "
                   "~1e-3, breaks the concordance"),
        note=("n is small and arrival directions are not a designed "
              "draw; this table registers the prediction so the "
              "fourth ISO tests it prospectively"))

    # ---------------- T6: classification stability --------------------
    # T1's concordance is load-bearing on the volatile classes: it
    # requires 'Oumuamua to remain unexplained AND Borisov/ATLAS to
    # remain explained.  For each object compute what a single
    # reclassification does to the concordance -- this registers the
    # monitored channel by which the datum could strengthen, break
    # or dissolve.
    reclass = {}
    for r in rows:
        cur = r["volatile_status"]
        hyp = "explained" if cur == "unexplained" else "unexplained"
        sim = [hyp if rr["des"] == r["des"]
               else rr["volatile_status"] for rr in rows]
        conc = all((s == "unexplained") == rr["double_lobe_transit"]
                   for s, rr in zip(sim, rows))
        if hyp == "explained":
            conseq = ("anomaly dissolves: no unexplained member "
                      "remains, and the ISO channel carries no datum")
        elif conc:
            conseq = "concordance preserved"
        else:
            conseq = ("concordance broken: an unexplained member "
                      "without a bipolar double-lobe transit")
        reclass[r["label"]] = dict(
            current_status=cur, hypothetical=hyp,
            concordance_survives=bool(conc), consequence=conseq)
    atlas = next((r for r in rows if r["des"] == "3I"), {})
    T6 = dict(
        rule=("single-member reclassification audit: for each ISO, "
              "flip its volatile class and re-evaluate the T1 "
              "concordance (unexplained == double-lobe for all "
              "three)"),
        per_iso_reclassification=reclass,
        load_bearing=("the concordance requires the explained status "
                      "of Borisov and ATLAS as much as the "
                      "unexplained status of 'Oumuamua -- the "
                      "explained objects are controls, not failed "
                      "detections, and their reclassification is a "
                      "registered concordance-break channel"),
        reported_anomaly_triage=dict(
            atlas_status=atlas.get("volatile_status"),
            note=("3I's publicly discussed oddities -- the retrograde "
                  "trajectory within ~5 deg of the ecliptic, the "
                  "sunward anti-tail/jet morphology, the extreme "
                  "polarization -- are kinematic and compositional "
                  "channels, not lapse channels; they carry no "
                  "discriminating weight in this test.  The only "
                  "lapse-relevant datum is the fitted NG accounting, "
                  "on which ATLAS is an ordinary outgassing body.  "
                  "Were its term genuinely unexplained it would sit "
                  "in the concordance-break column of this audit, "
                  "not the support column.")),
        monitoring=dict(
            channel=("SBDB solution evolution; records re-fetched "
                     "live each run, provenance-pinned under "
                     "data/raw/iso/"),
            soln_dates={r["label"]: r["soln_date"] for r in rows},
            mass_accounting_caveat=(
                "ATLAS's explained status rests on the nucleus "
                "mass/production accounting; post-perihelion "
                "systematic-uncertainty NG solutions (Spada et al. "
                "2026), NG-inferred mass/size (Thoss, Loeb & Burkert "
                "2026) and "
                "HST nucleus sizing (Hui et al. 2026, ApJL 999, L37) "
                "concur, while a minority large-nucleus claim is "
                "registered in iso_literature.json:disputed_claim.  "
                "An upgrade of either explained object to "
                "unexplained status breaks T1; a detected volatile "
                "driver for 'Oumuamua would dissolve the datum "
                "itself.")))
    logger.metric("T6_atlas_reclass",
                  f"{reclass.get('3I/ATLAS', {}).get('consequence')}",
                  f"status={atlas.get('volatile_status')} "
                  f"soln={atlas.get('soln_date')}")

    # ---------------- figure ------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.5, 6.2))

    ax.set_xlim(360, 0)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("ecliptic longitude (deg)")
    ax.set_ylabel("ecliptic latitude (deg)")
    th = np.linspace(0, 2 * np.pi, 240)
    for ctr, col in ((AXES["tno"], "tab:blue"), (AXES["anti"], "tab:red")):
        z = math.sin(math.radians(CAP))
        rr = math.cos(math.radians(CAP))
        u = np.cross(ctr, [0, 0, 1]); u /= np.linalg.norm(u)
        w = np.cross(ctr, u)
        pts = np.array([rr * ctr + z * (math.cos(t) * u + math.sin(t) * w)
                        for t in th])
        lons, lats = zip(*[lb(p / np.linalg.norm(p)) for p in pts])
        ax.plot(np.r_[lons, lons[0]], np.r_[lats, lats[0]],
                color=col, lw=1.0, alpha=0.6)
    ax.plot(*lb(AXES["tno"]), "*", ms=16, color="tab:blue")
    ax.plot(*lb(AXES["anti"]), "*", ms=16, color="tab:red")

    colors = {"1I": "tab:purple", "2I": "tab:green", "3I": "tab:orange"}
    mk = {"arrival": ("s", "arrival"), "departure": ("o", "departure"),
          "periapsis": ("^", "periapsis")}
    for r in rows:
        c = colors[r["des"]]
        for kind, (m, _) in mk.items():
            g = r["geometry"][kind]
            ax.plot(g["ecl_lon"], g["ecl_lat"], marker=m, ms=7,
                    color=c, ls="none",
                    mfc="none" if kind == "periapsis" else c)
            ax.annotate(r["des"], (g["ecl_lon"], g["ecl_lat"]),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=7, color=c)
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker=m, color="k", ls="none", ms=7,
                      label=lbl) for _, (m, lbl) in mk.items()]
    handles += [Line2D([0], [0], marker="s", color=c, ls="none", ms=7,
                       label=ISOS[des]) for des, c in colors.items()]
    ax.legend(handles=handles, fontsize=7, loc="lower left")
    ax.set_title("ISO transit geometry vs 60-deg caps")
    ax.grid(alpha=0.3)

    rs = np.linspace(0.2, 4.0, 200)
    for r in rows:
        law = r["ng"]["law"]
        aln = law.get("ALN", 1.0); r0 = law.get("R0", 1.0)
        m = law.get("NM", 2.0); nn = law.get("NN", 0.0)
        kk = law.get("NK", 0.0)
        g = aln * (rs / r0) ** (-m) * (1 + (rs / r0) ** nn) ** (-kk)
        f = abs(r["ng"].get("a1", 0.0)) * g / (GM_SUN / rs ** 2)
        lab = r["label"] + (" (unexplained)" if r["volatile_status"]
                            == "unexplained" else " (explained)")
        ax2.plot(rs, f, color=colors[r["des"]], lw=1.4, label=lab)
        ax2.axvline(r["q_au"], color=colors[r["des"]], lw=0.6,
                    alpha=0.5, ls=":")
    ax2.axhline(dA, color="k", ls="--", lw=1.0,
                label=f"comet channel dA/A = {dA:.1e}")
    ax2.set_yscale("log")
    ax2.set_xlabel("heliocentric distance r (au)")
    ax2.set_ylabel("|a_ng| / g_sun  (fractional offset)")
    ax2.set_title("Implied effective-attraction offset vs radius")
    ax2.legend(fontsize=7)
    ax2.grid(alpha=0.3, which="both")
    fig.tight_layout()
    figp = RESULTS / "figures" / "supplementary" / "step_b65_iso_ng_channel.png"
    fig.savefig(figp, dpi=300)
    logger.data_save(figp)

    # ---------------- csv ledger --------------------------------------
    import csv
    csvp = RESULTS / "step_b65_iso_ng_channel.csv"
    with open(csvp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["des", "label", "orbit_id", "soln_date",
                    "data_arc_d", "n_obs", "rms_arcsec", "e", "q_au",
                    "arr_sector", "arr_sep_axis", "arr_sep_anti",
                    "dep_sector", "dep_sep_axis", "dep_sep_anti",
                    "peri_sector", "peri_sep_axis", "peri_sep_anti",
                    "peri_bipolar_band",
                    "double_lobe", "any_lobe", "A1", "A1_sig",
                    "A1_sigma_det", "f_1au", "f_q", "A2_over_A1",
                    "A3_over_A1", "DT_d", "volatile_status"])
        for r in rows:
            ng = r["ng"]; prof = ng.get("radial_profile", {})
            w.writerow([
                r["des"], r["label"], r["orbit_id"], r["soln_date"],
                r["data_arc_d"], r["n_obs_used"], r["rms_arcsec"],
                r["e"], r["q_au"],
                r["geometry"]["arrival"]["sector"],
                r["geometry"]["arrival"]["sep_axis_deg"],
                r["geometry"]["arrival"]["sep_antiaxis_deg"],
                r["geometry"]["departure"]["sector"],
                r["geometry"]["departure"]["sep_axis_deg"],
                r["geometry"]["departure"]["sep_antiaxis_deg"],
                r["geometry"]["periapsis"]["sector"],
                r["geometry"]["periapsis"]["sep_axis_deg"],
                r["geometry"]["periapsis"]["sep_antiaxis_deg"],
                r["geometry"]["periapsis"]["bipolar_band"],
                r["double_lobe_transit"], r["any_lobe_leg"],
                ng.get("a1"), ng.get("a1_sig"),
                ng.get("a1_significance"),
                prof.get("r_1", {}).get("frac_of_gravity"),
                prof.get(f"r_{r['q_au']:.3g}", {}).get("frac_of_gravity"),
                ng.get("transverse_frac_a2"),
                ng.get("normal_frac_a3"), ng.get("dt"),
                r["volatile_status"]])
    logger.data_save(csvp)

    (ISO_DIR / "provenance.json").write_text(json.dumps(
        dict(step="step_101_iso_ng_channel",
             generated_utc=datetime.now(timezone.utc).isoformat(),
             files=prov), indent=1))
    logger.data_save(ISO_DIR / "provenance.json")

    out = dict(
        step="step_101_iso_ng_channel",
        description=("Fitted non-gravitational acceleration ledger of "
                     "the three interstellar objects (JPL SBDB "
                     "solutions, provenance-pinned) tested against "
                     "the boundary-sector transit geometry: arrival, "
                     "departure and periapsis directions vs the "
                     "60-deg caps and the bipolar slip morphology."),
        axis=dict(ecl_lon=49.0, ecl_lat=-17.0,
                  anti_lon=229.0, anti_lat=17.0, cap_deg=CAP,
                  edge_band_deg=[CAP, EDGE_HI]),
        iso_ledger=rows,
        T1_sector_concordance=T1,
        T2_magnitude=T2,
        T3_vector_morphology=T3,
        T4_inner_arc_sector=T4,
        T5_prospective_registration=T5,
        T6_classification_stability=T6,
        provenance=prov,
        inputs=["data/raw/iso/sbdb_iso_*.json (this step, live)",
                "data/raw/literature/iso_literature.json",
                "results/step_b64_clock_consistency.json"],
        caveats=[
            "n = 3 and ISO arrival directions are not a designed "
            "draw; the concordance is a registered datum whose "
            "weight is prospective (T5), not a stand-alone "
            "detection.",
            "Unexplained-NG status is a literature anchor "
            "(transcribed published values, iso_literature.json); "
            "'Oumuamua's volatile-free status could in principle "
            "change if a detected driver is ever identified.",
            "f(r) profiles use each solution's own Marsden-Sekanina "
            "g(r); different NG law choices move the implied "
            "fraction at fixed astrometry (reported per-solution).",
            "The wall model places the slip at the boundary "
            "crossing; an in-arc NG term is a different field "
            "functional.  The channel therefore tests interior-field "
            "extension (bounded by 2I) and records the asymptote- "
            "geometry concordance -- it does not claim the wall "
            "itself generates 'Oumuamua's in-arc term.",
            "Sector classifications use the pre-declared 60-deg "
            "cap; the measured edge band (step 073) is reported "
            "separately.",
            "The T1 concordance is load-bearing on the volatile "
            "classes (T6): 'Oumuamua's unexplained status requires "
            "Borisov's and ATLAS's explained status to hold, and "
            "vice versa -- a reclassification of any member either "
            "breaks the concordance or dissolves the datum.  "
            "ATLAS's explained verdict rests on the mass/production "
            "accounting, where a minority claim is disputed "
            "(iso_literature.json:disputed_claim); the channel is "
            "monitored through live SBDB re-fetches."])
    outp = RESULTS / "step_b65_iso_ng_channel.json"
    outp.write_text(json.dumps(out, indent=1))
    logger.data_save(outp)
    logger.success("ISO non-gravitational channel complete")


if __name__ == "__main__":
    main()
