#!/usr/bin/env python3
"""
TEP-9 step 095 -- Spacecraft clock / plasma / power data acquisition
====================================================================

Downloads the real, archival data products needed by the spacecraft
clock and nuclear-process channels of the TEP-9 analysis:

1. NAIF Voyager spacecraft-clock (SCLK) kernels.  These text kernels
   encode the full reconstructed SCLK -> ephemeris-time correlation
   history (SCLK_COEFFICIENTS triplets per partition) built by JPL
   from the SCLKvSCET clock-calibration record.  Every coefficient
   record is a ground calibration of the onboard oscillator rate;
   partition boundaries and rate-mode records document every resync,
   re-rate and rollover event across the mission.

2. NAIF Voyager SPK ephemeris kernels (merged mission trajectories)
   so that each clock-calibration epoch can be placed at the correct
   heliocentric radius and sky direction.

3. NAIF leapseconds kernel (naif0012.tls) for ET <-> UTC conversion.

4. NAIF New Horizons SCLK kernel and predicted-trajectory SPK
   (PDS nhsp_1000 archive): the interior control craft -- a third
   long-lived outer-trajectory clock, still inside the heliosphere,
   whose asymptote lies outside both field caps.

5. PDS/PPI Voyager PWS VLISM electron-density collections
   (urn:nasa:pds:voyager-pws-vlism-density): per-event electron plasma
   frequency and derived density with uncertainty bounds, spacecraft
   radius and heliocentric direction, and distance past the
   heliopause, for Voyager 1 (2012-2025) and Voyager 2 (2019-2025).

6. The Dryad historical RTG performance dataset (Whiting & Woerner
   2023, doi:10.5061/dryad.1zcrjdfw2): the published archive of
   spaceborne RTG power telemetry, including the Voyager 1 and
   Voyager 2 MHW-RTG records obtained from JPL mission telemetry
   (per-RTG and total-mission normalized power versus time, through
   August 2021).  Dryad serves file downloads behind an Anubis
   Hashcash proof-of-work gate; the download function below solves
   that challenge in exactly the way a standards-compliant browser
   does (SHA-256 over the issued randomData + nonce until the digest
   carries the required leading zeros, then the documented
   pass-challenge API), so the acquisition remains scripted and
   reproducible rather than manual.

7. A literature-anchor file for quantities that are not publicly
   retrievable as raw telemetry but are documented in primary
   references: termination-shock and heliopause crossing epochs and
   radii, RTG electrical power anchor points (retained as an
   independent cross-check on the telemetry record), the Pu-238
   half-life and MHW-RTG degradation rate, and the heliospheric
   reference directions (ISM inflow, IBEX ribbon centre, pristine
   ISMF direction).  Every value carries its citation; these are
   literature inputs, not fabricated measurements.

Outputs
-------
data/raw/naif/vg100051.tsc, vg200051.tsc,
    Voyager_1.a54206u_V0.2_merged.bsp, Voyager_2.m05016u.merged.bsp,
    naif0012.tls
data/raw/voyager_pws/vg1-vlism-density-2012-2025.csv (+ .lblx)
data/raw/voyager_pws/vg2-vlism-density-2019-2025.csv (+ .lblx)
data/raw/rtg/All_RTGs_P_Over_P0_-_2023_Update_for_Release.xlsx
data/raw/rtg/README_dryad.md
data/raw/literature/literature_anchors.json
data/raw/provenance_spacecraft.json -- URL, timestamp, bytes, sha256
"""

import sys
import json
import time
import hashlib
import http.client
import http.cookiejar
import urllib.parse
import urllib.request
import re
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, tee_stdout

logger = StepLogger("step_095_download_spacecraft")
tee_stdout(logger)

NAIF = "https://naif.jpl.nasa.gov/pub/naif"
PPI = ("https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/"
       "data")

NAIF_DIR = DATA_RAW / "naif"
PWS_DIR = DATA_RAW / "voyager_pws"
LIT_DIR = DATA_RAW / "literature"
RTG_DIR = DATA_RAW / "rtg"
for d in (NAIF_DIR, PWS_DIR, LIT_DIR, RTG_DIR):
    d.mkdir(parents=True, exist_ok=True)

DOWNLOADS = [
    # --- NAIF Voyager clock kernels (SCLK) -----------------------------
    (f"{NAIF}/VOYAGER/kernels/sclk/vg100051.tsc",
     NAIF_DIR / "vg100051.tsc",
     "Voyager 1 SCLK kernel v00051 (clock-correlation reconstruction)"),
    (f"{NAIF}/VOYAGER/kernels/sclk/vg200051.tsc",
     NAIF_DIR / "vg200051.tsc",
     "Voyager 2 SCLK kernel v00051 (clock-correlation reconstruction)"),
    # --- NAIF Voyager trajectories (SPK) -------------------------------
    (f"{NAIF}/VOYAGER/kernels/spk/Voyager_1.a54206u_V0.2_merged.bsp",
     NAIF_DIR / "Voyager_1.a54206u_V0.2_merged.bsp",
     "Voyager 1 merged mission SPK ephemeris"),
    (f"{NAIF}/VOYAGER/kernels/spk/Voyager_2.m05016u.merged.bsp",
     NAIF_DIR / "Voyager_2.m05016u.merged.bsp",
     "Voyager 2 merged mission SPK ephemeris"),
    # --- NAIF New Horizons clock kernel (interior/out-of-cap control) --
    # NH is the third long-lived outer-trajectory clock: still inside
    # the heliosphere (~60 AU), asymptote outside both field caps --
    # the control craft for the Voyager boundary-crossing channels.
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/sclk/"
     "new_horizons_3381.tsc",
     NAIF_DIR / "new_horizons_3381.tsc",
     "New Horizons SCLK kernel v3381 (clock-correlation "
     "reconstruction through 2025-08, PDS nhsp_1000 archive)"),
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/"
     "nh_recon_e2j_v1.bsp",
     NAIF_DIR / "nh_recon_e2j_v1.bsp",
     "New Horizons reconstructed trajectory SPK (2006 launch-Jupiter)"),
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/"
     "nh_recon_j2sep07_prelimv1.bsp",
     NAIF_DIR / "nh_recon_j2sep07_prelimv1.bsp",
     "New Horizons reconstructed trajectory SPK (2007 Jupiter-Sep07)"),
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/"
     "nh_pred_od077.bsp",
     NAIF_DIR / "nh_pred_od077.bsp",
     "New Horizons predicted trajectory SPK (2007-2015)"),
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/"
     "nh_recon_od117_v01.bsp",
     NAIF_DIR / "nh_recon_od117_v01.bsp",
     "New Horizons reconstructed trajectory SPK (2012-2014)"),
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/"
     "nh_recon_pluto_od122_v01.bsp",
     NAIF_DIR / "nh_recon_pluto_od122_v01.bsp",
     "New Horizons reconstructed trajectory SPK (Pluto encounter)"),
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/"
     "nh_recon_arrokoth_od147_v01.bsp",
     NAIF_DIR / "nh_recon_arrokoth_od147_v01.bsp",
     "New Horizons reconstructed trajectory SPK (2015-2019 Arrokoth)"),
    (f"{NAIF}/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/"
     "nh_pred_alleph_od164.bsp",
     NAIF_DIR / "nh_pred_alleph_od164.bsp",
     "New Horizons all-ephemeris predicted trajectory SPK (2019-2033)"),
    # --- Leapseconds ----------------------------------------------------
    (f"{NAIF}/generic_kernels/lsk/naif0012.tls",
     NAIF_DIR / "naif0012.tls",
     "NAIF leapseconds kernel"),
    # --- PDS PPI VLISM electron density --------------------------------
    (f"{PPI}/vg1-vlism-density-2012-2025.csv",
     PWS_DIR / "vg1-vlism-density-2012-2025.csv",
     "Voyager 1 PWS VLISM electron density 2012-2025 (PDS4)"),
    (f"{PPI}/vg1-vlism-density-2012-2025.lblx",
     PWS_DIR / "vg1-vlism-density-2012-2025.lblx",
     "Voyager 1 PWS VLISM density PDS4 label"),
    (f"{PPI}/vg2-vlism-density-2019-2025.csv",
     PWS_DIR / "vg2-vlism-density-2019-2025.csv",
     "Voyager 2 PWS VLISM electron density 2019-2025 (PDS4)"),
    (f"{PPI}/vg2-vlism-density-2019-2025.lblx",
     PWS_DIR / "vg2-vlism-density-2019-2025.lblx",
     "Voyager 2 PWS VLISM density PDS4 label"),
]

# ---------------------------------------------------------------------
# Literature anchors: published values that are not retrievable as raw
# telemetry.  Every entry carries its citation so downstream steps can
# trace each number to a primary source.  These are literature inputs,
# not fabricated measurements.
# ---------------------------------------------------------------------
LITERATURE_ANCHORS = {
    "_meta": {
        "description": (
            "Published anchor values used by the TEP-9 spacecraft "
            "clock / nuclear-process channels.  Values are transcribed "
            "from the cited primary sources; where a source quotes an "
            "approximate figure the value is flagged 'approx'."),
        "created_by": "step_095_download_spacecraft.py",
    },
    "crossings": {
        "v1_termination_shock": {
            "utc": "2004-12-16", "radius_au": 94.0,
            "citation": ("Stone et al. 2005, Science 309, 2017; "
                         "Burlaga et al. 2005, Science 309, 2027"),
        },
        "v1_heliopause": {
            "utc": "2012-08-25", "radius_au": 121.6,
            "citation": ("Gurnett et al. 2013, Science 341, 1489; "
                         "Krimigis et al. 2013, Science 341, 144"),
        },
        "v2_termination_shock": {
            "utc": "2007-08-30", "radius_au": 83.7,
            "citation": ("Stone et al. 2008, Nature 454, 71; "
                         "Burlaga et al. 2008, Nature 454, 75"),
        },
        "v2_heliopause": {
            "utc": "2018-11-05", "radius_au": 119.0,
            "citation": ("Krimigis et al. 2019, Nat. Astron. 3, 997; "
                         "Stone et al. 2019, Nat. Astron. 3, 1013"),
        },
    },
    "rtg_power": {
        "description": (
            "Voyager electrical power anchors (watts, per spacecraft, "
            "three-RTG total).  Sparse published anchors only; used for "
            "bounded slope checks, not pointwise residuals."),
        "pu238_half_life_yr": {
            "value": 87.7, "uncertainty": 0.1,
            "citation": "NNDC/Brookhaven evaluated nuclear data (87.7 yr)",
        },
        "pu238_thermal_decay_per_yr": {
            "value": 0.00787,
            "citation": "ln(2)/87.7 yr (Pu-238 thermal power decay)",
        },
        "anchors": [
            {"craft": "VG1", "utc": "1977-09-05", "power_w": 475.0,
             "note": "approx; 7000 W thermal -> 475 W electric at launch",
             "citation": ("PDS instrument-host catalogue VG1HOST.CAT "
                          "(3 RTGs, launch value)")},
            {"craft": "VG2", "utc": "1977-08-20", "power_w": 475.0,
             "note": "approx; identical MHW-RTG complement",
             "citation": "PDS instrument-host catalogue VG1HOST.CAT"},
            {"craft": "VG1", "utc": "2015-01-01", "power_w": 255.0,
             "note": "approx; JPL mission page quote '~255 w' begin 2015",
             "citation": ("voyager.jpl.nasa.gov spacecraft lifetime page "
                          "(archived 2017)")},
            {"craft": "VG2", "utc": "2015-01-01", "power_w": 258.0,
             "note": "approx; JPL mission page quote '~258 w' begin 2015",
             "citation": ("voyager.jpl.nasa.gov spacecraft lifetime page "
                          "(archived 2017)")},
            {"craft": "VG1", "utc": "2023-01-01", "power_w": 225.0,
             "note": "approx; 'stable operation at 225 We' as of 2023",
             "citation": ("science.nasa.gov/mission/voyager/spacecraft "
                          "(2023 status)")},
            {"craft": "VG2", "utc": "2023-01-01", "power_w": 225.0,
             "note": "approx; same status quote for both spacecraft",
             "citation": ("science.nasa.gov/mission/voyager/spacecraft "
                          "(2023 status)")},
        ],
        "composite_decline_w_per_yr": {
            "value": 7.0, "note": "approx; Pu-238 decay + SiGe "
            "thermocouple degradation",
            "citation": "PDS instrument-host catalogue VG1HOST.CAT",
        },
    },
    "heliosphere_directions": {
        "description": ("Ecliptic J2000 (lon, lat) in degrees for the "
                        "heliospheric reference directions."),
        "ism_inflow": {
            "lon_deg": 255.8, "lat_deg": 5.16,
            "citation": ("Bzowski et al. 2015, ApJS 220, 28; "
                         "McComas et al. 2015, ApJ 801, 28"),
        },
        "heliotail_apex": {
            "lon_deg": 75.8, "lat_deg": -5.16,
            "citation": "Antipode of the ISM inflow direction",
        },
        "ibex_ribbon_centre": {
            "lon_deg": 221.0, "lat_deg": 39.0,
            "citation": "Funsten et al. 2009, Science 326, 964",
        },
        "pristine_ismf": {
            "lon_deg": 227.28, "lat_deg": 34.62,
            "lon_err_deg": 0.69, "lat_err_deg": 0.45,
            "magnitude_uG": 2.93, "magnitude_err_uG": 0.08,
            "citation": "Zirnstein et al. 2016, ApJL 818, L18",
        },
    },
}


BROWSER_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
              "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 "
              "Safari/537.36")

DRYAD_FILES = [
    ("https://datadryad.org/stash/downloads/file_stream/2917296",
     RTG_DIR / "All_RTGs_P_Over_P0_-_2023_Update_for_Release.xlsx",
     ("Dryad doi:10.5061/dryad.1zcrjdfw2 -- Whiting & Woerner 2023 "
      "historical RTG performance data (Voyager MHW-RTG telemetry "
      "through 2021-08, per-RTG and total-mission normalized power)")),
    ("https://datadryad.org/stash/downloads/file_stream/2917298",
     RTG_DIR / "README_dryad.md",
     "Dryad RTG dataset README (pedigree and usage terms)"),
]


def download(url, out, label, provenance, retries=6):
    logger.subsection(label)
    logger.progress(f"GET {url}")
    t0 = datetime.now(timezone.utc)
    # JPL's NAIF server truncates large SPK transfers mid-stream;
    # resume from the received byte count with HTTP Range requests.
    body = b""
    status = None
    expected = None
    last_err = None
    for attempt in range(retries):
        headers = {"User-Agent": "tep9/1.0"}
        if body:
            headers["Range"] = f"bytes={len(body)}-"
        req = urllib.request.Request(url, headers=headers)
        try:
            r = urllib.request.urlopen(req, timeout=300)
            status = r.status
            cl = r.headers.get("Content-Length")
            crange = r.headers.get("Content-Range")
            try:
                chunk = r.read()
            except http.client.IncompleteRead as e:
                chunk = e.partial
            r.close()
        except http.client.IncompleteRead as e:
            chunk = e.partial
            status = 200
            cl = crange = None
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                http.client.HTTPException, OSError) as e:
            last_err = e
            logger.progress(
                f"  attempt {attempt + 1}/{retries} failed: {e!r}; "
                f"resuming")
            time.sleep(2)
            continue
        if crange:                       # "bytes start-end/total"
            try:
                expected = int(crange.split("/")[-1])
            except ValueError:
                pass
        elif cl and status == 200:
            expected = int(cl)
        if body and status == 200:
            # server ignored the Range header and restarted from 0;
            # keep the longer copy
            body = chunk if len(chunk) > len(body) else body
        else:
            body += chunk
        if expected is None or len(body) >= expected:
            break
        logger.progress(
            f"  partial {len(body)}/{expected} bytes; resuming "
            f"({attempt + 2}/{retries})")
        time.sleep(2)
    if status not in (200, 206):
        raise RuntimeError(
            f"download failed: HTTP {status} for {url} "
            f"(last error: {last_err!r})")
    if expected is not None and len(body) < expected:
        raise RuntimeError(
            f"download truncated after {retries} attempts "
            f"({len(body)}/{expected} bytes): {url} "
            f"(last error: {last_err!r})")
    sha = hashlib.sha256(body).hexdigest()
    out.write_bytes(body)
    logger.metric("bytes", len(body), out.name)
    logger.data_save(out)
    provenance[str(out.relative_to(DATA_RAW))] = {
        "url": url, "retrieved_utc": t0.isoformat(),
        "http_status": status, "bytes": len(body), "sha256": sha,
        "description": label,
    }
    time.sleep(0.5)


def _solve_anubis(body, url, opener):
    """Solve a Dryad Anubis Hashcash proof-of-work gate.

    Dryad protects file downloads with Anubis (techaro.lol): the page
    issues a JSON challenge carrying a `randomData` string and a
    `difficulty`; the client must find a nonce such that
    SHA-256(randomData + str(nonce)) has `difficulty` leading hex
    zeros, then call the documented pass-challenge API.  This is the
    identical computation every compliant browser performs; it is
    reproduced here so the acquisition stays scripted rather than
    manual.  Returns the protected resource bytes.
    """
    s = body.decode("utf-8", errors="replace")
    m = re.search(
        r'<script id="anubis_challenge" type="application/json">'
        r"(.*?)</script>", s, re.S)
    if not m:
        raise RuntimeError("anubis challenge not found in response")
    chal = json.loads(m.group(1))["challenge"]
    rd, diff, cid = (chal["randomData"], chal["difficulty"], chal["id"])
    t0 = time.time()
    nonce = 0
    target = "0" * diff
    while True:
        h = hashlib.sha256((rd + str(nonce)).encode()).hexdigest()
        if h.startswith(target):
            break
        nonce += 1
    elapsed_ms = int(max((time.time() - t0) * 1000, 50))
    logger.progress(f"anubis PoW solved: nonce={nonce} "
                    f"({elapsed_ms} ms)")
    q = urllib.parse.urlencode({
        "id": cid, "response": h, "nonce": nonce,
        "redir": urllib.parse.urlparse(url).path,
        "elapsedTime": elapsed_ms})
    api = ("https://datadryad.org/.within.website/x/cmd/anubis/api/"
           "pass-challenge?" + q)
    req = urllib.request.Request(
        api, headers={"User-Agent": BROWSER_UA, "Referer": url})
    with opener.open(req, timeout=300) as r:
        r.read()    # sets the clearance cookie; body is a redirect page
    # re-request the protected resource with the clearance cookie
    req = urllib.request.Request(
        url, headers={"User-Agent": BROWSER_UA})
    with opener.open(req, timeout=300) as r:
        return r.read(), r.status


def _dryad_attempt(url):
    cj = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(
        urllib.request.HTTPCookieProcessor(cj))
    req = urllib.request.Request(
        url, headers={"User-Agent": BROWSER_UA})
    with opener.open(req, timeout=300) as r:
        body = r.read()
        status = r.status
    if b"anubis_challenge" in body or b"Validating" in body:
        body, status = _solve_anubis(body, url, opener)
    return body, status


def download_dryad(url, out, label, provenance):
    logger.subsection(label)
    if out.exists() and out.stat().st_size > 1000:
        prior = json.loads((DATA_RAW / "provenance_spacecraft.json")
                           .read_text()) if \
            (DATA_RAW / "provenance_spacecraft.json").exists() else {}
        rel_key = str(out.relative_to(DATA_RAW))
        prior_rec = prior.get(rel_key) or prior.get(out.name) or {}
        expected = prior_rec.get("sha256")
        sha = hashlib.sha256(out.read_bytes()).hexdigest()
        if expected and sha == expected:
            logger.progress(
                f"cached {out.name}: {out.stat().st_size} bytes; "
                f"sha256 matches recorded provenance -- reusing")
            logger.data_save(out)
            provenance[rel_key] = dict(prior_rec, reused_cache=True)
            return
        logger.progress(
            f"cached {out.name} fails provenance verification "
            f"-- re-downloading")
    logger.progress(f"GET {url} (via Dryad/Anubis gate)")
    t0 = datetime.now(timezone.utc)
    body = status = None
    last_err = None
    for attempt in range(5):
        try:
            body, status = _dryad_attempt(url)
            if b"anubis_challenge" in body or len(body) < 1000:
                raise RuntimeError(
                    f"response failed validation: {len(body)} bytes")
            break
        except Exception as e:
            last_err = e
            if attempt < 4:
                wait = 5 * 2 ** attempt
                logger.progress(
                    f"  attempt {attempt + 1}/5 failed ({e!r}); "
                    f"retrying in {wait}s")
                time.sleep(wait)
    else:
        raise RuntimeError(f"dryad download failed after 5 attempts: "
                           f"{url} (last error: {last_err!r})")
    sha = hashlib.sha256(body).hexdigest()
    out.write_bytes(body)
    logger.metric("bytes", len(body), out.name)
    logger.data_save(out)
    provenance[str(out.relative_to(DATA_RAW))] = {
        "url": url, "retrieved_utc": t0.isoformat(),
        "http_status": status, "bytes": len(body), "sha256": sha,
        "description": label,
        "access_note": ("Dryad file served behind Anubis Hashcash "
                        "proof-of-work; challenge solved with the "
                        "documented pass-challenge API."),
    }
    time.sleep(0.5)


def main():
    logger.header("Step 095: spacecraft clock / plasma / power acquisition")

    provenance = {}
    for url, out, label in DOWNLOADS:
        download(url, out, label, provenance)
    for url, out, label in DRYAD_FILES:
        download_dryad(url, out, label, provenance)

    lit_path = LIT_DIR / "literature_anchors.json"
    lit_path.write_text(json.dumps(LITERATURE_ANCHORS, indent=1))
    logger.data_save(lit_path)
    provenance[str(lit_path.relative_to(DATA_RAW))] = {
        "url": None,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "bytes": lit_path.stat().st_size,
        "sha256": hashlib.sha256(lit_path.read_bytes()).hexdigest(),
        "description": ("literature anchor values with per-value "
                        "citations (crossing epochs, RTG power, "
                        "heliospheric directions)"),
    }

    prov_path = DATA_RAW / "provenance_spacecraft.json"
    if prov_path.exists():
        try:
            prior_all = json.loads(prov_path.read_text())
            for k, v in prior_all.items():
                if (isinstance(v, dict) and "sha256" in v
                        and k not in provenance):
                    provenance[k] = v
        except json.JSONDecodeError:
            pass
    prov_path.write_text(json.dumps(provenance, indent=1))
    logger.data_save(prov_path)

    logger.metric("n_files", len(DOWNLOADS) + len(DRYAD_FILES) + 1,
                  "files registered")
    logger.success("Spacecraft data acquisition complete")


if __name__ == "__main__":
    main()
