#!/usr/bin/env python3
"""
TEP-9 step 002 -- CDS catalogue download (OSSOS + Warsaw)
==========================================================

Downloads the two published survey catalogues used for bias
calibration and comet orbit reconstruction:

1. OSSOS characterized ensemble (Bannister et al. 2018, ApJS 236, 18;
   VizieR J/ApJS/236/18).  The t3char table is fetched as a VOTable
   via the VizieR votable endpoint, and the authoritative fixed-width
   t3char.dat + ReadMe are fetched from the CDS archive for provenance
   cross-checking.

2. Warsaw catalogue of cometary orbits (Krolikowska 2014, A&A 567,
   A126; VizieR J/A+A/567/A126): tables a1, b, b4, c, d plus the CDS
   ReadMe, from the CDS archive FTP.

3. One-apparition comet catalogue (Krolikowska 2014, A&A 571, A63;
   VizieR J/A+A/571/A63): osculating elements (tableb1), original/
   future 1/a with uncertainties and the Marsden-Williams 2008
   comparison values (tablec1), and original/future barycentric
   elements with per-element uncertainties (tabled1/tablee1).  These
   supply the uncertainty-normalized and third-lineage checks.

Outputs
-------
data/raw/ossos/ossos_t3char.vot
data/raw/ossos/t3char.dat, ReadMe         (provenance cross-check)
data/raw/warsaw/warsaw_{a1,b,b4,c,d}.dat, warsaw_ReadMe
data/raw/lpc/lpc_{osc,orig,fut,aaori}_2006_2010.vot
data/raw/{ossos,warsaw,lpc}/provenance.json
"""

import sys
import json
import time
import hashlib
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, tee_stdout

logger = StepLogger("step_002_download_cds")
tee_stdout(logger)

OSSOS_VOT_URL = ("https://vizier.cds.unistra.fr/viz-bin/votable"
                 "?-source=J/ApJS/236/18/t3char&-out.max=unlimited")
LPC_VOT_URL = ("https://vizier.cds.unistra.fr/viz-bin/votable"
               "?-source=J/A%2BA/571/A63/{tab}"
               "&-out.max=unlimited&-out.all")
CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp"

DOWNLOADS = [
    # (url, subdir, filename, note)
    (OSSOS_VOT_URL, "ossos", "ossos_t3char.vot",
     "OSSOS characterized ensemble, VOTable export (J/ApJS/236/18)"),
    (f"{CDS_FTP}/J/ApJS/236/18/t3char.dat", "ossos", "t3char.dat",
     "OSSOS t3char authoritative fixed-width table"),
    (f"{CDS_FTP}/J/ApJS/236/18/ReadMe", "ossos", "ReadMe",
     "OSSOS CDS ReadMe"),
    (f"{CDS_FTP}/J/A+A/567/A126/ReadMe", "warsaw", "warsaw_ReadMe",
     "Warsaw catalogue CDS ReadMe (J/A+A/567/A126)"),
    (f"{CDS_FTP}/J/A+A/567/A126/tablea1.dat", "warsaw", "warsaw_tablea1.dat",
     "Warsaw table A1: observational material / quality class"),
    (f"{CDS_FTP}/J/A+A/567/A126/tableb.dat", "warsaw", "warsaw_tableb.dat",
     "Warsaw table B: osculating heliocentric elements"),
    (f"{CDS_FTP}/J/A+A/567/A126/tableb4.dat", "warsaw", "warsaw_tableb4.dat",
     "Warsaw table B4: osculating elements, extended sample"),
    (f"{CDS_FTP}/J/A+A/567/A126/tablec.dat", "warsaw", "warsaw_tablec.dat",
     "Warsaw table C: original barycentric elements"),
    (f"{CDS_FTP}/J/A+A/567/A126/tabled.dat", "warsaw", "warsaw_tabled.dat",
     "Warsaw table D: future barycentric elements"),
    (LPC_VOT_URL.format(tab="tableb1"), "lpc", "lpc_osc_2006_2010.vot",
     "J/A+A/571/A63 tableb1: osculating elements, 1901-1950 cohort"),
    (LPC_VOT_URL.format(tab="tabled1"), "lpc", "lpc_orig_2006_2010.vot",
     "J/A+A/571/A63 tabled1: original barycentric elements + errors"),
    (LPC_VOT_URL.format(tab="tablee1"), "lpc", "lpc_fut_2006_2010.vot",
     "J/A+A/571/A63 tablee1: future barycentric elements + errors"),
    (LPC_VOT_URL.format(tab="tablec1"), "lpc", "lpc_aaori_2006_2010.vot",
     "J/A+A/571/A63 tablec1: orig/fut 1/a + errors + MW08 values"),
]


def fetch(url, dest, note, prov):
    logger.subsection(note)
    logger.progress(f"GET {url[:120]}")
    req = urllib.request.Request(url, headers={"User-Agent": "tep9/1.0"})
    t0 = datetime.now(timezone.utc)
    body = status = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                body = r.read()
                status = r.status
            if status == 200:
                break
            raise RuntimeError(f"download failed: HTTP {status} for {url}")
        except Exception as exc:
            if attempt == 5:
                raise
            wait = 5 * 2 ** (attempt - 1)
            logger.progress(f"attempt {attempt} failed ({exc}); "
                            f"retrying in {wait}s")
            time.sleep(wait)
    if status != 200:
        raise RuntimeError(f"download failed: HTTP {status} for {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(body)
    sha = hashlib.sha256(body).hexdigest()
    logger.metric("bytes", len(body), dest.name)
    logger.data_save(dest)
    prov[dest.name] = {"url": url, "retrieved_utc": t0.isoformat(),
                       "http_status": status, "bytes": len(body),
                       "sha256": sha, "note": note}
    # refresh flat compat symlink at data/raw root
    link = DATA_RAW / dest.name
    if link.is_symlink() or not link.exists():
        link.unlink(missing_ok=True)
        link.symlink_to(f"{dest.parent.name}/{dest.name}")
    return body


def main():
    logger.header("CDS catalogue download (OSSOS + Warsaw)")
    provs = {"ossos": {}, "warsaw": {}, "lpc": {}}
    for url, sub, fname, note in DOWNLOADS:
        body = fetch(url, DATA_RAW / sub / fname, note, provs[sub])
        if fname == "ossos_t3char.vot":
            n_rows = body.count(b"<TR>")
            logger.metric("votable_rows", n_rows,
                          "OSSOS characterized detections")
            if n_rows < 500:
                raise RuntimeError(
                    f"OSSOS VOTable looks truncated ({n_rows} rows)")
    for sub, prov in provs.items():
        prov_path = DATA_RAW / sub / "provenance.json"
        if prov_path.exists():
            try:
                prior = json.loads(prov_path.read_text())
                for k, v in prior.get("files", {}).items():
                    prov.setdefault(k, v)
            except json.JSONDecodeError:
                pass
        doc = {"step": "step_002_download_cds",
               "generated_utc": datetime.now(timezone.utc).isoformat(),
               "files": prov}
        prov_path.write_text(json.dumps(doc, indent=2))
        logger.data_save(prov_path)
    logger.progress("CDS download complete -> data/raw/{ossos,warsaw}/")


if __name__ == "__main__":
    main()
