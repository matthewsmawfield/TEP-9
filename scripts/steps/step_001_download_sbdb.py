#!/usr/bin/env python3
"""
TEP-9 step 001 -- JPL SBDB download
====================================

Downloads the JPL Small-Body Database (SBDB) orbital-element catalog
for the outer solar system via the public sbdb_query API
(https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html).

Three real-data samples are fetched:

1. Outer solar system asteroids: all bodies with heliocentric
   semimajor axis a >= 30 AU (classical Kuiper belt, resonant,
   scattered, detached and extreme trans-Neptunian objects).
2. Centaurs (orbit class CEN): the Neptune-crossing population,
   used for the cross-population coherence test.
3. Comets reaching the outer system (q >= 5 AU): the Jupiter-family
   and Halley-type population used for anti-axis directionality.

Outputs
-------
data/raw/sbdb/sbdb_outer_ss.json
data/raw/sbdb/sbdb_centaurs.json
data/raw/sbdb/sbdb_comets_ext.json
data/raw/sbdb/provenance.json  -- URL, timestamp, sha256, row counts
"""

import sys
import json
import time
import hashlib
import urllib.request
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout

logger = StepLogger("step_001_download_sbdb")
tee_stdout(logger)

API = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
OUT = DATA_RAW / "sbdb"

FIELDS = ["full_name", "a", "e", "i", "om", "w", "q", "ad", "H",
          "epoch", "class", "condition_code", "data_arc",
          "n_obs_used", "first_obs", "last_obs", "per"]


def retain_snapshot(label, outname, provenance, exc):
    """Live download exhausted its retries: keep the existing
    snapshot if one is present and valid, and say so explicitly in
    the log and provenance.  Raises when no usable snapshot exists --
    a first-time acquisition failure must abort the step rather than
    proceed silently on missing data."""
    out = OUT / outname
    if not out.exists():
        raise RuntimeError(
            f"{label}: live SBDB download failed ({exc}) and no "
            f"existing snapshot at {out}")
    try:
        d = json.loads(out.read_text())
        rows = len(d.get("data", []))
    except Exception as jexc:
        raise RuntimeError(
            f"{label}: live SBDB download failed ({exc}) and existing "
            f"snapshot {out} is unreadable ({jexc})")
    if rows == 0:
        raise RuntimeError(
            f"{label}: live SBDB download failed ({exc}) and existing "
            f"snapshot {out} contains no data rows")
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    logger.warning(
        f"{label}: live SBDB download failed ({exc}); retaining "
        f"existing snapshot {outname} ({rows} rows, "
        f"sha256={sha[:12]}...)")
    logger.progress("NOTE: catalogue content is UNCHANGED from the "
                    "previously downloaded snapshot; downstream "
                    "results remain fully traceable via the retained "
                    "sha256.")
    provenance[outname] = {
        "retained_existing_snapshot": True,
        "live_download_error": str(exc),
        "sha256": sha, "rows": rows,
        "bytes": out.stat().st_size,
        "api": "JPL SBDB sbdb_query.api",
    }
    return d


def sbdb_query(params, label, outname, provenance):
    url = API + "?" + urllib.parse.urlencode(params)
    logger.subsection(f"{label}")
    logger.progress(f"GET {url[:120]}...")
    req = urllib.request.Request(url, headers={"User-Agent": "tep9/1.0"})
    t0 = datetime.now(timezone.utc)
    body = status = None
    last_exc = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                body = r.read()
                status = r.status
            if status == 200:
                break
            raise RuntimeError(f"SBDB query failed: HTTP {status}")
        except Exception as exc:
            last_exc = exc
            if attempt == 5:
                break
            wait = 5 * 2 ** (attempt - 1)
            logger.progress(f"attempt {attempt} failed ({exc}); "
                            f"retrying in {wait}s")
            time.sleep(wait)
    if status != 200 or body is None:
        return retain_snapshot(label, outname, provenance, last_exc)
    d = json.loads(body.decode())
    n = d.get("count", len(d.get("data", [])))
    text = json.dumps(d, indent=1)
    out = OUT / outname
    out.write_text(text)
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    logger.metric("rows", len(d.get("data", [])), "objects returned")
    logger.metric("count", n, "catalogue count field")
    logger.data_save(out)
    provenance[outname] = {
        "url": url, "retrieved_utc": t0.isoformat(),
        "http_status": status, "bytes": out.stat().st_size,
        "wire_bytes": len(body),
        "wire_sha256": hashlib.sha256(body).hexdigest(),
        "sha256": sha,
        "rows": len(d.get("data", [])), "count_field": n,
        "api": "JPL SBDB sbdb_query.api",
    }
    # refresh flat compat symlink at data/raw root
    link = DATA_RAW / outname
    if link.is_symlink() or not link.exists():
        link.unlink(missing_ok=True)
        link.symlink_to(f"sbdb/{outname}")
    return d


def main():
    logger.header("JPL SBDB catalogue download")
    OUT.mkdir(parents=True, exist_ok=True)
    provenance = {"step": "step_001_download_sbdb",
                  "generated_utc": datetime.now(timezone.utc).isoformat(),
                  "fields": FIELDS, "files": {}}
    prov_files = {}

    # 1) Outer solar system: a >= 30 AU, asteroids only
    cdata = json.dumps({"AND": ["a|GE|30"]})
    sbdb_query({"fields": ",".join(FIELDS), "sb-kind": "a",
                "sb-cdata": cdata},
               "outer solar system asteroids (a>=30 AU)",
               "sbdb_outer_ss.json", prov_files)
    time.sleep(1)

    # 2) Centaurs (orbit class CEN)
    sbdb_query({"fields": ",".join(FIELDS), "sb-class": "CEN"},
               "Centaurs (class CEN)", "sbdb_centaurs.json", prov_files)
    time.sleep(1)

    # 3) Comets reaching the outer system: q >= 5 AU
    cdata_c = json.dumps({"AND": ["q|GE|5"]})
    sbdb_query({"fields": ",".join(FIELDS), "sb-kind": "c",
                "sb-cdata": cdata_c},
               "outer-reaching comets (q>=5 AU)",
               "sbdb_comets_ext.json", prov_files)

    prov_path = OUT / "provenance.json"
    if prov_path.exists():
        try:
            prior = json.loads(prov_path.read_text())
            for k, v in prior.get("files", {}).items():
                prov_files.setdefault(k, v)
        except json.JSONDecodeError:
            pass
    provenance["files"] = prov_files
    prov_path.write_text(json.dumps(provenance, indent=2))
    logger.data_save(prov_path)
    logger.progress("SBDB download complete -> data/raw/sbdb/")


if __name__ == "__main__":
    main()
