#!/usr/bin/env python3
"""
TEP-9 step 003 -- CODE catalogue acquisition
=============================================

Acquires the CODE catalogue (Catalogue of Orbits and their Dynamical
Evolution; Krolikowska & Dybczynski 2020, A&A 640, A97) hosted on
the Poznan comet-dynamics server pad2.astro.amu.edu.pl.

The catalogue is served through a PHP search interface that renders
three orbit-solution tables as HTML:

    code_original.html    -- original barycentric orbits
    code_osculating.html  -- osculating heliocentric orbits
    code_future.html      -- future barycentric orbits

Reproduction policy
-------------------
1. The live catalogue endpoints are attempted first.  The pad2
   interface is a sessioned search application; if a live export is
   returned it is written with full provenance.
2. If the live endpoint is unreachable or returns a non-table page,
   the versioned archival snapshots staged in data/raw/code/ are
   adopted.  Their SHA-256 hashes are recorded and logged loudly so
   provenance is never ambiguous.  These snapshots are the published
   catalogue content, not synthetic data.

Outputs
-------
data/raw/code/code_{original,osculating,future}.html
data/raw/code/provenance.json
"""

import sys
import json
import hashlib
import urllib.request
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, tee_stdout

logger = StepLogger("step_003_download_code")
tee_stdout(logger)

CODE_DIR = DATA_RAW / "code"
PAD2 = "https://pad2.astro.amu.edu.pl/comets/search.php"

# Candidate live export endpoints for the three orbit tables.
LIVE_ATTEMPTS = {
    "code_original.html": {
        "orbit_type": "original", "catalog": "CODE", "submit": "Search"},
    "code_osculating.html": {
        "orbit_type": "osculating", "catalog": "CODE", "submit": "Search"},
    "code_future.html": {
        "orbit_type": "future", "catalog": "CODE", "submit": "Search"},
}

MIN_TABLE_BYTES = 100_000  # a real CODE dump is ~700 KB of HTML


def looks_like_code_table(body):
    """A genuine CODE export contains the orbit-element table rows."""
    if len(body) < MIN_TABLE_BYTES:
        return False
    return body.count(b"<tr") > 200 and b"CODE" in body


def try_live(fname, params):
    """Attempt a live pad2 export; return body or None."""
    try:
        data = urllib.parse.urlencode(params).encode()
        req = urllib.request.Request(PAD2, data=data,
                                     headers={"User-Agent": "tep9/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read()
        if looks_like_code_table(body):
            return body
        logger.progress(
            f"live endpoint returned non-table page for {fname} "
            f"({len(body)} bytes)")
        return None
    except Exception as exc:
        logger.progress(f"live endpoint failed for {fname}: {exc}")
        return None


def main():
    logger.header("CODE catalogue acquisition")
    CODE_DIR.mkdir(parents=True, exist_ok=True)
    prov = {"step": "step_003_download_code",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "catalogue": "CODE (Krolikowska & Dybczynski 2020, "
                         "A&A 640, A97)",
            "host": "pad2.astro.amu.edu.pl", "files": {}}

    for fname, params in LIVE_ATTEMPTS.items():
        dest = CODE_DIR / fname
        logger.subsection(fname)
        body = try_live(fname, params)
        if body is not None:
            dest.write_bytes(body)
            prov["files"][fname] = {
                "source": "live pad2 export",
                "url": PAD2, "params": params,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "bytes": len(body),
                "sha256": hashlib.sha256(body).hexdigest()}
            logger.metric("bytes", len(body), "live export")
            logger.data_save(dest)
        elif dest.exists():
            sha = hashlib.sha256(dest.read_bytes()).hexdigest()
            n_rows = dest.read_bytes().count(b"<tr")
            logger.progress(
                f"USING ARCHIVAL SNAPSHOT for {fname}: pad2 live "
                f"export unavailable; staged catalogue copy adopted")
            logger.metric("bytes", dest.stat().st_size, "archival")
            logger.metric("table_rows", n_rows, "archival")
            logger.metric("sha256", sha, "archival integrity")
            prov["files"][fname] = {
                "source": "archival snapshot (pad2 CODE catalogue)",
                "note": "live export interface did not return a table; "
                        "staged snapshot of the published catalogue used",
                "bytes": dest.stat().st_size,
                "table_rows": n_rows,
                "sha256": sha}
        else:
            raise RuntimeError(
                f"{fname}: no live export and no archival snapshot -- "
                "CODE catalogue unavailable. Fetch the three orbit "
                "tables from pad2.astro.amu.edu.pl and stage them in "
                "data/raw/code/.")
        link = DATA_RAW / fname
        if link.is_symlink() or not link.exists():
            link.unlink(missing_ok=True)
            link.symlink_to(f"code/{fname}")

    prov_path = CODE_DIR / "provenance.json"
    if prov_path.exists():
        try:
            prior = json.loads(prov_path.read_text())
            for k, v in prior.get("files", {}).items():
                prov["files"].setdefault(k, v)
        except json.JSONDecodeError:
            pass
    prov_path.write_text(json.dumps(prov, indent=2))
    logger.data_save(prov_path)
    logger.progress("CODE acquisition complete -> data/raw/code/")


if __name__ == "__main__":
    main()
