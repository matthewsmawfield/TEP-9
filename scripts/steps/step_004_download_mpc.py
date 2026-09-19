#!/usr/bin/env python3
"""
TEP-9 step 004 -- MPC / JPL element files
==========================================

Downloads the two public element files used by the injection-chain
and catalogue-concordance steps:

1. MPCORB.DAT.gz -- the Minor Planet Center's full orbit catalogue
   (all numbered and multi-opposition minor planets).  Used to
   recompute the TNO cluster on MPC orbit solutions.

2. CometEls.txt -- the JPL comet elements file used for the
   injected-population (JFC) depth chain.

3. de440s.bsp -- the JPL DE440s planetary ephemeris (NAIF generic
   kernel, coverage 1849-2150).  Used by step_042 to draw
   solar-system barycentric planet states at each comet's
   perihelion epoch.

Outputs
-------
data/raw/mpc/MPCORB.DAT.gz
data/raw/mpc/CometEls.txt
data/raw/mpc/provenance.json
data/raw/spice/de440s.bsp
data/raw/spice/provenance.json
"""

import sys
import json
import hashlib
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, tee_stdout

logger = StepLogger("step_004_download_mpc")
tee_stdout(logger)

OUT = DATA_RAW / "mpc"
SPICE_OUT = DATA_RAW / "spice"

DOWNLOADS = [
    ("https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz",
     "MPCORB.DAT.gz", "MPC full orbit catalogue (gzipped fixed-width)"),
    ("https://minorplanetcenter.net/iau/MPCORB/CometEls.txt",
     "CometEls.txt", "MPC comet elements file (MPCAT format)"),
]

SPICE_DOWNLOADS = [
    ("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
     "de440s.bsp", "JPL DE440s ephemeris, coverage 1849-2150"),
]


def fetch(url, dest, note, prov):
    logger.subsection(note)
    logger.progress(f"GET {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    t0 = datetime.now(timezone.utc)
    body = None
    buf = bytearray()
    for attempt in range(1, 6):
        try:
            headers = {"User-Agent": "tep9/1.0"}
            if buf:
                headers["Range"] = f"bytes={len(buf)}-"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=600) as r:
                if r.status not in (200, 206):
                    raise RuntimeError(f"HTTP {r.status}")
                if r.status == 200 and buf:
                    # server ignored the Range request -> restart
                    buf = bytearray()
                expected = r.headers.get("Content-Length")
                expected = int(expected) if expected else None
                if r.status == 206 and expected is not None:
                    expected += len(buf)
                while True:
                    block = r.read(1 << 20)
                    if not block:
                        break
                    buf.extend(block)
            if expected is None or len(buf) >= expected:
                body = bytes(buf)
                break
            logger.progress(
                f"attempt {attempt}: truncated at {len(buf)} of "
                f"{expected} bytes; resuming")
        except Exception as exc:
            logger.progress(f"attempt {attempt} failed at {len(buf)} bytes: {exc}")
    if body is not None and dest.suffix == ".gz":
        # integrity: the gzip stream must reach its end marker
        import zlib
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            d.decompress(body)
            d.flush()
            ok = d.eof
        except zlib.error:
            ok = False
        if not ok:
            logger.progress(
                f"gzip integrity check FAILED for {dest.name} "
                f"({len(body)} bytes); treating as incomplete")
            body = None
    if body is not None:
        dest.write_bytes(body)
    elif dest.exists():
        sha = hashlib.sha256(dest.read_bytes()).hexdigest()
        logger.progress(
            f"USING ARCHIVAL SNAPSHOT for {dest.name}: live download "
            f"incomplete; staged copy adopted (sha256={sha[:16]}...)")
        prov[dest.name] = {"url": url,
                           "retrieved_utc": t0.isoformat(),
                           "source": "archival snapshot (live download "
                                     "truncated by remote server)",
                           "bytes": dest.stat().st_size, "sha256": sha,
                           "note": note}
        link = DATA_RAW / dest.name
        if link.is_symlink() or not link.exists():
            link.unlink(missing_ok=True)
            link.symlink_to(f"{dest.parent.name}/{dest.name}")
        return
    else:
        raise RuntimeError(f"download failed for {url}: all attempts "
                           "truncated and no archival snapshot present")
    sha = hashlib.sha256(body).hexdigest()
    logger.metric("bytes", len(body), dest.name)
    logger.data_save(dest)
    prov[dest.name] = {"url": url, "retrieved_utc": t0.isoformat(),
                       "bytes": len(body),
                       "sha256": sha, "note": note}
    link = DATA_RAW / dest.name
    if link.is_symlink() or not link.exists():
        link.unlink(missing_ok=True)
        link.symlink_to(f"{dest.parent.name}/{dest.name}")


def main():
    logger.header("MPC / JPL element-file download")
    OUT.mkdir(parents=True, exist_ok=True)
    prov = {}
    for url, fname, note in DOWNLOADS:
        fetch(url, OUT / fname, note, prov)
    doc = {"step": "step_004_download_mpc",
           "generated_utc": datetime.now(timezone.utc).isoformat(),
           "files": prov}
    prov_path = OUT / "provenance.json"
    prov_path.write_text(json.dumps(doc, indent=2))
    logger.data_save(prov_path)
    logger.progress("MPC download complete -> data/raw/mpc/")

    SPICE_OUT.mkdir(parents=True, exist_ok=True)
    sprov = {}
    for url, fname, note in SPICE_DOWNLOADS:
        fetch(url, SPICE_OUT / fname, note, sprov)
    sdoc = {"step": "step_004_download_mpc",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "files": sprov}
    sprov_path = SPICE_OUT / "provenance.json"
    sprov_path.write_text(json.dumps(sdoc, indent=2))
    logger.data_save(sprov_path)
    logger.progress("SPICE kernel download complete -> data/raw/spice/")


if __name__ == "__main__":
    main()
