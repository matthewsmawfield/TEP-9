"""Step 162 -- NIMA ephemeris-residual pilot on the detached-resident
population.

The resident channel of the paper rests on discovery-geometry
clustering of the detached extreme TNOs.  A logically distinct check
is whether the *astrometric residuals* of those same residents carry
boundary-organized structure.  The Lucky Star / NIMA programme
maintains independent orbit solutions and publishes per-observation
O-C files (omc_ast_NIMAvN.txt) for its tracked object list, with the
underlying astrometry drawn from MPC plus dedicated WFI/ESO runs --
so independence from the comet catalogues is real but partial.

This pilot intersects the step_b77 resident ledger with the NIMA
object list (parsed from the nima.php page source), downloads the
highest-version O-C file per object, parses the fixed-format rows,
and builds per-object residual statistics: level, RMS, chi
distribution, and a linear drift of O-C versus epoch (an unmodelled
ephemeris offset surfaces as a secular slope).  Results are joined
to the ledger's axis separation and cap membership.

Registered tests

  T1  coverage ledger: how many detached-resident objects NIMA
      tracks with usable O-C files, split by cap membership.
  T2  per-object residual level and drift: |O-C| median/RMS and the
      fitted O-C-vs-epoch slope with its uncertainty, per object.
  T3  in-cap vs out-of-cap contrast on residual level and on |drift|
      -- reported descriptively; with ~7 in-cap objects and ~1
      control the channel is underpowered by construction and is
      registered as a pilot, not a detection.

Outputs
  results/step_b128_nima_pilot.json    per-object summary + tests
  results/step_b128_nima_pilot.csv     per-observation O-C records
  data/raw/nima/<obj>/omc_ast_NIMAvN.txt   downloaded O-C files
  data/raw/nima/nima.php.html              object-list page source
  data/raw/nima/provenance.json            URL + sha256 ledger
"""

import csv
import hashlib
import json
import math
import re
import sys as _sys
import urllib.request
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

import numpy as np

from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import tee_stdout
from scripts.utils.tep9_common import DATA_RAW, RESULTS

logger = StepLogger("step_162_nima_pilot")
tee_stdout(logger)

NIMA_DIR = DATA_RAW / "nima"
NIMA_BASE = "https://lesia.obspm.fr/lucky-star/"
NIMA_DATA = NIMA_BASE + "data/nima/"
PROV_PATH = NIMA_DIR / "provenance.json"

logger.header("NIMA ephemeris-residual pilot on detached residents")


def fetch(url, dest, required=True):
    """Download url -> dest with a sha256 provenance record.
    Returns None on 404 when required=False."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        blob = dest.read_bytes()
    else:
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                blob = r.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 404 and not required:
                return None
            raise
        dest.write_bytes(blob)
    prov = {}
    if PROV_PATH.exists():
        prov = json.loads(PROV_PATH.read_text())
    prov[str(dest.relative_to(DATA_RAW))] = {
        "url": url, "sha256": hashlib.sha256(blob).hexdigest(),
        "bytes": len(blob)}
    PROV_PATH.write_text(json.dumps(prov, indent=1))
    return blob


# ------------------------------------------------------------------
# resident ledger x NIMA object list
# ------------------------------------------------------------------

html = fetch(NIMA_BASE + "nima.php", NIMA_DIR / "nima.php.html").decode()
m = re.search(r"var Mtable = (\[.*?\]);", html, re.S)
nima_tab = json.loads(m.group(1))
nima_names = {row[1] for row in nima_tab}
logger.info(f"NIMA tracked objects: {len(nima_names)}")

ledger = list(csv.DictReader(
    open(RESULTS / "step_b77_resident_ledger.csv")))


def key(n):
    return n.replace(" ", "").replace("(", "").replace(")", "")


nima_keys = {key(n) for n in nima_names}
# Sedna's NIMA directory uses the unnumbered name; TG387 carries O-C
# files although absent from the Mtable row block -- the directory
# listing is authoritative, so membership is probed per object.
candidates = {}
for r in ledger:
    k = key(r["name"])
    if k not in candidates:
        candidates[k] = {"name": r["name"], "in_cap": r["in_cap"] == "1",
                         "d_axis": float(r["d_axis"]),
                         "a": float(r["a"]), "q": float(r["q"]),
                         "in_mtable": k in nima_keys}

# ------------------------------------------------------------------
# per-object O-C files: probe directory listing, take newest version
# ------------------------------------------------------------------

# negative cache: directories that 404'd once stay absent -- the
# ledger probes ~70 objects, most of which NIMA does not track, and
# re-fetching every dead directory on each run is wasted traffic
MISS_PATH = NIMA_DIR / "missing_dirs.json"
missing = set(json.loads(MISS_PATH.read_text())
              if MISS_PATH.exists() else [])

objects = {}
for k, meta in candidates.items():
    if k in missing:
        meta["nima_versions"] = []
        meta["nima_dir"] = False
        continue
    blob_idx = fetch(NIMA_DATA + k + "/", NIMA_DIR / k / "index.html",
                     required=False)
    if blob_idx is None:
        missing.add(k)
        meta["nima_versions"] = []
        meta["nima_dir"] = False
        continue
    idx = blob_idx.decode()
    vers = sorted(set(re.findall(r"omc_ast_NIMAv(\d+)\.txt", idx)),
                  key=int)
    if not vers:
        meta["nima_versions"] = []
        continue
    v = vers[-1]
    fn = f"omc_ast_NIMAv{v}.txt"
    blob = fetch(NIMA_DATA + k + "/" + fn, NIMA_DIR / k / fn)
    meta["nima_versions"] = vers
    meta["omc_file"] = fn
    objects[k] = (meta, blob.decode().splitlines())
    logger.info(f"{meta['name']:>12s}  in_cap={meta['in_cap']}  "
                f"versions={vers}  -> {fn}")

MISS_PATH.write_text(json.dumps(sorted(missing), indent=1))

n_in = sum(1 for m, _ in objects.values() if m["in_cap"])
logger.info(f"objects with O-C files: {len(objects)} "
            f"({n_in} in-cap, {len(objects) - n_in} control)")

# ------------------------------------------------------------------
# parse the O-C records
# ------------------------------------------------------------------
# row:  O  yyyy  mm  dd.frac  ra_deg  dec_deg  obs  bias_ra  bias_dec
#       w_ra  w_dec  n  catflag  flg2  i1  i2  oc_ra  oc_dec  chi

records = []
for k, (meta, lines) in objects.items():
    for ln in lines:
        f = ln.split()
        if len(f) < 18 or f[0] != "O":
            continue
        # 19-field rows carry a catalogue flag at [12]; 18-field rows
        # omit it, shifting the trailing O-C block left by one
        sh = 0 if len(f) >= 19 else 1
        try:
            yr = int(f[1]); mo = int(f[2]); dy = float(f[3])
            ra, dec = float(f[4]), float(f[5])
            bias_ra, bias_dec = float(f[7]), float(f[8])
            w_ra, w_dec = float(f[9]), float(f[10])
            oc_ra, oc_dec = float(f[16 - sh]), float(f[17 - sh])
            chi = float(f[18 - sh])
        except ValueError:
            continue
        records.append({
            "object": meta["name"], "key": k,
            "in_cap": meta["in_cap"], "d_axis_deg": meta["d_axis"],
            "a_au": meta["a"], "q_au": meta["q"],
            "yr": yr + (mo - 1 + dy / 31.0) / 12.0,
            "ra_deg": ra, "dec_deg": dec, "obs_code": f[6],
            # bias cols carry NIMA's catalogue correction (mas); the
            # w cols are inverse-variance weights, not mas
            "bias_ra_mas": bias_ra, "bias_dec_mas": bias_dec,
            "w_ra": w_ra, "w_dec": w_dec,
            "cat_flag": f[12] if sh == 0 else "",
            "flag_b": f[13 - sh],
            "src_i1": f[14 - sh], "src_i2": f[15 - sh],
            "oc_ra": oc_ra, "oc_dec": oc_dec, "chi": chi,
            "nima_file": meta["omc_file"]})

logger.info(f"parsed {len(records)} O-C records")

# ------------------------------------------------------------------
# T1/T2: coverage + per-object level & drift
# ------------------------------------------------------------------

per_obj = []
for k, (meta, _) in objects.items():
    rr = [r for r in records if r["key"] == k]
    if not rr:
        continue
    t = np.array([r["yr"] for r in rr])
    ocra = np.array([r["oc_ra"] for r in rr])
    ocdec = np.array([r["oc_dec"] for r in rr])
    amp = np.hypot(ocra, ocdec)
    # linear drift of each component vs epoch, with OLS sigma
    def drift(y):
        X = np.column_stack([np.ones(len(t)), t - t.mean()])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        dof = max(len(t) - 2, 1)
        s2 = float(np.sum((y - X @ coef) ** 2) / dof)
        se = math.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
        return float(coef[1]), float(se)
    s_ra, se_ra = drift(ocra)
    s_de, se_de = drift(ocdec)
    per_obj.append({
        "object": meta["name"], "in_cap": meta["in_cap"],
        "d_axis_deg": meta["d_axis"], "a_au": meta["a"],
        "q_au": meta["q"], "nima_file": meta["omc_file"],
        "n_obs": len(rr),
        "yr0": round(float(t.min()), 2), "yr1": round(float(t.max()), 2),
        "med_amp": float(np.median(amp)),
        "rms_amp": float(np.sqrt(np.mean(amp ** 2))),
        "med_chi": float(np.median([r["chi"] for r in rr])),
        "slope_ra_per_yr": s_ra, "slope_ra_se": se_ra,
        "slope_dec_per_yr": s_de, "slope_dec_se": se_de,
        "n_src_flags": len({(r["cat_flag"], r["flag_b"],
                            r["src_i1"]) for r in rr}),
        "obs_codes": sorted({r["obs_code"] for r in rr})})

res = {
    "step": "step_162_nima_pilot",
    "description": ("NIMA O-C pilot on detached-resident objects: "
                    "per-observation residuals from the Lucky Star "
                    "ephemeris files joined to the step_b77 cap "
                    "ledger"),
    "T1_coverage": {
        "n_ledger": len({r['name'] for r in ledger}),
        "n_in_cap_ledger": len({r['name'] for r in ledger
                                if r['in_cap'] == '1'}),
        "n_with_oc": len(per_obj),
        "n_in_cap_oc": sum(1 for p in per_obj if p["in_cap"]),
        "n_control_oc": sum(1 for p in per_obj if not p["in_cap"]),
        "n_records": len(records)},
    "T2_per_object": per_obj}

# T3: descriptive in/out contrast (underpowered by construction)
inn = [p for p in per_obj if p["in_cap"]]
outc = [p for p in per_obj if not p["in_cap"]]
t3 = {"n_in": len(inn), "n_out": len(outc)}
if inn:
    t3["in_med_amp"] = float(np.median([p["med_amp"] for p in inn]))
    t3["in_med_slope_abs"] = float(np.median(
        [abs(p["slope_ra_per_yr"]) + abs(p["slope_dec_per_yr"])
         for p in inn]))
if outc:
    t3["out_med_amp"] = float(np.median([p["med_amp"] for p in outc]))
    t3["out_med_slope_abs"] = float(np.median(
        [abs(p["slope_ra_per_yr"]) + abs(p["slope_dec_per_yr"])
         for p in outc]))
# pooled significance-weighted drift across in-cap objects
sig = [(abs(p["slope_ra_per_yr"]) / max(p["slope_ra_se"], 1e-9),
        abs(p["slope_dec_per_yr"]) / max(p["slope_dec_se"], 1e-9))
       for p in inn]
t3["in_slope_sig_pairs"] = [[round(a, 2), round(b, 2)]
                            for a, b in sig]
res["T3_contrast"] = t3

res["caveats"] = [
    "NIMA residuals are computed against NIMA's own orbit solution, "
    "which absorbed the same astrometry -- this is a "
    "post-fit-residual channel, not an absolute ephemeris error.",
    "NIMA astrometry descends largely from MPC data plus dedicated "
    "runs, so independence from the comet lineages is partial.",
    "The in-cap detached-resident intersection is ~7 objects with a "
    "single out-of-cap resident control; the channel is registered "
    "as a pilot and cannot support a detection claim at this N.",
    "No occultation-anchor rows are distinguished in the omc_ast "
    "format; occultation ephemeris ties enter through the separate "
    "SOSB results database and are not folded in here."]

with open(RESULTS / "step_b128_nima_pilot.json", "w") as f:
    json.dump(res, f, indent=1)
logger.data_save(RESULTS / "step_b128_nima_pilot.json")

keys = sorted({k for r in records for k in r})
with open(RESULTS / "step_b128_nima_pilot.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in records:
        w.writerow({k: r.get(k) for k in keys})
logger.data_save(RESULTS / "step_b128_nima_pilot.csv")

logger.info(json.dumps(res, indent=1))
logger.info("done")
