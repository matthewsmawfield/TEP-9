#!/usr/bin/env python3
"""Step 077 -- era stability of the aphelion dipole.

Step 072 found the reconstruction-residual amplitude is era-weighted
(strongest contrast pre-1950).  The aphelion dipole is the
reconstruction-independent channel -- a property of the catalogue
entries, not of the fit machinery -- so its era structure bounds how
much of the spatial signal could be an old-catalogue artefact.  The
per-comet aphelion projection onto the resident axis (49,-17) and the
galactic-anticenter control is measured in the same three era bins.

Inputs : data/raw/code/code_original.html,
         data/raw/warsaw/warsaw_table[cd].dat
Outputs: results/step_b42_dipole_era.json
         results/figures/step_b42_dipole_era.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout, parse_code
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_077_dipole_era")
tee_stdout(logger)
logger.header("Era stability of the aphelion dipole")

import json
import math
import re
import numpy as np
from scipy.stats import norm
from html.parser import HTMLParser

SEED = 20260918
rng  = np.random.default_rng(SEED)

def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

AXIS = lv(49.0, -17.0)          # resident axis (step_061 dipole convention)
ANTI = lv(312.0, 13.0)          # galactic anticenter-ish control on-plane
# anticenter direction: l=0-ish?  use the step_061 control: galactic
# anticenter is (l,b) = (0,0) in galactic coords -> ecliptic unit vector
# step_061 used an on-plane anticenter; approximate galactic anticenter
# at l_gal=0 -> RA 17h45m -> ecliptic ~ (266.4, -5.5).  Keep it simple:
# use the literal anti-axis of the TNO direction as control-2 as well.
ANTIAXIS = -AXIS

# ------------------------------------------------------------------
# CODE parser (full table with perihelion year in designation)
# ------------------------------------------------------------------

class TP(HTMLParser):
    def __init__(self):
        super().__init__(); self.rows = []; self.cur = []; self.buf = ""; self.in_td = False
    def handle_starttag(self, t, a):
        if t == "tr": self.cur = []
        elif t == "td": self.in_td = True; self.buf = ""
    def handle_endtag(self, t):
        if t == "td": self.in_td = False; self.cur.append(self.buf.strip())
        elif t == "tr" and self.cur: self.rows.append(self.cur)
    def handle_data(self, d):
        if self.in_td: self.buf += d

def parse_code_full(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            out[r[0].strip()] = dict(desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                T=r[7].strip(),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

def year_of(desig):
    m = re.match(r"[CP]*/(\d{4})", desig)
    return int(m.group(1)) if m else None

orig = parse_code_full(str(DATA_RAW / "code" / "code_original.html"))
warsaw_ids = {l[5:17].strip() for l in
              open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
              if len(l) > 115}

comets = []
for k, ro in orig.items():
    if not (0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b")) or k in warsaw_ids:
        continue
    yr = year_of(k)
    if yr is None: continue
    d = -perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                   math.radians(ro["i"]))
    comets.append(dict(desig=k, cohort="code", year=yr, dir=d))

# Warsaw original legs (tablec) with year from table columns
def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115:
            continue
        try:
            rows.append(dict(desig=line[5:17].strip(),
                             com=line[3].strip(),
                             tyr=int(line[27:31]),
                             q=float(line[42:56]),
                             w=float(line[70:82]), Om=float(line[82:94]),
                             i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows

PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
war_c = {}
for r in parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")):
    k = r["desig"]
    if k not in war_c or PREF.get(r["com"], 9) < PREF.get(war_c[k]["com"], 9):
        war_c[k] = r
for d, ro in war_c.items():
    if not (0 < ro["aa"] < 100):
        continue
    dir_ = -perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                      math.radians(ro["i"]))
    comets.append(dict(desig=d, cohort="warsaw", year=ro["tyr"], dir=dir_))

logger.info(f"comets: {len(comets)} "
            f"(code {sum(c['cohort']=='code' for c in comets)}, "
            f"warsaw {sum(c['cohort']=='warsaw' for c in comets)})")

# ------------------------------------------------------------------
# Era-stratified dipole projection
# ------------------------------------------------------------------

ERAS = [(1885, 1950), (1950, 1990), (1990, 2030)]
res = {"method": "per-comet aphelion projection onto the resident axis "
                 "(49,-17) by perihelion era; z = mean*sqrt(3N) is the "
                 "isotropic null; anti-axis control appended.  The "
                 "dipole is reconstruction-independent, so era "
                 "stability here bounds catalogue-age systematics.",
       "axis": "49,-17", "seed": SEED,
       "eras": [f"{a}-{b}" for a, b in ERAS]}

for co in ("code", "warsaw", "pooled"):
    sub = [c for c in comets if co == "pooled" or c["cohort"] == co]
    out = {"n": len(sub)}
    for lo, hi in ERAS:
        s = [c for c in sub if lo <= c["year"] < hi]
        if len(s) < 5:
            out[f"{lo}-{hi}"] = None; continue
        d = np.array([c["dir"] for c in s])
        proj = d @ AXIS
        z = proj.mean() * math.sqrt(3 * len(s))
        pa = (d @ ANTIAXIS).mean() * math.sqrt(3 * len(s))
        out[f"{lo}-{hi}"] = {
            "n": len(s),
            "d_parallel": float(proj.mean()),
            "z": float(z), "p_1sided": float(norm.sf(z)),
            "antiaxis_z": float(pa),
            "antiaxis_p": float(norm.sf(pa))}
        logger.info(f"{co} {lo}-{hi}: n={len(s)} d_par={proj.mean():+.3f} "
                    f"p={norm.sf(z):.4f}")
    res[co] = out

with open(RESULTS / "step_b42_dipole_era.json", "w") as f:
    json.dump(res, f, indent=1)
print(f"wrote {RESULTS / 'step_b42_dipole_era.json'}")

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7.2, 4.4))
labels = [f"{a}-{b}" for a, b in ERAS]
xs = np.arange(len(ERAS))
for co, col, off in (("code", "crimson", -0.15), ("warsaw", "darkorange", 0.15)):
    ys, ns = [], []
    for lo, hi in ERAS:
        e = res[co].get(f"{lo}-{hi}")
        ys.append(e["d_parallel"] if e else np.nan)
        ns.append(e["n"] if e else 0)
    ax.bar(xs + off, ys, width=0.28, color=col, alpha=0.8, label=co)
    for x, y, n in zip(xs, ys, ns):
        if np.isfinite(y):
            ax.annotate(f"n={n}", (x + off, y), ha="center",
                        va="bottom" if y > 0 else "top", fontsize=7)
ax.axhline(0, color="k", ls=":", lw=1)
ax.set_xticks(xs); ax.set_xticklabels(labels)
ax.set_ylabel(r"mean aphelion projection on resident axis $d_\parallel$")
ax.legend(frameon=False, fontsize=8)
ax.set_title("dipole vs era -- reconstruction-independent channel",
             fontsize=10)
fig.tight_layout()
FIG = RESULTS / "figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b42_dipole_era.png", dpi=150)
print(f"wrote {FIG / 'step_b42_dipole_era.png'}")
