"""step_13: axis convergence and influence audit on the CODE-only sample.

Does the independent comet sample recover the TNO axis on its own?
A free sky scan finds the direction of maximum cap-contrast in d_of;
if the anomaly is a real spatial structure, the recovered axis should
land on the pre-declared TNO direction -- cross-population coherence
measured without ever fitting to it.  Leave-one-out checks that no
single comet drives the result.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_036_axis_convergence")
tee_stdout(logger)
logger.header("Independent axis recovery")
import json, math, os, re
import numpy as np
from scipy.stats import mannwhitneyu
from html.parser import HTMLParser


def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def lv(l, b):
    l, b = math.radians(l), math.radians(b)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

TNO = lv(34, -13)

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

def parse(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            out[r[0].strip()] = dict(desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

orig = parse(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip() for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")) if len(l) > 115}

co = []
for k, ro in orig.items():
    if k not in fut: continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    if 0 < ro["aa"] < 100 and ro["q"] < 3.1 and ro["cls"] in ("1a", "1a+", "1b") and k not in warsaw:
        co.append(dict(desig=k, aph=-po, d_of=sep(-po, -pf)))

# free axis scan
grid = [(p, lam, b) for lam in np.arange(0, 360, 10) for b in np.arange(-60, 61, 10)
        for p in [mannwhitneyu(
            np.array([r["d_of"] for r in co])[np.array([sep(r["aph"], lv(lam, b)) for r in co]) < 60],
            np.array([r["d_of"] for r in co])[np.array([sep(r["aph"], lv(lam, b)) for r in co]) >= 60],
            alternative="greater").pvalue]
        if (np.array([sep(r["aph"], lv(lam, b)) for r in co]) < 60).sum() >= 4
        and (np.array([sep(r["aph"], lv(lam, b)) for r in co]) < 60).sum() <= len(co) - 4]
grid.sort()
bp, blam, bb = grid[0]
recovered = {"lam": float(blam), "beta": float(bb), "p": float(bp),
             "sep_from_tno": float(sep(lv(blam, bb), TNO)),
             "top10": [{"lam": float(l), "beta": float(b), "p": float(p)} for p, l, b in grid[:10]]}

# LOO
loo = []
for i, r in enumerate(co):
    sub = [x for j, x in enumerate(co) if j != i]
    th = np.array([sep(x["aph"], TNO) for x in sub]); v = np.array([x["d_of"] for x in sub])
    loo.append({"dropped": r["desig"],
        "p": float(mannwhitneyu(v[th < 60], v[th >= 60], alternative="greater").pvalue)})
loo_p = sorted(x["p"] for x in loo)

res = {"n_code_only": len(co), "recovered_axis": recovered,
       "loo": {"min": loo_p[0], "max": loo_p[-1],
               "weakest": sorted(loo, key=lambda x: -x["p"])[:5]}}
out = str(RESULTS / "step_13_axis_convergence.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)