"""step_17: the planetary baseline -- the test the whole analysis needed.

Every earlier test compared raw d_of against the axis.  The correct
TEP observable is the residual AFTER the standard planetary baseline.
This step builds that baseline:

  For each comet: node heliocentric distances r_asc, r_desc (analytic
  from q, e, w) and Jupiter's ecliptic longitude at the perihelion
  epoch (mean longitude 34.35 deg at J2000, rate 30.349 deg/yr --
  an analytic approximation adequate for a covariance baseline).
  Encounter proxy jproxy = min over nodes of the combined
  node-distance/node-timing separation (a comet is perturbed
  strongly only when a node sits near 5.2 AU AND Jupiter is there).

  Baseline regression: log(d_of) ~ 1 + jproxy + q over the full
  sample -> residuals.  Then the axis tests run on the residuals:
  cap contrast and the continuous Spearman test, no cap, no cuts.

The baseline works: d_of correlates with jproxy globally
(rho = -0.47).  The in-cap comets carry NO geometric advantage
(jproxy, node_min, i, q, |daa| all flat vs out-cap).  And the
residual axis preference STRENGTHENS: p = 0.0018 (matched),
p = 0.0036 cap / rho = -0.28, p = 0.0012 continuous (all class-1).
Standard planetary geometry is now explicitly subtracted; what
remains is the non-planetary residual, and it still points at
the TNO axis.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_040_planetary_baseline")
tee_stdout(logger)
logger.header("Planetary baseline residual")
import json, math, os, re
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
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

def parse_code(path):
    p = TP(); p.feed(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for r in p.rows:
        if len(r) < 14: continue
        try:
            out[r[0].strip()] = dict(desig=r[0].strip(),
                cls=re.sub(r"^\d", "", r[3].strip()),
                tyr=int(r[7].split()[0]) if r[7].split() else 0,
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError, IndexError):
            continue
    return out

orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
fut  = parse_code(str(DATA_RAW / "code" / "code_future.html"))
warsaw = {l[5:17].strip() for l in open(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat")) if len(l) > 115}

def jupiter_proxy(ro):
    e, q, w = ro["e"], ro["q"], math.radians(ro["w"])
    p_par = q * (1 + e)
    r_asc = p_par / (1 + e * math.cos(-w)); r_desc = p_par / (1 + e * math.cos(math.pi - w))
    lamJ = (34.35 + 30.349 * (ro["tyr"] - 2000)) % 360
    dasc = abs((lamJ - ro["Om"] + 180) % 360 - 180)
    ddesc = abs((lamJ - (ro["Om"] + 180) % 360 + 180) % 360 - 180)
    return min(math.sqrt((r_asc - 5.2) ** 2 + (dasc / 36.) ** 2),
               math.sqrt((r_desc - 5.2) ** 2 + (ddesc / 36.) ** 2))

rows = []
for k, ro in orig.items():
    if k not in fut: continue
    if not (0 < ro["aa"] < 100 and ro["cls"] in ("1a", "1a+", "1b")) or k in warsaw: continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    rows.append(dict(desig=k, cls=ro["cls"], q=ro["q"], i=ro["i"],
        theta=sep(-po, TNO), d_of=sep(-po, -pf), jproxy=jupiter_proxy(ro)))

res = {"n_all_c1": len(rows)}

def run(sub, label):
    th = np.array([r["theta"] for r in sub]); v = np.array([r["d_of"] for r in sub])
    J = np.array([r["jproxy"] for r in sub]); Q = np.array([r["q"] for r in sub])
    inc = th < 60
    out = {"n": len(sub), "n_in": int(inc.sum())}
    # geometry distributions flat?
    for key, x in [("jproxy", J), ("q", Q), ("i", np.array([r["i"] for r in sub]))]:
        u = mannwhitneyu(x[inc], x[~inc], alternative="two-sided")
        out[f"{key}_p_2s"] = float(u.pvalue)
    rho, p = spearmanr(J, v)
    out["spearman_dof_jproxy"] = {"rho": float(rho), "p": float(p)}
    # residual after baseline
    y = np.log(v); X = np.column_stack([np.ones(len(y)), J, Q])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    u = mannwhitneyu(resid[inc], resid[~inc], alternative="greater")
    rho2, p2 = spearmanr(th, resid)
    out["residual_cap"] = {"med_in": float(np.median(resid[inc])),
                           "med_out": float(np.median(resid[~inc])), "p": float(u.pvalue)}
    out["residual_continuous"] = {"rho": float(rho2), "p_2sided": float(p2)}
    return out

res["matched"] = run([r for r in rows if r["q"] < 3.1], "matched")
res["all_c1"]  = run(rows, "all_c1")

out = str(RESULTS / "step_17_planetary_baseline.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
print("wrote", out)
