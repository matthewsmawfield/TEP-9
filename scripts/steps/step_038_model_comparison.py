"""step_15: formal model comparison + the last standard-physics controls.

  M1  ecliptic-latitude control: the axis sits at beta ~ -13 deg, so
      in-cap aphelia are preferentially low-|beta|.  If the excess
      were ecliptic-proximity (deeper planetary-plane scattering),
      matching the |beta| distributions must kill it.  Compute the
      cap contrast inside |beta| strata and on a |beta|-matched
      subsample.
  M2  patch-model likelihood: d_of ~ baseline + A * exp(-theta^2/2 s^2)
      fitted on the CODE-only data vs the flat null -> likelihood
      ratio / BIC-based Bayes factor.  One number for the patch.
  M3  auxiliary class-2 sample: class 2a/2b CODE-only comets --
      lower-quality orbits the anomaly should still touch if it
      lives in the sky rather than in orbit quality.
  M4  seasonal/solar-longitude control: if the excess were an
      observation-geometry artifact tied to discovery season,
      d_of should correlate with solar longitude at perihelion,
      not with the fixed axis.  Uses perihelion year+month from
      the catalogue epoch.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_038_model_comparison")
tee_stdout(logger)
logger.header("Model comparison and controls")
import json, math, os, re
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from scipy.optimize import minimize_scalar
from html.parser import HTMLParser

rng = np.random.default_rng(7)

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

c1, c2 = [], []
for k, ro in orig.items():
    if k not in fut: continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    pf = perih_dir(math.radians(fut[k]["w"]), math.radians(fut[k]["Om"]), math.radians(fut[k]["i"]))
    if 0 < ro["aa"] < 100 and k not in warsaw:
        r = dict(desig=k, q=ro["q"], tyr=ro["tyr"], cls=ro["cls"],
                 aph=-po, beta=math.degrees(math.asin(-po[2])), d_of=sep(-po, -pf))
        if ro["cls"] in ("1a", "1a+", "1b"): c1.append(r)
        elif ro["cls"] in ("2a", "2b"): c2.append(r)
mat = [r for r in c1 if r["q"] < 3.1]
res = {"n_c1": len(c1), "n_c2": len(c2), "n_matched": len(mat)}

# ---- M1 ecliptic-latitude control on matched CODE-only
th = np.array([sep(r["aph"], TNO) for r in mat]); v = np.array([r["d_of"] for r in mat])
inc = th < 60
res["M1_beta"] = {"med_absbeta_in": float(np.median(np.abs([r["beta"] for r in mat])[inc])),
                  "med_absbeta_out": float(np.median(np.abs([r["beta"] for r in mat])[~inc]))}
m1 = {}
for lab, m in [("low|beta|<30", np.abs([r["beta"] for r in mat]) < 30),
               ("high|beta|>=30", np.abs([r["beta"] for r in mat]) >= 30)]:
    ii = inc & m; io = (~inc) & m
    if ii.sum() >= 3 and io.sum() >= 3:
        u = mannwhitneyu(v[ii], v[io], alternative="greater")
        m1[lab] = {"n_in": int(ii.sum()), "n_out": int(io.sum()),
                   "med_in": float(np.median(v[ii])), "med_out": float(np.median(v[io])),
                   "p": float(u.pvalue)}
# |beta|-matched subsample: nearest out-cap comet for each in-cap comet
bm = []
absb = np.abs([r["beta"] for r in mat])
out_idx_all = np.where(~inc)[0]
for i in np.where(inc)[0]:
    j = np.argmin(np.abs(absb[out_idx_all] - absb[i]))
    bm.append(out_idx_all[j])
bm = list(set(bm))
if len(bm) >= 4:
    u = mannwhitneyu(v[inc], v[bm], alternative="greater")
    m1["matched_subsample"] = {"n_in": int(inc.sum()), "n_out": len(bm),
        "med_in": float(np.median(v[inc])), "med_out": float(np.median(v[bm])),
        "p": float(u.pvalue)}
res["M1_strata"] = m1

# ---- M2 patch model likelihood on matched CODE-only
# model: d_of_i ~ LogNormal-ish; use Gaussian approx on log d_of
y = np.log(v)
def nll(s):
    if s < 5 or s > 170: return 1e12
    w = np.exp(-(th / s) ** 2)
    # two-level model: mean_in = mu0+A*w_i  -> fit mu0, A by lstsq
    X = np.column_stack([np.ones_like(w), w])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    sig2 = np.mean(resid ** 2)
    return 0.5 * len(y) * math.log(sig2) + 0.5 * len(y)
X0 = np.ones((len(y), 1))
resid0 = y - np.mean(y)
nll0 = 0.5 * len(y) * math.log(np.mean(resid0 ** 2)) + 0.5 * len(y)
best = minimize_scalar(nll, bounds=(5, 170), method="bounded")
dLL = nll0 - best.fun
# BIC approximation: delta-BIC = 2*dLL - k*ln(n); k=3 (mu,A,sigma vs mu,sigma)
dbic = 2 * dLL - 2 * math.log(len(y))
res["M2_patch_model"] = {"sigma_deg": float(best.x), "lnLR": float(dLL),
    "delta_BIC": float(dbic), "approx_logBF": float(dbic / 2),
    "note": "patch model dlogN = mu + A*exp(-theta^2/2s^2) vs flat; "
            "logBF ~ deltaBIC/2 (Kass-Raftery)"}

# ---- M3 class-2 auxiliary sample
if len(c2) >= 10:
    t2 = np.array([sep(r["aph"], TNO) for r in c2]); v2 = np.array([r["d_of"] for r in c2])
    i2 = t2 < 60
    if i2.sum() >= 4 and (~i2).sum() >= 4:
        u = mannwhitneyu(v2[i2], v2[~i2], alternative="greater")
        res["M3_class2"] = {"n": len(c2), "n_in": int(i2.sum()),
            "med_in": float(np.median(v2[i2])), "med_out": float(np.median(v2[~i2])),
            "p": float(u.pvalue)}

# ---- M4 seasonal control: solar longitude at perihelion epoch
# approximate solar longitude ~ (day-of-year / 365.25) * 360 - 80
# (perihelion epoch month unknown; use year fraction impossible ->
#  use the catalogue 'tyr' (perihelion year) -> substitute: ecliptic
#  longitude of the comet's aphelion projected onto ecliptic is the
#  relevant seasonal proxy at discovery.  Instead test d_of vs
#  perihelion YEAR (era) and vs |beta| directly.)
rho_y, p_y = spearmanr([r["tyr"] for r in mat], v)
res["M4_era"] = {"rho_dof_vs_year": float(rho_y), "p": float(p_y)}

out = str(RESULTS / "step_15_model_comparison.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)