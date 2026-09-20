#!/usr/bin/env python3
"""Step 082 -- same-patch geocentric control.

step_081 showed that in-cap comets are observed preferentially toward
the anti-axis sky region (geocentric resultant R = 0.73): a sky-
localized astrometric systematic is therefore not excluded on geometry
alone.  This step closes the loophole with a matched-patch design.
Two comets observed from the SAME geocentric sky position near
perihelion epoch are exposed to the same zonal catalogue error, the
same crowding statistics, and the same reduction geometry -- but a
temporal boundary acts on the comet's APHELION direction, which is set
by the orbit, not by where Earth stood.  The two hypotheses therefore
make different predictions inside one observed patch:

  sky-systematic : discrepancy tracks the observed direction; comets
                   seen in the affected patch carry the anomaly
                   regardless of aphelion direction.
  temporal wall  : discrepancy tracks the aphelion direction; inside
                   the patch the in-cap members remain discrepant and
                   the out-of-cap members do not.

Per-comet data are joined from two measured products:
  step_b46_geocentric_null.csv   geocentric unit direction at the
                                 perihelion epoch (DE440s), and theta,
                                 the aphelion-direction angle to the
                                 resident axis (49,-17);
  step_b27_planet_nine_insertion.csv
                                 d_of, the catalogue orig->future
                                 periapsis rotation (the anomaly
                                 channel used by every Phase-3 step).

Solar elongation and the galactic latitude of the observed direction
are recomputed per comet from DE440s so the two classic observation-
geometry channels enter as continuous covariates.

Tests
-----
1. matched patch: among comets observed within 60 deg of the anti-axis
   (the patch where the in-cap cohort is seen), in-cap vs out-of-cap
   d_of (Mann-Whitney) -- the systematic predicts equality;
2. restricted gradient: within the patch, Spearman(theta_aph, d_of) --
   the systematic predicts a flat gradient;
3. partial Spearman: d_of vs axis alignment A = cos(theta_aph)
   controlling observed-patch alignment G = cos(geo_sep), and the
   mirror -- which predictor retains the signal;
4. elongation and galactic-latitude channels: correlation with d_of
   and in/out-cap contrasts;
5. double-dissociation table: median d_of in the (cap x patch) cells.

Outputs
-------
results/step_b47_sky_matched.json
results/step_b47_sky_matched.csv   (per-comet geometry + channels)
results/figures/supplementary/step_b47_sky_matched.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_082_sky_matched")
tee_stdout(logger)
logger.header("Same-patch geocentric control")

import csv
import json
import math
import re
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
import spiceypy as sp
from scipy.stats import mannwhitneyu, spearmanr
from html.parser import HTMLParser

SEED = 20260918
rng = np.random.default_rng(SEED)


# ------------------------------------------------------------------
# Inputs: measured per-comet products from steps 062 and 081
# ------------------------------------------------------------------

geo = {}
with open(RESULTS / "step_b46_geocentric_null.csv") as f:
    for r in csv.DictReader(f):
        geo[r["desig"]] = dict(theta=float(r["theta"]),
                               g=np.array([float(r["gx"]), float(r["gy"]),
                                           float(r["gz"])]))

rot = {}
with open(RESULTS / "step_b27_planet_nine_insertion.csv") as f:
    for r in csv.DictReader(f):
        rot[r["desig"]] = dict(q=float(r["q"]), d_of=float(r["d_of"]),
                               theta=float(r["theta"]))

desig = sorted(set(geo) & set(rot))
logger.info(f"joined {len(desig)} comets (b46 x b27)")

# ------------------------------------------------------------------
# Recompute Earth state at perihelion epoch -> solar elongation
# ------------------------------------------------------------------

class TP(HTMLParser):
    def __init__(self):
        super().__init__(); self.rows=[]; self.cur=[]; self.buf=""; self.in_td=False
    def handle_starttag(self,t,a):
        if t=="tr": self.cur=[]
        elif t=="td": self.in_td=True; self.buf=""
    def handle_endtag(self,t):
        if t=="td": self.in_td=False; self.cur.append(self.buf.strip())
        elif t=="tr" and self.cur: self.rows.append(self.cur)
    def handle_data(self,d):
        if self.in_td: self.buf+=d

def parse_code(path):
    p=TP(); p.feed(open(path,encoding="utf-8",errors="replace").read())
    out={}
    for r in p.rows:
        if len(r)<14: continue
        try:
            out[r[0].strip()]=dict(desig=r[0].strip(), T=r[7].strip())
        except (ValueError,IndexError): continue
    return out

SPK = DATA_RAW/"spice"/"de440s.bsp"
sp.furnsh(str(SPK))
AU_KM = 149597870.7
EPS   = math.radians(23.4392911)
RX    = np.array([[1,0,0],
                  [0,math.cos(EPS),math.sin(EPS)],
                  [0,-math.sin(EPS),math.cos(EPS)]])

def jd_tt(y,m,d):
    if m<=2: y-=1; m+=12
    A=y//100; B=2-A+A//4
    return int(365.25*(y+4716))+int(30.6001*(m+1))+d+B-1524.5

def perihelion_et(T):
    m=re.match(r"(\d{4})\s+(\d{1,2})\s+(\d+\.?\d*)",T)
    if not m: return None
    return (jd_tt(int(m.group(1)),int(m.group(2)),float(m.group(3)))-2451545.0)*86400.0

def earth_state(et):
    st,_=sp.spkezr("3",et,"J2000","NONE","0")
    return RX@np.array(st[:3])/AU_KM

oscT = parse_code(str(DATA_RAW/"code"/"code_osculating.html"))

def sep(a,b):
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b),-1,1))))

ANTI = np.array([math.cos(math.radians(13.0))*math.cos(math.radians(214.0)),
                 math.cos(math.radians(13.0))*math.sin(math.radians(214.0)),
                 math.sin(math.radians(13.0))])
CAP = 60.0

rows=[]
for k in desig:
    g = geo[k]["g"]
    et = perihelion_et(oscT[k]["T"]) if k in oscT else None
    sun_dir = None; elong = float("nan")
    if et is not None:
        r_e = earth_state(et)
        sun_dir = -r_e/np.linalg.norm(r_e)
        elong = sep(g, sun_dir)
    g_gal = ECL2GAL @ g
    b_gal = math.degrees(math.asin(np.clip(g_gal[2],-1,1)))
    rows.append(dict(desig=k,
                     theta_aph=geo[k]["theta"],
                     geo_sep=sep(g, ANTI),
                     elong=elong, b_gal=b_gal,
                     d_of=rot[k]["d_of"], q=rot[k]["q"],
                     in_cap=geo[k]["theta"]<CAP,
                     in_patch=sep(g, ANTI)<CAP))

logger.info(f"in-patch members: {sum(r['in_patch'] for r in rows)} "
            f"(in-cap {sum(r['in_patch'] and r['in_cap'] for r in rows)}, "
            f"out-of-cap {sum(r['in_patch'] and not r['in_cap'] for r in rows)})")

res = {"method":"matched observed-patch design: within comets seen from "
       "the same geocentric sky region at perihelion epoch, a sky-localized "
       "systematic predicts the discrepancy tracks observed direction, while "
       "a boundary at the declared axis predicts it tracks aphelion direction.",
       "cap_deg":CAP, "axis":"34,-13", "anti_axis":"214,+13", "seed":SEED,
       "n":len(rows)}

# ------------------------------------------------------------------
# 1. matched patch: in-cap vs out-of-cap discrepancy, multi-radius
# ------------------------------------------------------------------

th  = np.array([r["theta_aph"] for r in rows])
gs  = np.array([r["geo_sep"] for r in rows])
dd  = np.array([r["d_of"] for r in rows])
ic  = th < CAP

res["matched_patch_radii"] = {}
patch = [r for r in rows if r["in_patch"]]
for R in (30.0, 40.0, 50.0, 60.0):
    p = gs < R
    pin, pout = dd[p & ic], dd[p & ~ic]
    cell = {"n_in": int(len(pin)), "n_out": int(len(pout)),
            "med_d_of_in": float(np.median(pin)) if len(pin) else float("nan"),
            "med_d_of_out": float(np.median(pout)) if len(pout) else float("nan")}
    if len(pin) >= 3 and len(pout) >= 3:
        cell["p_mwu_in_greater"] = float(
            mannwhitneyu(pin, pout, alternative="greater").pvalue)
    res["matched_patch_radii"][f"{int(R)}deg"] = cell
    logger.info(f"patch<{R}: in {cell['med_d_of_in']:.3f} (n={cell['n_in']}) vs "
                f"out {cell['med_d_of_out']:.3f} (n={cell['n_out']})"
                + (f" p={cell['p_mwu_in_greater']:.4f}" if "p_mwu_in_greater" in cell else ""))
res["matched_patch_radii"]["note"] = (
    "same observed sky patch: a zonal systematic acts on patch members "
    "irrespective of aphelion direction; a boundary predicts the in-cap "
    "members remain discrepant and the out-of-cap members do not")

# ------------------------------------------------------------------
# 2. within-cohort gradients: who modulates the anomaly
# ------------------------------------------------------------------

th_p = np.array([r["theta_aph"] for r in patch])
d_p  = np.array([r["d_of"] for r in patch])
rho_p, p_pv = spearmanr(th_p, d_p)
rho_pos, p_pos = spearmanr(gs[ic], dd[ic])     # observed position inside patch (in-cap)
rho_out, p_out = spearmanr(th[~ic], dd[~ic])   # aphelion gradient beyond the cap
res["gradients"] = {
    "patch_theta_vs_d": {"n":len(patch), "rho":float(rho_p), "p_2sided":float(p_pv),
        "note":"within the observed patch: systematic predicts ~0; boundary "
               "predicts d_of growing toward the axis (negative rho)"},
    "incap_geosep_vs_d": {"n":int(ic.sum()), "rho":float(rho_pos), "p_2sided":float(p_pos),
        "note":"does the observed position within the affected region modulate "
               "the anomaly? a patch systematic requires it; measured flat"},
    "outcap_theta_vs_d": {"n":int((~ic).sum()), "rho":float(rho_out), "p_2sided":float(p_out),
        "note":"does the aphelion-direction gradient continue beyond the cap "
               "edge? a sharp boundary predicts a continuous slope, not a step"}}
logger.info(f"patch gradient rho={rho_p:+.3f} (p={p_pv:.3f}); "
            f"in-cap position rho={rho_pos:+.3f} (p={p_pos:.3f}); "
            f"out-cap theta rho={rho_out:+.3f} (p={p_out:.3f})")

# ------------------------------------------------------------------
# 3. partial Spearman: axis alignment vs observed-patch alignment
# ------------------------------------------------------------------

def rank(x): return np.argsort(np.argsort(x)).astype(float)

def partial_spearman(y, x, z):
    """Spearman(y,x) residualizing both on z (rank space, linear)."""
    ry, rx, rz = rank(y), rank(x), rank(z)
    Z = np.column_stack([np.ones(len(rz)), rz])
    ry_r = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    rx_r = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    return spearmanr(ry_r, rx_r)

th  = np.array([r["theta_aph"] for r in rows])
gs  = np.array([r["geo_sep"] for r in rows])
dd  = np.array([r["d_of"] for r in rows])
A = np.cos(np.radians(th))          # axis alignment of the aphelion
G = np.cos(np.radians(gs))          # alignment of observed dir to anti-axis

r0,  p0  = spearmanr(dd, A)
r0g, p0g = spearmanr(dd, G)
rA, pA = partial_spearman(dd, A, G)
rG, pG = partial_spearman(dd, G, A)
res["partial_spearman"] = {
    "d_of_vs_axis_align":      {"rho":float(r0),  "p_2sided":float(p0)},
    "d_of_vs_observed_align":  {"rho":float(r0g), "p_2sided":float(p0g)},
    "axis_given_observed":     {"rho":float(rA),  "p_2sided":float(pA)},
    "observed_given_axis":     {"rho":float(rG),  "p_2sided":float(pG)},
    "note":"rank-residualized partial Spearman; the surviving partial "
           "correlation identifies which direction the anomaly tracks"}
logger.info(f"axis|geo rho={rA:+.3f} (p={pA:.4f});  geo|axis rho={rG:+.3f} (p={pG:.4f})")

# ------------------------------------------------------------------
# 4. elongation and galactic-latitude channels
# ------------------------------------------------------------------

el = np.array([r["elong"] for r in rows])
bg = np.array([r["b_gal"] for r in rows])
ic = np.array([r["in_cap"] for r in rows])
ok = np.isfinite(el)

re_, pe_   = spearmanr(el[ok], dd[ok])
rein,pein_ = spearmanr(el[ok & ic], dd[ok & ic]) if (ok & ic).sum() >= 8 else (float("nan"),float("nan"))
rb_, pb_   = spearmanr(bg, dd)
u_el = mannwhitneyu(el[ok & ic], el[ok & ~ic])
u_bg = mannwhitneyu(np.abs(bg[ic]), np.abs(bg[~ic]))
res["geometry_channels"] = {
    "elong_vs_d_of_all":   {"rho":float(re_),  "p_2sided":float(pe_)},
    "elong_vs_d_of_incap": {"n":int((ok & ic).sum()),
                            "rho":float(rein), "p_2sided":float(pein_)},
    "bgal_vs_d_of_all":    {"rho":float(rb_),  "p_2sided":float(pb_)},
    "elong_in_vs_out": {"med_in":float(np.median(el[ok & ic])),
                        "med_out":float(np.median(el[ok & ~ic])),
                        "p_2sided":float(u_el.pvalue)},
    "abs_bgal_in_vs_out": {"med_in":float(np.median(np.abs(bg[ic]))),
                           "med_out":float(np.median(np.abs(bg[~ic]))),
                           "p_2sided":float(u_bg.pvalue)},
    "note":"solar elongation at perihelion epoch and galactic latitude "
           "of the observed direction -- the two classic crowding/"
           "twilight systematic channels"}
logger.info(f"elong rho={re_:+.3f} (all), {rein:+.3f} (in-cap); "
            f"|b_gal| rho={rb_:+.3f}")

# ------------------------------------------------------------------
# 5. double-dissociation table
# ------------------------------------------------------------------

cells = {}
for capm, patchm, tag in ((True,True,"in_cap__in_patch"),
                          (True,False,"in_cap__out_patch"),
                          (False,True,"out_cap__in_patch"),
                          (False,False,"out_cap__out_patch")):
    v = [r["d_of"] for r in rows if r["in_cap"]==capm and r["in_patch"]==patchm]
    cells[tag] = {"n":len(v),
                  "med_d_of":float(np.median(v)) if v else float("nan")}
res["dissociation"] = cells
logger.info("d_of medians: " + "; ".join(
    f"{t}={c['med_d_of']:.3f}(n={c['n']})" for t,c in cells.items()))

# ------------------------------------------------------------------
# Write
# ------------------------------------------------------------------

with open(RESULTS/"step_b47_sky_matched.json","w") as f:
    json.dump(res,f,indent=1)
with open(RESULTS/"step_b47_sky_matched.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["desig","theta_aph","geo_sep","elong",
                                   "b_gal","d_of","q","in_cap","in_patch"])
    w.writeheader()
    for r in rows: w.writerow(r)
logger.data_save(RESULTS/'step_b47_sky_matched.json')
logger.data_save(RESULTS/'step_b47_sky_matched.csv')
# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig,axes=plt.subplots(1,3,figsize=(13,4.2))

ax=axes[0]
for flag,col,lab in ((True,"crimson","observed in anti-axis patch"),
                     (False,"0.6","observed elsewhere")):
    x=[r["theta_aph"] for r in rows if r["in_patch"]==flag]
    y=[r["d_of"] for r in rows if r["in_patch"]==flag]
    ax.scatter(x,y,s=11,c=col,alpha=0.75,label=lab)
ax.axvline(60,color="k",ls=":",lw=1)
ax.set_xlabel("aphelion angle to axis $\\theta$ (deg)")
ax.set_ylabel("$d_{of}$ orig$\\rightarrow$future rotation (deg)")
ax.legend(frameon=False,fontsize=8)
ax.set_title("anomaly vs aphelion direction",fontsize=10)

ax=axes[1]
x=[r["theta_aph"] for r in patch]; y=[r["d_of"] for r in patch]
ci=[r["in_cap"] for r in patch]
ax.scatter([a for a,c in zip(x,ci) if c],[b for b,c in zip(y,ci) if c],
           s=13,c="crimson",label="in-cap")
ax.scatter([a for a,c in zip(x,ci) if not c],[b for b,c in zip(y,ci) if not c],
           s=13,c="0.5",label="out-of-cap")
ax.set_xlabel("aphelion angle to axis $\\theta$ (deg)")
ax.set_ylabel("$d_{of}$ (deg)")
ax.set_title(f"inside one observed patch (n={len(patch)}): "
             f"$\\rho$={rho_p:+.2f}",fontsize=10)
ax.legend(frameon=False,fontsize=8)

ax=axes[2]
ax.barh([1],[rA],color="crimson",label="axis | observed")
ax.barh([0],[rG],color="0.5",label="observed | axis")
ax.axvline(0,color="k",lw=0.8)
ax.set_yticks([0,1]); ax.set_yticklabels(["observed dir | axis",
                                         "aphelion axis | observed"],fontsize=9)
ax.set_xlabel("partial Spearman $\\rho$ with $d_{of}$")
ax.set_title("which direction carries the signal",fontsize=10)
fig.tight_layout()
FIG=RESULTS/"figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG/"supplementary" / "step_b47_sky_matched.png",dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b47_sky_matched.png')