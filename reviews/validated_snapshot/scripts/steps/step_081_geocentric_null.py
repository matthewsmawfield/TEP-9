#!/usr/bin/env python3
"""Step 081 -- geocentric observation-direction null.

Any astrometric or observational systematic localized on the sky
(zonal catalogue errors, galactic-plane crowding, twilight geometry)
must act on the direction from which the comet is actually seen --
its geocentric direction near perihelion.  For Oort-spike comets that
direction is dominated by Earth's position at the perihelion epoch,
not by the comet's aphelion: in-cap membership and observed-sky
position are geometrically decoupled.  This step computes each
class-1 CODE comet's geocentric direction at perihelion epoch from
DE440s and tests whether the in-cap cohort's observed directions
cluster -- the necessary condition for a sky-position systematic to
manufacture the cap-correlated anomaly.

Inputs : data/raw/code/code_*.html, data/raw/warsaw/warsaw_tablec.dat,
         data/raw/spice/de440s.bsp
Outputs: results/step_b46_geocentric_null.json
         results/step_b46_geocentric_null.csv
         results/figures/step_b46_geocentric_null.png
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
from scripts.utils.step_logger import StepLogger

logger = StepLogger("step_081_geocentric_null")
tee_stdout(logger)
logger.header("Geocentric observation-direction null")

import csv
import json
import math
import re
import numpy as np
import spiceypy as sp
from html.parser import HTMLParser

SEED = 20260918
rng  = np.random.default_rng(SEED)

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
            out[r[0].strip()]=dict(desig=r[0].strip(),
                cls=re.sub(r"^\d","",r[3].strip()), T=r[7].strip(),
                q=float(r[8]), e=float(r[9]), w=float(r[10]),
                Om=float(r[11]), i=float(r[12]), aa=float(r[13]))
        except (ValueError,IndexError): continue
    return out

SPK = DATA_RAW/"spice"/"de440s.bsp"
if not SPK.exists():
    raise FileNotFoundError(f"JPL ephemeris missing: {SPK}")
sp.furnsh(str(SPK))

AU_KM  = 149597870.7
EPS    = math.radians(23.4392911)
RX     = np.array([[1,0,0],
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

def perih_dir(om,Om,inc):
    co,so,cO,sO,ci,si=np.cos(om),np.sin(om),np.cos(Om),np.sin(Om),np.cos(inc),np.sin(inc)
    return np.array([cO*co-sO*so*ci, sO*co+cO*so*ci, so*si])

def lv(l,b):
    l,b=math.radians(l),math.radians(b)
    return np.array([math.cos(b)*math.cos(l),math.cos(b)*math.sin(l),math.sin(b)])

def sep(a,b):
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b),-1,1))))

TNO=lv(49.0,-17.0)
CAP=60.0

osc  = parse_code(str(DATA_RAW/"code"/"code_osculating.html"))
orig = parse_code(str(DATA_RAW/"code"/"code_original.html"))
fut  = parse_code(str(DATA_RAW/"code"/"code_future.html"))
warsaw={l[5:17].strip() for l in open(str(DATA_RAW/"warsaw"/"warsaw_tablec.dat")) if len(l)>115}

rows=[]
for k,ro in orig.items():
    if k not in fut or k not in osc or k in warsaw: continue
    if not (0<ro["aa"]<100 and ro["cls"] in ("1a","1a+","1b")): continue
    et=perihelion_et(osc[k]["T"])
    if et is None: continue
    d=perih_dir(math.radians(ro["w"]),math.radians(ro["Om"]),math.radians(ro["i"]))
    r_c=ro["q"]*d                       # comet heliocentric pos at perihelion
    r_e=earth_state(et)
    g=r_c-r_e
    g/=np.linalg.norm(g)                 # geocentric direction
    rows.append(dict(desig=k,theta=sep(-d,TNO),
                     gx=g[0],gy=g[1],gz=g[2]))
logger.info(f"comets: {len(rows)}")

# ------------------------------------------------------------------
# Clustering of geocentric directions, in-cap vs out-of-cap
# ------------------------------------------------------------------

def resultant(dirs):
    return float(np.linalg.norm(np.mean(np.array(dirs),axis=0)))

res={"method":"geocentric comet direction at perihelion epoch "
       "(heliocentric q*periapsis_dir minus Earth position, DE440s); "
       "resultant length and isotropic p-value per cohort.  A sky-"
       "localized systematic requires clustered observed directions; "
       "Earth-orbit smearing scatters them.",
     "cap_deg":CAP,"axis":"49,-17","seed":SEED}

def cohort_stats(sub,tag):
    dirs=np.array([[r["gx"],r["gy"],r["gz"]] for r in sub])
    R=resultant(dirs)
    n=len(dirs)
    p=math.exp(-n*R*R)                     # Rayleigh on the sphere
    # projection onto the axis and anti-axis
    pr_axis=float(np.mean(dirs@TNO))
    pr_anti=float(np.mean(dirs@(-TNO)))
    res[tag]={"n":n,"R":R,"p_rayleigh":p,
              "mean_proj_axis":pr_axis,"mean_proj_antiaxis":pr_anti}
    logger.info(f"{tag}: n={n} R={R:.3f} p={p:.4f} "
                f"proj(axis)={pr_axis:+.3f} proj(anti)={pr_anti:+.3f}")

inn=[r for r in rows if r["theta"]<CAP]
outt=[r for r in rows if r["theta"]>=CAP]
cohort_stats(inn,"in_cap"); cohort_stats(outt,"out_cap"); cohort_stats(rows,"all")

# Monte-Carlo: is in-cap R unusual given the out-of-cap distribution?
g_all=np.array([[r["gx"],r["gy"],r["gz"]] for r in inn])
nin=len(inn)
Robs=resultant(g_all)
null=np.empty(20000)
for j in range(20000):
    idx=rng.choice(len(rows),nin,replace=False)
    null[j]=resultant(np.array([[rows[i]["gx"],rows[i]["gy"],rows[i]["gz"]] for i in idx]))
res["perm"]={"R_obs":Robs,
             "p_ge":float((int((null>=Robs).sum())+1)/(len(null)+1)),
             "null_med":float(np.median(null))}
logger.info(f"in-cap R={Robs:.3f} vs resampled null p={res['perm']['p_ge']:.4f}")

with open(RESULTS/"step_b46_geocentric_null.json","w") as f:
    json.dump(res,f,indent=1)
with open(RESULTS/"step_b46_geocentric_null.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["desig","theta","gx","gy","gz"])
    w.writeheader()
    for r in rows: w.writerow(r)
print(f"wrote {RESULTS/'step_b46_geocentric_null.json'}")
print(f"wrote {RESULTS/'step_b46_geocentric_null.csv'}")

# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig,axes=plt.subplots(1,2,figsize=(9,4.2))
ax=axes[0]
for sub,col,lab in ((inn,"crimson","in-cap"),(outt,"0.6","out-of-cap")):
    lam=[math.degrees(math.atan2(r["gy"],r["gx"]))%360 for r in sub]
    bet=[math.degrees(math.asin(r["gz"])) for r in sub]
    ax.scatter(lam,bet,s=10,c=col,alpha=0.7,label=lab)
ax.scatter([49],[-17],marker="*",s=160,c="crimson",edgecolors="k",label="axis")
ax.scatter([229],[17],marker="*",s=160,facecolors="none",edgecolors="k",label="anti-axis")
ax.set_xlim(0,360); ax.set_xlabel("geocentric ecliptic longitude (deg)")
ax.set_ylabel("geocentric ecliptic latitude (deg)")
ax.legend(frameon=False,fontsize=8)
ax.set_title("observed sky positions at perihelion epoch",fontsize=10)
ax=axes[1]
ax.hist(null,bins=50,color="0.7",label="resampled null")
ax.axvline(Robs,color="crimson",lw=1.5,label=f"in-cap R={Robs:.2f}")
ax.set_xlabel("resultant length of geocentric directions")
ax.set_ylabel("count"); ax.legend(frameon=False,fontsize=8)
ax.set_title("in-cap observed-direction clustering",fontsize=10)
fig.tight_layout()
FIG=RESULTS/"figures"; FIG.mkdir(exist_ok=True)
fig.savefig(FIG/"step_b46_geocentric_null.png",dpi=150)
print(f"wrote {FIG/'step_b46_geocentric_null.png'}")
