"""step_08: the clock channel -- osculating-to-original orbit discrepancy
vs direction to the TNO boundary axis.

TEP nuance missed by the directional-excess test (step_07): under TEP a
proper-time boundary does not primarily change WHICH direction comets
arrive from -- it changes the comet's CLOCK during passage.  The Warsaw
catalogue's "original" orbits are back-integrated under standard
dynamics; if a comet crossed a temporal boundary, that reconstruction is
biased.  The observable is therefore:

  - the DISCREPANCY between the measured (osculating) orbit and the
    back-integrated (original) orbit -- direction-selective if and only
    if the discrepancy source sits on the boundary axis;
  - direction-selective structure in the original ENERGY (1/a_orig);
  - direction-selective need for non-gravitational (NG) solutions --
    unexplained momentum exchange is the standard patch for exactly the
    kind of anomaly a temporal gradient produces.

Data: J/A+A/567/A126 (Krolikowska 2014) tables a1, b, b4, c -- all real.
Nulls: permutation of directions (preserves aa/element distributions),
plus control axes (galactic pole, apex, anti-axis) -- the signature must
be axis-specific, not generic.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_031_clock_channel")
tee_stdout(logger)
logger.header("Orbit-reconstruction clock channel")
import json, math, os
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu

RNG = np.random.default_rng(77)


def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def parse_orbit_table(path):
    rows = []
    for line in open(path):
        if len(line) < 115: continue
        try:
            rows.append(dict(sample=line[0:2].strip(), com=line[3].strip(),
                desig=line[5:17].strip(),
                q=float(line[42:56]), e=float(line[56:70]),
                w=float(line[70:82]), Om=float(line[82:94]),
                i=float(line[94:106]), aa=float(line[106:115])))
        except ValueError:
            continue
    return rows

def parse_a1(path):
    out = {}
    for line in open(path):
        if len(line) < 160: continue
        desig = line[3:15].strip()
        if not desig: continue
        out.setdefault(desig, []).append(dict(
            qosc=line[40:46].strip(), model=line[132:140].strip(),
            qnew=line[156:159].strip(), datat=line[122:132].strip()))
    return out

# ---- load ----
orig_rows = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
osc_rows  = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tableb.dat"))
meta      = parse_a1(str(DATA_RAW / "warsaw" / "warsaw_tablea1.dat"))

# preferred original orbit per comet: 'a'/'h' (entire-data) first
PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
orig = {}
for r in orig_rows:
    k = r["desig"]
    if k not in orig or PREF.get(r["com"], 9) < PREF.get(orig[k]["com"], 9):
        orig[k] = r
# preferred osculating per comet: 'a' or 'g' (entire data), then 'd'
PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}
osc = {}
osc_alt = {}   # keep secondary solutions for the duplicate-solution test
for r in osc_rows:
    k = r["desig"]
    if k not in osc:
        osc[k] = r
    else:
        osc_alt.setdefault(k, []).append(r)
        if PREF_OSC.get(r["com"], 9) < PREF_OSC.get(osc[k]["com"], 9):
            osc[k] = r

ax = math.radians(34.0); bx = math.radians(-13.0)
TNO = np.array([math.cos(bx)*math.cos(ax), math.cos(bx)*math.sin(ax), math.sin(bx)])
CTRL = {
    "gal_pole": np.array([0, 0, -1.0]) @ np.linalg.inv(ECL2GAL) * -1,  # N gal pole in ecl
    "anti_tno": -TNO,
}
# north galactic pole in ecliptic J2000 = (192.86, +27.13) ecl -> unit vec
lgp, bgp = math.radians(192.86), math.radians(27.13)
CTRL["gal_pole"] = np.array([math.cos(bgp)*math.cos(lgp), math.cos(bgp)*math.sin(lgp), math.sin(bgp)])

comets = []
for k, ro in orig.items():
    if k not in osc: continue
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    ps = perih_dir(math.radians(osc[k]["w"]), math.radians(osc[k]["Om"]), math.radians(osc[k]["i"]))
    m = meta.get(k, [{}])[0]
    comets.append(dict(desig=k, sample=ro["sample"], com=ro["com"],
        aph_orig=-po, aph_osc=-ps,
        aa_orig=ro["aa"], aa_osc=osc[k]["aa"],
        ddir=sep(-po, -ps), daa=osc[k]["aa"] - ro["aa"],
        e_orig=ro["e"], q=ro["q"],
        model=m.get("model", "?"), qnew=m.get("qnew", "?"),
        theta=sep(-po, TNO)))
print(f"{len(comets)} comets matched osc+orig")

spike = [c for c in comets if 0 < c["aa_orig"] < 100]
new   = [c for c in comets if 0 < c["aa_orig"] < 35]
allb  = [c for c in comets if c["aa_orig"] > -50]

def perm_p(vals, thetas, stat_fn, n_perm=20000, greater=True):
    obs = stat_fn(vals, thetas)
    cnt = 0
    for _ in range(n_perm):
        t = RNG.permutation(thetas)
        s = stat_fn(vals, t)
        if (s >= obs) == greater: cnt += 1
    return obs, (cnt + 1) / (n_perm + 1)

res = {"n_matched": len(comets)}
for name, samp in [("spike", spike), ("new", new), ("all_boundish", allb)]:
    if len(samp) < 10:
        res[name] = {"n": len(samp)}; continue
    th = np.array([c["theta"] for c in samp])
    s = {"n": len(samp)}
    # 1. does the osc->orig correction magnitude depend on axis distance?
    dd = np.array([c["ddir"] for c in samp])
    rho, p0 = spearmanr(th, dd)
    obs, p = perm_p(dd, th, lambda v, t: spearmanr(t, v)[0])
    s["corr_ddir_vs_theta"] = {"rho": rho, "p_perm": p}
    # signed: is the correction toward or away from the axis?
    # 2. energy: does |Daa| or signed Daa correlate with theta?
    da = np.array([c["daa"] for c in samp])
    rho2, _ = spearmanr(th, np.abs(da))
    obs2, p2 = perm_p(np.abs(da), th, lambda v, t: spearmanr(t, v)[0])
    s["corr_absDaa_vs_theta"] = {"rho": rho2, "p_perm": p2}
    rho3, _ = spearmanr(th, da)
    obs3, p3 = perm_p(da, th, lambda v, t: spearmanr(t, v)[0])
    s["corr_signedDaa_vs_theta"] = {"rho": rho3, "p_perm": p3}
    # 3. original energy vs theta
    aa = np.array([c["aa_orig"] for c in samp])
    rho4, _ = spearmanr(th, aa)
    obs4, p4 = perm_p(aa, th, lambda v, t: spearmanr(t, v)[0])
    s["corr_aaorig_vs_theta"] = {"rho": rho4, "p_perm": p4}
    # 4. inside/outside 60 deg cap
    inc = th < 60
    if inc.sum() >= 5 and (~inc).sum() >= 5:
        u = mannwhitneyu(dd[inc], dd[~inc], alternative="greater")
        s["ddir_in_vs_out_60"] = {"n_in": int(inc.sum()),
            "median_in": float(np.median(dd[inc])),
            "median_out": float(np.median(dd[~inc])), "p_greater": u.pvalue}
        u2 = mannwhitneyu(np.abs(da)[inc], np.abs(da)[~inc], alternative="greater")
        s["absDaa_in_vs_out_60"] = {"median_in": float(np.median(np.abs(da)[inc])),
            "median_out": float(np.median(np.abs(da)[~inc])), "p_greater": u2.pvalue}
    # 5. NG-need by direction (NG/NGun full-arc solutions = orbit needs non-GR;
    #    NGun entries carry published NG parameters in tableb4)
    ng = np.array([1 if c["model"] in ("NG", "NGun") else 0 for c in samp])
    if inc.sum() >= 5 and (~inc).sum() >= 5:
        tab = [[int(ng[inc].sum()), int(inc.sum() - ng[inc].sum())],
               [int(ng[~inc].sum()), int((~inc).sum() - ng[~inc].sum())]]
        orr, pf = fisher_exact(tab, alternative="greater")
        s["ng_frac_in_vs_out_60"] = {"in": tab[0], "out": tab[1],
            "p_fisher_greater": pf}
    res[name] = s

# control: same correlations vs control axes on spike sample
thg = np.array([sep(c["aph_orig"], CTRL["gal_pole"]) for c in spike])
tha = np.array([sep(c["aph_orig"], CTRL["anti_tno"]) for c in spike])
dd = np.array([c["ddir"] for c in spike])
res["controls_spike"] = {
    "gal_pole_corr_ddir": float(spearmanr(thg, dd)[0]),
    "gal_pole_p": float(spearmanr(thg, dd)[1]),
    "antitno_corr_ddir": float(spearmanr(tha, dd)[0]),
    "antitno_p": float(spearmanr(tha, dd)[1]),
}

# duplicate-solution comets: |aa_a - aa_e| discrepancy vs direction
dups = []
seen2 = {}
for r in orig_rows:
    seen2.setdefault(r["desig"], []).append(r)
for k, rs in seen2.items():
    if len(rs) > 1 and k in orig:
        po = perih_dir(math.radians(orig[k]["w"]), math.radians(orig[k]["Om"]), math.radians(orig[k]["i"]))
        dups.append(dict(desig=k, theta=sep(-po, TNO),
            daa=abs(rs[0]["aa"] - rs[1]["aa"]),
            ddir=sep(-po, -perih_dir(math.radians(rs[1]["w"]), math.radians(rs[1]["Om"]), math.radians(rs[1]["i"])))))
res["duplicate_solutions"] = dups

out = str(RESULTS / "step_08_clock_channel.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print("RESULT PAYLOAD:\n" + json.dumps(res, indent=1, default=float))
logger.data_save(out)