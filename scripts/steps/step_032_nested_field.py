"""step_09: nested-field comet test -- layered depth, element decomposition,
multi-axis structure, and orbital-dynamics discriminators.

The proper-time field is layered: Solar-System domain nested inside the
galactic domain, nested inside cosmological structure (TEP-C0: the same
suppression operator propagates across environments).  A comet inbound
from the Oort cloud crosses multiple boundaries -- the accumulated
reconstruction discrepancy should therefore:

  L1. scale with traversal depth (1/a_orig = which shell it came from),
  L2. decompose into specific elements -- a temporal field couples
      primarily to ENERGY (1/a) and orientation, not uniformly,
  L3. appear on each physically-motivated axis separately (TNO boundary
      axis, heliosphere/ISM inflow axis, CMB dipole apex = field rest
      frame, galactic axes as controls),
  L4. survive the orbital-dynamics confound: retrograde comets feel
      weaker planetary impulse (faster relative velocity), so a cap
      excess among retrograde comets is non-planetary.

Data: Warsaw catalogue tables a1/b/c (Krolikowska 2014), all real.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_032_nested_field")
tee_stdout(logger)
logger.header("Nested-field multi-axis test")
import json, math, os
import numpy as np
from scripts.utils.coordinates import ECL2GAL, GAL2ECL
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact

RNG = np.random.default_rng(9)


def perih_dir(om, Om, inc):
    co, so, cO, sO, ci, si = np.cos(om), np.sin(om), np.cos(Om), np.sin(Om), np.cos(inc), np.sin(inc)
    return np.array([cO*co - sO*so*ci, sO*co + cO*so*ci, so*si])

def sep(a, b):
    return math.degrees(math.acos(np.clip(np.dot(a, b), -1, 1)))

def gal_to_ecl(l_deg, b_deg):
    l, b = math.radians(l_deg), math.radians(b_deg)
    g = np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])
    return GAL2ECL @ g

def lonlat_to_vec(l_deg, b_deg):
    l, b = math.radians(l_deg), math.radians(b_deg)
    return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])

AXES = {
    "tno":      lonlat_to_vec(34.0, -13.0),    # confined-TNO cluster axis
    "anti_tno": lonlat_to_vec(214.0, 13.0),
    "ism":      lonlat_to_vec(255.8, 5.16),    # ISM inflow (Bzowski+2015, ecliptic)
    "cmb_apex": gal_to_ecl(264.02, 48.25),     # CMB dipole apex (field rest frame)
    "cmb_anti": gal_to_ecl(84.02, -48.25),
    "gal_ctr":  gal_to_ecl(0.0, 0.0),
    "gal_pole": gal_to_ecl(0.0, 90.0),
    "ecl_pole": np.array([0.0, 0.0, 1.0]),
}

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
        d = line[3:15].strip()
        if d:
            out.setdefault(d, []).append(dict(
                model=line[132:140].strip(), qnew=line[156:159].strip(),
                datat=line[122:132].strip()))
    return out

orig_rows = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tablec.dat"))
osc_rows  = parse_orbit_table(str(DATA_RAW / "warsaw" / "warsaw_tableb.dat"))
meta      = parse_a1(str(DATA_RAW / "warsaw" / "warsaw_tablea1.dat"))

PREF = {"a": 0, "h": 0, "e": 1, "b": 2}
orig = {}
for r in orig_rows:
    k = r["desig"]
    if k not in orig or PREF.get(r["com"], 9) < PREF.get(orig[k]["com"], 9):
        orig[k] = r
PREF_OSC = {"a": 0, "g": 0, "d": 1, "e": 2, "f": 2, "b": 3, "c": 3}
osc = {}
for r in osc_rows:
    k = r["desig"]
    if k not in osc or PREF_OSC.get(r["com"], 9) < PREF_OSC.get(osc[k]["com"], 9):
        osc[k] = r

def ang_elem_diff(w1, O1, w2, O2):
    """circular-aware angular differences (deg)."""
    dw = (w2 - w1 + 180) % 360 - 180
    dO = (O2 - O1 + 180) % 360 - 180
    return dw, dO

comets = []
for k, ro in orig.items():
    if k not in osc: continue
    rs = osc[k]
    po = perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]), math.radians(ro["i"]))
    ps = perih_dir(math.radians(rs["w"]), math.radians(rs["Om"]), math.radians(rs["i"]))
    m = meta.get(k, [{}])[0]
    dw, dO = ang_elem_diff(ro["w"], ro["Om"], rs["w"], rs["Om"])
    thetas = {name: sep(-po, v) for name, v in AXES.items()}
    comets.append(dict(desig=k, sample=ro["sample"], aph=-po,
        aa_orig=ro["aa"], aa_osc=rs["aa"],
        ddir=sep(-po, -ps), daa=rs["aa"] - ro["aa"],
        dq=rs["q"] - ro["q"], de=rs["e"] - ro["e"],
        dw=dw, dOm=dO, di=rs["i"] - ro["i"],
        i_orig=ro["i"], q_orig=ro["q"],
        model=m.get("model", "?"), datat=m.get("datat", "?"),
        qnew=m.get("qnew", "?"), thetas=thetas))
spike = [c for c in comets if 0 < c["aa_orig"] < 100]
full  = [c for c in spike if c["datat"] == "full"]
print(f"{len(comets)} matched; spike n={len(spike)}, full-arc n={len(full)}")

res = {"n_matched": len(comets), "n_spike": len(spike), "n_full": len(full)}

# ---- L3: multi-axis cap scan on ddir (full-arc spike) ----
scan = {}
for name, v in AXES.items():
    th = np.array([sep(c["aph"], v) for c in full])
    d = np.array([c["ddir"] for c in full])
    inc = th < 60
    if inc.sum() >= 5 and (~inc).sum() >= 5:
        u = mannwhitneyu(d[inc], d[~inc], alternative="greater")
        scan[name] = {"n_in": int(inc.sum()), "med_in": float(np.median(d[inc])),
                      "med_out": float(np.median(d[~inc])), "p": float(u.pvalue)}
res["axis_scan_ddir_full"] = scan

# ---- L1: depth -- correction vs 1/a within and outside cap ----
depth = {}
for lab, sub in [("in60", [c for c in spike if c["thetas"]["tno"] < 60]),
                 ("out60", [c for c in spike if c["thetas"]["tno"] >= 60])]:
    aa = np.array([c["aa_orig"] for c in sub])
    dd = np.array([c["ddir"] for c in sub])
    depth[lab] = {"n": len(sub),
        "spearman_ddir_vs_aa": {"rho": float(spearmanr(aa, dd)[0]),
                                 "p": float(spearmanr(aa, dd)[1])}}
res["depth_L1"] = depth
# global depth trend
aa = np.array([c["aa_orig"] for c in spike]); dd = np.array([c["ddir"] for c in spike])
res["depth_L1"]["global"] = {"rho": float(spearmanr(aa, dd)[0]),
                              "p": float(spearmanr(aa, dd)[1])}

# ---- L2: element decomposition of the cap excess ----
elems = {}
for key, lab in [("ddir", "dir"), ("dw", "omega"), ("dOm", "Omega"),
                 ("di", "incl"), ("daa", "inv_a"), ("de", "ecc"), ("dq", "q")]:
    v = np.abs(np.array([c[key] for c in spike]))
    th = np.array([c["thetas"]["tno"] for c in spike])
    inc = th < 60
    u = mannwhitneyu(v[inc], v[~inc], alternative="greater")
    rho, p = spearmanr(th, v)
    elems[lab] = {"med_in": float(np.median(v[inc])), "med_out": float(np.median(v[~inc])),
                  "p_cap_greater": float(u.pvalue),
                  "rho_vs_theta": float(rho), "p_spearman": float(p)}
res["element_decomp_L2"] = elems

# ---- L4: orbital dynamics -- retrograde split ----
dyn = {}
for lab, sub in [("retro_i>90", [c for c in full if c["i_orig"] > 90]),
                 ("prog_i<90", [c for c in full if c["i_orig"] <= 90])]:
    th = np.array([c["thetas"]["tno"] for c in sub])
    d = np.array([c["ddir"] for c in sub])
    inc = th < 60
    if inc.sum() >= 4 and (~inc).sum() >= 4:
        u = mannwhitneyu(d[inc], d[~inc], alternative="greater")
        dyn[lab] = {"n_in": int(inc.sum()), "n_out": int((~inc).sum()),
            "med_in": float(np.median(d[inc])), "med_out": float(np.median(d[~inc])),
            "p": float(u.pvalue)}
res["dynamics_L4"] = dyn

# ---- L1b: layered shell test -- bin by 1/a, cap test per shell ----
shells = {}
for lo, hi in [(0, 20), (20, 50), (50, 100)]:
    sub = [c for c in full if lo < c["aa_orig"] <= hi]
    if len(sub) < 8: continue
    th = np.array([c["thetas"]["tno"] for c in sub])
    d = np.array([c["ddir"] for c in sub])
    inc = th < 60
    if inc.sum() >= 4 and (~inc).sum() >= 4:
        u = mannwhitneyu(d[inc], d[~inc], alternative="greater")
        shells[f"aa_{lo}_{hi}"] = {"n": len(sub), "n_in": int(inc.sum()),
            "med_in": float(np.median(d[inc])), "med_out": float(np.median(d[~inc])),
            "p": float(u.pvalue)}
res["shells_L1b"] = shells

# ---- nested axes: does ISM or CMB axis carry signal too? ----
nested = {}
for name in ["ism", "cmb_apex", "cmb_anti"]:
    th = np.array([sep(c["aph"], AXES[name]) for c in spike])
    d = np.array([c["ddir"] for c in spike])
    for cap in (45, 60):
        inc = th < cap
        if inc.sum() >= 5 and (~inc).sum() >= 5:
            u = mannwhitneyu(d[inc], d[~inc], alternative="greater")
            nested[f"{name}_{cap}"] = {"n_in": int(inc.sum()),
                "med_in": float(np.median(d[inc])), "med_out": float(np.median(d[~inc])),
                "p": float(u.pvalue)}
res["nested_axes"] = nested

# ---- epoch trend inside cap (dynamic boundary?) ----
sub = [c for c in spike if c["thetas"]["tno"] < 60]
if len(sub) > 10:
    yrs = []
    for c in sub:
        pass  # perihelion year not parsed; skip -- noted below
res["note_epoch"] = "perihelion-year trend not implemented (year not parsed in this step)"

out = str(RESULTS / "step_09_nested_field.json")
json.dump(res, open(out, "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
print("wrote", out)
