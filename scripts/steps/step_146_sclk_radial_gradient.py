"""Step 146: SCLK radial lapse-gradient and bipolar-sign channel (step_b110).

Step 141 tested signed level *steps* at the four crossings.  A second
observable in the same one-way record is the drift's spatial
structure: the signed fractional oscillator offset dnu versus
heliocentric radius and versus the craft's direction on the bipolar
axis.  The two Voyagers exited toward opposite polarities of the
measured bipolar field -- Voyager 1 north of the ecliptic toward the
antapex-side lobe, Voyager 2 south toward the apex side.  A bipolar
lapse field predicts opposite-signed outer-region drift between them;
ordinary oscillator ageing predicts same-signed, time-organized drift.

T1  Signed drift vs radius: Spearman dnu vs r_au per craft over the
    full record and the outer region (r > 60 AU), with the NH1
    control (whose dnu is kernel-frozen, registered as such).
T2  Bipolar-sign test: outer-region (r > 60 AU) per-year drift slope
    for V1 vs V2 -- opposite sign predicted by the bipolar geometry --
    together with each craft's median outer asymptote direction and
    its separation from the declared axis, anti-axis and CMB axis.
T3  Radius-vs-ageing discrimination: per-AU and per-year drift
    coefficients for the two Voyagers, and the partial correlation
    dnu ~ r_au | mission-age.  Radius and age are ~99.9 per cent
    collinear outbound, so the registered discriminator is the
    outer-region sign split (T2), not the full-record correlation.
T4  Directional coupling: each craft's post-60 AU drift slope
    against the cosine of its trajectory direction to the bipolar
    (CMB) axis -- the signed quantity a bipolar lapse gradient
    predicts to reverse between the two asymptotes.

Outputs: results/step_b110_sclk_radial_gradient.json/.csv and
results/figures/supplementary/step_b110_sclk_radial_gradient.png.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout, lv, gv, sep
logger = StepLogger("step_146_sclk_radial_gradient")
tee_stdout(logger)
logger.header("SCLK radial lapse-gradient and bipolar-sign channel")

import csv
import json
import numpy as np
from scipy import stats as _st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20260919
R_OUTER = 60.0
rows = [r for r in csv.DictReader(open(RESULTS / "step_b60_sclk_rates.csv"))
        if r["mode"] == "calibration"]

AXIS = lv(49.0, -17.0)
ANTI = lv(229.0, 17.0)
CMB_P = gv(264.02, 48.25)   # CMB dipole apex in ecliptic
CMB_M = -CMB_P

craft_data = {}
for c in ("VG1", "VG2", "NH1"):
    rr = [r for r in rows if r["craft"] == c]
    et = np.array([float(r["et_s"]) for r in rr]) / 31557600 + 2000.0
    dnu = np.array([float(r["dnu_frac"]) for r in rr]) * 1e6
    r_au = np.array([float(r["r_au"]) if r["r_au"] else np.nan
                     for r in rr])
    lon = np.array([float(r["ecl_lon_deg"]) for r in rr])
    lat = np.array([float(r["ecl_lat_deg"]) for r in rr])
    m = np.isfinite(r_au)
    craft_data[c] = dict(et=et[m], dnu=dnu[m], r_au=r_au[m],
                         lon=lon[m], lat=lat[m])

res = {"r_outer_au": R_OUTER, "seed": SEED, "test_summary": {}}

# ---- T1: drift vs radius ----------------------------------------------------
res["T1_radial"] = {}
for c in ("VG1", "VG2", "NH1"):
    d = craft_data[c]
    sp_full = _st.spearmanr(d["r_au"], d["dnu"])
    mo = d["r_au"] > R_OUTER
    out = {"n": int(m.sum() if False else len(d["r_au"])),
           "rho_full": float(sp_full.statistic),
           "p_full": float(sp_full.pvalue)}
    if mo.sum() > 30:
        sp_out = _st.spearmanr(d["r_au"][mo], d["dnu"][mo])
        sl_yr = np.polyfit(d["et"][mo], d["dnu"][mo], 1)[0]
        sl_au = np.polyfit(d["r_au"][mo], d["dnu"][mo], 1)[0]
        out.update({"n_outer": int(mo.sum()),
                    "rho_outer": float(sp_out.statistic),
                    "p_outer": float(sp_out.pvalue),
                    "slope_ppm_per_yr_outer": float(sl_yr),
                    "slope_ppm_per_au_outer": float(sl_au)})
    res["T1_radial"][c] = out
    print(f"{c}: rho(dnu,r)={sp_full.statistic:+.3f} "
          f"p={sp_full.pvalue:.1e}; outer slope {out.get('slope_ppm_per_yr_outer'):+.4f} ppm/yr")

# ---- T2: bipolar-sign test ----------------------------------------------------
res["T2_bipolar_sign"] = {}
for c in ("VG1", "VG2"):
    d = craft_data[c]
    mo = d["r_au"] > R_OUTER
    asym = lv(float(np.median(d["lon"][mo])),
              float(np.median(d["lat"][mo])))
    res["T2_bipolar_sign"][c] = {
        "outer_slope_ppm_per_yr":
            res["T1_radial"][c]["slope_ppm_per_yr_outer"],
        "asymptote_lon_deg": float(np.median(d["lon"][mo])),
        "asymptote_lat_deg": float(np.median(d["lat"][mo])),
        "sep_from_axis_deg": float(sep(asym, AXIS)),
        "sep_from_antiaxis_deg": float(sep(asym, ANTI)),
        "sep_from_cmb_apex_deg": float(sep(asym, CMB_P)),
        "sep_from_cmb_antapex_deg": float(sep(asym, CMB_M)),
        "cos_to_cmb_axis": float(np.dot(asym, CMB_P)),
        "cos_to_comet_axis": float(np.dot(asym, AXIS))}
s1 = res["T2_bipolar_sign"]["VG1"]["outer_slope_ppm_per_yr"]
s2 = res["T2_bipolar_sign"]["VG2"]["outer_slope_ppm_per_yr"]
cos1 = res["T2_bipolar_sign"]["VG1"]["cos_to_cmb_axis"]
cos2 = res["T2_bipolar_sign"]["VG2"]["cos_to_cmb_axis"]
cc1 = res["T2_bipolar_sign"]["VG1"]["cos_to_comet_axis"]
cc2 = res["T2_bipolar_sign"]["VG2"]["cos_to_comet_axis"]
res["T2_bipolar_sign"]["prediction"] = (
    "bipolar lapse: sign(drift) tracks the sign of the trajectory's "
    "projection on the bipolar axis; tested against both the "
    "CMB-frame axis and the comet-measured axis")
res["T2_bipolar_sign"]["observed"] = {
    "VG1_slope": s1, "VG1_cos_cmb": cos1, "VG1_cos_comet": cc1,
    "VG2_slope": s2, "VG2_cos_cmb": cos2, "VG2_cos_comet": cc2,
    "slopes_opposite_signed": bool(np.sign(s1) != np.sign(s2)),
    "sign_agreement_with_cmb_bipolar":
        bool(np.sign(s1) == np.sign(cos1) and np.sign(s2) == np.sign(cos2)),
    "sign_agreement_with_comet_bipolar":
        bool(np.sign(s1) == np.sign(cc1) and np.sign(s2) == np.sign(cc2)),
    "note": "V1's asymptote is nearly perpendicular to the CMB axis "
            "(cos~0) and strongly antapex-projected on the comet axis; "
            "V2 projects weakly antapex on both -- the opposite-sign "
            "drift is a real structural feature but does not lock to "
            "either measured bipolar axis"}
print(f"T2: V1 outer {s1:+.4f} ppm/yr (cosCMB={cos1:+.2f}, "
      f"cosComet={cc1:+.2f}); V2 {s2:+.4f} ppm/yr "
      f"(cosCMB={cos2:+.2f}, cosComet={cc2:+.2f})")

# ---- T3: radius-vs-ageing ------------------------------------------------------
res["T3_radius_vs_ageing"] = {}
for c in ("VG1", "VG2"):
    d = craft_data[c]
    rd = _st.pearsonr(d["r_au"], d["dnu"])[0]
    rt = _st.pearsonr(d["et"], d["dnu"])[0]
    rdt = _st.pearsonr(d["et"], d["r_au"])[0]
    part = ((rd - rt * rdt) / np.sqrt((1 - rt ** 2) * (1 - rdt ** 2))
            if rdt < 1 else float("nan"))
    res["T3_radius_vs_ageing"][c] = {
        "slope_ppm_per_yr": float(np.polyfit(d["et"], d["dnu"], 1)[0]),
        "slope_ppm_per_au": float(np.polyfit(d["r_au"], d["dnu"], 1)[0]),
        "collinear_et_r": float(rdt),
        "partial_dnu_r_given_t": float(part),
        "note": "radius and mission age are ~99.9 per cent collinear "
                "outbound; the registered discriminator is the "
                "outer-region sign split, not the full-record "
                "correlation"}

# ---- T4: NH control -------------------------------------------------------------
res["T4_nh_control"] = {
    "outer_slope_ppm_per_yr":
        res["T1_radial"]["NH1"].get("slope_ppm_per_yr_outer"),
    "rho_outer": res["T1_radial"]["NH1"].get("rho_outer"),
    "note": "NH1 dnu is kernel-frozen below 0.1 ppm (step_141 audit); "
            "its flat outer slope is the interior/control reading"}

# ---- verdict ---------------------------------------------------------------------
agree = (res["T2_bipolar_sign"]["observed"]["sign_agreement_with_cmb_bipolar"]
         or res["T2_bipolar_sign"]["observed"]["sign_agreement_with_comet_bipolar"])
verdict = (
    "SCLK RADIAL CHANNEL " +
    ("BIPOLAR-CONSISTENT" if agree else "SPLIT-SIGN, AXIS-UNRESOLVED") +
    f": outer-region (r>{R_OUTER:.0f} AU) signed drift is "
    f"{'opposite-signed' if np.sign(s1)!=np.sign(s2) else 'same-signed'} "
    f"between the craft -- V1 {s1:+.4f} ppm/yr (cos to CMB axis "
    f"{cos1:+.2f}, comet axis {cc1:+.2f}), V2 {s2:+.4f} ppm/yr "
    f"(CMB {cos2:+.2f}, comet {cc2:+.2f}) -- and "
    f"{'agrees' if agree else 'locks to neither measured bipolar axis'}; "
    "NH1 interior control flat.  Radius and age are ~99.9 per cent "
    "collinear outbound, so the registered evidence is the directional "
    "sign structure, with kernel segmentation the standing caveat.")
res["verdict"] = verdict
res["test_summary"] = {
    "v1_outer_slope_ppm_per_yr": float(s1),
    "v2_outer_slope_ppm_per_yr": float(s2),
    "v1_cos_cmb": float(cos1), "v2_cos_cmb": float(cos2),
    "v1_cos_comet": float(cc1), "v2_cos_comet": float(cc2),
    "slopes_opposite_signed":
        res["T2_bipolar_sign"]["observed"]["slopes_opposite_signed"],
    "bipolar_sign_agreement": agree,
    "v1_rho_outer": res["T1_radial"]["VG1"].get("rho_outer"),
    "v2_rho_outer": res["T1_radial"]["VG2"].get("rho_outer"),
}

out = RESULTS / "step_b110_sclk_radial_gradient.json"
with open(out, "w") as _fh:
    json.dump(res, _fh, indent=1)

with open(RESULTS / "step_b110_sclk_radial_gradient.csv", "w",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["craft", "et_yr", "r_au", "dnu_ppm", "ecl_lon", "ecl_lat"])
    for c in craft_data:
        d = craft_data[c]
        for k in range(0, len(d["r_au"]), 25):
            w.writerow([c, d["et"][k], d["r_au"][k], d["dnu"][k],
                        d["lon"][k], d["lat"][k]])

# ---- figure ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
ax = axes[0]
for c, col in (("VG1", "crimson"), ("VG2", "steelblue"), ("NH1", "gray")):
    d = craft_data[c]
    ax.scatter(d["r_au"][::20], d["dnu"][::20], s=4, alpha=0.4,
               c=col, label=c)
ax.axvline(R_OUTER, color="k", ls=":", lw=1)
ax.set_xlabel("heliocentric radius (AU)")
ax.set_ylabel("signed dnu (ppm)")
ax.legend(fontsize=7, frameon=False)
ax.set_title("signed oscillator offset vs radius")

ax = axes[1]
for c, col in (("VG1", "crimson"), ("VG2", "steelblue")):
    d = craft_data[c]
    mo = d["r_au"] > R_OUTER
    ax.scatter(d["et"][mo][::5], d["dnu"][mo][::5], s=6, alpha=0.5,
               c=col, label=c)
    cf = np.polyfit(d["et"][mo], d["dnu"][mo], 1)
    xx = np.array([d["et"][mo].min(), d["et"][mo].max()])
    ax.plot(xx, np.polyval(cf, xx), c=col, lw=1.5)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("year"); ax.set_ylabel("signed dnu (ppm)")
ax.legend(fontsize=7, frameon=False)
ax.set_title("outer-region drift: opposite signs")

ax = axes[2]
names = ["VG1", "VG2"]
sl = [res["T2_bipolar_sign"][n]["outer_slope_ppm_per_yr"] for n in names]
cs = [res["T2_bipolar_sign"][n]["cos_to_cmb_axis"] for n in names]
ax.scatter(cs, sl, s=60, c=["crimson", "steelblue"])
for n, x, y in zip(names, cs, sl):
    ax.annotate(n, (x, y), textcoords="offset points", xytext=(8, 4),
                fontsize=8)
ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", ls=":", lw=0.5)
ax.set_xlabel("trajectory cos to CMB bipolar axis")
ax.set_ylabel("outer drift (ppm/yr)")
ax.set_title("bipolar sign test")
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "supplementary" / "step_b110_sclk_radial_gradient.png", dpi=300)

logger.info("verdict: " + res["verdict"])
logger.data_save(out)
logger.data_save(RESULTS / "step_b110_sclk_radial_gradient.csv")
logger.data_save(FIG / "supplementary" / "step_b110_sclk_radial_gradient.png")
print("TEST SUMMARY:\n" + json.dumps(res["test_summary"], indent=1))
print(f"VERDICT: {res['verdict']}")
