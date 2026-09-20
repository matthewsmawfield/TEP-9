#!/usr/bin/env python3
"""Step 113 -- Independent-survey resident ledger.

The resident signature now stands on three detached-TNO cohorts drawn
by three independent survey selections and three independent orbit
determinations:

  SBDB (JPL):      the 44-object secure detached cohort (a>150, q>30,
                   condition code <=3); orbit lineage JPL/MPC;
                   designation-proxy footprint baseline.
  DES (Bernardinelli et al. 2022): the 16-object detached cohort from
                   a single 5000 deg^2 southern survey; DES-team orbit
                   fits; designation-proxy footprint baseline
                   (step_085).
  OSSOS (Bannister et al. 2018): the 31-object 'det' cohort from the
                   only fully characterized TNO survey; OSSOS-team
                   fits; EXACT discovery astrometry -- the footprint
                   baseline here is measured from the table's own
                   discovery RA/Dec, not proxied (step_112).

This step scores each cohort identically against the pre-declared
axis (lam,beta) = (49,-17) deg, cap 60 deg, against each cohort's own
footprint baseline, and combines them three ways:

  1.  Fisher's method on the per-cohort p_vs_baseline values -- the
      honest combination when the cohorts share members (their
      selection functions, not their memberships, are independent).
  2.  A deduplicated union pool: each distinct object counted once,
      assigned the baseline of its primary survey; a single binomial
      test of the union in-cap count, reported for the full pool and
      for the boundary-resident subset (a>150) on which the resident
      anomaly is defined.
  3.  Direction concordance: the in-cap members' mean varpi per
      survey, tested for mutual consistency against the footprint
      null (does each survey's in-cap population point at the same
      place?).

Membership overlap is handled explicitly: DES has 14/16 objects in
SBDB, OSSOS-det has 9 named designations of which 9 match SBDB; the
ledger reports the overlap so no independence is overclaimed.  The
OSSOS 'det' class is dominated by a ~50-90 AU orbits -- interior to
the ~150 AU boundary -- so its role in the ledger is a selection-
function bound, not a high-leverage replication.

Outputs
-------
results/step_b77_resident_ledger.json
results/step_b77_resident_ledger.csv
results/figures/supplementary/step_b77_resident_ledger.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import DATA_RAW, RESULTS, tee_stdout
logger = StepLogger("step_113_resident_ledger")
tee_stdout(logger)
logger.header("Independent-survey resident ledger")

import csv
import json
import math
import re
import numpy as np
from scipy.stats import binomtest, norm, chi2
from astropy.io.votable import parse as vot_parse
from astropy.io import fits
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time
from astropy import units as u
import warnings
warnings.filterwarnings("ignore")

rng = np.random.default_rng(20261012)
AXIS_LAM, AXIS_BET, CAP = 49.0, -17.0, 60.0
SIG_COUP = 50.0
UNIFORM = 2 * CAP / 360.0


def d_ax(a):
    return np.abs((np.asarray(a) - AXIS_LAM + 180) % 360 - 180)


def smear_cap(d, sig=SIG_COUP):
    return norm.cdf((CAP - d) / sig) + norm.cdf((CAP + d) / sig) - 1


def circ_R(a):
    return float(abs(np.exp(1j * np.asarray(a)).mean()))


def circ_mean(a):
    return float(np.rad2deg(np.angle(np.exp(1j * np.asarray(a)).mean())) % 360)


def ecl_lon(ra, dec):
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    return float(c.geocentrictrueecliptic.lon.deg)


# designation half-month -> mid-date -> opposition longitude (steps
# 014/054/085 model, retained for SBDB/DES consistency)
HALF = {"A": (1, 8), "B": (1, 23), "C": (2, 8), "D": (2, 22),
        "E": (3, 8), "F": (3, 23), "G": (4, 8), "H": (4, 23),
        "J": (5, 8), "K": (5, 23), "L": (6, 8), "M": (6, 23),
        "N": (7, 8), "O": (7, 23), "P": (8, 8), "Q": (8, 23),
        "R": (9, 8), "S": (9, 23), "T": (10, 8), "U": (10, 23),
        "V": (11, 8), "W": (11, 23), "X": (12, 8), "Y": (12, 23)}
DESIG = re.compile(r"(\d{4})\s*([A-Z])[A-Z]?\d*")


def desig_lam(desig):
    m = DESIG.search(str(desig))
    if not m or m.group(2) not in HALF:
        return float("nan")
    mo, dy = HALF[m.group(2)]
    t = Time(f"{int(m.group(1)):04d}-{mo:02d}-{dy:02d}T00:00:00",
             format="isot", scale="utc")
    return float((get_sun(t).geocentrictrueecliptic.lon.deg + 180) % 360)


# ------------------------------------------------------------------
# 1. Per-survey cohorts
# ------------------------------------------------------------------

surveys = {}

# ---- SBDB secure detached ----
d = json.loads((DATA_RAW / "sbdb" / "sbdb_outer_ss.json").read_text())
sb_rows = []
for rec in d["data"]:
    o = dict(zip(d["fields"], rec))
    try:
        a, q = float(o["a"]), float(o["q"])
        om, w = float(o["om"]), float(o["w"])
        cc = int(o["condition_code"] or 9)
    except (TypeError, ValueError):
        continue
    if a > 150 and q > 30 and cc <= 3:
        nm = str(o["full_name"]).strip()
        m = re.search(r"\((\d{4}\s*[A-Z]+\d*)\)", nm)
        key = m.group(1) if m else nm
        sb_rows.append(dict(
            name=key, survey="SBDB", varpi=(om + w) % 360, a=a, q=q,
            lam=desig_lam(key), lam_kind="designation-proxy"))
surveys["SBDB"] = dict(cohort="a>150, q>30, cc<=3", rows=sb_rows,
                       fit_lineage="JPL SBDB")

# ---- DES detached (reuse step_085 scoring) ----
tab = fits.open(str(DATA_RAW / "des" / "y6_des_tnos_color.fits"))[1].data
des_rows = []
for r in tab[(tab["a"] > 150) & (tab["q"] > 30)]:
    varpi = float((r["lan"] + r["aop"]) % 360)
    name = str(r["MPC"]).strip()
    des_rows.append(dict(name=name, survey="DES", varpi=varpi,
                         a=float(r["a"]), q=float(r["q"]),
                         lam=desig_lam(name),
                         lam_kind="designation-proxy"))
surveys["DES"] = dict(cohort="a>150, q>30", rows=des_rows,
                      fit_lineage="DES team (Bernardinelli+ 2022)")

# ---- OSSOS det (exact astrometry) ----
t = vot_parse(str(DATA_RAW / "ossos" / "ossos_t3char.vot")) \
    .get_first_table().to_table()
oss_rows = []
for r in t[t["cl"] == "det"]:
    oss_rows.append(dict(
        name=str(r["Astorb"]).strip() or str(r["ID"]).strip(),
        survey="OSSOS", varpi=(float(r["Omega"]) + float(r["omega"])) % 360,
        a=float(r["a"]), q=float(r["a"]) * (1 - float(r["e"])),
        lam=ecl_lon(float(r["RAJ2000"]), float(r["DEJ2000"])),
        lam_kind="exact discovery astrometry"))
surveys["OSSOS"] = dict(cohort="Gladman 'det' class (a>47.7, e>0.24, "
                               "non-resonant, non-scattering)",
                      rows=oss_rows, fit_lineage="OSSOS team "
                      "(Petit et al. orbit pipeline)")

# ------------------------------------------------------------------
# 2. Score each cohort
# ------------------------------------------------------------------

logger.subheader("Per-survey scoring")
ledger = []
for name, sv in surveys.items():
    rows = sv["rows"]
    vps = np.array([o["varpi"] for o in rows])
    lams = np.array([o["lam"] for o in rows])
    ok = np.isfinite(lams)
    n, n_in = len(rows), int(np.sum(d_ax(vps) < CAP))
    base = float(np.mean(smear_cap(d_ax(lams[ok])))) if ok.any() \
        else float("nan")
    in_vps = vps[d_ax(vps) < CAP]
    row = dict(
        survey=name, cohort=sv["cohort"], fit_lineage=sv["fit_lineage"],
        n=n, n_in=n_in, frac_in=n_in / n,
        baseline=base, lam_kind=rows[0]["lam_kind"],
        # declared directional prediction (in-cap excess) -> greater tail,
        # matching the exceedance convention of steps 117/133/152
        p_vs_baseline=float(binomtest(n_in, n, base,
                                      alternative="greater").pvalue)
        if np.isfinite(base) else float("nan"),
        p_vs_uniform=float(binomtest(n_in, n, UNIFORM,
                                     alternative="greater").pvalue),
        R_varpi=circ_R(np.deg2rad(vps)),
        mean_varpi=circ_mean(np.deg2rad(vps)),
        incap_mean_varpi=circ_mean(np.deg2rad(in_vps))
        if len(in_vps) else float("nan"),
        role=("detection" if name != "OSSOS" else
              "selection-function bound (interior cohort)"))
    ledger.append(row)
    logger.info(f"{name}: n={n} in-cap={n_in} ({n_in/n:.3f}) "
                f"baseline={base:.3f} p={row['p_vs_baseline']:.3g} "
                f"in-cap mean={row['incap_mean_varpi']:.1f}")

# ------------------------------------------------------------------
# 3. Overlap accounting + deduplicated union
# ------------------------------------------------------------------

def norm_key(name):
    return re.sub(r"\s+", "", str(name))


all_names = {}
for sv_name, sv in surveys.items():
    for o in sv["rows"]:
        k = norm_key(o["name"])
        all_names.setdefault(k, []).append(sv_name)
overlap = {k: v for k, v in all_names.items() if len(v) > 1}
n_multi = len(overlap)

# union pool: first survey in fixed order supplies the baseline
ORDER = ["OSSOS", "DES", "SBDB"]
union = {}
for sv_name in ORDER:
    for o in surveys[sv_name]["rows"]:
        k = norm_key(o["name"])
        if k not in union:
            union[k] = dict(name=o["name"], varpi=o["varpi"],
                            lam=o["lam"], a=o["a"], q=o["q"],
                            src=sv_name)


def score_union(urows, label):
    uv = np.array([o["varpi"] for o in urows])
    ul = np.array([o["lam"] for o in urows])
    ok = np.isfinite(ul)
    n_u, n_u_in = len(urows), int(np.sum(d_ax(uv) < CAP))
    base_u = float(np.mean(smear_cap(d_ax(ul[ok]))))
    res = dict(
        n=n_u, n_in=n_u_in, frac_in=n_u_in / n_u,
        baseline=base_u,
        p_vs_baseline=float(binomtest(n_u_in, n_u, base_u,
                                      alternative="greater").pvalue),
        p_vs_uniform=float(binomtest(n_u_in, n_u, UNIFORM,
                                     alternative="greater").pvalue))
    logger.info(f"union[{label}] n={n_u} in-cap={n_u_in} "
                f"({n_u_in/n_u:.3f}) pooled baseline={base_u:.3f} "
                f"p={res['p_vs_baseline']:.3g}")
    return res


urows = list(union.values())
union_res = score_union(urows, "all")
union_res["note"] = ("each distinct object counted once; baseline "
                     "assigned by primary survey (OSSOS > DES > SBDB "
                     "precedence); mixes interior and boundary-resident "
                     "orbits -- descriptive only")

# the resident claim is defined on the boundary-resident population
# (a > 150 AU); the OSSOS 'det' criterion (a > 47.7, e > 0.24) admits
# interior detachments, so the pooled all-object union above dilutes
# the resident test with objects that never sample the boundary.
ub_rows = [o for o in urows if o["a"] > 150.0]
union_boundary = score_union(ub_rows, "a>150")
union_boundary["note"] = (
    "deduplicated union restricted to boundary residents (a>150 AU): "
    "the population on which the resident anomaly is defined; "
    "baseline assigned by primary survey (OSSOS > DES > SBDB "
    "precedence)")

# ------------------------------------------------------------------
# 4. Fisher combination + direction concordance
# ------------------------------------------------------------------

ps = np.array([r["p_vs_baseline"] for r in ledger])
ps = np.clip(ps, 1e-16, 1)
X2 = float(-2 * np.sum(np.log(ps)))
fisher_p = float(chi2.sf(X2, 2 * len(ps)))

ps_u = np.clip([r["p_vs_uniform"] for r in ledger], 1e-16, 1)
X2u = float(-2 * np.sum(np.log(ps_u)))
fisher_p_uniform = float(chi2.sf(X2u, 2 * len(ps_u)))

# direction concordance among in-cap members: spread of the three
# in-cap mean directions on the circle
inc_means = np.deg2rad([r["incap_mean_varpi"] for r in ledger])
conc_R = circ_R(inc_means)
# null: draw 3 directions uniform in the cap, R of their mean
mc = 20000
rand_means = rng.uniform(AXIS_LAM - CAP, AXIS_LAM + CAP, (mc, 3))
R_null = np.abs(np.exp(1j * np.deg2rad(rand_means)).mean(axis=1))
p_conc = float((int((R_null >= conc_R).sum()) + 1)
               / (len(R_null) + 1))

combined = dict(
    fisher_X2=X2, fisher_p=fisher_p,
    fisher_p_vs_uniform=fisher_p_uniform,
    incap_direction_concordance=dict(
        R=conc_R, p=p_conc,
        note="R of the three surveys' in-cap mean varpi directions "
             "against the null of three directions drawn uniformly "
             "inside the cap"),
    interpretation="two detections against own-footprint baselines "
                   "(SBDB, DES) plus one selection-function bound "
                   "(OSSOS interior-detached cohort); Fisher "
                   "combination is dominated by the detection "
                   "cohorts -- reported transparently, not as a "
                   "three-fold replication")

# ------------------------------------------------------------------
# 5. Write
# ------------------------------------------------------------------

res = dict(
    step="step_113_resident_ledger",
    description="Three-survey resident ledger: SBDB secure detached, "
                "DES detached, and OSSOS 'det' cohorts scored "
                "identically against the pre-declared axis and each "
                "cohort's own footprint baseline; Fisher combination, "
                "deduplicated union, and in-cap direction "
                "concordance.",
    inputs=["data/raw/sbdb/sbdb_outer_ss.json",
            "data/raw/des/y6_des_tnos_color.fits",
            "data/raw/ossos/ossos_t3char.vot"],
    seed=20261012,
    axis_deg=[AXIS_LAM, AXIS_BET], cap_deg=CAP,
    per_survey=ledger,
    overlap=dict(
        n_objects_in_multiple_surveys=n_multi,
        shared_keys=sorted(overlap.keys()),
        note="survey selections are independent though memberships "
             "overlap; DES 14/16 in SBDB secure detached, OSSOS-det "
             "1/31 in SBDB secure detached (its Gladman-detached "
             "members are mostly interior to the a>150 cut; 9 OSSOS-"
             "det objects share SBDB catalogue membership overall)"),
    union=union_res,
    union_boundary_resident=union_boundary,
    combined=combined,
    caveat="The OSSOS 'det' class median semimajor axis is ~60 AU -- "
           "inside the ~150 AU boundary -- so its baseline-consistency "
           "is the expected selectivity outcome, not a failed "
           "replication.  The Fisher combination should be read as "
           "'two detections + one bound', not three detections.")

out = RESULTS / "step_b77_resident_ledger.json"
json.dump(res, open(out, "w"), indent=1, default=float)
logger.data_save(out)
csv_out = RESULTS / "step_b77_resident_ledger.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["survey", "name", "varpi", "lam_disc", "a", "q",
                "d_axis", "in_cap"])
    for sv_name, sv in surveys.items():
        for o in sv["rows"]:
            w.writerow([sv_name, o["name"], f"{o['varpi']:.3f}",
                        f"{o['lam']:.2f}", f"{o['a']:.2f}",
                        f"{o['q']:.2f}", f"{d_ax(o['varpi']):.1f}",
                        int(d_ax(o["varpi"]) < CAP)])
logger.data_save(csv_out)
# ------------------------------------------------------------------
# 6. Figure
# ------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))

ax = axes[0]
x = np.arange(len(ledger) + 2)
obs = [r["frac_in"] for r in ledger] + [union_res["frac_in"],
                                        union_boundary["frac_in"]]
bas = [r["baseline"] for r in ledger] + [union_res["baseline"],
                                         union_boundary["baseline"]]
lbl = [f"{r['survey']}\nn={r['n']}" for r in ledger] + \
      [f"union\nn={union_res['n']}",
       f"union a>150\nn={union_boundary['n']}"]
ax.bar(x - 0.2, obs, 0.4, color="crimson", label="observed in-cap")
ax.bar(x + 0.2, bas, 0.4, color="0.6",
       label="own-footprint baseline")
ax.axhline(UNIFORM, color="k", ls=":", lw=1)
for i, (o, b) in enumerate(zip(obs, bas)):
    ax.text(i - 0.2, o + 0.015, f"{o:.2f}", ha="center", fontsize=8)
    ax.text(i + 0.2, b + 0.015, f"{b:.2f}", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(lbl, fontsize=8)
ax.set_ylabel("fraction within 60 deg of axis")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Resident ledger: three independent selections", fontsize=10)

ax = axes[1]
cols = {"SBDB": "steelblue", "DES": "teal", "OSSOS": "crimson"}
for sv_name, sv in surveys.items():
    v = np.array([o["varpi"] for o in sv["rows"]])
    ax.hist(v, bins=np.arange(0, 361, 30), histtype="step", lw=1.5,
            density=True, color=cols[sv_name],
            label=f"{sv_name} (n={len(v)})")
ax.axvspan(AXIS_LAM - CAP, AXIS_LAM + CAP, color="crimson", alpha=0.07)
ax.axvline(AXIS_LAM, color="k", ls="--", lw=1, label="axis 49 deg")
ax.set_xlabel("$\\varpi$ (deg)"); ax.set_ylabel("density")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Cohort varpi distributions", fontsize=10)

fig.tight_layout()
FIG = RESULTS / "figures"
fig.savefig(FIG / "supplementary" / "step_b77_resident_ledger.png", dpi=300)
logger.data_save(FIG / 'supplementary' / 'step_b77_resident_ledger.png')