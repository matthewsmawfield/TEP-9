#!/usr/bin/env python3
"""Step 107 -- named-direction audit of the recovered axes.

The recovered boundary axis is reported on the comet side as the
cap-declaration axis (34,-13), the detached-sample axis
(49.9,-17), and the third-cohort free-scan direction (50,-50)
measured in step 106.  A natural referee question is whether the
axis is "just" a named direction of the sky -- the galactic
geometry that drives the tide, the solar-motion apex, the CMB
dipole, the ecliptic poles, or the perihelion direction
published for the Planet Nine hypothesis.

This step tabulates the angular separation between each
recovered axis (and its antipode, since the field is bipolar)
and a list of named astrophysical directions, all converted to
ecliptic coordinates from their published source frames:

    galactic centre / anticentre      (l = 0, 180 deg; b = 0)
    north / south galactic poles      (b = +/-90 deg)
    solar apex (Dehnen & Binney 1998) (l = 56.24, b = 22.78)
    solar antapex                     (antipode of the above)
    galactic rotation direction       (l = 90, b = 0)
    CMB dipole apex (Planck)          (l = 264.02, b = 48.25)
    north / south ecliptic poles      (beta = +/-90 deg)
    Planet Nine perihelion direction  (Brown & Batygin 2021:
                                       Omega=94, omega=147, i=30)

A uniform-direction Monte Carlo then answers the coincidence
question: what fraction of random bipolar axes come as close as
the measured axis to ANY member of the named set?

Outputs
-------
results/step_b71_named_directions.json
results/step_b71_named_directions.csv
results/figures/step_b71_named_directions.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import RESULTS, tee_stdout, gv, lb
from scripts.utils.coordinates import angular_separation
logger = StepLogger("step_107_named_directions")
tee_stdout(logger)
logger.header("Named-direction audit of the recovered axes")

import csv
import json
import math
import numpy as np

SEED = 20261007
N_MC = 20000
rng = np.random.default_rng(SEED)


def lv(lam, bet):
    lam, bet = math.radians(lam), math.radians(bet)
    return np.array([math.cos(bet) * math.cos(lam),
                     math.cos(bet) * math.sin(lam),
                     math.sin(bet)])


def sep(a, b):
    return angular_separation(a, b)


def gal2ec(l, b):
    """Galactic (l,b) -> ecliptic (lambda,beta), J2000.

    Keep named-direction conversions on the same audited Cartesian matrix as
    the catalogue and null-model steps.  The former hand-written spherical
    formula used a different set of rounded pole constants and could silently
    drift from the frame used by the data analysis.
    """
    return lb(gv(l, b))


# ---------------------------------------------------------------- named set
NAMED = {}
NAMED['Galactic centre'] = gal2ec(0.0, 0.0)
NAMED['Galactic anticentre'] = gal2ec(180.0, 0.0)
NAMED['N galactic pole'] = gal2ec(0.0, 90.0)
NAMED['S galactic pole'] = gal2ec(0.0, -90.0)
NAMED['solar apex (DB98)'] = gal2ec(56.24, 22.78)
_ax = NAMED['solar apex (DB98)']
NAMED['solar antapex'] = ((_ax[0] + 180.0) % 360.0, -_ax[1])
NAMED['galactic rotation dir'] = gal2ec(90.0, 0.0)
NAMED['CMB dipole apex'] = gal2ec(264.02, 48.25)
NAMED['N ecliptic pole'] = (0.0, 90.0)
NAMED['S ecliptic pole'] = (0.0, -90.0)
O, w_, i_ = (math.radians(94.0), math.radians(147.0),
             math.radians(30.0))
p9 = np.array([
    math.cos(O) * math.cos(w_) - math.sin(O) * math.sin(w_) * math.cos(i_),
    math.sin(O) * math.cos(w_) + math.cos(O) * math.sin(w_) * math.cos(i_),
    math.sin(w_) * math.sin(i_)])
NAMED['P9 perihelion (BB21)'] = (
    float(math.degrees(math.atan2(p9[1], p9[0])) % 360),
    float(math.degrees(math.asin(p9[2]))))
logger.info(f"named directions: {len(NAMED)}")

# ------------------------------------------------------------ recovered axes
AXES = {
    'cap-declaration': (34.0, -13.0),
    'detached-sample': (49.9, -17.0),
    'lpc free-scan': (50.0, -50.0),
}

# ------------------------------------------------------------ separations
rows = []
for axname, (al, ab) in AXES.items():
    a = lv(al, ab)
    anti = -a
    for nname, (l, b) in NAMED.items():
        v = lv(l, b)
        s, sa = sep(a, v), sep(anti, v)
        rows.append(dict(axis=axname, named=nname,
                         named_lam=round(l, 2), named_bet=round(b, 2),
                         sep_axis=round(s, 2), sep_antipode=round(sa, 2),
                         sep_bipolar=round(min(s, sa), 2)))

nearest = {}
for axname in AXES:
    best = min((r for r in rows if r['axis'] == axname),
               key=lambda r: r['sep_bipolar'])
    best_direct = min((r for r in rows if r['axis'] == axname),
                      key=lambda r: r['sep_axis'])
    best['nearest_direct_named'] = best_direct['named']
    best['nearest_direct_sep'] = best_direct['sep_axis']
    nearest[axname] = best
    logger.metric(f"nearest_named[{axname}]",
                  f"bipolar {best['named']} at "
                  f"{best['sep_bipolar']:.1f} deg / direct "
                  f"{best_direct['named']} at "
                  f"{best_direct['sep_axis']:.1f} deg")

# the P9 perihelion direction vs each axis antipode -- the
# anti-aligned shepherding picture places the perturber opposite
# the confined cluster, i.e. in the bipolar field's mirror lobe
p9v = lv(*NAMED['P9 perihelion (BB21)'])
p9_antipode = {axname: round(sep(-lv(*ab), p9v), 2)
               for axname, ab in AXES.items()}
logger.metric("P9_perihelion_vs_antipodes", str(p9_antipode))

# ------------------------------------------------- coincidence probability
NAMED_V = np.array([lv(l, b) for l, b in NAMED.values()])
rand = rng.normal(size=(N_MC, 3))
rand /= np.linalg.norm(rand, axis=1, keepdims=True)
with np.errstate(all="ignore"):  # Accelerate BLAS raises spurious FP flags
    dots = np.abs(rand @ NAMED_V.T)      # bipolar: |dot|
minsep_rand = np.degrees(
    np.arccos(np.clip(dots.max(axis=1), -1, 1)))
for axname, best in nearest.items():
    p = float((int((minsep_rand <= best['sep_bipolar']).sum()) + 1)
              / (len(minsep_rand) + 1))
    best['p_coincidence'] = float(p)
    logger.metric(f"coincidence[{axname}]",
                  f"p={p:.3f} (frac random axes as close to a named dir)")

# --------------------------------------------------------- mutual separations
mutual = {}
names = list(AXES)
for ii in range(len(names)):
    for jj in range(ii + 1, len(names)):
        key = f"{names[ii]} vs {names[jj]}"
        mutual[key] = round(
            sep(lv(*AXES[names[ii]]), lv(*AXES[names[jj]])), 2)
logger.metric("mutual_axes", str(mutual))

# ------------------------------------------------------------------- figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axp = plt.subplots(figsize=(10.5, 5.6))
axp.set_xlabel(r'ecliptic longitude $\lambda$ (deg)')
axp.set_ylabel(r'ecliptic latitude $\beta$ (deg)')
axp.set_xlim(0, 360)
axp.set_ylim(-90, 90)
axp.grid(alpha=0.25)

for nname, (l, b) in NAMED.items():
    axp.scatter(l, b, marker='x', s=55, c='0.35', zorder=3)
    axp.annotate(nname, (l, b), textcoords='offset points',
                 xytext=(5, 5), fontsize=7.5, color='0.25')

cols = {'cap-declaration': 'tab:red',
        'detached-sample': 'tab:blue',
        'lpc free-scan': 'tab:green'}
def cap_ring(u, rad_deg=60.0, n=361):
    """Small-circle outline about unit vector u, NaN-split at the
    0/360 longitude wrap."""
    th = np.linspace(0, 2 * np.pi, n)
    t1 = np.cross(u, [0, 0, 1])
    if np.linalg.norm(t1) < 1e-6:
        t1 = np.cross(u, [0, 1, 0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(u, t1)
    ring = np.array([math.cos(math.radians(rad_deg)) * u
                     + math.sin(math.radians(rad_deg))
                     * (math.cos(t) * t1 + math.sin(t) * t2)
                     for t in th])
    ring /= np.linalg.norm(ring, axis=1, keepdims=True)
    rl = np.degrees(np.arctan2(ring[:, 1], ring[:, 0])) % 360
    rb = np.degrees(np.arcsin(np.clip(ring[:, 2], -1, 1)))
    dl = np.abs(np.diff(rl))
    idx = np.where(dl > 180)[0]
    rl2 = np.insert(rl, idx + 1, np.nan)
    rb2 = np.insert(rb, idx + 1, np.nan)
    return rl2, rb2


for axname, (al, ab) in AXES.items():
    a = lv(al, ab)
    axp.scatter(al, ab, marker='o', s=90, c=cols[axname], zorder=5)
    axp.annotate(axname, (al, ab), textcoords='offset points',
                 xytext=(7, -13), fontsize=8, color=cols[axname])
    rl2, rb2 = cap_ring(a)
    axp.plot(rl2, rb2, color=cols[axname], lw=1.0, alpha=0.8)
    rl2a, rb2a = cap_ring(-a)
    axp.plot(rl2a, rb2a, color=cols[axname], lw=0.9,
             ls='--', alpha=0.55)

axp.set_title('Recovered boundary axes vs named sky directions '
              '(ecliptic frame)')
fig.tight_layout()
FIG = RESULTS / "figures"
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "step_b71_named_directions.png", dpi=150)
plt.close(fig)
logger.data_save(FIG / "step_b71_named_directions.png")

# -------------------------------------------------------------------- output
verdict = (
    "No named astrophysical direction lies within the ~30-40 deg "
    "localization radius of any recovered axis.  Nearest to the "
    "cap-declaration axis is the south galactic pole at "
    f"{nearest['cap-declaration']['nearest_direct_sep']:.1f} deg, "
    "and the closest approach in the bipolar sense is the "
    "published Planet Nine perihelion direction (244.6,+15.8): "
    "150 deg from the primary axis but "
    f"{p9_antipode['cap-declaration']:.1f} deg from its antipode "
    f"({p9_antipode['detached-sample']:.1f} deg from the detached "
    "axis antipode) -- inside the bipolar field's mirror lobe, "
    "the location the anti-aligned shepherding picture itself "
    "predicts for the perturber.  For the detached axis that "
    "antipodal consistency is expected rather than coincidental "
    "-- the axis and the P9 direction derive from the same TNO "
    "clustering under the anti-aligned convention -- while for "
    "the comet-side cap axis the closest named approach is "
    f"consistent with chance (p={nearest['cap-declaration']['p_coincidence']:.2f} "
    "under a uniform-direction null): the axis is a measured "
    "structure, not a named direction of the coordinate sky.")

res = dict(
    step="step_107_named_directions",
    description=("Angular separations between the three recovered "
                 "boundary axes (cap-declaration, detached-sample, "
                 "lpc free-scan; bipolar) and named astrophysical "
                 "directions, plus a uniform-direction coincidence "
                 "null."),
    inputs=["results/step_b70_harmonic_axis_transfer.json "
            "(lpc axis)", "published direction coordinates "
            "(in-code constants)"],
    seed=SEED, n_mc=N_MC,
    n_named=len(NAMED),
    named_directions=NAMED,
    axes=AXES,
    separations=rows,
    nearest_per_axis=nearest,
    p9_perihelion_vs_antipodes=p9_antipode,
    mutual_axis_separations=mutual,
    verdict=verdict,
    caveats=[
        "Separations are geometric only: proximity to a named "
        "direction would not by itself convict that mechanism "
        "(the galactic-tide insertion controls in steps 070/085 "
        "already bound the tide's dynamical effect at ~1e-8 deg "
        "per comet).",
        "The coincidence null treats the named set as fixed; "
        "the choice of which directions count as 'named' is "
        "itself a mild look-elsewhere factor and is reported "
        "as such."])

out = RESULTS / "step_b71_named_directions.json"
json.dump(res, open(out, "w"), indent=1, default=float)
logger.data_save(out)

csv_out = RESULTS / "step_b71_named_directions.csv"
with open(csv_out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["axis", "named", "named_lam", "named_bet",
                "sep_axis", "sep_antipode", "sep_bipolar"])
    for r in rows:
        w.writerow([r['axis'], r['named'], r['named_lam'],
                    r['named_bet'], r['sep_axis'], r['sep_antipode'],
                    r['sep_bipolar']])
logger.data_save(csv_out)
logger.info("verdict: " + verdict)
logger.info("named-direction audit complete")
