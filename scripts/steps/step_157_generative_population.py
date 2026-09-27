"""step_157: conventional generative population model for the arrival
anisotropy (b123).

Step 149 (b113) showed that Newtonian dynamics *transmits* an
incoming-direction anisotropy into the periapsis-rotation field, and
step 152 (b116) priced the conditioned nulls on the *observed* arrival
directions.  Neither is generative: both condition on the catalogued
directions.  This step closes the loop -- it GENERATES an incoming
population under conventional physics, passes it through a simulated
discovery-selection function, and asks whether the resulting synthetic
arrival-direction distribution reproduces the observed cap excess.

Conventional chain
  1. Source: isotropic Oort-spike aphelion directions u on the sphere.
  2. Galactic-tide injection: for each u, the per-orbit secular drift
     of perihelion distance is computed by integrating the tidal torque
     r x (T . r) over the unperturbed eccentric orbit (impulse/secular
     approximation; the same tide tensor as step_070:  diag(+OmG^2,
     -OmG^2, -4 pi G rho) in the Galactic frame).  The orbit's
     angular-momentum direction h sweeps a uniform circle perpendicular
     to u.  Two injection weights are evaluated:
       w_flux(u)  = <max(0, -Dq)> over h orientations -- the standard
                    loss-cone filling-flux proxy (per-orbit downward
                    drift across the observability boundary);
       w_hard(u)  = fraction of h draws with q0 + Dq < q_loss for
                    q0 = 8 AU, q_loss = 5 AU.
  3. Stellar-impulse diffusion (model M2): the injected density is
     mixed with an isotropic component -- passing stars randomise
     the source on ~Gyr timescales, so a fraction f_iso of the arrival
     density is isotropic.
  4. Discovery selection: rejection sampling against the empirical
     ecliptic-latitude detection kernel, kappa(beta) = smoothed
     observed latitude density / uniform density -- the first-order
     pointing model, exactly the conditioning step_152 applies but
     used here as a generative filter.
  5. Planetary diffusion to observable q is treated as direction-
     independent (inclination-averaged); stated as an assumption.

Models
  M0  isotropic source x selection            (selection-only control)
  M1  isotropic source x tide  x selection    (canonical conventional)
  M2  M1 + isotropic stellar remixing f_iso=0.5
  M3  tide injection alone, no selection      (raw tidal signature)

Tests (per model, matched n=54 and all-class-1 n=131 samples)
  T1  predicted in-cap fraction distribution vs the observed
      22/54 = 40.7% and 48/131 = 36.6% (60-deg cap on the resident
      axis).  p_model = fraction of synthetic samples with f >= f_obs.
  T2  the model density's own preferred direction (resultant of the
      predicted sky density) and its separation from the resident
      axis -- does conventional generation have any reason to point
      at (49, -17)?
  T3  per-sample resultant direction vs the axis: rank of the observed
      sample's axis-proximity among synthetic samples.

Inputs
  results/step_b28_bidirectional_rotation.csv   (arrival cohort)
  data/raw/code/code_original.html              (periapsis elements)
  Galactic-tide tensor as in step_070

Outputs
  results/step_b123_generative_population.json
  results/step_b123_generative_population.csv
  results/figures/supplementary/step_b123_generative_population.png
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent.parent))
from scripts.utils.step_logger import StepLogger
from scripts.utils.tep9_common import (
    DATA_RAW, RESULTS, lv, lb, sep, perih_dir, parse_code,
    tee_stdout)
logger = StepLogger("step_157_generative_population")
tee_stdout(logger)
logger.header("Conventional generative population model -- "
              "Oort source + Galactic tide + stellar remixing + "
              "discovery selection")

import csv
import json
import math
import numpy as np
from scipy.stats import binomtest

SEED = 20261101
CAP = 60.0
N_DIR = 4096            # source-direction grid (Fibonacci sphere)
K_PLANES = 24           # orbit-plane orientations per direction
N_REAL = 2000           # synthetic sample realisations per model
rng = np.random.default_rng(SEED)

TNO = lv(49.0, -17.0)          # resident cluster axis
CAX = lv(34.0, -13.0)          # declared comet transit axis

# ------------------------------------------------------------------
# frames -- identical to step_152
# ------------------------------------------------------------------
_EPS = math.radians(23.4392911)

def _radec_vec(ra, dec):
    ra, dec = math.radians(ra), math.radians(dec)
    return np.array([math.cos(ra) * math.cos(dec),
                     math.sin(ra) * math.cos(dec), math.sin(dec)])

_GX = _radec_vec(266.405, -28.936)
_GZ = _radec_vec(192.859508, 27.128336)
_GY = np.cross(_GZ, _GX); _GY /= np.linalg.norm(_GY)
_GX = np.cross(_GY, _GZ)
_EQ2GAL = np.stack([_GX, _GY, _GZ])

def _ecl2eq(v):
    return np.array([v[0],
                     math.cos(_EPS) * v[1] - math.sin(_EPS) * v[2],
                     math.sin(_EPS) * v[1] + math.cos(_EPS) * v[2]])

def ecl2gal(v):
    return _EQ2GAL @ _ecl2eq(v)

def gal2ecl(v):
    # inverse of ecl2gal
    eq = _EQ2GAL.T @ v
    x, y, z = eq[0], math.cos(_EPS) * eq[1] + math.sin(_EPS) * eq[2], \
        -math.sin(_EPS) * eq[1] + math.cos(_EPS) * eq[2]
    return np.array([x, y, z])

# ------------------------------------------------------------------
# Galactic tide tensor -- identical coefficients to step_070
# ------------------------------------------------------------------
G_AU   = 4 * math.pi ** 2                        # G in AU^3/(yr^2 Msun)
KMS_AUYR = 0.210805
AU_KPC = 206264.806 * 1000.0
OMG    = 26.0 * KMS_AUYR / AU_KPC                # Omega_G [1/yr]
RHO_G  = 0.1 / (AU_KPC / 1000.0) ** 3            # Msun/AU^3
KZ     = -4.0 * math.pi * G_AU * RHO_G           # -4 pi G rho [1/yr^2]
KX, KY = OMG ** 2, -OMG ** 2
TIDE = np.diag([KX, KY, KZ])
MU = 4 * math.pi ** 2
logger.info(f"tide tensor (gal frame): Kx={KX:.3e} Ky={KY:.3e} "
            f"Kz={KZ:.3e} yr^-2")

# ------------------------------------------------------------------
# observed cohort -- identical join to steps 149/152
# ------------------------------------------------------------------
orig = parse_code(str(DATA_RAW / "code" / "code_original.html"))
rows = []
for r in csv.DictReader(open(RESULTS / "step_b28_bidirectional_rotation.csv")):
    k = r["desig"]
    if k not in orig:
        continue
    ro = orig[k]
    u = -perih_dir(math.radians(ro["w"]), math.radians(ro["Om"]),
                   math.radians(ro["i"]))
    beta = math.degrees(math.asin(np.clip(u[2], -1, 1)))
    rows.append(dict(desig=k, q=float(r["q"]), u=u, beta=beta))
logger.info(f"joined {len(rows)} comets (step_b28 x CODE original)")

subsets = {"matched": [r for r in rows if r["q"] < 3.1],
           "all_c1": rows}
f_obs = {nm: np.mean([sep(r["u"], TNO) < CAP for r in rr])
         for nm, rr in subsets.items()}
n_obs = {nm: len(rr) for nm, rr in subsets.items()}
logger.info(f"observed in-cap fractions: "
            f"matched {f_obs['matched']:.3f} (n={n_obs['matched']}), "
            f"all_c1 {f_obs['all_c1']:.3f} (n={n_obs['all_c1']})")

# observed arrival resultant (per subset) and its axis separation
obs_res = {}
for nm, rr in subsets.items():
    U = np.array([r["u"] for r in rr])
    R = U.sum(axis=0); R /= np.linalg.norm(R)
    obs_res[nm] = dict(u=R.tolist(), sep_TNO=sep(R, TNO),
                       sep_CAX=sep(R, CAX))
    logger.info(f"  {nm}: observed resultant {np.round(lb(R),1)}, "
                f"sep to resident axis {obs_res[nm]['sep_TNO']:.1f} deg")

# ------------------------------------------------------------------
# discovery-selection kernel: empirical ecliptic-latitude density
# kappa(beta) = rho_obs(beta)/rho_uniform(beta), normalised to <=1
# ------------------------------------------------------------------
def build_kappa(rr, nbin=36):
    betas = np.array([math.radians(r["beta"]) for r in rr])
    sb = np.sin(betas)
    edges = np.linspace(-1, 1, nbin + 1)
    h, _ = np.histogram(sb, bins=edges)
    h = h.astype(float)
    # uniform in sin(beta) is the isotropic expectation, so kappa
    # is just the smoothed normalised histogram
    k = np.convolve(h, np.ones(3) / 3, mode="same")
    k[0] = (2 * k[0] + k[1]) / 3; k[-1] = (2 * k[-1] + k[-2]) / 3
    if k.max() > 0:
        k /= k.max()
    return edges, k

def kappa_accept(u, edges, k):
    sb = np.clip(np.asarray(u)[..., 2], -0.9999, 0.9999)
    idx = np.clip(np.searchsorted(edges, sb) - 1, 0, len(k) - 1)
    return k[idx]

kappa_tab = {nm: build_kappa(rr) for nm, rr in subsets.items()}

# ------------------------------------------------------------------
# source grid: Fibonacci sphere (uniform aphelion directions)
# ------------------------------------------------------------------
def fib_sphere(n):
    i = np.arange(n)
    z = 1 - 2 * (i + 0.5) / n
    phi = math.pi * (1 + 5 ** 0.5) * i
    r = np.sqrt(np.clip(1 - z * z, 0, 1))
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)

U_SRC = fib_sphere(N_DIR)          # ecliptic-frame aphelion dirs
U_SRC_G = np.array([ecl2gal(u) for u in U_SRC])   # galactic frame

# ------------------------------------------------------------------
# secular tidal Dq per orbit via orbit-integrated torque
# Dh = oint r x (T . r) dt  over the unperturbed orbit
# Expand r(E) = a[(cosE-e) ep + s sinE eQ]:
#   Dh = C1 (ep x T ep) + C2 (ep x T eQ + eQ x T ep) + C3 (eQ x T eQ)
# with Ci the quadrature over eccentric anomaly.
# ------------------------------------------------------------------
A_SEMI = 25000.0                    # Oort-spike semi-major axis [AU]
Q0_HARD, Q_LOSS = 8.0, 5.0          # hard-loss-cone variant [AU]
q_ref = 10.0                        # reference q for the orbit family
e_orb = 1 - q_ref / A_SEMI          # near-parabolic (e ~ 0.9996)
s_e = math.sqrt(1 - e_orb ** 2)
N_E = 512
E = np.linspace(0, 2 * math.pi, N_E, endpoint=False)
dE = E[1] - E[0]
n_mean = math.sqrt(MU / A_SEMI ** 3)
# dt = (1 - e cosE)/n dE
wdt = (1 - e_orb * np.cos(E)) / n_mean * dE
c1 = np.sum((np.cos(E) - e_orb) ** 2 * wdt)
c2 = np.sum((np.cos(E) - e_orb) * s_e * np.sin(E) * wdt)
c3 = np.sum(s_e ** 2 * np.sin(E) ** 2 * wdt)
C = np.array([c1, c2, c3]) * A_SEMI ** 2
h_orb = math.sqrt(MU * A_SEMI * (1 - e_orb ** 2))
logger.info(f"orbit family: a={A_SEMI:.0f} AU, q_ref={q_ref}, "
            f"e={e_orb:.5f}, h={h_orb:.3f} AU^2/yr")

def tide_dh(u_g, psi):
    """secular angular-momentum change per orbit for arrival dir u
    (galactic frame) with orbit-plane orientation parameter psi."""
    ep = -u_g                                   # periapsis dir = -u
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(ep, ref)) > 0.95:
        ref = np.array([1.0, 0.0, 0.0])
    eQ0 = np.cross(ref, ep); eQ0 /= np.linalg.norm(eQ0)
    h0 = np.cross(ep, eQ0)
    eQ = np.cos(psi)[:, None] * eQ0[None, :] + \
        np.sin(psi)[:, None] * h0[None, :]
    hh = np.cross(np.broadcast_to(ep, eQ.shape), eQ)
    pep = ep[None, :]
    Tep = TIDE @ ep; TeQ = (TIDE @ eQ.T).T
    a1 = np.cross(pep, np.broadcast_to(Tep, eQ.shape))
    a2 = np.cross(pep, TeQ) + np.cross(eQ, np.broadcast_to(Tep, eQ.shape))
    a3 = np.cross(eQ, TeQ)
    Dh = C[0] * a1 + C[1] * a2 + C[2] * a3      # [K,3]
    return np.einsum("ki,ki->k", Dh, hh)        # Dh . h

PSI = np.linspace(0, 2 * math.pi, K_PLANES, endpoint=False)

w_flux = np.zeros(N_DIR)
w_hard = np.zeros(N_DIR)
for i, ug in enumerate(U_SRC_G):
    dh = tide_dh(ug, PSI)                       # [K]
    dq = h_orb * dh / MU                        # Dq per orbit [AU]
    w_flux[i] = np.mean(np.maximum(0.0, -dq))
    w_hard[i] = np.mean(q_ref + dq < Q_LOSS)
for w, nm in ((w_flux, "w_flux"), (w_hard, "w_hard")):
    logger.info(f"  {nm}: max={w.max():.4g} mean={w.mean():.4g} "
                f"nonzero-frac={np.mean(w > 0):.3f}")

# model density resultant (preferred direction of the injected map)
for w, nm in ((w_flux, "flux"), (w_hard, "hard")):
    if w.max() > 0:
        Rv = (w[:, None] * U_SRC).sum(axis=0)
        rn = np.linalg.norm(Rv)
        if rn > 0:
            logger.info(f"  model {nm}: density resultant "
                        f"{np.round(lb(Rv / rn), 1)} deg, "
                        f"sep to resident axis {sep(Rv / rn, TNO):.1f} deg")

# ------------------------------------------------------------------
# generative models -> synthetic arrival samples
# ------------------------------------------------------------------
def sample_model(w_src, edges_k, k_k, n, n_real):
    """Draw n_real samples of n arrivals: direction density propto
    w_src on the source grid, then ecliptic-latitude selection."""
    p = w_src / w_src.sum() if w_src.sum() > 0 else \
        np.full(N_DIR, 1.0 / N_DIR)
    acc = kappa_accept(U_SRC, edges_k, k_k)
    # joint draw: pick grid dir ~ p, accept with kappa; equivalent
    # density p_i * acc_i -- draw directly from reweighted grid
    pj = p * acc
    pj /= pj.sum()
    out = np.empty((n_real, n))
    dirs = np.empty((n_real, n, 3))
    for ir in range(n_real):
        idx = rng.choice(N_DIR, size=n, p=pj)
        out[ir] = np.array([sep(U_SRC[j], TNO) < CAP for j in idx])
        dirs[ir] = U_SRC[idx]
    return out, dirs

MODELS = {}
iso = np.ones(N_DIR)
# M0: isotropic x selection
MODELS["M0_iso_sel"] = iso
# M1: tide(flux) x selection
MODELS["M1_tide_flux_sel"] = w_flux
# M2: tide + 50% isotropic stellar remixing, x selection
MODELS["M2_tide_remix_sel"] = 0.5 * w_flux + 0.5 * w_flux.mean() * iso
# M3: tide alone (no selection -> acceptance = 1 handled via kappa=1)
MODELS["M3_tide_nosel"] = w_flux

results = {"configs": {}, "models": {}}
results["configs"] = dict(
    a_semi=A_SEMI, q_ref=q_ref, e_orb=e_orb, q0_hard=Q0_HARD,
    q_loss=Q_LOSS, k_planes=K_PLANES, n_dir=N_DIR, n_real=N_REAL,
    seed=SEED,
    cap=CAP, axis="resident (49,-17)",
    tide=dict(KX=KX, KY=KY, KZ=KZ, OMG=OMG, RHO_G=RHO_G),
    observed={nm: dict(n=n_obs[nm], f_cap=f_obs[nm],
                       resultant=obs_res[nm]) for nm in subsets})

csv_rows = []
for nm, rr in subsets.items():
    edges_k, k_k = kappa_tab[nm]
    for mnm, w_src in MODELS.items():
        if mnm == "M3_tide_nosel":
            edges_use, k_use = np.linspace(-1, 1, 37), np.ones(36)
        else:
            edges_use, k_use = edges_k, k_k
        frac, dirs = sample_model(w_src, edges_use, k_use,
                                  n_obs[nm], N_REAL)
        # T1: in-cap fraction distribution
        fcap = frac.mean(axis=1)
        p_model = (np.sum(fcap >= f_obs[nm] - 1e-12) + 1) / (N_REAL + 1)
        # T2/T3: per-sample resultant separation to resident axis
        seps = np.array([sep(d.mean(axis=0) / np.linalg.norm(
            d.mean(axis=0)), TNO) for d in dirs])
        p_dir = (np.sum(seps <= obs_res[nm]["sep_TNO"] + 1e-12) + 1) / \
            (N_REAL + 1)
        key = f"{mnm}|{nm}"
        results["models"][key] = dict(
            model=mnm, subset=nm, n=n_obs[nm],
            f_obs=f_obs[nm],
            f_cap_med=float(np.median(fcap)),
            f_cap_p95=float(np.quantile(fcap, 0.95)),
            f_cap_p99=float(np.quantile(fcap, 0.99)),
            p_cap=float(p_model),
            sep_med=float(np.median(seps)),
            sep_p05=float(np.quantile(seps, 0.05)),
            obs_sep=obs_res[nm]["sep_TNO"],
            p_dir=float(p_dir))
        for ir in range(N_REAL):
            csv_rows.append(dict(model=mnm, subset=nm, real=ir,
                                 f_cap=float(fcap[ir]),
                                 res_sep=float(seps[ir])))
        logger.info(f"  {key}: f_cap med={np.median(fcap):.3f} "
                    f"p95={np.quantile(fcap,0.95):.3f} "
                    f"(obs {f_obs[nm]:.3f}, p={p_model:.4f}); "
                    f"res-sep med={np.median(seps):.1f} "
                    f"(obs {obs_res[nm]['sep_TNO']:.1f}, "
                    f"p={p_dir:.4f})")

# model sky-density diagnostics (the tidal signature itself)
wmaps = {}
for mnm, w in (("tide_flux", w_flux), ("tide_hard", w_hard)):
    g = U_SRC_G
    sg = np.sin(np.arctan2(g[:, 2], np.hypot(g[:, 0], g[:, 1])))
    # galactic-latitude dependence of the injected density
    bins = np.linspace(-1, 1, 19)
    prof = [float(w[(sg >= bins[i]) & (sg < bins[i + 1])].mean())
            for i in range(len(bins) - 1)]
    Rv = (w[:, None] * U_SRC).sum(axis=0)
    rn = np.linalg.norm(Rv)
    wmaps[mnm] = dict(bins=bins.tolist(), gal_lat_profile=prof,
                      density_resultant=list(lb(Rv / rn)),
                      density_resultant_sep=float(sep(Rv / rn, TNO)))
results["wmaps"] = wmaps

# ------------------------------------------------------------------
# outputs
# ------------------------------------------------------------------
out_json = RESULTS / "step_b123_generative_population.json"
out_csv = RESULTS / "step_b123_generative_population.csv"
out_fig = RESULTS / "figures/supplementary/step_b123_generative_population.png"

json.dump(dict(step="step_157_generative_population", seed=SEED,
               results=results), open(out_json, "w"), indent=1)
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["model", "subset", "real",
                                      "f_cap", "res_sep"])
    w.writeheader(); w.writerows(csv_rows)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
ax = axes[0, 0]
g = U_SRC_G
glon = np.degrees(np.arctan2(g[:, 1], g[:, 0]))
glat = np.degrees(np.arcsin(np.clip(g[:, 2], -1, 1)))
sc = ax.scatter(glon, glat, c=w_flux, s=4, cmap="viridis")
tno_g = ecl2gal(TNO)
tg_lon = math.degrees(math.atan2(tno_g[1], tno_g[0]))
tg_lat = math.degrees(math.asin(np.clip(tno_g[2], -1, 1)))
ax.plot(tg_lon, tg_lat, "*", ms=16, color="red",
        label="resident axis")
ax.set_xlabel("Galactic longitude (deg)"); ax.set_ylabel("Galactic latitude (deg)")
ax.set_title("tidal injection weight w_flux(u)")
ax.legend(loc="lower right", fontsize=8)
plt.colorbar(sc, ax=ax)

ax = axes[0, 1]
b = np.array(wmaps["tide_flux"]["bins"])
prof = wmaps["tide_flux"]["gal_lat_profile"]
ax.plot(np.degrees(np.arcsin(0.5 * (b[:-1] + b[1:]))), prof, "o-")
ax.set_xlabel("Galactic latitude (deg)")
ax.set_ylabel("mean injection weight")
ax.set_title("tidal signature: Galactic-latitude band")

for col, (nm, rr) in enumerate(subsets.items()):
    ax = axes[1, col]
    for mnm in ("M0_iso_sel", "M1_tide_flux_sel", "M2_tide_remix_sel"):
        key = f"{mnm}|{nm}"
        d = results["models"][key]
        vals = [r["f_cap"] for r in csv_rows
                if r["model"] == mnm and r["subset"] == nm]
        ax.hist(vals, bins=40, histtype="step", density=True,
                label=mnm.replace("_sel", ""))
    ax.axvline(f_obs[nm], color="k", lw=2, label=f"obs {f_obs[nm]:.3f}")
    ax.set_xlabel("synthetic in-cap fraction")
    ax.set_title(f"{nm} (n={n_obs[nm]})")
    ax.legend(fontsize=7)

fig.suptitle("step_157: conventional generative population model "
             "vs observed cap excess", y=0.99)
fig.tight_layout()
fig.savefig(out_fig, dpi=140)
logger.data_save(out_json)
logger.data_save(out_csv)
logger.data_save(out_fig)
logger.info(f"saved {out_json.name}, {out_csv.name}, {out_fig.name}")

# headline verdict
m1 = results["models"]["M1_tide_flux_sel|all_c1"]
m0 = results["models"]["M0_iso_sel|all_c1"]
logger.info(f"VERDICT: tide+selection predicts median f_cap "
            f"{m1['f_cap_med']:.3f} (selection-only {m0['f_cap_med']:.3f}) "
            f"vs observed {f_obs['all_c1']:.3f} -- model p_cap "
            f"{m1['p_cap']:.4f}, p_dir {m1['p_dir']:.4f}")
