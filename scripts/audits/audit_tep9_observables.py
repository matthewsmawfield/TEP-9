#!/usr/bin/env python3
"""Focused, non-mutating audit of TEP-9 geometry and observable definitions.

Run from any directory with the repository Python environment. This writes only
its own audit JSON; it does not regenerate or alter the research pipeline.
The alternative calculations are diagnostics, not new confirmatory tests.
"""

import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import astropy
from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord
from astropy import units as u
import numpy as np
from scipy.stats import chi2, mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.utils.tep9_common import ECL2GAL as CORRECTED_ECL2GAL, lv, perih_dir

# The historical scripts applied this equatorial-to-Galactic matrix directly
# to ecliptic vectors.  Keep it here only as an explicit legacy comparator so
# the audit continues to measure the frame bug after the production code has
# been repaired.
LEGACY_EQ2GAL = np.array([
    [-0.0548755604, -0.8734370902, -0.4838350155],
    [ 0.4941094279, -0.4448296300,  0.7469822445],
    [-0.8676661490, -0.1980763734,  0.4559837762],
])
LEGACY_ECL2GAL = LEGACY_EQ2GAL

N_MC = 20000
SEED = 20261030
SOURCE_124 = ROOT / "scripts/steps/step_124_lineage_aphelion.py"
INPUTS = [
    SOURCE_124,
    ROOT / "scripts/steps/step_123_leg_coherence.py",
    ROOT / "scripts/steps/step_063_bidirectional_rotation.py",
    ROOT / "scripts/steps/step_065_proper_time_slip.py",
    ROOT / "scripts/utils/tep9_common.py",
    ROOT / "data/raw/sbdb/sbdb_comets_all.json",
    ROOT / "results/step_b84_pre2018_sbdb.csv",
    ROOT / "results/step_b28_bidirectional_rotation.csv",
    ROOT / "results/step_b88_lineage_aphelion.json",
]


def read_csv(relative):
    with (ROOT / relative).open() as stream:
        return list(csv.DictReader(stream))


def lonlat(vector):
    vector = vector / np.linalg.norm(vector)
    return [float(np.degrees(np.arctan2(vector[1], vector[0])) % 360),
            float(np.degrees(np.arcsin(vector[2])))]


def separation(a, b):
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b)), a @ b)))


def frame_audit():
    # Each transformed Cartesian basis vector is a column of the rotation.
    basis = SkyCoord(lon=[0, 90, 0] * u.deg, lat=[0, 0, 90] * u.deg,
                     frame=BarycentricMeanEcliptic(equinox="J2000"))
    corrected = basis.galactic.cartesian.xyz.value
    eps = np.radians(23.4392911)
    ecl_to_eq = np.array([[1, 0, 0], [0, np.cos(eps), -np.sin(eps)],
                         [0, np.sin(eps), np.cos(eps)]])
    composed = CORRECTED_ECL2GAL
    # Small differences reflect IAU frame-bias and obliquity conventions.
    assert np.max(np.abs(corrected - composed)) < 1e-6
    assert np.allclose(corrected @ corrected.T, np.eye(3), atol=1e-12)
    targets = {}
    for label, gal in [("galactic_center", np.array([1., 0., 0.])),
                       ("galactic_north", np.array([0., 0., 1.]))]:
        old = np.linalg.solve(LEGACY_ECL2GAL, gal)
        new = np.linalg.solve(corrected, gal)
        targets[label] = dict(current_ecliptic_deg=lonlat(old),
                              corrected_ecliptic_deg=lonlat(new),
                              displacement_deg=separation(old, new))
    definitions, consumers = [], []
    for path in sorted((ROOT / "scripts").rglob("*.py")):
        if path.resolve() == Path(__file__).resolve():
            continue
        source = path.read_text()
        if "-0.0548755604" in source:
            definitions.append(str(path.relative_to(ROOT)))
        if "ECL2GAL" in source or "GAL2ECL" in source:
            consumers.append(str(path.relative_to(ROOT)))
    return corrected, dict(
        reference=f"Astropy {astropy.__version__}, BarycentricMeanEcliptic J2000",
        corrected_ecl_to_gal=corrected.tolist(),
        explicit_obliquity_composition_max_difference=float(np.max(np.abs(corrected-composed))),
        targets=targets, duplicate_definitions=definitions,
        direct_matrix_consumers=consumers,
        limitation="The consumer list omits indirect dependencies and cached downstream results.")


def load_eras():
    with INPUTS[5].open() as stream:
        data = json.load(stream)
    eras = {"post2017": [], "pre2018": []}
    for values in data["data"]:
        row = dict(zip(data["fields"], values))
        name = str(row["full_name"]).strip()
        match = re.match(r"C/(\d{4})", name)
        if not match or re.match(r"C/\d{4}\s+\S+-\w", name):
            continue
        try:
            e, q, w, om, inc, arc, nobs = [float(row[k]) for k in
                ("e", "q", "w", "om", "i", "data_arc", "n_obs_used")]
        except (TypeError, ValueError):
            continue
        if not np.all(np.isfinite([e, q, w, om, inc, arc, nobs])):
            continue
        if not (0.95 <= e < 1.5 and q >= 0.1 and arc >= 30 and nobs >= 20):
            continue
        aph = -perih_dir(*np.radians([w, om, inc]))
        eras["post2017" if int(match[1]) >= 2018 else "pre2018"].append(aph)
    assert [len(eras[k]) for k in eras] == [287, 880]
    return {k: np.array(v) for k, v in eras.items()}


def original_null_function(matrix):
    # Execute only these three pure function definitions, avoiding every
    # top-level pipeline action, logger, result write, and plotting operation.
    tree = ast.parse(SOURCE_124.read_text())
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef)
             and node.name in {"gal_of", "ecl_lb", "dipole_null"}]
    assert len(nodes) == 3
    namespace = dict(np=np, math=math, N_MC=N_MC, ECL2GAL=matrix,
                     GAL2ECL=np.linalg.inv(matrix), rng=np.random.default_rng(SEED))
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE_124), "exec"), namespace)
    return namespace["dipole_null"]


def era_audit(corrected):
    eras = load_eras()
    output = {}
    for label, matrix in [("legacy_frame_full_precision", LEGACY_ECL2GAL),
                          ("corrected_frame_full_precision", corrected)]:
        print(f"Evaluating {label}...", flush=True)
        function = original_null_function(matrix)
        axes = [lv(34, -13), lv(49.9, -17), lv(214, 13),
                np.linalg.solve(matrix, np.array([-1., 0., 0.])), lv(120, -40)]
        out = {}
        # Preserve exactly the original random-call order of the first two
        # dipole_block calls, including their non-primary direction checks.
        for era, vectors in eras.items():
            values = [function(vectors, axis, "gal") for axis in axes]
            ecliptic = function(vectors, axes[0], "ecl")
            p = values[0][1]
            out[era] = dict(n=len(vectors), d_parallel=values[0][0],
                            exceedances=int(round(p*(N_MC+1)-1)),
                            p_galactic=p, p_ecliptic=ecliptic[1])
        ps = [out[era]["p_galactic"] for era in eras]
        fisher = -2 * sum(math.log(p) for p in ps)
        out["fisher"] = dict(statistic=fisher, df=4,
                              p_nominal=float(chi2.sf(fisher, 4)),
                              inputs=ps)
        output[label] = out
    with INPUTS[-1].open() as stream:
        saved = json.load(stream)["results"]
    for era, key in [("post2017", "post2017_near_parabolic"),
                     ("pre2018", "pre2018_near_parabolic")]:
        repaired = output["corrected_frame_full_precision"][era]
        saved_p = saved[key]["T1_tide_aware"]["cap_34_-13"]["p_tide_null"]
        # The result file is regenerated after the frame repair; compare it to
        # the corrected pure-function replay, not to the retired legacy value.
        assert math.isclose(repaired["p_galactic"], saved_p, rel_tol=1e-12,
                           abs_tol=1e-12)
    output["saved_era_combination"] = saved["T6_era_combination"]
    output["mc_draws_per_test"] = N_MC
    output["mc_p_floor"] = 1/(N_MC+1)
    output["interpretation"] = (
        "Fisher p-values here are nominal under independent, valid component nulls. "
        "Disjoint objects do not remove shared catalogue selection or null-model "
        "misspecification. A zero exceedance count is Monte Carlo resolution, "
        "not a measurement of an arbitrarily small true tail probability. "
        "The original |latitude| resampling and hemisphere randomization are "
        "retained; neither is a complete discovery-selection function.")
    return output


def twist_audit():
    rows = read_csv("results/step_b84_pre2018_sbdb.csv")
    deg = np.array([[float(row[k]) for k in ("d_in", "d_out", "drot")] for row in rows])
    a, b, c = np.radians(deg).T
    sa, sb, sc = np.sin(a/2)**2, np.sin(b/2)**2, np.sin(c/2)**2
    denominator = np.sin(a)*np.sin(b)
    cosine = (2*sa + 2*sb - 4*sa*sb - 2*sc)/denominator
    assert np.all(np.isfinite(cosine)) and np.all(np.abs(cosine) <= 1+1e-6)
    # A long-double independent evaluation verifies the arithmetic on the
    # stored angle triples without making the audit depend on optional
    # ``mpmath``.  This checks numerical cancellation only; it cannot supply
    # the observational covariance absent from those triples.
    ld = np.longdouble
    aa, bb, cc = np.deg2rad(deg.astype(np.longdouble)).T
    reference = ((np.cos(cc) - np.cos(aa) * np.cos(bb)) /
                 (np.sin(aa) * np.sin(bb)))
    error = float(np.max(np.abs(cosine.astype(np.longdouble)-reference)))
    # The stored angles are decimal-degree values, so the independent
    # long-double replay is limited by input quantisation at roughly 1e-9 in
    # the near-singular regime.
    assert error < 1e-8
    inc = np.array([float(row["theta"]) < 60 for row in rows])
    cells = {}
    for label, flag in [("overlap", "True"), ("non_overlap", "False")]:
        mask = np.array([row["in_code"] == flag for row in rows])
        x, y = cosine[mask & inc], cosine[mask & ~inc]
        cells[label] = dict(n=int(mask.sum()), n_in=len(x), n_out=len(y),
                            median_cos_in=float(np.median(x)),
                            median_cos_out=float(np.median(y)),
                            p_less=float(mannwhitneyu(x, y, alternative="less").pvalue))
    return dict(n=len(rows), old_retained=int(np.sum(denominator >= 0.01)),
                numerically_recovered=len(rows),
                denominator_min=float(denominator.min()),
                denominator_max=float(denominator.max()),
                max_error_against_longdouble_arithmetic=error, cohorts=cells,
                limitation="Numerical recoverability does not establish precision of the fitted orbital geometry; orbit covariance is still required.")


def residual_audit():
    rows = read_csv("results/step_b28_bidirectional_rotation.csv")
    output = {}
    for label, limit in [("q_lt_3p1", 3.1), ("all_class1", np.inf)]:
        pool = [row for row in rows if float(row["q"]) < limit]
        def column(key):
            return np.array([float(row[key]) for row in pool])
        inc = column("theta") < 60
        cat, sim = column("drot_cat"), column("drot_sim")
        delta = cat-sim
        x = np.column_stack([np.ones(len(pool)), np.log10(np.abs(column("daa_sim"))+1),
                             np.log10(column("denc")), column("q"), column("i")])
        y = np.log(sim)
        proxy = y-x@np.linalg.lstsq(x, y, rcond=None)[0]
        output[label] = dict(
            n=len(pool), cat_sim_rho=float(spearmanr(cat, sim).statistic),
            median_absolute_cat_minus_sim_deg=float(np.median(np.abs(delta))),
            direct_delta_median_in_deg=float(np.median(delta[inc])),
            direct_delta_median_out_deg=float(np.median(delta[~inc])),
            direct_delta_p_greater=float(mannwhitneyu(delta[inc], delta[~inc], alternative="greater").pvalue),
            proxy_residual_p_greater=float(mannwhitneyu(proxy[inc], proxy[~inc], alternative="greater").pvalue))
    output["interpretation"] = (
        "The regression outcome is log(drot_sim), generated by the standard "
        "dynamics model. Its residual is not an observed-minus-dynamical-model "
        "discrepancy. The direct difference is an integrator/catalogue consistency "
        "check, not an independent TEP exclusion: catalogue boundary legs also "
        "depend on standard-dynamics propagation. Population anisotropy needs a "
        "separately calibrated population and selection model.")
    return output


def time_conversion_counterexample():
    # Analytic Kepler states at two different orbital phases. The eccentricity
    # vector (periapsis direction) is unchanged by a pure phase advance.
    mu, q, eccentricity = 4*np.pi**2, 1.0, 0.99
    p = q*(1+eccentricity)
    vectors = []
    for f in (0.4, 0.7):
        radius = p/(1+eccentricity*np.cos(f))
        r = radius*np.array([np.cos(f), np.sin(f), 0.])
        v = np.sqrt(mu/p)*np.array([-np.sin(f), eccentricity+np.cos(f), 0.])
        vectors.append(np.cross(v, np.cross(r, v))/mu-r/np.linalg.norm(r))
    rotation = separation(*vectors)
    assert rotation < 1e-10
    return dict(phase_advance_deg=float(np.degrees(0.3)),
                periapsis_direction_rotation_deg=rotation,
                statement="h/r^2 is the true-anomaly sweep rate. Dividing a periapsis-direction rotation by it defines an equivalent time scale, not a derived proper-time measurement. A TEP field and observation/fit response are needed to relate them.")


def main():
    corrected, frames = frame_audit()
    result = dict(
        scope="Focused audit; original pipeline and manuscript are unmodified.",
        provenance={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS},
        frames=frames, twist=twist_audit(), residual_definition=residual_audit(),
        time_conversion=time_conversion_counterexample(),
        log_units=dict(source="steps 117-123 residual models use np.log10",
                       statement="Residual amplitudes are now expressed in true dex. The earlier natural-log labels were a unit error; rank p-values are unchanged by this unit correction.",
                       published_0p3532556112920151_in_dex=0.3532556112920151/np.log(10)),
        eras=era_audit(corrected))
    output = ROOT / "results/audits/audit_tep9_observables.json"
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"output": str(output), "twist": result["twist"],
                      "residual_definition": result["residual_definition"],
                      "eras": result["eras"]}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
