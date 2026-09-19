"""Statistical primitives with explicit probability and precision contracts."""
import numpy as np
from scipy.stats import chi2


def monte_carlo_p(exceedances, draws):
    """Finite-simulation p-value; observed statistic included in the rank."""
    if int(draws) != draws or draws < 1 or int(exceedances) != exceedances or not 0 <= exceedances <= draws:
        raise ValueError("Require integer 0 <= exceedances <= positive draws")
    return (int(exceedances)+1)/(int(draws)+1)


def fisher_combination(p_values):
    """Nominal Fisher combination. Dependence/selection require calibration."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not len(p) or not np.all(np.isfinite(p)) or np.any((p <= 0) | (p > 1)):
        raise ValueError("Fisher inputs must be nonempty, finite and in (0, 1]")
    statistic = float(-2*np.log(p).sum())
    return dict(statistic=statistic, df=2*len(p), p=float(chi2.sf(statistic, 2*len(p))),
                inputs=p.tolist(), interpretation="Nominal under independent valid null p-values; not a TEP posterior probability.")


def monte_carlo_tail(null_statistics, observed, alternative="greater"):
    """Rank against finite simulations without ever returning a zero p-value."""
    values = np.asarray(null_statistics, dtype=float)
    if values.ndim != 1 or not len(values) or not np.all(np.isfinite(values)) or not np.isfinite(observed):
        raise ValueError("Require finite observed statistic and nonempty finite null")
    if alternative == "greater":
        count = np.count_nonzero(values >= observed)
    elif alternative == "less":
        count = np.count_nonzero(values <= observed)
    elif alternative == "two-sided":
        count = np.count_nonzero(np.abs(values) >= abs(observed))
    else:
        raise ValueError("Unknown alternative")
    return monte_carlo_p(count, len(values))


def nearest_cone_probability(folded_angle_deg, targets_deg):
    """Isotropic tail for distance to the nearest of several folded cones.

    The null event is the UNION of bands about every candidate cone. Selecting
    the nearest cone and integrating only that band undercounts the tail.
    Folded polar angle is in [0, 90] with density sin(theta).
    """
    targets = np.asarray(targets_deg, dtype=float)
    if targets.ndim != 1 or not len(targets) or not np.all(np.isfinite(targets)) or np.any((targets < 0) | (targets > 90)):
        raise ValueError("Cone targets must be finite and in [0, 90]")
    if not np.isfinite(folded_angle_deg) or not 0 <= folded_angle_deg <= 90:
        raise ValueError("Folded angle must be in [0, 90]")
    delta = float(np.min(abs(targets-folded_angle_deg)))
    intervals = sorted((max(0., t-delta), min(90., t+delta)) for t in targets)
    merged = []
    for lo, hi in intervals:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return float(sum(np.cos(np.deg2rad(lo))-np.cos(np.deg2rad(hi)) for lo, hi in merged))


def rayleigh_3d(vectors):
    """3-D directional Rayleigh statistic: 3*n*Rbar^2 ~ chi-square(3).

This is the large-sample isotropic reference, not the 2-D circular
exp(-n*Rbar^2) law. Direction-dependent selection needs a separate null.
"""
    v = np.asarray(vectors, dtype=float)
    if v.ndim != 2 or v.shape[1] != 3 or not len(v) or not np.all(np.isfinite(v)):
        raise ValueError("Expected a nonempty array of finite three-vectors")
    norm = np.linalg.norm(v, axis=1)
    if np.any(norm == 0):
        raise ValueError("Directions cannot have zero length")
    rbar = float(np.linalg.norm((v/norm[:, None]).mean(axis=0)))
    return rbar, float(chi2.sf(3*len(v)*rbar*rbar, 3))
