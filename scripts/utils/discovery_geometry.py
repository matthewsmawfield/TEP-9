"""Keplerian geometry and integration measures for discovery simulations."""
import numpy as np


def visibility_eligible(q, r_vis, dmax):
    """Necessary reach of the stated heliocentric/geocentric cuts (AU)."""
    return q < r_vis and q < dmax + 1.0


def orientation_cdf(weights):
    """Normalize each fixed object's orientation choices independently."""
    w = np.asarray(weights, dtype=float)
    if w.ndim != 2 or np.any(w < 0) or not np.all(np.isfinite(w)):
        raise ValueError('Require finite nonnegative per-object orientation weights')
    total = w.sum(axis=1)
    if np.any(total <= 0):
        raise ValueError('An object has no visible orientations in this model')
    cdf = np.cumsum(w / total[:, None], axis=1)
    cdf[:, -1] = 1.0
    return cdf


def haar_rotations(rng, n):
    """Haar rotations on SO(3), with determinant +1 (no reflections)."""
    q, r = np.linalg.qr(rng.normal(size=(n, 3, 3)))
    signs = np.sign(np.diagonal(r, axis1=1, axis2=2))
    signs[signs == 0] = 1
    q *= signs[:, None, :]
    q[:, :, 2] *= np.linalg.det(q)[:, None]
    return q


def conic_window(q, e, mu, n, outer_radius=4.0):
    """r, true anomaly, time from perihelion, and dt quadrature weights.

    Sample only the physically accessible bright window. Elliptic and
    hyperbolic anomalies preserve the supplied eccentricity; the exactly
    parabolic case uses Barker's equation. Units follow q and mu.
    """
    if q <= 0 or e < 0 or mu <= 0 or n < 3:
        raise ValueError("Invalid conic/window parameters")
    p = q * (1 + e)
    edge = max(outer_radius, q * (1 + 1e-8))
    limit = (np.arccos(np.clip((p / edge - 1) / e, -1, 1))
             if e > 0 else np.pi)
    nu = np.linspace(-limit, limit, n)
    radius = p / (1 + e * np.cos(nu))
    if e < 1:
        a = q / (1 - e)
        anomaly = 2 * np.arctan2(np.sqrt(1-e) * np.sin(nu/2),
                                 np.sqrt(1+e) * np.cos(nu/2))
        time = (anomaly - e * np.sin(anomaly)) * np.sqrt(a**3 / mu)
    elif e > 1:
        a_abs = q / (e - 1)
        anomaly = 2 * np.arctanh(np.sqrt((e-1)/(e+1)) * np.tan(nu/2))
        time = (e * np.sinh(anomaly) - anomaly) * np.sqrt(a_abs**3 / mu)
    else:
        d = np.tan(nu/2)
        time = np.sqrt(2*q**3/mu) * (d + d**3/3)
    # Kepler's area law, with trapezoid endpoint weights in true anomaly.
    dt = radius**2 / np.sqrt(mu*p) * (nu[1] - nu[0])
    dt[[0, -1]] *= 0.5
    return radius, nu, time, dt
