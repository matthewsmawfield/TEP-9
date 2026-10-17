"""Numerical controls for the boundary-transport benchmark (AU, years)."""
import numpy as np


def radial_crossing_speed(inverse_a_micro, angular_momentum, radius_au,
                          mu=4 * np.pi**2):
    """Radial speed from vis-viva and h/r, with no silent imputation.

    inverse_a_micro is in 1e-6 AU^-1; h is in AU^2/year. This is the
    two-body boundary approximation, not an independently refitted orbit.
    """
    aa, h = np.broadcast_arrays(np.asarray(inverse_a_micro, dtype=float),
                               np.asarray(angular_momentum, dtype=float))
    if (not np.isfinite(radius_au) or radius_au <= 0
            or not np.all(np.isfinite(aa)) or not np.all(np.isfinite(h))
            or np.any(h < 0)):
        raise ValueError("Require finite orbital inputs and positive radius")
    radial_squared = mu * (2 / radius_au - aa * 1e-6) - (h / radius_au)**2
    if np.any(radial_squared <= 0):
        raise ValueError("Orbit has no nonzero radial crossing at this radius")
    return np.sqrt(radial_squared)


def angular_harmonics(values, theta_deg):
    """Fit monopole, cos(theta), cos(2 theta) jointly.

    A finite uniform-theta grid does not make the intercept orthogonal to
    cos(2 theta); separate dot products leak a constant into the quadrupole.
    """
    t = np.radians(np.asarray(theta_deg, dtype=float))
    y = np.asarray(values, dtype=float)
    design = np.column_stack([np.ones_like(t), np.cos(t), np.cos(2*t)])
    coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return dict(monopole=float(coeff[0]), cos1=float(coeff[1]), cos2=float(coeff[2]))


def llr_planck_statistics(result):
    """Keep observed-sky localization and the fixed-axis null test separate."""
    rank = result['planck_dipole']['rank_by_delta_aic'] / result['n_directions']
    p = result['look_elsewhere_correction']['p_planck_axis_exceeds_observed']
    if not 0 < p <= 1 or not 0 < rank <= 1:
        raise ValueError("Invalid LLR rank or null probability")
    return dict(fixed_axis_null_p=float(p), observed_sky_rank=float(rank))


def conic_boundary_direction(q, e, periapsis, orbit_normal, radius_au, outbound=False):
    """Spatial crossing direction; nearly parabolic crossings are nearby.

    Incoming velocity and incoming position have opposite radial signs;
    the spatial field must be sampled at the position, not the velocity.
    """
    p = q*(1+e)
    cosine = (p/radius_au-1)/e
    if not -1 <= cosine <= 1:
        raise ValueError('Orbit does not cross the requested sphere')
    transverse = np.cross(orbit_normal, periapsis)
    transverse /= np.linalg.norm(transverse)
    sine = np.sqrt(max(0., 1-cosine*cosine)) * (1 if outbound else -1)
    direction = cosine*np.asarray(periapsis) + sine*transverse
    return direction / np.linalg.norm(direction)
