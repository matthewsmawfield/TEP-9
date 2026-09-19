"""Fixed J2000 frames used by orbital elements and DE440s states.

Ecliptic J2000 is SPICE ECLIPJ2000 (IAU 1976 mean obliquity).
The ICRS/equatorial-to-Galactic rotation must follow the ecliptic-to-
equatorial rotation; it cannot be applied directly to orbital vectors.
Matrices act on column vectors. Row arrays use the transpose.
"""
import numpy as np

OBLIQUITY_J2000_RAD = np.deg2rad(84381.448 / 3600.0)
_c, _s = np.cos(OBLIQUITY_J2000_RAD), np.sin(OBLIQUITY_J2000_RAD)
ECL2EQ = np.array([[1., 0., 0.], [0., _c, -_s], [0., _s, _c]])
EQ2ECL = ECL2EQ.T
EQ2GAL = np.array([
    [-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
    [0.4941094278755837, -0.4448296299600112, 0.7469822444972189],
    [-0.8676661490190047, -0.1980763734312015, 0.4559837761750669],
])
ECL2GAL = EQ2GAL @ ECL2EQ
GAL2ECL = ECL2GAL.T


def angular_separation(a, b):
    """Stable angle in degrees; vector lengths need not be exactly one."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != (3,) or b.shape != (3,):
        raise ValueError("angular_separation requires two three-vectors")
    if not np.all(np.isfinite([a, b])) or min(np.linalg.norm(a), np.linalg.norm(b)) == 0:
        raise ValueError("angular_separation requires finite, nonzero vectors")
    return float(np.rad2deg(np.arctan2(np.linalg.norm(np.cross(a, b)), a @ b)))


def leg_twist_cosine(d_in, d_out, d_between):
    """Relative leg azimuth from degree-valued spherical triangle sides.

The half-angle identity avoids cancellation of cosines near one. Only
numerically singular or impossible triangles return NaN; observational
uncertainty must be propagated separately using the orbital covariance.
"""
    a, b, c = np.deg2rad([d_in, d_out, d_between])
    if not np.all(np.isfinite([a, b, c])) or min(a, b, c) < 0 or max(a, b, c) > np.pi:
        return float("nan")
    denominator = np.sin(a) * np.sin(b)
    if denominator <= 64 * np.finfo(float).eps:
        return float("nan")
    sa, sb, sc = np.sin(np.array([a, b, c]) / 2)**2
    value = (2*sa + 2*sb - 4*sa*sb - 2*sc) / denominator
    tolerance = 64*np.finfo(float).eps / denominator
    if abs(value) > 1 + tolerance:
        return float("nan")
    return float(np.clip(value, -1, 1))
