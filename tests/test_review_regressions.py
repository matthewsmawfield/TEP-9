"""Mathematical controls, not empirical TEP evidence."""
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import chi2

from scripts.utils.statistics import (
    fisher_combination, monte_carlo_p, monte_carlo_tail,
    nearest_cone_probability)
from scripts.utils.coordinates import ECL2EQ, EQ2ECL, ECL2GAL, leg_twist_cosine
from scripts.utils import mpc_refit


def test_finite_mc_floor_and_all_exceedances():
    assert monte_carlo_p(0, 20000) == 1/20001
    assert monte_carlo_tail([0, 1, 2], 3) == 0.25
    assert monte_carlo_tail([0, 1, 2], 0) == 1
    assert monte_carlo_tail([-2, 0, 2], 2, 'two-sided') == 0.75


@pytest.mark.parametrize('count,draws', [(-1, 10), (11, 10), (1, 0), (1.2, 10), (1, 3.2)])
def test_invalid_mc_counts_fail(count, draws):
    with pytest.raises(ValueError):
        monte_carlo_p(count, draws)


def test_fisher_does_not_drop_probability_one_or_use_product_as_p():
    result = fisher_combination([0.01, 1])
    assert result['df'] == 4
    assert result['p'] == pytest.approx(chi2.sf(-2*math.log(0.01), 4))
    assert result['p'] > 0.01
    with pytest.raises(ValueError):
        fisher_combination([0, 0.01])


def test_multiple_cone_null_uses_union_including_overlap():
    rng = np.random.default_rng(1284)
    folded = np.degrees(np.arccos(rng.uniform(0, 1, 300000)))
    for obs, targets in [(30, [25, 60]), (48, [45, 50]), (80, [20, 60])]:
        delta = min(abs(obs-t) for t in targets)
        reference = (np.min(abs(folded[:, None]-np.array(targets)), axis=1) <= delta).mean()
        assert nearest_cone_probability(obs, targets) == pytest.approx(reference, abs=0.003)
    assert nearest_cone_probability(45, [45, 45]) == 0


def test_ecliptic_transform_against_spice():
    import spiceypy as sp
    assert EQ2ECL == pytest.approx(sp.pxform('J2000', 'ECLIPJ2000', 0), abs=1e-12)
    assert ECL2EQ @ EQ2ECL == pytest.approx(np.eye(3), abs=1e-14)
    assert ECL2GAL @ ECL2GAL.T == pytest.approx(np.eye(3), abs=1e-14)


def test_small_angle_triangle_retains_leg_twist():
    a = b = 0.001
    c = np.degrees(np.arccos(np.cos(np.radians(a))*np.cos(np.radians(b))))
    assert leg_twist_cosine(a, b, c) == pytest.approx(0, abs=1e-6)
    assert np.isnan(leg_twist_cosine(0, b, b))


def test_propagation_returns_requested_order(monkeypatch):
    class Sim:
        def __init__(self):
            self.particles = [SimpleNamespace(x=0., y=0., z=0., vx=0., vy=0., vz=0.)]
        @property
        def N(self):
            return len(self.particles)
        def add(self, **state):
            self.particles.append(SimpleNamespace(**state))
        def integrate(self, t, **kwargs):
            self.particles[-1].x = t
    monkeypatch.setattr(mpc_refit, 'init_sim', lambda _: Sim())
    year = 86400*mpc_refit.DAY_YR
    result = mpc_refit.propagate_states(np.zeros(3), np.zeros(3), 0, np.array([2., -1., 0.])*year)
    assert result[:, 0] == pytest.approx([2, -1, 0])


def test_failed_damping_keeps_last_accepted_orbit(monkeypatch):
    seed = np.zeros(3)
    def artificial_residual(r, v, epoch, obs):
        # Positive finite-difference derivative proposes a negative step;
        # every negative trial is deliberately worse than the initial orbit.
        delta = float(np.sum(r) + np.sum(v))
        value = 1 + delta if delta >= 0 else 2 - delta
        return np.full((len(obs), 2), value)
    monkeypatch.setattr(mpc_refit, 'residuals', artificial_residual)
    r, v, rms, keep, _ = mpc_refit.fit(seed, seed, 0, [{}]*8)
    assert r == pytest.approx(seed)
    assert v == pytest.approx(seed)
    assert rms == pytest.approx(1)


def test_phase_advance_does_not_rotate_periapsis():
    eccentricity, q, mu = .8, 1., 4*np.pi**2
    p = q*(1+eccentricity)
    ev = []
    for anomaly in (.2, .9):
        r = p/(1+eccentricity*np.cos(anomaly))*np.array([np.cos(anomaly), np.sin(anomaly), 0])
        v = np.sqrt(mu/p)*np.array([-np.sin(anomaly), eccentricity+np.cos(anomaly), 0])
        ev.append(np.cross(v, np.cross(r,v))/mu-r/np.linalg.norm(r))
    assert ev[0] == pytest.approx(ev[1], abs=1e-14)


def test_pdf_css_is_css_not_nested_style_markup():
    from scripts.utils.html_to_pdf import HTMLToPDFConverter
    css = HTMLToPDFConverter()._get_css_for_pdf({'custom_css': 'p { orphans: 3; }'})
    assert '<style' not in css
    assert 'page-break-after: avoid' in css
    assert 'orphans: 3' in css
