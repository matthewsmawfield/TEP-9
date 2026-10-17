"""Mathematical controls, not empirical TEP evidence."""
import json
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
    def artificial_residual(r, v, epoch, obs, geom=None):
        # Positive finite-difference derivative proposes a negative step;
        # every negative trial is deliberately worse than the initial orbit.
        delta = float(np.sum(r) + np.sum(v))
        value = 1 + delta if delta >= 0 else 2 - delta
        return np.full((len(obs), 2), value)
    fake_geom = {"ets": np.zeros(8), "robs": np.zeros((8, 3)),
                 "bad": np.zeros(8, bool), "ra": np.zeros(8),
                 "dec": np.zeros(8), "cosdec": np.ones(8)}
    monkeypatch.setattr(mpc_refit, 'residuals', artificial_residual)
    monkeypatch.setattr(mpc_refit, 'observer_geom', lambda obs: fake_geom)
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


def test_disformal_branch_preserves_sign_separately_from_unsigned_observable():
    result_path = Path(__file__).resolve().parents[1] / 'results' / 'step_b131_disformal_slip.json'
    result = json.loads(result_path.read_text())
    comparison = result['branch_comparison']
    admissible = comparison['admissible_signed_transport_contrast_yr']
    excluded = comparison['excluded_signed_transport_contrast_yr']
    amplitude = comparison['observable_amplitude_contrast_yr']

    assert result['constants']['B0'] > 0
    assert admissible == pytest.approx(-excluded)
    assert abs(admissible) == pytest.approx(amplitude)
    assert amplitude == pytest.approx(result['headline']['measured_contrast_yr'], abs=1e-10)


def test_boundary_speed_uses_radius_and_radial_component():
    from scripts.utils.boundary_transport import radial_crossing_speed
    mu, q, radius = 4*np.pi**2, 2., 250.
    h = np.sqrt(2*mu*q)  # parabolic orbit
    radial = radial_crossing_speed(0., h, radius)
    assert radial**2 + (h/radius)**2 == pytest.approx(2*mu/radius)
    assert radial < np.sqrt(2*mu/100.)
    with pytest.raises(ValueError):
        radial_crossing_speed(1e6/50., h, radius)
    with pytest.raises(ValueError):
        radial_crossing_speed(np.nan, h, radius)


def test_harmonic_fit_does_not_turn_constant_into_quadrupole():
    from scripts.utils.boundary_transport import angular_harmonics
    theta = np.linspace(0, 180, 721)
    flat = angular_harmonics(np.ones_like(theta)*7, theta)
    assert flat['cos1'] == pytest.approx(0, abs=1e-14)
    assert flat['cos2'] == pytest.approx(0, abs=1e-14)
    mixed = 7 + 3*np.cos(np.radians(theta)) - 2*np.cos(2*np.radians(theta))
    assert angular_harmonics(mixed, theta) == pytest.approx(
        dict(monopole=7, cos1=3, cos2=-2))


def test_llr_null_probability_is_not_observed_sky_rank():
    from scripts.utils.boundary_transport import llr_planck_statistics
    result = {'n_directions': 2664,
              'planck_dipole': {'rank_by_delta_aic': 226},
              'look_elsewhere_correction': {
                  'p_planck_axis_exceeds_observed': 1/1001}}
    stats = llr_planck_statistics(result)
    assert stats['fixed_axis_null_p'] == 1/1001
    assert stats['observed_sky_rank'] == 226/2664
    assert stats['fixed_axis_null_p'] != stats['observed_sky_rank']


def test_discovery_rotations_are_proper_and_isotropic():
    from scripts.utils.discovery_geometry import haar_rotations
    rotations = haar_rotations(np.random.default_rng(1981), 20000)
    assert np.linalg.det(rotations) == pytest.approx(np.ones(20000), abs=1e-12)
    assert np.max(abs(rotations[:, :, 0].mean(0))) < .02
    assert (rotations[:, :, 0]**2).mean(0) == pytest.approx(np.ones(3)/3, abs=.01)


@pytest.mark.parametrize('e', [.95, 1., 1.1])
def test_discovery_time_measure_matches_kepler_elapsed_time(e):
    from scripts.utils.discovery_geometry import conic_window
    radius, anomaly, time, dt = conic_window(1., e, 4*np.pi**2, 20001)
    assert np.all(radius > 0)
    assert np.all(np.diff(time) > 0)
    assert radius[[0, -1]] == pytest.approx([4., 4.])
    assert dt.sum() == pytest.approx(time[-1]-time[0], rel=1e-7)
    # Independent area-law differential check, away from endpoint stencils.
    derivative = np.gradient(time, anomaly)
    assert derivative[1:-1] == pytest.approx(
        radius[1:-1]**2 / np.sqrt(4*np.pi**2*(1+e)), rel=1e-6)


def test_raw_input_inventory_keeps_leaf_paths():
    from scripts.audits.audit_publication import (literal_raw_dependencies, input_record,
                                                declared_input_strings)
    source = "x = DATA_RAW / 'code' / 'original.html'\ny = DATA_RAW / 'sbdb' / name"
    assert literal_raw_dependencies(source) == ['code/original.html', 'sbdb']
    assert list(declared_input_strings({'catalogue': 'data/raw/catalogue.csv'})) == ['data/raw/catalogue.csv']
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        record = input_record(Path(directory), 'absent.csv', 0)
        assert record['missing'] is True
        assert input_record(Path(directory), 'catalogue', 0)['annotated'] is True


def test_transport_quadrature_against_analytic_small_field_limit():
    # Load only pure numerical definitions, never the step's data reads or scan.
    import ast
    root = Path(__file__).resolve().parents[1]
    tree = ast.parse((root/'scripts/steps/step_165_disformal_slip.py').read_text())
    definitions = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                   and n.name in {'env_b', '_layer', 'K_matrix'}]
    namespace = {'np': np, 'SIG_B': 1., 'N_S': 8001}
    exec(compile(ast.Module(body=definitions, type_ignores=[]), '<quadrature>', 'exec'), namespace)
    u, width = 1e-6, 8*1.496e11
    g = np.array([0., .5, 1., 2.])
    computed = namespace['K_matrix'](g, u, width)
    # Change variable x=logistic(s): integral x^3(1-x) dx / width.
    lo, hi = 1/(1+np.exp(6)), 1/(1+np.exp(-6))
    primitive = lambda x: x**4/4 - x**5/5
    expected = (u*g)**4/width*(primitive(hi)-primitive(lo))
    assert computed == pytest.approx(expected, rel=1e-7, abs=1e-60)


def test_discovery_brightness_is_integrated_over_time():
    import ast
    root = Path(__file__).resolve().parents[1]
    tree = ast.parse((root/'scripts/steps/step_163_footprint_forward.py').read_text())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'weights')
    namespace = {'np': np}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), '<weights>', 'exec'), namespace)
    geometry = dict(rr=np.array([1., 2.]), gd=np.ones((1, 2)),
                    el=np.full((1, 2), 100.), dec=np.zeros((1, 2)),
                    dt_weights=np.array([1., 4.]))
    weighted = namespace['weights'](geometry, 3., 4., 60., -40.)
    assert weighted == pytest.approx([2.])  # 1/1^2*1 + 1/2^2*4


def test_code_exports_reject_duplicate_solutions_and_mismatched_legs():
    from scripts.utils.catalogue_exports import code_export_models, require_matching_code_models
    def row(name, model):
        return '<tr>' + ''.join('<td>'+v+'</td>' for v in [name, model]+['0']*12) + '</tr>'
    body = ('<table>'+row('C/1950 K1', 'a1')+'</table>').encode()
    assert code_export_models(body) == {'C/1950 K1': 'a1'}
    with pytest.raises(ValueError, match='multiple solutions'):
        code_export_models((row('C/1950 K1', 'a1')+row('C/1950 K1', 'a2')).encode())
    with pytest.raises(ValueError, match='disagree'):
        require_matching_code_models({'original': {'C/1950 K1': 'a1'},
                                      'future': {'C/1950 K1': 'a2'}})


def test_refit_provenance_hashes_stored_cache_separately_from_wire(monkeypatch, tmp_path):
    import hashlib
    cache = tmp_path/'obs'/'comet.json'
    cache.parent.mkdir()
    cache.write_bytes(b'[1, 2]')
    monkeypatch.setattr(mpc_refit, 'PROV_PATH', tmp_path/'provenance.json')
    wire = b'{"ADES_DF": [1, 2]}'
    mpc_refit._prov_add('obs/comet.json', 'https://example.invalid', wire, 'control')
    record = json.loads((tmp_path/'provenance.json').read_text())['files']['obs/comet.json']
    assert record['sha256'] == hashlib.sha256(cache.read_bytes()).hexdigest()
    assert record['wire_sha256'] == hashlib.sha256(wire).hexdigest()


def test_visibility_scope_and_per_object_weights_preserve_fixed_cohort():
    from scripts.utils.discovery_geometry import visibility_eligible, orientation_cdf
    assert visibility_eligible(2., 3., 4.)
    assert not visibility_eligible(14., 3., 4.)
    assert not visibility_eligible(3., 3., 4.)
    # Absolute brightness must not change an object's weight in a null
    # conditioned on that object's already-observed shape and epoch.
    cdf = orientation_cdf([[1., 3.], [100., 300.]])
    assert cdf[0] == pytest.approx(cdf[1])
    with pytest.raises(ValueError, match='no visible'):
        orientation_cdf([[0., 0.]])


def test_parabolic_spatial_crossings_are_not_antipodal():
    from scripts.utils.boundary_transport import conic_boundary_direction
    peri, normal = np.array([1., 0., 0.]), np.array([0., 0., 1.])
    incoming = conic_boundary_direction(2., 1., peri, normal, 250.)
    outgoing = conic_boundary_direction(2., 1., peri, normal, 250., True)
    angle = np.degrees(np.arccos(incoming @ outgoing))
    assert angle == pytest.approx(20.526191797245076)
    assert incoming[0] < 0 and outgoing[0] < 0
    assert incoming[1] == pytest.approx(-outgoing[1])


def test_even_transport_shape_requires_the_field_parity_assumption():
    from scripts.utils.boundary_transport import angular_harmonics
    theta = np.linspace(0, 180, 721)
    signed = np.cos(np.radians(theta))**4
    offset = (1+np.cos(np.radians(theta)))**4
    assert angular_harmonics(signed, theta)['cos1'] == pytest.approx(0, abs=1e-12)
    assert angular_harmonics(offset, theta)['cos1'] > 1


def test_nasa_daily_archive_validation_rejects_access_pages_and_wrong_cadence():
    from scripts.utils.nasa_archive import validate_daily_ascii
    body = b'1963 1 0 ' + b'1 '*11 + b'\n1963 2 0 ' + b'2 '*11 + b'\n'
    assert validate_daily_ascii(body, 14, min_rows=2) == body
    with pytest.raises(ValueError):
        validate_daily_ascii(b'<html>access gate</html>', 14, min_rows=2)
    with pytest.raises(ValueError, match='time fields'):
        validate_daily_ascii(body.replace(b'1 0 ', b'1 12 '), 14, min_rows=2)


def test_discovery_null_draws_one_orientation_per_fixed_object():
    import ast
    import math
    from scripts.utils.discovery_geometry import orientation_cdf
    root = Path(__file__).resolve().parents[1]
    tree = ast.parse((root/'scripts/steps/step_163_footprint_forward.py').read_text())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'null_draws')
    namespace = dict(np=np, math=math, orientation_cdf=orientation_cdf,
                     rng=np.random.default_rng(123), N_MC=100, CAP=60,
                     TNO=np.array([1., 0., 0.]),
                     sep=lambda a, b: np.degrees(np.arccos(np.clip(a @ b, -1, 1))))
    exec(compile(ast.Module(body=[fn], type_ignores=[]), '<null_draws>', 'exec'), namespace)
    # Every orientation of object 1 points in x, every orientation of
    # object 2 in y. A pooled resample could duplicate either object;
    # the conditional null must always retain exactly one of each.
    directions = np.array([[1., 0., 0.], [1., 0., 0.],
                           [0., 1., 0.], [0., 1., 0.]])
    rv, _, cap, _ = namespace['null_draws'](
        directions, np.array([1., 0., 0., 1.]), 2, np.array([1., 0., 0.]))
    assert rv == pytest.approx(np.full(100, np.sqrt(.5)))
    assert cap == pytest.approx(np.full(100, .5))
