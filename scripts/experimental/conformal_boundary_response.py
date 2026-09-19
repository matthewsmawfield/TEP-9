#!/usr/bin/env python3
"""Synthetic control for a prescribed, static conformal boundary.

This checks the leading nonrelativistic response of the canonical B=0
matter metric. It is not a field-equation solution, Solar-System fit,
or observational evidence. Coordinates and time are dimensionless:
z=z_physical/L, t=t_physical*v_normal_initial/L.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

C = 299792458.0


def profile(z, kind):
    tanh = np.tanh(z)
    sech2 = 1.0 - tanh*tanh
    if kind == 'step':
        return 0.5*(1.0+tanh), 0.5*sech2
    if kind == 'bump':
        return sech2, -2.0*sech2*tanh
    raise ValueError(kind)


def crossing(epsilon, kind, rtol=2e-11):
    vn0, vp = 3000.0, 1500.0
    kappa = C*C*epsilon/(vn0*vn0)
    start, stop = -12.0, 12.0

    def rhs(t, y):
        _, slope = profile(y[0], kind)
        return [y[1], -kappa*slope]

    def exit_boundary(t, y):
        return y[0]-stop

    exit_boundary.terminal = True
    exit_boundary.direction = 1
    sol = solve_ivp(rhs, (0., 100.), [start, 1.],
                    events=exit_boundary, rtol=rtol, atol=rtol/100,
                    max_step=0.15, method='DOP853')
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError('The chosen synthetic trajectory did not transmit')
    initial_f = profile(start, kind)[0]
    final_f = profile(stop, kind)[0]
    expected_u = np.sqrt(1.-2.*kappa*(final_f-initial_f))
    energy = 0.5*sol.y[1]**2 + kappa*profile(sol.y[0], kind)[0]
    result = dict(profile=kind, prescribed_delta_lnA_amplitude=epsilon,
        initial_normal_velocity_m_s=vn0, parallel_velocity_m_s=vp,
        final_normal_velocity_m_s=float(sol.y[1,-1]*vn0),
        analytic_final_normal_velocity_m_s=float(expected_u*vn0),
        analytic_agreement_m_s=float(abs(sol.y[1,-1]-expected_u)*vn0),
        max_energy_error_normalized=float(np.max(abs(energy-energy[0]))),
        velocity_direction_change_deg=float(np.degrees(
            np.arctan2(vp, sol.y[1,-1]*vn0)-np.arctan2(vp,vn0))),
        travel_time_change_in_L_over_vn0=float(sol.t[-1]-(stop-start)),
        endpoint_delta_lnA=float(epsilon*(final_f-initial_f)))
    if result['analytic_agreement_m_s'] > 1e-5:
        raise AssertionError(result)
    if result['max_energy_error_normalized'] > 1e-8:
        raise AssertionError(result)
    return result


def main():
    root = Path(__file__).resolve().parents[2]
    rows = [crossing(eps, kind) for kind in ['step', 'bump']
            for eps in [-1e-11, -1e-12, 0., 1e-12, 1e-11]]
    for row in rows:
        if row['profile'] == 'bump' or row['prescribed_delta_lnA_amplitude'] == 0:
            assert abs(row['velocity_direction_change_deg']) < 1e-7
    refinement = crossing(1e-11, 'step', rtol=2e-13)
    base = rows[4]
    convergence = abs(base['final_normal_velocity_m_s']-
                      refinement['final_normal_velocity_m_s'])
    assert convergence < 1e-5
    output = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
        evidence_status='Synthetic mathematical response control; no empirical fit or detection',
        assumptions=['B=0', 'prescribed static planar field',
                     'leading weak-field nonrelativistic test-particle dynamics',
                     'no solar gravity, backreaction, screening or observation model'],
        equation='a_extra = -c^2 grad(ln A)',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        canonical_definition_sha256=hashlib.sha256((root/'core/definitions.md').read_bytes()).hexdigest(),
        n_cases=len(rows), cases=rows,
        tolerance_refinement_delta_m_s=convergence,
        interpretation='A step can deflect while a transmitting symmetric planar bump returns the asymptotic velocity. Neither velocity angle is automatically a heliocentric periapsis rotation.')
    out = root/'results/experiments/conformal_boundary_response.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(rows),
        max_analytic_error_m_s=max(r['analytic_agreement_m_s'] for r in rows),
        tolerance_refinement_delta_m_s=convergence, output=str(out)), indent=2))


if __name__ == '__main__':
    main()
