"""Replay retained statistics and compare the October 14 regression snapshot.

This does not substitute for replaying the missing raw catalogue inputs.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'reviews/regression_investigation_20261014'


def read(path):
    return json.loads(path.read_text())


def main():
    before = OUT/'before'
    results = ROOT/'results'
    old = read(before/'step_b131_disformal_slip.json')
    now = read(results/'step_b131_disformal_slip.json')
    o99 = read(before/'step_b99_cross_channel.json')
    n99 = read(results/'step_b99_cross_channel.json')
    global_result = read(results/'step_b90_global_synthesis.json')
    landscape = np.load(results/'step_b90_s_landscape.npz')
    checks = {}
    for label in ['eval_at_tno_axis', 'eval_at_comet_axis']:
        a = global_result[label]
        rank = float((1+np.count_nonzero(landscape['S'] >= a['S_stat']))
                     / (len(landscape['S'])+1))
        score = float(-2*np.log(a['channel_p']).sum())
        checks[label] = dict(recorded_rank=a['p_global'], replayed_rank=rank,
                             recorded_score=a['S_stat'], recomputed_score=score,
                             matches=bool(rank == a['p_global'] and
                                          np.isclose(score, a['S_stat'])))
    report = dict(
        scope='All five publication-review findings and additional input, simulation and audit defects; not a full pipeline rerun.',
        global_saved_landscape_replay=checks,
        fisher=dict(before=o99['T3_fisher'], after=n99['T3_fisher']),
        boundary=dict(before_headline=old['headline'], after_headline=now['headline'],
                      velocity_provenance=now['velocity_provenance'],
                      approximation_diagnostics=now['admissibility_ledger']['approximation_diagnostics']),
        energy_step145_replay_identical=(read(before/'step_b109_transfer_function.json') ==
                                         read(results/'step_b109_transfer_function.json')),
        blocked_replays={
            'step_034': 'Missing original CODE and Warsaw raw catalogues.',
            'step_126': 'Missing original SBDB, MPC and CODE raw catalogues; stored landscape replayed.',
            'step_163': 'Attempted run stops at missing CODE original catalogue; numerical controls tested independently.',
            'step_164': 'Missing original CODE/SBDB catalogues; source reviewed; saved results retained.'},
        not_regressions=[
            'The fixed-catalogue directional rank is not a null-data global discovery probability.',
            'Replacing encounter distance our_denc with energy change our_daa repaired a wrong-variable defect.',
            'Fitted source modulation is not an independently derived field profile.'])
    OUT.mkdir(exist_ok=True)
    (OUT/'results_comparison.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    assert all(x['matches'] for x in checks.values())
    assert report['energy_step145_replay_identical']
    print('Saved-landscape ranks and scores reproduce; corrected energy result reproduces unchanged.')
    print('Nominal Fisher p:', report['fisher']['before']['fisher_p_nominal'],
          '->', report['fisher']['after']['fisher_p_nominal'])


if __name__ == '__main__':
    main()
