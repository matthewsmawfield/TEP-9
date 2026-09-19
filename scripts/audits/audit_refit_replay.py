#!/usr/bin/env python3
"""Deterministic stratified replay of cached real comets after fitter repairs.

Samples low/median/high observation counts in each of three cohorts. This is
regression validation, not a new detection test and not a full cohort refit.
"""
from pathlib import Path
import hashlib
import json
import multiprocessing as mp
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.utils import mpc_refit as R


def initialize():
    m=ROOT/'data/raw/mpc'
    R.configure(m/'obs',m/'sbdb_fp',m/'obscodes.json',m/'provenance.json',
                ROOT/'data/raw/spice/de440s.bsp',ROOT/'data/raw/naif/naif0012.tls','review_refit_replay')


def replay(item):
    cohort, old = item
    start=time.monotonic()
    new,err=R.fit_comet(old['des'])
    out=dict(cohort=cohort, designation=old['des'], n_obs=old['n_obs'],
             elapsed_seconds=time.monotonic()-start,error=err)
    if new:
        out.update(old_rms=old['fit_all_rms'],new_rms=new['fit_all_rms'],
                   old_drot=old['our_drot'],new_drot=new.get('drot'))
        out['unchanged_to_1e_9']=bool(np.isclose(out['old_rms'],out['new_rms'],rtol=0,atol=1e-9)
            and out['new_drot'] is not None and np.isclose(out['old_drot'],out['new_drot'],rtol=0,atol=1e-9))
    return out


def main():
    chosen=[]
    for cohort in ['step_b91_refit','step_b92_refit','step_b98_lpc_refit']:
        rows=[json.loads(s) for s in (ROOT/'results'/f'{cohort}.jsonl').read_text().splitlines()]
        rows=sorted([r for r in rows if 'our_drot' in r and 'fit_all_rms' in r],key=lambda r:r['n_obs'])
        # Quantiles avoid the exceptionally long-arc outliers while testing
        # different observation-count regimes in each cohort.
        for quantile in [.1,.5,.9]:
            chosen.append((cohort,rows[int(quantile*(len(rows)-1))]))
    with mp.get_context('spawn').Pool(3,initializer=initialize) as pool:
        result=[]
        for item in pool.imap_unordered(replay,chosen):
            result.append(item)
            print(json.dumps(item),flush=True)
    report=dict(scope='Nine stratified real-data regression replays; not a full cohort regeneration',
        source_sha256=hashlib.sha256(Path(R.__file__).read_bytes()).hexdigest(),
        n=len(result),n_unchanged=sum(r.get('unchanged_to_1e_9',False) for r in result),results=result)
    (ROOT/'results/audits/REFIT_REPLAY_20260919.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return any(r.get('error') for r in result)
if __name__=='__main__':
    sys.exit(main())
