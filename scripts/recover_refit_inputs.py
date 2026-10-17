"""Recover raw MPC observations and SBDB seeds for the archived refit cohort.

Downloads inputs only. The registered pipeline performs the fresh fits.
"""
import concurrent.futures
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.utils import mpc_refit as refit


def main():
    out = ROOT/'reviews/completion_run_20261014'
    retired = out/'retired_checkpoints'
    names = set()
    for filename in ('step_b91_refit.jsonl', 'step_b92_refit.jsonl', 'step_b98_lpc_refit.jsonl'):
        for line in (retired/filename).read_text().splitlines():
            record = json.loads(line)
            if not record.get('failed'):
                names.add(record.get('des') or record['desig'].split(' (')[0].strip())
    mpc = ROOT/'data/raw/mpc'
    refit.OBS_DIR = mpc/'obs'
    refit.SBDB_DIR = mpc/'sbdb_fp'
    refit.PROV_PATH = mpc/'provenance.json'
    for folder in (refit.OBS_DIR, refit.SBDB_DIR):
        folder.mkdir(exist_ok=True)
    start = time.monotonic()
    errors = []
    def recover(name):
        try:
            obs = refit.get_obs(name)
            if not obs:
                raise ValueError('No raw astrometry returned')
            refit.get_sbdb(name)
            return name, len(obs), None
        except Exception as exc:
            return name, None, str(exc)
    print(f'Recovering raw inputs for {len(names)} comets', flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for count, (name, nobs, error) in enumerate(pool.map(recover, sorted(names)), 1):
            if error:
                errors.append(dict(designation=name, error=error))
                print('FAILED', name, error, flush=True)
            if count % 25 == 0 or count == len(names):
                status = dict(completed=count, total=len(names), errors=errors,
                              elapsed_seconds=time.monotonic()-start)
                (out/'refit_input_recovery.json').write_text(json.dumps(status, indent=2)+'\n')
                print(f'{count}/{len(names)} recovered; {len(errors)} errors; {status["elapsed_seconds"]:.1f}s', flush=True)
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
