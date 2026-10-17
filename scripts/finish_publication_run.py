"""Finish one authorized pipeline run with checks and local publication builds.

This is a one-off continuation, with no deployment or recurring schedule.
The generated PDF still requires visual review before release certification.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reviews/completion_run_20261014'


def status(stage, **details):
    record = dict(stage=stage, updated_utc=datetime.datetime.now(
        datetime.timezone.utc).isoformat(), **details)
    (OUT/'publication_completion_status.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wait-pid', type=int, required=True)
    parser.add_argument('--started-utc', required=True)
    args = parser.parse_args()
    status('waiting_for_pipeline', pipeline_pid=args.wait_pid)
    while True:
        try:
            os.kill(args.wait_pid, 0)
        except ProcessLookupError:
            break
        time.sleep(10)
    audits = []
    for path in (ROOT/'results/audits').glob('RESEARCH_AUDIT_*.json'):
        audit = json.loads(path.read_text())
        if audit.get('telemetry', {}).get('timestamp_utc', '') >= args.started_utc:
            audits.append((path, audit))
    if not audits:
        status('failed', reason='Pipeline ended without a new execution audit')
        return 1
    audits.sort(key=lambda pair: pair[1]['telemetry']['timestamp_utc'])
    path, audit = audits[-1]
    by_position = {}
    for source_path, source_audit in audits:
        for entry in source_audit.get('pipeline_results', []):
            by_position[entry['step']] = dict(entry, source_audit=str(source_path))
    entries = list(by_position.values())
    from scripts.audits.audit_publication import registered_steps
    expected = len(registered_steps())
    failed = [x for x in entries if x['status'] != 'SUCCESS']
    if failed or len(entries) != expected:
        status('pipeline_requires_repair', audit=str(path), failed=failed,
               attempted=len(entries), expected=expected)
        return 1
    (OUT/'completion_run_execution.json').write_text(json.dumps(
        dict(expected=expected, completed=len(entries), pipeline_results=entries), indent=2)+'\n')
    commands = [
        ('tests', ['python3', '-m', 'pytest', '-q', 'tests']),
        ('publication_bindings', ['python3', 'scripts/audits/audit_publication.py']),
        ('static_site', ['node', 'site/build.js']),
        ('markdown', ['node', 'site/html-to-markdown.js']),
        ('pdf', ['python3', 'scripts/generate_site_pdf.py', '--skip-build', '--quality', 'high', '--wait-time', '10']),
    ]
    for stage, command in commands:
        status(stage, command=command)
        result = subprocess.run(command, cwd=ROOT)
        if result.returncode:
            status('requires_repair', failed_stage=stage, returncode=result.returncode)
            return result.returncode
    status('artifacts_built_visual_review_pending', pipeline_audit=str(path))
    return 0


if __name__ == '__main__':
    import sys
    sys.path.insert(0, str(ROOT))
    raise SystemExit(main())
