#!/usr/bin/env python3
"""Exact publication bindings and repository inventory; no scientific certification.

Unlike numeric-pool matching, a binding names a specific result and JSON path.
The inventory records nonfinite legacy leaves and missing provenance explicitly.
It does not infer producers from arbitrary filenames mentioned by consumers.
"""
import ast
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nonfinite_leaves(value, path=''):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from nonfinite_leaves(item, f'{path}/{key}')
    elif isinstance(value, list):
        for key, item in enumerate(value):
            yield from nonfinite_leaves(item, f'{path}/{key}')
    elif isinstance(value, float) and not math.isfinite(value):
        yield path


def registered_steps():
    tree = ast.parse((ROOT/'scripts/run_all.py').read_text())
    for node in tree.body:
        if isinstance(node, ast.AnnAssign) and getattr(node.target, 'id', None) == 'CORE_STEPS':
            return ast.literal_eval(node.value)
    raise ValueError('CORE_STEPS missing')


def audit():
    errors = []
    registered = registered_steps()
    script_inventory = []
    for path in sorted([*(ROOT/'scripts').rglob('*.py'), *(ROOT/'core').rglob('*.py')]):
        try:
            ast.parse(path.read_text())
        except SyntaxError as exc:
            errors.append(f'{path.relative_to(ROOT)}: {exc}')
        script_inventory.append({'file': str(path.relative_to(ROOT)), 'sha256': digest(path)})
    for name, _ in registered:
        if not (ROOT/'scripts/steps'/name).exists():
            errors.append(f'Registered step missing: {name}')
    if len({name for name, _ in registered}) != len(registered):
        errors.append('Duplicate registered steps')
    registry = json.loads((ROOT/'site/claims.json').read_text())
    used, bindings = set(), []
    for comp in sorted((ROOT/'site/components').glob('*.html')):
        source = comp.read_text()
        for ident in re.findall(r'\{\{claim:([a-z0-9_]+)\}\}', source):
            used.add(ident)
            try:
                spec = registry[ident]
                result = (ROOT/'results'/spec['file']).resolve()
                if not result.is_relative_to(ROOT/'results'):
                    raise ValueError('result outside results directory')
                value = json.loads(result.read_text())
                for key in spec['path']:
                    value = value[int(key)] if isinstance(value, list) else value[key]
                if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
                    raise ValueError(f'invalid numeric value {value}')
                if spec.get('probability') and not 0 < value <= 1:
                    raise ValueError(f'invalid probability {value}')
                bindings.append(dict(claim=ident, component=comp.name,
                                     result=spec['file'], path=spec['path'],
                                     value=value, source_sha256=digest(result)))
            except (KeyError, ValueError, IndexError, OSError, TypeError) as exc:
                errors.append(f'{comp.name}: {ident}: {exc}')
        for name in re.findall(r'src="figures/([^"?]+)"', source):
            if not (ROOT/'results/figures'/name).exists():
                errors.append(f'Missing figure: {comp.name}: {name}')
    # Exercise the actual JS renderer, including strict JSON handling.
    check = subprocess.run(['node', '-e', r'''const fs=require('fs');
const {renderClaims}=require('./site/claims');
for (const f of fs.readdirSync('site/components').filter(f=>f.endsWith('.html'))) {
 const text=renderClaims(fs.readFileSync('site/components/'+f,'utf8'));
 if (/\{\{/.test(text)) throw Error('Unresolved placeholder '+f);
}'''], cwd=ROOT, text=True, capture_output=True)
    if check.returncode:
        errors.append('Publication renderer: '+check.stderr[-1800:])
    result_inventory = []
    for path in sorted((ROOT/'results').glob('step_*.json')):
        if path.name == 'step_b78_claims_trace.json':
            continue  # never hash the audit's own superseded report
        try:
            result = json.loads(path.read_text())
            nonfinite = list(nonfinite_leaves(result))
            producer = result.get('step') if isinstance(result, dict) else None
            source = ROOT/'scripts/steps'/f'{producer}.py' if producer else None
            result_inventory.append(dict(file=path.name, sha256=digest(path),
                nonfinite_count=len(nonfinite), nonfinite_paths=nonfinite,
                declared_producer=producer,
                older_than_declared_producer=(path.stat().st_mtime < source.stat().st_mtime)
                    if source and source.exists() else None))
        except (ValueError, OSError) as exc:
            result_inventory.append(dict(file=path.name, parse_error=str(exc)))
    packages = {}
    for package in ['numpy','scipy','astropy','matplotlib','rebound','spiceypy','requests','pytest','PyYAML','playwright']:
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    return dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
        scope='Repository syntax/inventory and exact publication bindings. Not a complete numerical rerun or a proof of TEP.',
        summary=dict(registered_steps=len(registered), python_sources=len(script_inventory),
            bound_claim_occurrences=len(bindings), unique_bound_claims=len(used),
            result_files=len(result_inventory), hard_errors=len(errors),
            legacy_files_with_nonfinite_values=sum(x.get('nonfinite_count',0)>0 for x in result_inventory)),
        errors=errors, unused_claims=sorted(set(registry)-used), bindings=bindings,
        script_inventory=script_inventory, result_inventory=result_inventory,
        runtime=dict(python=sys.version, executable=sys.executable, packages=packages),
        provenance_limit='mtimes are hints only; the producer field is self-reported. Utilities, data and external dependencies require hashes captured when each result is produced.')


def main():
    report = audit()
    target = ROOT/'results/step_b78_claims_trace.json'
    target.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    # Do not retain a stale old CSV with numeric-pool matching described as provenance.
    import csv
    with (ROOT/'results/step_b78_claims_trace.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['claim','component','result','path','value','source_sha256'])
        writer.writeheader()
        for row in report['bindings']:
            writer.writerow({**row, 'path': '/'.join(row['path'])})
    print(json.dumps(report['summary'], indent=2))
    for error in report['errors']:
        print(error)
    return bool(report['errors'])


if __name__ == '__main__':
    sys.exit(main())
