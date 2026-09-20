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


_hash_cache = {}

def digest(path):
    if path not in _hash_cache:
        _hash_cache[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    return _hash_cache[path]


def input_record(root, rel, result_mtime):
    """Fingerprint one declared input and cross-check sibling provenance.json."""
    if '*' in rel or '(' in rel:
        return dict(input=rel, annotated=True)  # descriptive entry, not a literal path
    ip = (root / rel).resolve()
    if not ip.is_relative_to(root) or not ip.exists():
        return dict(input=rel, missing=not ip.exists())
    rec = dict(input=rel, sha256=digest(ip))
    if ip.stat().st_mtime > result_mtime:
        rec['newer_than_result'] = True
    prov = ip.parent / 'provenance.json'
    if prov.exists():
        try:
            entry = json.loads(prov.read_text()).get('files', {}).get(ip.name, {})
            rec['retrieved_utc'] = entry.get('retrieved_utc')
            rec['url'] = entry.get('url')
            if entry.get('sha256'):
                rec['matches_download_hash'] = entry['sha256'] == rec['sha256']
        except (ValueError, OSError):
            pass
    return rec


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
    # Caption numbering must be strictly sequential with no gaps or
    # duplicates, and every in-text Figure/Table reference must resolve.
    # Lettered sub-captions (Table 4a, 4b) occupy one numeric slot and must
    # carry consecutive suffixes beginning at 'a'.
    caption_counts = {'Figure': {}, 'Table': {}}
    ref_targets = {'Figure': set(), 'Table': set()}
    for comp in sorted((ROOT/'site/components').glob('*.html')):
        source = comp.read_text()
        for kind, num in re.findall(
                r'<figcaption>\s*(?:<[^>]+>)*\s*(Figure|Table)\s+'
                r'(\d+[a-z]?)', source):
            caption_counts[kind].setdefault(num, []).append(comp.name)
        for num in re.findall(
                r'<caption>\s*(?:<[^>]+>)*\s*Table\s+(\d+[a-z]?)', source):
            caption_counts['Table'].setdefault(num, []).append(comp.name)
        stripped = re.sub(r'<figcaption>.*?</figcaption>', '', source,
                          flags=re.S)
        stripped = re.sub(r'<caption>.*?</caption>', '', stripped,
                          flags=re.S)
        for kind, num in re.findall(r'\b(Figure|Table)s?\s+(\d+[a-z]?)',
                                    stripped):
            ref_targets[kind].add(num)
        # capture range endpoints: "Tables 2–7", "Figures 4-6"
        for kind, num in re.findall(r'\b(Figure|Table)s\s+\d+[\u2013-](\d+)',
                                    stripped):
            ref_targets[kind].add(num)
    for kind, counts in caption_counts.items():
        for num, files in counts.items():
            if len(files) > 1:
                errors.append(f'Duplicate {kind} {num} caption: {files}')
        by_base = {}
        for label in counts:
            m = re.fullmatch(r'(\d+)([a-z]?)', label)
            by_base.setdefault(int(m.group(1)), []).append(m.group(2))
        bases = sorted(by_base)
        if bases and bases != list(range(1, len(bases) + 1)):
            errors.append(f'{kind} captions not sequential 1..N: {bases}')
        for base, suffixes in by_base.items():
            want = [''] if len(suffixes) == 1 else \
                [chr(97 + i) for i in range(len(suffixes))]
            if sorted(suffixes) != want:
                errors.append(
                    f'{kind} {base} sub-captions not sequential: {suffixes}')
        for label in sorted(ref_targets[kind]):
            if label in counts:
                continue
            # a bare numeric ref ("Table 4") resolves to a sub-captioned base
            if re.fullmatch(r'\d+', label) and int(label) in by_base:
                continue
            errors.append(f'{kind} {label} referenced but no caption')
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
    stale_inputs = []
    provenance_mismatches = []
    for path in sorted((ROOT/'results').glob('step_*.json')):
        if path.name == 'step_b78_claims_trace.json':
            continue  # never hash the audit's own superseded report
        try:
            result = json.loads(path.read_text())
            nonfinite = list(nonfinite_leaves(result))
            producer = result.get('step') if isinstance(result, dict) else None
            source = ROOT/'scripts/steps'/f'{producer}.py' if producer else None
            inputs = [input_record(ROOT, rel, path.stat().st_mtime)
                      for rel in (result.get('inputs') or [])
                      if isinstance(rel, str)] if isinstance(result, dict) else []
            for rec in inputs:
                if rec.get('newer_than_result'):
                    stale_inputs.append(dict(result=path.name, input=rec['input']))
                if rec.get('matches_download_hash') is False:
                    provenance_mismatches.append(dict(result=path.name, input=rec['input']))
            result_inventory.append(dict(file=path.name, sha256=digest(path),
                nonfinite_count=len(nonfinite), nonfinite_paths=nonfinite,
                declared_producer=producer,
                inputs=inputs,
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
            legacy_files_with_nonfinite_values=sum(x.get('nonfinite_count',0)>0 for x in result_inventory),
            results_with_newer_inputs=len(stale_inputs),
            provenance_hash_mismatches=len(provenance_mismatches)),
        errors=errors, unused_claims=sorted(set(registry)-used), bindings=bindings,
        stale_inputs=stale_inputs, provenance_mismatches=provenance_mismatches,
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
    for warning in report['stale_inputs']:
        print('STALE-INPUT WARNING:', warning['result'], '<-', warning['input'])
    for mismatch in report['provenance_mismatches']:
        print('PROVENANCE MISMATCH:', mismatch['result'], '<-', mismatch['input'])
    for error in report['errors']:
        print(error)
    return bool(report['errors'])


if __name__ == '__main__':
    sys.exit(main())
