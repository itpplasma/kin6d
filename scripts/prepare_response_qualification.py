#!/usr/bin/env python3
"""Optional file adapter for the Fortran-owned TC24 qualification contract.

Read the existing frozen ladder and raw references. Never execute or submit
cases. Existing single-ring *_res.dat can be inspected, but cannot acquire
coupled/calibrated status. Coupled records supply explicit diagnostics via
<case>/qualification.json; all scientific gates run in response_qualify.
"""
import argparse
import json
import math
from pathlib import Path

AXES = ('radial', 'core', 'grading', 'harmonic', 'poloidal', 'extent',
        'resistivity', 'frequency')
CUTS = ('c0995', 'c0997', 'full')
PHASES = ('000', '010')
TARGETS = (3, 10 / 3)


def number(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('nonfinite numeric input')
    return f'{value:.17e}'


def complex_value(value):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError('complex values must be [real, imaginary]')
    return f'({number(value[0])},{number(value[1])})'


def raw_gpec(path, q):
    candidates = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 7:
            continue
        try:
            values = [float(x.replace('D', 'E')) for x in fields]
        except ValueError:
            continue
        if abs(values[0] - q) < 0.001:
            if abs(values[2] - 5e-4) > 1e-12:
                raise ValueError(f'{path}: unexpected GPEC spot')
            candidates.append(values[3:5])
    if len(candidates) != 1:
        raise ValueError(f'{path}: ambiguous/missing q={q}')
    return candidates[0]


def single_ring(path, phase, q):
    if not path.exists():
        return None
    matches = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        values = [float(x.replace('D', 'E')) for x in line.split()]
        if len(values) < 13:
            raise ValueError(f'{path}: malformed response record')
        if int(values[0]) == phase and abs(abs(values[1]) - q) < 1e-4:
            # Legacy L4 columns: iph,q,psi,m,area,Phi_vac(2),b_vac(2),
            # lambda(2),b_shield(2). No jump fit is imported here.
            matches.append({'b': values[11:13], 'b_vac': values[7:9]})
    if len(matches) > 1:
        raise ValueError(f'{path}: duplicate response record')
    return matches[0] if matches else None


def prepare(root, results, output, decks):
    ladder = json.loads((root / 'ladder.json').read_text())
    cases = ladder['cases']
    identities = {(c['cut'], c['variant']) for c in cases}
    expected = {(c, f'{a}_{i}') for c in CUTS for a in AXES for i in range(3)}
    if identities != expected or len(cases) != len(expected):
        raise ValueError('input does not describe the frozen 72-case experiment')
    lines = ['&qualification', f' nrows={len(cases) * 4}']
    reference_summary = []
    for c, cut in enumerate(CUTS, 1):
        for p, ph in enumerate(PHASES, 1):
            path = root / 'references' / f'gpec_{cut}_{ph}.out'
            for t, q in enumerate(TARGETS, 1):
                ref = raw_gpec(path, q)
                lines.append(f' gpec({t},{p},{c})={complex_value(ref)}')
                reference_summary.append({'cut': cut, 'phase': ph, 'q': q,
                                          'singflx_T': ref, 'source': str(path)})
    mars = json.loads((root / 'references/mars_comparison.json').read_text())
    row = 0
    seen = set()
    for case in cases:
        name = case['name']
        cut = CUTS.index(case['cut']) + 1
        axis, level = case['variant'].rsplit('_', 1)
        axis = AXES.index(axis) + 1
        level = int(level) + 1
        diagnostics = results / name / 'qualification.json'
        records = []
        if diagnostics.exists():
            payload = json.loads(diagnostics.read_text())
            if payload['controls'] != case['controls']:
                raise ValueError(f'{name}: controls differ from frozen ladder')
            if payload['cut'] != case['cut'] or payload['spot'] != 5e-4:
                raise ValueError(f'{name}: changed cut or spot')
            records = payload['records']
        for p, ph in enumerate(PHASES, 1):
            for t, q in enumerate(TARGETS, 1):
                row += 1
                identity = (cut, axis, level, p, t)
                if identity in seen:
                    raise ValueError('duplicate ladder identity')
                seen.add(identity)
                for key, value in zip(('cut', 'axis', 'level', 'phase', 'target'), identity):
                    lines.append(f' {key}({row})={value}')
                qlabel = '3' if t == 1 else '10/3'
                key = f'{case["cut"]}_phi{ph}_q{qlabel}'
                comparison = mars.get('comparison', {}).get(key, {}).get('ideal', {})
                amplitude = comparison.get('singflx_cyl_estimate_T')
                if amplitude is not None:
                    lines.append(f' mars_amplitude({t},{p},{cut})={number(amplitude)}')
                # Admission is always held until the controller supplies
                # independent MARS jump/unit-sheet calibration evidence.
                matches = [r for r in records if r['phase'] == ph and
                           abs(r['q'] - q) < 1e-6]
                if len(matches) > 1:
                    raise ValueError(f'{name}: duplicate coupled record')
                if matches:
                    record = matches[0]
                    for key in ('b', 'b_vac', 'total_flux', 'current'):
                        lines.append(f' {key}({row})={complex_value(record[key])}')
                    for key in ('residual', 'balance'):
                        lines.append(f' {key}({row})={number(record[key])}')
                    width = record['extent']
                    if len(width) != 2:
                        raise ValueError('extent must give actual inner/outer distances')
                    lines.append(f' extent(:,{row})={number(width[0])},{number(width[1])}')
                    for key in ('available', 'coupled', 'calibrated', 'full_pressure',
                                'all_rationals', 'interface_test', 'whole_edge'):
                        if type(record[key]) is not bool:
                            raise ValueError(f'{name}: {key} must be boolean')
                        lines.append(f' {key}({row})=.{str(record[key]).lower()}.')
                else:
                    legacy = single_ring(results / name / f'{name}_res.dat', p, q)
                    if legacy:
                        lines.append(f' available({row})=.true.')
                        for key in ('b', 'b_vac'):
                            lines.append(f' {key}({row})={complex_value(legacy[key])}')
                if decks:
                    deck = decks / f'{name}_q{t}.nml'
                    # This validates the unit-sheet readback at the requested
                    # outer mesh. It is not an invented coupled solver deck.
                    if p == 1:
                        ctrl = case['controls']
                        eq = (root / case['equilibrium']).resolve()
                        deck.write_text('&spot_calibration\n'
                                        f" eqdsk='{str(eq).replace(chr(39), chr(39)*2)}'\n"
                                        f' n=3, qabs={number(q)}, spot=5.0e-4\n'
                                        f' nrad={ctrl["nrad"]}, nth_max={ctrl["nth_max"]}\n'
                                        f' nfine={ctrl["np_sample"]}, order=2\n/\n')
    lines.append('/')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('\n'.join(lines) + '\n')
    return reference_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--results', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--calibration-decks', type=Path)
    args = parser.parse_args()
    if args.calibration_decks:
        args.calibration_decks.mkdir(parents=True, exist_ok=True)
    refs = prepare(args.root, args.results or args.root / 'runs', args.output,
                   args.calibration_decks)
    args.output.with_suffix('.references.json').write_text(json.dumps(refs, indent=2) + '\n')
    print(f'Prepared frozen ladder input: {args.output}; launch readiness remains held')


if __name__ == '__main__':
    main()
