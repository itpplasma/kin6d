"""Behavioral external-file adapter/Fortran CLI tests; no fitted calibration."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'adapter', ROOT / 'scripts/prepare_response_qualification.py')
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
EXE = Path(sys.argv.pop(1)).resolve()


class AdapterTest(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.root = Path(self.scratch.name)
        (self.root / 'references').mkdir()
        for path in (ROOT / 'test/fixtures/response_gpec').glob('*.out'):
            shutil.copy(path, self.root / 'references' / path.name)
        (self.root / 'references/mars_comparison.json').write_text('{}')
        self.cases = []
        for c, cut in enumerate(adapter.CUTS, 1):
            for axis in adapter.AXES:
                for l in range(3):
                    name = f'{cut}_{axis}_{l}'
                    self.cases.append({'name': name, 'cut': cut,
                                       'variant': f'{axis}_{l}',
                                       'equilibrium': 'inputs/equilibrium',
                                       'controls': {'nrad': 8, 'nth_max': 16,
                                                    'np_sample': 64}})
                    folder = self.root / 'runs' / name
                    folder.mkdir(parents=True)
                    text = '# synthetic legacy single-ring physics data\n'
                    for p in (1, 2):
                        for t, q in enumerate(adapter.TARGETS, 1):
                            # signed q, area, peak vacuum flux/field, lambda,
                            # shielding peak field, L3 vacuum check columns.
                            values = [p, -q, .9, -round(3*q), 1, 1, 0, 1, 0, 0, 0,
                                      .001*c, .002*p + .0003*t, 1, 0]
                            text += ' '.join(map(str, values)) + '\n'
                    (folder / f'{name}_res.dat').write_text(text)
        (self.root / 'ladder.json').write_text(json.dumps({'cases': self.cases}))

    def tearDown(self):
        self.scratch.cleanup()

    def test_all_six_raw_references_and_single_ring_no_verdict(self):
        output = self.root / 'qualification.nml'
        adapter.prepare(self.root, self.root / 'runs', output, None)
        result = subprocess.run([str(EXE), '--config', str(output)], cwd=self.root,
                                text=True, capture_output=True, check=True)
        observed = {}
        for line in result.stdout.splitlines():
            fields = line.split()
            if len(fields) == 7 and fields[0].isdigit():
                c, p, t, gates = map(int, fields[:4])
                observed[c, p, t] = float(fields[4])
                self.assertEqual(gates, 0)
                self.assertEqual(int(fields[6]), 0)
        for c, cut in enumerate(adapter.CUTS, 1):
            for p, ph in enumerate(adapter.PHASES, 1):
                path = self.root / 'references' / f'gpec_{cut}_{ph}.out'
                # Independent fixed-format reference reader, separate from adapter.
                physical = [line.split() for line in path.read_text().splitlines()
                            if line.startswith('  3.000') or line.startswith('  3.333')]
                for t, fields in enumerate(physical, 1):
                    reference = complex(float(fields[3]), float(fields[4]))
                    b = complex(.001*c, .002*p + .0003*t)
                    self.assertAlmostEqual(observed[c, p, t] / (abs(b)/abs(reference)),
                                           1.0, delta=1e-5)
        self.assertIn('STOP_NO_TIEBREAK: 12/12', result.stdout)

    def test_duplicate_reference_refused(self):
        path = self.root / 'references/gpec_full_000.out'
        row = next(line for line in path.read_text().splitlines()
                   if line.startswith('  3.000'))
        path.write_text(path.read_text() + row + '\n')
        with self.assertRaisesRegex(ValueError, 'ambiguous/missing'):
            adapter.prepare(self.root, self.root / 'runs', self.root / 'input.nml', None)

    def test_changed_controls_refused(self):
        case = self.cases[0]
        path = self.root / 'runs' / case['name'] / 'qualification.json'
        path.write_text(json.dumps({'controls': {'nrad': 99}}))
        with self.assertRaisesRegex(ValueError, 'controls differ'):
            adapter.prepare(self.root, self.root / 'runs', self.root / 'input.nml', None)

    def test_fortran_rejects_duplicate_ladder_member(self):
        path = self.root / 'duplicate.nml'
        path.write_text('&qualification nrows=2, cut=1,1, phase=1,1, target=1,1, '
                        'axis=1,1, level=1,1 /\n')
        result = subprocess.run([str(EXE), '--config', str(path)], cwd=self.root,
                                text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('duplicate row', result.stderr)


if __name__ == '__main__':
    unittest.main()
