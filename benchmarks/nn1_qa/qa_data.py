#!/usr/bin/env python3
"""Generate the data files of the NN1 precise-QA vacuum benchmark from the
published simsopt files (input.LandremanPaul2021_QA, biot_savart_opt.json).

Usage: qa_data.py INPUT_QA BIOT_SAVART_JSON OUTDIR [npts_per_coil]
Writes OUTDIR/qa_boundary.dat (nfp, RBC/ZBS table) and qa_coils.dat (16 closed
polylines with current in A). The wall map is built in Fortran from the table.
Standard library + numpy only.
"""
import json, re, sys
import numpy as np


def read_boundary(path):
    nfp = None
    rows = {}
    for line in open(path):
        line = line.split('!')[0]
        m = re.search(r'NFP\s*=\s*(\d+)', line, re.I)
        if m:
            nfp = int(m.group(1))
        for mm in re.finditer(r'(RBC|ZBS)\(\s*(-?\d+),\s*(-?\d+)\)\s*=\s*([-+0-9.eEdD]+)', line):
            kind, n, m_, v = mm.groups()
            rows.setdefault((int(n), int(m_)), {})[kind] = float(v.replace('d', 'e').replace('D', 'e'))
    return nfp, rows


def fourier_curve(x, order, t):
    # simsopt CurveXYZFourier dof layout per coordinate: c0, s1, c1, s2, c2, ...
    x = np.asarray(x).reshape(3, 2 * order + 1)
    g = np.zeros((len(t), 3))
    for c in range(3):
        g[:, c] = x[c, 0]
        for k in range(1, order + 1):
            g[:, c] += x[c, 2 * k - 1] * np.sin(2 * np.pi * k * t) + x[c, 2 * k] * np.cos(2 * np.pi * k * t)
    return g


def coils(path, npts):
    d = json.load(open(path))['simsopt_objs']
    t = np.arange(npts) / npts
    base = {}
    for k, v in d.items():
        if isinstance(v, dict) and v.get('@class') == 'CurveXYZFourier':
            dofs = d[v['dofs']['value']]
            names = dofs['names']
            assert names[0].endswith('xc(0)') or 'xc(0)' in names[0], names[:3]
            assert 'xs(1)' in names[1] and 'xc(1)' in names[2], names[:3]
            base[v['@name']] = fourier_curve(dofs['x']['data'], v['order'], t)
    cur = {v['@name']: v['current'] for v in d.values()
           if isinstance(v, dict) and v.get('@class') == 'Current'}
    out = []
    for i in range(1, 5):  # base coils, as listed in the file: Coil1-4, then rotated copies
        pass
    names = sorted((v['@name'] for v in d.values() if isinstance(v, dict) and v.get('@class') == 'Coil'),
                   key=lambda s: int(s[4:]))
    for nm in names:
        coil = next(v for v in d.values() if isinstance(v, dict) and v.get('@name') == nm)
        cv = d[coil['curve']['value']] if coil['curve']['value'] in d else next(
            w for w in d.values() if isinstance(w, dict) and w.get('@name') == coil['curve']['value'])
        cr = coil['current']['value']
        crd = next(w for w in d.values() if isinstance(w, dict) and w.get('@name') == cr)
        # current
        sign = 1.0
        if crd['@class'] == 'ScaledCurrent':
            sign = crd['scale']
            crd = next(w for w in d.values() if isinstance(w, dict) and w.get('@name') == crd['current_to_scale']['value'])
        i0 = crd['current'] * sign
        # curve
        if cv['@class'] == 'RotatedCurve':
            g = base[cv['curve']['value']].copy()
            if cv['flip']:
                g[:, 1] *= -1; g[:, 2] *= -1
            p = cv['phi']
            r = np.array([[np.cos(p), -np.sin(p), 0], [np.sin(p), np.cos(p), 0], [0, 0, 1]])
            g = g @ r.T
        else:
            g = base[cv['@name']]
        out.append((i0, g))
    return out


if __name__ == '__main__':
    inp, js, outdir = sys.argv[1:4]
    npts = int(sys.argv[4]) if len(sys.argv) > 4 else 192
    nfp, rows = read_boundary(inp)
    with open(f'{outdir}/qa_boundary.dat', 'w') as f:
        f.write(f'{nfp} {len(rows)}\n')
        for (n, m), v in sorted(rows.items(), key=lambda kv: (kv[0][1], kv[0][0])):
            f.write(f'{n:4d} {m:3d} {v.get("RBC", 0.0):.15e} {v.get("ZBS", 0.0):.15e}\n')
    cs = coils(js, npts)
    with open(f'{outdir}/qa_coils.dat', 'w') as f:
        f.write(f'{len(cs)} {npts}\n')
        for i0, g in cs:
            f.write(f'{i0:.12e}\n')
            for p in g:
                f.write(f'{p[0]:.15e} {p[1]:.15e} {p[2]:.15e}\n')
