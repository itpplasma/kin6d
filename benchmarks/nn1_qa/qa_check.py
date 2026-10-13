#!/usr/bin/env python3
"""Independent numpy checks of the NN1 QA benchmark data (run with python3 -I).
Usage: qa_check.py DIR [cutoff]. Prints: toroidal flux constancy, vacuum residual
B.n/|B| of the published surface, minimum coil-surface distance, map Jacobian
and boundary truncation error for the coefficient cutoff.
"""
import sys
import numpy as np
MU0 = 4e-7 * np.pi


def load(d):
    with open(f'{d}/qa_boundary.dat') as f:
        nfp, nr = map(int, f.readline().split())
        tab = np.loadtxt(f)
    with open(f'{d}/qa_coils.dat') as f:
        nc, npt = map(int, f.readline().split())
        coils = []
        for _ in range(nc):
            i0 = float(f.readline())
            coils.append((i0, np.loadtxt([f.readline() for _ in range(npt)])))
    return nfp, tab, coils


def bfield(coils, x):
    """Biot-Savart of closed polygons at points x (N,3), exact per segment."""
    b = np.zeros_like(x, dtype=float)
    for i0, g in coils:
        a = g; c = np.roll(g, -1, axis=0)
        dl = c - a
        for s in range(len(a)):
            r1 = x - a[s]; r2 = x - c[s]
            n1 = np.linalg.norm(r1, axis=1); n2 = np.linalg.norm(r2, axis=1)
            cr = np.cross(dl[s], r1)
            den = n1 * n2 * (n1 * n2 + np.einsum('ij,ij->i', r1, r2))
            L = np.linalg.norm(dl[s])
            # field of a straight segment: mu0 I/(4pi) (n1+n2) (dl x r1)/(n1 n2 (n1 n2 + r1.r2))
            b += MU0 * i0 / (4 * np.pi) * ((n1 + n2) / den)[:, None] * cr
    return b


def surface(nfp, tab, th, ph, cutoff=0.0, nmax=99):
    R = np.zeros_like(th); Z = np.zeros_like(th); dR = np.zeros_like(th)
    dth = {'R': np.zeros_like(th), 'Z': np.zeros_like(th)}
    dph = {'R': np.zeros_like(th), 'Z': np.zeros_like(th)}
    for n, m, rb, zb in tab:
        n = int(n); m = int(m)
        if max(abs(rb), abs(zb)) < cutoff or abs(n) > nmax:
            continue
        a = m * th - nfp * n * ph
        R += rb * np.cos(a); Z += zb * np.sin(a)
        dth['R'] += -rb * m * np.sin(a); dth['Z'] += zb * m * np.cos(a)
        dph['R'] += rb * nfp * n * np.sin(a); dph['Z'] += -zb * nfp * n * np.cos(a)
    return R, Z, dth, dph


if __name__ == '__main__':
    d = sys.argv[1]
    cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    nfp, tab, coils = load(d)
    # toroidal flux F = (1/2pi) int R B_phi dphi at two points
    ph = np.linspace(0, 2 * np.pi, 96, endpoint=False)
    for (r, z) in [(1.0, 0.0), (1.08, 0.12)]:
        p = np.stack([r * np.cos(ph), r * np.sin(ph), z * np.ones_like(ph)], 1)
        b = bfield(coils, p)
        bphi = -np.sin(ph) * b[:, 0] + np.cos(ph) * b[:, 1]
        print(f'F at R={r} Z={z}: {np.mean(r * bphi):.8e}')
    # vacuum residual on the surface (full coefficient set)
    th, ph = np.meshgrid(np.linspace(0, 2 * np.pi, 48, endpoint=False),
                         np.linspace(0, 2 * np.pi / nfp, 48, endpoint=False), indexing='ij')
    R, Z, dth, dph = surface(nfp, tab, th, ph)
    x = np.stack([R * np.cos(ph), R * np.sin(ph), Z], -1)
    tth = np.stack([dth['R'] * np.cos(ph), dth['R'] * np.sin(ph), dth['Z']], -1)
    tph = np.stack([dph['R'] * np.cos(ph) - R * np.sin(ph), dph['R'] * np.sin(ph) + R * np.cos(ph), dph['Z']], -1)
    nrm = np.cross(tth, tph); area = np.linalg.norm(nrm, axis=-1); nrm /= area[..., None]
    b = bfield(coils, x.reshape(-1, 3)).reshape(x.shape)
    bn = np.einsum('...i,...i', b, nrm); bm = np.linalg.norm(b, axis=-1)
    print(f'surface: |B| mean {np.sum(bm*area)/np.sum(area):.4e} T; '
          f'<|B.n|>/<|B|> = {np.sum(abs(bn)*area)/np.sum(bm*area):.3e}; '
          f'rms(B.n)/rms(B) = {np.sqrt(np.sum(bn**2*area)/np.sum(bm**2*area)):.3e}; '
          f'max|B.n|/|B| = {np.max(abs(bn)/bm):.3e}')
    # min coil-surface distance
    xs = x.reshape(-1, 3)[::7]
    dmin = min(np.min(np.linalg.norm(xs[:, None, :] - g[None, ::2, :], axis=-1)) for _, g in coils)
    print(f'min coil-surface distance (polyline vertices) {dmin:.4f} m')
    # wall map of test_nonnested_qa.f90: harmonic extension of the VMEC boundary with the
    # fold-removing term -beta (1-|z|^2) m (ZBS Re w, RBC Im w); z in the unit disk
    beta = 0.1
    def mapRZ(z, p):
        R = 0; Z = 0
        for n, m, rb, zb in tab:
            w = z ** int(m) * np.exp(-1j * nfp * n * p); damp = 1 - np.abs(z) ** 2
            R = R + (rb - beta * damp * m * zb) * w.real
            Z = Z + (zb - beta * damp * m * rb) * w.imag
        return R, Z
    for b_ in (0.0, beta):
        beta_save = beta; beta = b_
        jmin = 1e9; jmax = 0.0
        uu, tt = np.meshgrid(np.linspace(0.0, 1, 61), np.linspace(0, 2 * np.pi, 121, endpoint=False), indexing='ij')
        z0 = uu * np.exp(1j * tt); e = 1e-6
        for p in np.linspace(0, 2 * np.pi / nfp, 17):
            R1, Z1 = mapRZ(z0 + e, p); R2, Z2 = mapRZ(z0 - e, p)
            R3, Z3 = mapRZ(z0 + 1j * e, p); R4, Z4 = mapRZ(z0 - 1j * e, p)
            j = ((R1 - R2) * (Z3 - Z4) - (R3 - R4) * (Z1 - Z2)) / (4 * e * e)
            jmin = min(jmin, j.min()); jmax = max(jmax, j.max())
        print(f'map beta={b_}: det d(R,Z)/d(x,y) on the unit disk in [{jmin:.4f}, {jmax:.4f}]')
        beta = beta_save
    # boundary truncation
    for nmax in (1, 2, 3, 4, 5):
        R0, Z0, _, _ = surface(nfp, tab, th, ph)
        R1, Z1, _, _ = surface(nfp, tab, th, ph, cutoff, nmax)
        print(f'n<={nmax}: max boundary deviation {np.max(np.hypot(R0-R1, Z0-Z1)):.3e} m')
