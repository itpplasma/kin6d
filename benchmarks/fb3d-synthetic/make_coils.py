"""Synthetic public nfp = 3 modular coils for the KIN6D 3D free-boundary case.

Planar circular coils of radius RC centred on the helical axis
    R(phi) = R0 + EPS cos(NFP phi),  Z(phi) = EPS sin(NFP phi)
in the plane normal to the axis tangent, NCOIL coils per field period at
phi_j = 2 pi (j + 1/2) / (NFP NCOIL) (stellarator symmetric). All coils form
one circuit (MAKEGRID group 1) carrying CURRENT ampere-turns each, oriented so
that the vacuum field points along +phi. Written in MAKEGRID format with the
winding current in column 4, read by VMEC++ (makegrid) and by KIN6D through
libneo's neo_makegrid_coils.
"""
import argparse
import numpy as np

R0, EPS, NFP, NCOIL, RC, NPOINT = 1.0, 0.15, 3, 8, 0.30, 96
B0 = 1.0
MU0 = 4e-7 * np.pi


def axis(phi):
    r = R0 + EPS * np.cos(NFP * phi)
    z = EPS * np.sin(NFP * phi)
    x = np.array([r * np.cos(phi), r * np.sin(phi), z])
    dr = -NFP * EPS * np.sin(NFP * phi)
    dz = NFP * EPS * np.cos(NFP * phi)
    t = np.array([dr * np.cos(phi) - r * np.sin(phi), dr * np.sin(phi) + r * np.cos(phi), dz])
    return x, t / np.linalg.norm(t)


def coils(ncoil=NCOIL):
    # ncoil > NCOIL places the same winding on more coils with the same total
    # current (vacuum_reference.py --subdivide); the written file uses NCOIL.
    ntotal = NFP * ncoil
    current = 2 * np.pi * R0 * B0 / (MU0 * ntotal)
    out = []
    for j in range(ntotal):
        phi = 2 * np.pi * (j + 0.5) / ntotal
        c, t = axis(phi)
        ez = np.array([0.0, 0.0, 1.0])
        n = ez - ez.dot(t) * t
        n /= np.linalg.norm(n)
        b = np.cross(t, n)
        alpha = 2 * np.pi * np.arange(NPOINT) / NPOINT
        # (n, b, t) right-handed: increasing alpha circulates right-handed about t.
        pts = c[:, None] + RC * (np.cos(alpha)[None, :] * n[:, None]
                                 + np.sin(alpha)[None, :] * b[:, None])
        assert np.allclose(np.cross(n, b), t)
        out.append(pts.T)
    return out, current


def write(path, polys, current):
    with open(path, "w") as f:
        f.write("periods %d\nbegin filament\nmirror NIL\n" % NFP)
        for pts in polys:
            for p in pts:
                f.write("% .15e % .15e % .15e % .15e\n" % (p[0], p[1], p[2], current))
            p = pts[0]
            f.write("% .15e % .15e % .15e % .15e 1 modular\n" % (p[0], p[1], p[2], 0.0))
        f.write("end\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("output")
    a = ap.parse_args()
    polys, current = coils()
    write(a.output, polys, current)
    print("coils", len(polys), "points", NPOINT, "current per coil %.6e A" % current)
