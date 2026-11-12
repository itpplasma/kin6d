"""Synthetic coil-box case for test_coil_box: 2 field periods, 16 planar
circular coils (radius 1.2 m) about R0 = 3 m, tilted alternately by +-0.25 rad
about the radial direction so the vacuum field has period-2 helical content,
and a rotating-ellipse wall (axes 0.5 and 0.35 m) well inside the coils.
Writes coils.box (MAKEGRID format, one circuit) and wall.txt (n m rc zs for
cos/sin(m theta - n P phi)). Usage: python3 make_case.py"""
import numpy as np

r0, rc, ncoil, npt, tilt, period = 3.0, 1.2, 16, 48, 0.25, 2
with open("coils.box", "w") as f:
    f.write("periods 2\nbegin filament\nmirror NIL\n")
    for k in range(ncoil):
        ph = 2 * np.pi * (k + 0.5) / ncoil
        a = tilt * (-1) ** k
        t = 2 * np.pi * np.arange(npt + 1) / npt
        # circle in the (e_R, e_Z) plane, tilted about e_R by a
        er = np.array([np.cos(ph), np.sin(ph), 0.0]); ez = np.array([0.0, 0.0, 1.0])
        ep = np.array([-np.sin(ph), np.cos(ph), 0.0])
        ezt = np.cos(a) * ez + np.sin(a) * ep
        pts = r0 * er[None] + rc * (np.cos(t)[:, None] * er[None] + np.sin(t)[:, None] * ezt[None])
        for i, p in enumerate(pts):
            if i < npt:
                f.write(f"{p[0]: .14e} {p[1]: .14e} {p[2]: .14e} {1.0: .6e}\n")
            else:
                f.write(f"{p[0]: .14e} {p[1]: .14e} {p[2]: .14e} {0.0: .6e} 1 COIL\n")
    f.write("end\n")
# rotating ellipse: R = r0 + a cos(th) cos(P phi/2)... written exactly in (m, n)
# R - r0 + i Z = e^{i th} (a+b)/2 + e^{-i th + i P phi} (a-b)/2
a, b = 0.5, 0.35
rows = [(0, 0, r0, 0.0), (0, 1, (a + b) / 2, (a + b) / 2)]
# second term: (a-b)/2 [cos(-th + P phi) + i sin(-th + P phi)]: m = 1, n = -1 in cos(m th - n P phi)
rows.append((-1, 1, (a - b) / 2, -(a - b) / 2))
with open("wall.txt", "w") as f:
    for n, m, r, z in rows:
        f.write(f"{n:4d} {m:4d} {r: .12e} {z: .12e}\n")
