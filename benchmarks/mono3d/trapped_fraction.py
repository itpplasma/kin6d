"""Trapped fraction f_t = 1 - (3/4) <b^2> int_0^{1/b_max} lam dlam/<sqrt(1 - lam b)>
of a surface grid written by surface_grid.py (Boozer weight 1/b^2)."""
import sys
import numpy as np
lines = open(sys.argv[1]).read().split('\n')
nt, nz, _ = map(int, lines[0].split())
b = np.array([float(x) for x in lines[2:2 + nt*nz]])
w = 1/b**2
avg = lambda f: (w*f).sum()/w.sum()
x, wq = np.polynomial.legendre.leggauss(400)
t, wt = (x + 1)/2, wq/2
lam = (1 - t**2)/b.max()
integral = sum(wt[k]*2*t[k]/b.max()*lam[k]/avg(np.sqrt(np.maximum(1 - lam[k]*b, 0)))
               for k in range(len(t)))
print(1 - 0.75*avg(b**2)*integral)
