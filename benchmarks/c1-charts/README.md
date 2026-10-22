# C1: optional field-aligned charts from a nested solution (IH2 pilot)

[IH2](../../docs/island-equilibrium.md) allows curve-following charts built
from a nested solution as an optional accelerator. The unaligned axis-regular
path stays complete and gated. This pilot measures whether such a chart pays
off for anisotropic transport on
[FortFEM](https://github.com/lazy-fortran/fortfem)'s MF1 matrix-free tensor
operator (lazy-fortran/fortfem#110, pinned at its head 52cd3be).

- `fit_wout.py` fits a stellarator-symmetric VMEC/VMEC++ wout to an
  axis-regular map: each mode becomes rho^m P(s), with P a polynomial in s.
  The fitted map, not the wout, is the test geometry.
- `nested_chart.f90` evaluates that map (jets, inverse, b, Laplacian of s) and
  builds MF1 mesh arrays (vertices, isoparametric geometry, outer-face flags)
  for a five-block O-grid times zeta over one field period. In the aligned
  chart the ring blocks follow the nested surfaces, while in the unaligned
  chart they follow scaled copies of the boundary. Both charts use the same
  core Coons patch, so the axis is an ordinary interior point.
- `c1_charts.f90` solves -div(K grad u) = f with
  K = chi_perp I + (chi_par - chi_perp) b b^T and u = cos(pi s/2). Since
  b . grad s = 0, f does not depend on chi_par and every chi_par dependence of
  the error is pollution.

Build and test (standalone, KIN6D's main FortFEM pin is untouched):

    cmake -S benchmarks/c1-charts -B build-c1 -DCMAKE_BUILD_TYPE=Release \
        -DKIN6D_C1_MAP=/path/map.txt
    cmake --build build-c1 -j4 && ctest --test-dir build-c1 -L c1
    build-c1/c1_charts map.txt aligned P N NPHI CHI_RATIO [TOL] [RHO_CORE]

Remote runs use `c1.sbatch`.
