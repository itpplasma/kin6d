# Matrix-free Newton-MR for the mapped NN1 lane

Status: candidate, opt-in. Item 1 of the
[relaxation comparators](relaxation-comparators.md#algorithms-to-adopt), first
slice. Sources `src/core/newton_mr.f90` (`kin6d_newton_mr`) and
`src/mhd/nonnested_newton_mr.f90` (`kin6d_nonnested_newton_mr`); tests
`test/test_newton_mr.f90` and `test/test_nonnested_newton_mr.f90`. The default
path is unchanged.

## Scope of this slice

The target is the energy minimisation of the mapped non-nested FEEC lane
([nn1-mapped-3d.md](nn1-mapped-3d.md)): with `h = F/R e_phi` and fixed traces,

    W(A) = int |h + curl A|^2 / (2 mu0),    S = C^T M2(K) C,   K = DF^T DF / J.

`nn_mapped_solve` assembles `S` with mode-dense element blocks
(`O((nl (2N+1))^2)` per element) and solves by FGMRES(100), right-preconditioned
by a factor of the assembled phi-averaged operator. `nn_mapped_solve_mr` takes
the same arguments plus `nn_mr_options_t`, returns the same
`nn_mapped_report_t` and never assembles `S`:

- **Hessian action.** FortFEM's sum-factorised `meridian_fourier_curl_curl_apply`
  (FFT in phi per meridional quadrature point, OpenMP over elements). It equals
  the assembled product to `8e-16` relative (test `agree`).
- **Newton-MR** (Liu and Roosta 2022). Each step solves `H d = -g` by
  preconditioned MINRES to the forcing `eta = min(1/2, (|g|/|g0|)^(1/2))`, the
  rule of the map lane's `krylov_step`, floored at what the stopping test needs.
  MINRES stops at the first residual of nonpositive curvature (the
  `c_(k-1) gamma1_k >= 0` test) and returns `M r`; such a step is
  forward-tracked, a solution step backtracked under the Armijo rule. Inside the
  roundoff band of `W` steps are ranked by the trapezoidal change, as in
  `kin6d_ideal_map`. The singular `n = 0` block (kernel `grad V0`) needs no
  special treatment: the right-hand side is consistent and MINRES stays in the
  range. `B = C A` is unique.
- **Preconditioners** (`options%preconditioner`): `NN_MR_PC_NONE`,
  `NN_MR_PC_JACOBI` (diagonal of the lumped operator) and `NN_MR_PC_LUMPED`
  (default), described next.

Callbacks are type-bound procedures of caller-owned objects, so the new code
needs no executable-stack trampolines.

## Metric lumping on this lane

MRX lumps the metric along each logical coordinate of a tensor-product spline
space, so each block becomes a Kronecker product or sum of 1-D matrices with a
dense block at the polar axis. Here only phi is a tensor direction:

- **phi (applies).** `K` is replaced by its phi mean `Kbar(R, Z)` at each
  meridional quadrature point (`lumped_metric_t`). The Fourier V2 mass of `Kbar`
  is exactly `blockdiag(M0, M0/2, M0/2, ...)`, so the lumped operator is the
  Kronecker sum of the meridional forms with the diagonal 1-D phi mass and the
  stiffness `n^2`. The fast diagonalisation in phi is the Fourier basis itself:
  mode `n` is one meridional system `C_n^T (M0/2 (+) M0/2) C_n` on its (cos, sin)
  pair, with `C_n = [[a, n b], [-n b, a]]`, built from FortFEM's per-n curl
  pieces and factored once by `nn_factor`. `n = 0` adds `eps M1` on the Nedelec
  block, as the default preconditioner does.
- **(R, Z) (does not apply).** The meridional mesh is unstructured triangles
  (P_k Lagrange/Nedelec/RT), so the in-plane forms are not Kronecker products
  even on the structured polar ring topology. The per-mode sparse factor is the
  exact inverse of the lumped model instead of 1-D fast diagonalisation.
  Triangles have no singular logical axis, so no dense axis block is needed.
- **Multigrid (does not apply).** `flux_mesh_build` meshes are not nested under
  refinement (new rings re-triangulate each band), so a geometric V-cycle on
  the curl-curl stiffness has no hierarchy. The multilevel route for this
  complex is the Hiptmair-Xu auxiliary-space preconditioner of FortFEM (open
  PR #102 for the per-mode NN1 solve), which could replace each per-mode factor.

The lumped preconditioner is spectrally equivalent to `S` with constants set by
`sup K / Kbar` and `sup Kbar / K`, independent of `h` and `N`; the identity map
makes it exact (one MINRES iteration).

## Evidence (this VM, 4 shared cores, `OMP_NUM_THREADS=2`, Release)

Rotating-ellipse vacuum oracle of `test_nonnested_mapped` (P2, `nphi = 4N + 4`).

Agreement with the factored path (`agree 4 4` is the CTest; `agree 2 3` shown):

| case | max\|b_MR - b_ref\|/max\|b\| | energy rel. diff | oracle error ref / MR | residual ref / MR |
|---|---:|---:|---:|---:|
| nr 2, N 3 | 2.6e-12 | 1.8e-16 | 1.047e-1 / 1.047e-1 | 4.1e-14 / 1.5e-13 |
| identity map vs per-n NN1 | 8.3e-14 | | | |

MINRES iterations of one solve to relative residual `1e-8` (`ladder pc N nr...`):

| N | nr (ndof) | lumped | Jacobi | none |
|---|---|---:|---:|---:|
| 4 | 2 (1345) | 11 | 656 | 1070 |
| 4 | 4 (4945) | 11 | 1307 | 2697 |
| 4 | 8 (18161) | 12 | 2404 | > 5000 (1.6e-6) |
| 4 | 16 (68657) | 12 | | |
| 4 | 32 (277937) | 12 | | |
| 2 / 8 / 12 | 4 | 10 / 12 / 12 | | |

The lumped preconditioner is mesh- and mode-independent; Jacobi and no
preconditioner grow like `1/h`. The full Newton-MR solve to `1e-12` takes 6-7
Newton steps and 16-19 MINRES iterations in total (FGMRES: 17-18), with no
NPC step on this convex problem.

Wall time and peak RSS of one solve (`path factored|mr nr N`, fresh process):

| nr, N | ndof | factored s / MiB | Newton-MR s / MiB | MR split: operator / PC / MINRES / evaluate s |
|---|---:|---:|---:|---|
| 4, 4 | 4945 | 0.72 / 126 | 0.53 / 42 | 0.02 / 0.15 / 0.25 / 0.07 |
| 8, 4 | 18161 | 3.29 / 408 | 1.41 / 82 | 0.09 / 0.33 / 0.51 / 0.22 |
| 16, 4 | 68657 | 17.2 / 1346 | 4.51 / 214 | 0.28 / 0.79 / 1.82 / 0.90 |
| 8, 10 | 41393 | 23.8 / 1875 | 4.88 / 123 | 0.17 / 0.55 / 2.28 / 1.22 |
| 16, 10 | 156401 | 101 / 6888 | 19.4 / 373 | 0.64 / 1.94 / 9.24 / 4.73 |

Oracle field errors agree to all printed digits in every row. Timings were
taken on a host shared with two other agents and are indicative only.

Unit checks (`test_newton_mr`): MINRES against a dense solve for SPD
(cond 1e3), diagonally preconditioned SPD (cond 1e9; recurred and explicit
residuals differ by about cond eps there, so that gate is 1e-10), indefinite
and singular consistent systems (minimum-norm solution to 4e-12); the NPC exit
returns a direction with negative curvature that descends. Newton-MR from a
saddle of a nonconvex quartic takes 5 NPC steps and ends at a minimiser with a
positive definite Hessian.

## Limits and open items

- The objective is quadratic, so this slice exercises the Newton-MR machinery
  (forcing, line search, NPC exit) only in its convex regime on the physics
  lane; the nonconvex path is covered by the unit test. The nonlinear lanes are
  the next consumers: the map lane's Newton (`ideal_map_solve`, whose
  `krylov_step` CG with a retained factor and LU fallbacks this is meant to
  replace) has an indefinite Hessian but is assembled per element without a
  matrix-free action; NN2's coupled Newton is nonsymmetric, so MINRES does not
  apply there (GMRES/FGMRES remain); the Beltrami closure is symmetric
  indefinite and is a direct MINRES candidate.
- MINRES belongs in FortNum (`fortnum_krylov` has CG/GMRES/FGMRES only) and
  should move there; no Lanczos reorthogonalisation is implemented (MRX keeps
  the basis; not needed at these iteration counts).
- The FortFEM apply allocates scratch per call and the curl products use the
  serial `csc_matvec`; both are upstream hot-path items.
- Claims are candidate: residuals and sampled oracle errors only, no
  certificate. Derivatives are unchanged in meaning: the converged state solves
  the same discrete equation, whose tangent/adjoint is the symmetric solve with
  `S`.

Reproduce:

```sh
cmake --build build --target test_newton_mr test_nonnested_newton_mr
ctest --test-dir build -R newton_mr
OMP_NUM_THREADS=2 build/test_nonnested_newton_mr ladder jacobi 4 2 4 8
OMP_NUM_THREADS=2 build/test_nonnested_newton_mr path factored 16 10
OMP_NUM_THREADS=2 build/test_nonnested_newton_mr path mr 16 10
```
