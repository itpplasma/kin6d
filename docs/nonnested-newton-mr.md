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

# Taylor/Beltrami closure: MINRES on the saddle system

Status: candidate, opt-in. Second slice: the symmetric indefinite Taylor
closure of [equilibrium-beltrami.md](equilibrium-beltrami.md). Source
`src/mhd/nonnested_beltrami_mr.f90` (submodule of `kin6d_nonnested_beltrami`),
tests `test/test_nonnested_beltrami_mr.f90`, the `mr` variants of
`test_nonnested_beltrami` and `test_newton_mr saddle`. Passing
`mr=beltrami_mr_options_t()` to `beltrami_prepare` routes every later
`beltrami_solve` / `beltrami_constraint` on that system to this path; without
it nothing changes.

## Equations

After lifting the trace, a fixed-mu solve is

    [ S - mu H_s   G ] [A]   [b(mu)]     S = C^T M2(K) C,  H_s = (P C + (P C)^T)/2,
    [ G^T          0 ] [p] = [c    ],    G = M1 grad_0 (Coulomb multiplier, n = 0).

`S` is applied by `meridian_fourier_curl_curl_apply` and never assembled; `C`,
`P`, `G` are mode-block-diagonal and stay sparse. The lifted trace enters
through `P C` itself, the free-free block through `H_s`, exactly as the
factored path normalises its assembled matrix. MINRES (`npc_exit` off) solves
to `|K x - b| <= tol max(1, |b|)` (default `tol = 1e-10`); the result is
accepted with the factored path's gate (relative residual `<= 1e-9`), and a
MINRES failure is reported as `BELTRAMI_RESONANT` without publishing B. The
`dmu` tangent is a second MINRES solve with the same operator.

**Helicity constraint (nonlinear in mu).** `H(A) = H0` is the stationary
point of the Taylor Lagrangian `L = W(A) + p.G^T A - mu (H(A) - H0)/(4 pi)`,
whose Hessian in `(A, p, mu)` is the symmetric bordered saddle
`[[K, G, -g], [G^T, 0, 0], [-g^T, 0, 0]]`, `g = H_s A + P h`. It is solved by
Newton-MR from `mu0` = mid-bracket and the fixed-`mu0` state, with the line
search on `|grad L|^2` (`MR_MERIT_GRADIENT`, new in `kin6d_newton_mr`; the
original Newton-MR of Roosta et al. for non-minimising stationary points: every
MINRES iterate descends `|g|^2` because MINRES residuals do not increase).
The converged mu must lie in the bracket; the state is re-solved at that mu
and published only if `|H - H0| <= tol max(1, |H0|)`. Closed walls only, as
the factored constraint.

## Preconditioner

Block-diagonal SPD `blockdiag((Sbar + tau M1)^-1, tau L0^-1)`,
`L0 = grad_0^T M1 grad_0`, `Sbar` the phi-lumped curl-curl of the first slice
(one meridional factor per Fourier mode, `nn_factor`):

- **Multiplier block exact.** `C grad = 0` gives `(S + tau M1) grad q =
  tau M1 grad q`, so the Schur complement `G^T (S + tau M1)^-1 G = L0/tau`
  exactly, for every `tau > 0` (Greif and Schoetzau 2007 for the mixed
  Maxwell system). On `grad V0 x V0` the preconditioned saddle has the
  eigenvalues `+-1`.
- **A block.** On the Coulomb complement `|a.H_s a| <= |a||curl a| <=
  C_P |curl a|^2` (discrete Poincare, uniform in `h`), so `S^-1 K = I -
  mu S^-1 H_s` is the identity plus a compact, `h`-uniformly bounded
  perturbation: its spectrum clusters at 1, and its negative part is the
  finitely many Beltrami eigenvalues below `|mu|` (the inertia of `K`). The
  MINRES count is therefore bounded in `h` and grows with that number and with
  the inverse distance of mu to the nearest resonance, as for SPD-preconditioned
  indefinite Helmholtz problems. Mapped, `Sbar ~ S` adds the lumping constants.
- **tau.** `tau` only regularises the `n = 0` gradient kernel, so it is small
  and mu-independent, `tau = 0.01/a^2` (`a` the meridional half-width); the
  factors are reused while Newton moves mu. The H(curl) choice `tau = mu^2`
  (bounded from above by 3/2, but it squeezes the low curl modes towards zero)
  was measured to need 2-5x more iterations (mapped, `mu = 1.1/4.5/8`:
  26/213/1227 vs 13/62/227 with the small `tau`). The mapped 1-form mass
  (phi mean of `K^-1`) in place of `M1` was also tried and was worse.
- No inertia is available on this path; resonance proximity shows up only as
  iteration growth or a MINRES failure.

## Evidence (this VM, `OMP_NUM_THREADS=2`, Release, P2 unless stated)

Agreement with the factored LDLT path on the mapped helical case (rotating
ellipse `eps = 0.1`, nfp 2, R0 = 6, flux 1, closed wall; `agree 2 2 2`, MINRES
tolerance `1e-11`). `mu = 4.5` lies above the first Taylor eigenvalue: the
factored inertia has 9 negative A-block eigenvalues besides the 49 multiplier
ones, so MINRES runs on a genuinely indefinite `K`:

| mu | rel. diff b | energy | helicity | dA/dmu | residual factored / MINRES | MINRES its (solve + tangent) |
|---|---:|---:|---:|---:|---|---:|
| 1.1 | 8.3e-14 | 2.5e-16 | 6.3e-16 | 5.4e-13 | 2.3e-15 / 8.7e-13 | 33 |
| 4.5 | 6.7e-13 | 2.4e-15 | 4.6e-15 | 1.7e-12 | 4.5e-13 / 9.5e-12 | 114 |

Helicity constraint (target `H(1.1)`, bracket `[0.5, 1.9]`): factored
bisection `mu = 1.09999999993`, Newton-MR from `mu0 = 1.2`
`mu = 1.10000000000` in 8 Newton steps and 77 MINRES iterations in total,
`|H - H0|/|H0| = 8.8e-16`.

Analytic/manufactured ladders of `test_nonnested_beltrami` through the MINRES
path (`... mr`, tolerance `1e-12`) reproduce the factored sample and volume L2
errors to all printed digits, hence the same observed orders: manufactured P2
1.940/1.879, P3 2.971/2.923; Lundquist P2 1.766/1.834, P3 2.653/2.832
(sample/L2). The mu = 0 NN1 reduction differs by 2.1e-11 (P2) and 1.5e-10 (P3)
(gate 1e-9), the mu tangent check passes (4.1e-8 Lundquist P2, gate 1e-7), and
the constraints case recovers mu = 1.1 by Newton-MR. Whole-test wall times
factored / MINRES: Lundquist P2 12.3 / 0.8 s, manufactured P2 33.4 / 5.0 s,
Lundquist P3 3.3 / 1.7 s, manufactured P3 5.0 / 5.8 s (unmapped: here the
factored system is mode-block-diagonal and the MINRES advantage comes from
avoiding the saddle LDLT).

MINRES iterations of one fixed-mu solve to `1e-8` (`ladder pc 2 N mu eps 0 nr...`):

| case | nr = 2 (1285) | 4 (4245) | 8 (14645) | 16 (53365) |
|---|---:|---:|---:|---:|
| block, unmapped, mu = 0 | 1 | 1 | 1 | 1 |
| block, unmapped, mu = 1.1 | 10 | 10 | 10 | 12 |
| block, unmapped, mu = 4.5 | 27 | 33 | 33 | 34 |
| block, mapped, mu = 0 | 6 | 6 | 6 | 6 |
| block, mapped, mu = 1.1 | 12 | 13 | 13 | 13 |
| block, mapped, mu = 4.5 | 48 | 62 | 65 | 65 |
| Jacobi, mapped, mu = 1.1 | 1483 | 6560 | > 20000 | |
| Jacobi, mapped, mu = 4.5 | 3444 | 15788 | | |
| none, mapped, mu = 1.1 | 3618 | > 20000 (4.8e-8) | | |
| none, mapped, mu = 4.5 | 7908 | > 20000 (7.2e-3) | | |

(V1 dofs in parentheses, N = 2.) Unmapped, the right-hand side excites only
`n = 0`, which is why its counts are lower than the mapped ones at the same mu.
Against the mode number (mapped, nr = 4): mu = 1.1 needs 13/14/14/14 iterations
for N = 2/4/6/8 (bounded); mu = 4.5 needs 62/93/122/138, tracking the factored
count of negative A-block eigenvalues 9/13/17/21 (each added toroidal mode
brings its Beltrami eigenvalues below mu, up to `n ~ mu R0`).

Wall time and peak RSS of one prepare + solve (`path factored|mr 2 nr N mu`,
fresh process, mapped case):

| nr, N, mu | V1 dofs | factored prepare + solve s / MiB | MINRES prepare + solve s / MiB | MINRES its |
|---|---:|---:|---:|---:|
| 4, 2, 1.1 | 4245 | 0.15 + 0.78 / 32 | 0.08 + 0.30 / 27 | 16 |
| 8, 2, 1.1 | 14645 | 0.43 + 24.4 / 149 | 0.29 + 1.14 / 76 | 16 |
| 16, 2, 1.1 | 53365 | 1.70 + 505 / 801 | 1.08 + 4.05 / 206 | 16 |
| 4, 6, 1.1 | 11037 | 0.51 + 4.14 / 150 | 0.17 + 0.64 / 47 | 16 |
| 8, 6, 1.1 | 38077 | 1.76 + 134 / 496 | 0.62 + 2.32 / 117 | 16 |
| 16, 6, 1.1 | 138749 | > 900 s (stopped) | 2.46 + 8.11 / 375 | 16 |
| 8, 6, 4.5 | 38077 | 1.80 + 141 / 496 | 0.66 + 15.4 / 117 | 137 |
| 8, 10, 1.1 | 61509 | 9.83 + 447 / 1190 | 1.05 + 2.65 / 154 | 16 |

MINRES tolerance `1e-10`; energies agree to the 7 printed digits in every row.
The factored time is the unchanged `beltrami_solve` (assembly, native
symmetric-indefinite LDLT with inertia and pivot checks); it was not profiled
here and grows much faster than the SPD factor of the first slice. Timings were taken on a host shared with other
agents and are indicative only.

Unit checks (`test_newton_mr saddle`): MINRES with the exact block-diagonal
preconditioner on a dense KKT matrix converges in 3 iterations (three
eigenvalues, Murphy-Golub-Wathen); gradient-merit Newton-MR on an
equality-constrained convex quartic (a Lagrangian saddle with 4 negative
Hessian eigenvalues) converges in 9 steps to `|g|/|g0| = 1e-13`, where the
objective merit stalls. The Lanczos step of `minres_solve` now subtracts the
previous vector before forming `alpha` (the Paige-Saunders order); counts of
the existing unit checks drop slightly (SPD 83 -> 79, indefinite 129 -> 119)
and the mapped NN1 tests are unchanged.

## Limits

- No inertia and no near-resonance rejection by pivots: the factored path
  stays the default and the reference for resonance diagnostics.
- Iterations grow with the number of Beltrami eigenvalues below mu; an SPD
  preconditioner cannot remove that (an indefinite per-mode factor would need
  a non-MINRES Krylov method or an `L |D| L^T` factor, which FortSparse does
  not expose).
- The Newton-MR helicity solve finds the root nearest to mid-bracket and checks
  only that it lies in the bracket; it does not enforce the branch/inertia
  condition of the bisection.

Reproduce:

```sh
cmake --build build --target test_newton_mr test_nonnested_beltrami test_nonnested_beltrami_mr
ctest --test-dir build -R 'newton_mr|beltrami.*_mr'
OMP_NUM_THREADS=2 build/test_nonnested_beltrami_mr ladder block 2 2 4.5 0.1 0 2 4 8 16
OMP_NUM_THREADS=2 build/test_nonnested_beltrami_mr path factored 2 8 6 1.1
OMP_NUM_THREADS=2 build/test_nonnested_beltrami_mr path mr 2 8 6 1.1
```

