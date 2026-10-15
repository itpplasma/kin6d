# Nonlinear 3D ideal equilibrium: frozen-in energy map (N1)

Status: first verified stage of step N1 of the 3D TC24 design (nonlinear lane,
DESIGN section 7). Fixed boundary, ideal, nested-compatible, gamma = 0.
Module `kin6d_ideal_map` (`src/mhd/ideal_map.f90`); generated leaves
`src/generated/kin6d_ideal_map_leaves.f90` from `gen/gen_ideal_map.f90`.

## Formulation

The unknown is a map X from abstract coordinates s = (s1, s2, s3) into
Cartesian space: (s1, s2) are the master coordinates of each triangle of the
KIN6D 2D reference mesh (curved P2/P3) and s3 = phi-hat. The frozen reference
flux density b (a closed 2-form, contravariant density in s) and pressure p are
carried by the reference. With J = sigma det DX > 0,

    B = DX b / J,      W[X] = int ( |DX b|^2 / (2 mu0 J) - p J ) ds

per field period (Kruskal--Kulsrud frozen-in principle, the VMEC energy without
flux coordinates). B is divergence-free by construction for every map: within
a cell it is the push-forward of a closed form, and face fluxes equal reference
fluxes because X is continuous. The reference data are a label chi (poloidal
flux), the toroidal density b3 and p(chi), all given in reference coordinates
and pulled back to master coordinates exactly (`ideal_map_freeze`).

At N = 0 with the cylindrical chart this is the inverse (prescribed toroidal
flux, i.e. prescribed q) Grad--Shafranov problem; the generator proves the
pointwise reduction w = R B^2/(2 mu0) - p R at the identity map.

## Charts

The energy, quadrature, leaves and solver see the geometry only through
X = H(q) and its q-derivatives (`ideal_chart_t%jet`: H, dH, d2H, d3H). The map
unknowns are q1, q2 (Lagrange P_k in the plane times a real Fourier series in
phi-hat, n in nfp*Z, |n| <= nfp*N) and omega with q3 = phi-hat + omega.
`cylindrical_chart_t` is H = (q1 cos q3, q1 sin q3, q2). A curve-following
frame (GVEC hmap_axisNB/hmap_frenet type, H = Gamma(q3) + q1 N(q3) + q2 B(q3)
for a closed axis curve Gamma with nfp periodicity) is a further extension of
`ideal_chart_t` supplying the same jet; the reference mesh then lives in the
normal plane. Nothing else changes. The chart types and their generated jets
are the equilibrium core's (`src/equilibrium/equilibrium_chart.f90`,
`gen/gen_equilibrium_core.f90`); a new chart is added there.

## Discretisation

- Energy density, DX-gradient (9) and DX-Hessian (45) are fortsym leaves; the
  chain rule through the chart is assembled in `local_derivatives`.
- Planes: N_phi = 3N + 1 by default, trapezoidal in phi-hat; Gauss order 7 in
  the plane. The integrand is rational, so neither is exact.
- Relabelling gauge (DESIGN "known difficulties"): each node moves in the
  poloidal plane only along a fixed transversal d_a = grad chi / |grad chi|
  (a VMEC-like fixed-angle gauge, omega in the role of lambda), plus a global
  axis expansion T(chi) sum_j sigma_j P_j(y), T = 1 at the O-point and 0 at the
  edge, P_j the monomials of degree <= k - 1 about the O-point node. Without
  the axis expansion fixed transversals cannot translate small surfaces and the
  map becomes cone-like at the axis; the observed order then drops to 2-2.7
  for P2.
- omega gauge: the n = 0 rigid rotation is pinned at one node. omega = f(chi)
  remains a continuum null direction at n = 0 and is discretely near-null. A
  small term (gauge_weight = 1e-4 of the largest physical omega stiffness)
  (w/2) int R^2 (det G grad chi . grad omega_0)^2 / det G penalises the
  across-surface derivative; the physics never sees the surface function it
  selects. The pullback ladders are unchanged by it to four digits. Hard
  orthogonality to hat functions of chi biased the solution (order 1.8,
  measured) and is not used. The N4 remedy (J = boundary nodes/(8 pi) hat
  band constraints C^T du = 0 on omega_0, bordered KKT system, merit on the
  gradient minus its part along C) was tried on the A3 P2 seed ladder: Newton
  steps drop (16/25/52 against 27/40/118) but the unprojected residual stalls
  at 4e-6..9e-6 of the seed gradient and the end points are biased (4.8e-5,
  2.0e-5 against 3.4e-5, 5.4e-6; order 1.3). At the discrete level omega =
  f(chi) is not a gauge: the discrete minimiser uses it to absorb the P_k
  error of B_phi, so freezing it per band costs order.
  `gauge_hessian_only` keeps the form in the Hessian (regularisation of the
  near-null block) and removes it from energy and gradient, so the end point
  is the stationary point of W itself. A3 P2: weights 1e-4, 1e-6, 1e-2 give
  2.77626362e-4, 3.36274479e-5, 5.4509006e-6 (1e-4 and 1e-6 agree to 1e-11,
  3e-9, 9e-9; 1e-2 to 2e-8, 3e-7, 3e-7), and the penalised default differs
  from it by 3e-5, 2e-4, 8e-4 relative. In the P2 crawl (finest mesh,
  steps 110) the Newton step is 98.5 % omega while 1 % of the gradient is
  omega: the crawl is the near-null omega block under Levenberg--Marquardt
  damping 1e-8.
  `omega_elimination` (nested solve, ntor = 0 only): W is exactly quadratic
  in omega at fixed transversal dofs (J does not depend on omega for the
  cylindrical chart at ntor = 0), so omega is eliminated by one undamped
  solve of the omega block at the start and at every trial point, and Newton
  runs on W_red(nu) = min_omega W (no gauge form). With g_omega = 0 the
  full-system solve with right-hand side -g is the Schur-complement step;
  damping and trust region act on nu only. `ideal_map_eliminate_omega`
  exposes W_red; the `reduced` test checks g_omega after elimination
  (1e-14), the W_red gradient against differences (3e-8) and the Schur
  complement (1e-8).
  Gauge directions are excluded from the unknowns, not damped:
  (1) Nested axis basis (`nested_axis_basis`): of the 2 npoly axis fields
  only the 2K tangential harmonics (shift, elongation, triangularity) are
  kept; the near-normal fields (duplicated by nu; tangential fractions
  <= 6e-3) and the field-line slides g b_pol, g = T P (deg P <= K - 2;
  cos^2 >= 0.98 with the dropped directions) are removed.
  (2) Omega gauge (`nested_omega_gauge`): at n = 0 every omega = f(chi) is
  a continuum-exact null (b_pol . grad f(chi) = 0); discretely the smooth f
  are near-nulls of H_oo (curvature ~ h^(2k)), along which omega* jittered
  with the roundoff of g_omega and set a nu-gradient floor of 1e-11 (P2) to
  5e-9 (P3) x g_seed on the finest meshes. Omega is constrained by its L2
  moments against the C1 partition of unity cos^2(pi/2 (L rho - l)),
  rho = (1 - T)^1/2, L = one per element layer, entered as Lagrange rows of
  the omega and Schur solves (the omega pin, a point condition, is
  released; the moments exclude the constant). The moments are smooth
  functionals, so a smooth f restores them for any continuum minimiser.
  The gradient reported for omega is projected onto ker C^T.
  `gauge <k>` (`ideal_map_nested_slide_null_p{2,3}`) measures the W_red
  curvature at the exact A3 map relative to an axis shift (Richardson
  second differences; roundoff floor 1.2e-7):

  | probe | P2 32/64 | P3 24/48 | P3 order |
  |---|---|---|---|
  | slide g = T | 6.1e-4/1.3e-5 | -3.6e-5/-2.0e-6 | 4.2 |
  | slide g = T u1 | 2.5e-3/3.9e-4 | 6.0e-5/1.1e-6 | 5.8 |
  | slide g = T u2 | 2.9e-3/1.8e-4 | 9.2e-5/-8.0e-7 | 6.8 |
  | r^2-translation, tangential | n/a | 7.8e-4/1.8e-5 | 5.4 |
  | translation, tangential | 1.9/2.4 | 4.6/3.3 | (physical) |
  | omega f = T (fixed nu) | 6.2e-5/3.2e-6 (order 4.3) | 5.5e-8/floor | |
  | omega f = sin(2 pi T) | 0.98/0.14 | 0.15/0.029 | 2.4 |

  Outer step (`nested_line_search`): damped Schur Newton direction d
  (g.d >= 0 raises the damping), line search with Armijo 1e-4 and
  g(t).d <= |g.d|/2, cubic interpolation in [0.1 t, 0.9 t]; inside the
  roundoff band of W the change is the trapezoid (g.d + g(t).d) t/2. There
  is no stall exit: convergence is the gradient tolerance only.
Local seed ladders (SuperLU, OMP 4, tol 1e-11 x g_seed per mesh, all
  converged to 6e-14 to 3.4e-12; damping 1e-4 and 1e-2 reach the same end
  points to all printed digits):

  | case | Newton steps (1e-4; 1e-2) | d3D_exact | orders |
  |---|---|---|---|
  | A3 P2 32/64/128 | 9/4/5; 9/5/6 | 2.7686e-4, 3.4298e-5, 5.3369e-6 | 3.01, 2.68 |
  | A3 P3 24/48/96 | 14/6/8; 15/7/9 | 1.0106e-4, 6.3603e-6, 4.1899e-7 | 3.99, 3.92 |
  | shaped P2 32/64/128 | 13/6/6; 13/7/7 | 2.1376e-3, 2.0351e-4, 2.0125e-5 | 3.39, 3.34 |
  | shaped P3 32/64/128 | 26/12/15; 24/13/16 | 2.3029e-4, 1.5962e-5, 9.9914e-7 | 3.85, 4.00 |

  The constrained omega space is a different discrete gauge, so the end
  points differ from the pinned ones at the discretisation level (A3 P3
  finest 4.190e-7 against 4.155e-7) with unchanged orders. Non-nested
  paths are bitwise unchanged (`seed 2`, `solovev 2` end points).
  sCluster (MUMPS, job 3225970) reproduces all eight ladders to the printed
  digits (|g|/g_seed 6e-14 to 8.1e-12; steps within 8 of SuperLU); A3 P3
  with damping 1e-6 takes 14/6/7 steps to the same end points.
  An inner Newton for ntor > 0 converged only linearly in a probe; ntor > 0
  is rejected.
- Fixed boundary: nu = 0 on boundary nodes for all modes.

## Solver

Damped Newton (Levenberg--Marquardt, Hessian diagonal metric) on the reduced
(nu, omega, sigma) coefficients with the exact assembled Hessian (all Fourier
modes coupled), factorised by threaded MUMPS LDL^T (nested dissection) when
KIN6D is built with MUMPS, serial SuperLU LU otherwise. The CSR pattern is
built once per map. Steps must keep J > 0
at every quadrature point. A step is accepted if W decreases clearly; inside
the roundoff band of W it is accepted only if the trapezoidal energy change
(g(x) + g(x + s)).s/2 is not positive (a plain gradient decrease test kept the
damping at 1e-6..1e-9 and 300-600 steps of linear crawl). A trust
region caps the nodal displacement at 5% of the reference extent and omega at
1 rad per step, because the near-null directions carry large Newton components.
Converged when the reduced free gradient falls by 1e-6 (tests; seed ladders
stop at 1e-11 x the gradient of X0 = id). Element energy, gradient and
Hessian blocks are OpenMP-parallel over cells. The report gives wall time for
energy, Hessian and factor+solve.

Shaped Solovev P3 ladder on sCluster (one exclusive node, MUMPS, 433c1a5,
cached CSR pattern), finest mesh (5.1e4 dofs, 35-41 Newton steps), seconds:

| threads | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| energy/gradient | 2.43 | 1.08 | 0.65 | 0.36 | 0.23 | 0.17 |
| Hessian assembly | 9.33 | 5.71 | 5.13 | 4.00 | 3.70 | 3.61 |
| factor + solve | 28.0 | 21.5 | 20.9 | 19.0 | 19.1 | 19.8 |
| whole ladder | 48.4 | 36.3 | 33.8 | 30.4 | 30.2 | 31.0 |

Before caching the pattern the Hessian took 24.1 s (1 thread) and 16.4 s (32);
serial SuperLU took 627 s for the ladder locally (4 threads). Element work
scales (energy 14x at 32 threads). MUMPS factorisation of 5e4 unknowns
is now the floor and stops scaling beyond 4 threads. Its symbolic analysis
(METIS ordering, repeated every Newton step) is 14.4 of 27.7 s at 1 thread and
12.2-12.7 of 19-20 s at 4-32 threads (aaae9b1), so reusing it across steps
would remove about 60% of the factor time; fortsparse has no refactor entry
point yet (numerical JOB=2 on a kept analysis). The per-n block
preconditioner (L4) is the other lever.

## Nested seed (`src/mhd/ideal_seed.f90`)

`nested_seed_t` builds the frozen data from a boundary b(t), q(s), p(s) in
normalised poloidal flux s and the total toroidal flux Phi_edge (the KIN6D 2D
prescribed-q inputs). The boundary is extended harmonically into the unit disk
(Fourier coefficients of b, kept above 1e-13 of the size), optionally with the
seed axis moved to a chosen centre; labels chi = psi_axis (1 - |zeta|^2) and
b3 = 4 pi |psi_axis| q(s)/(R K(s)), K(s) = int det DY/R dtheta (toroidal flux
like 1/R on each surface), are then smooth at the seed axis. A ray map
y = c + rho (b(t) - c) is only C1 at c for a non-circular boundary and left an
error floor (measured). `ideal_map_prolong` initialises a finer map from a
solved coarser one (grid sequencing, VMEC's ns ladder).

## Verification (`test/test_ideal_map_axisymmetric.f90`)

Exact states: Solovev A3 (R0 = 6.2, a = R0/3, FF' = 0) and a shaped Solovev
(kappa 1.61, delta 0.18, FF' /= 0). The labels are the exact state pulled back
through an interior shift Y (Y = id on the LCFS, axis offset about 0.2 m), so
the seed X0 = id is not an equilibrium and the solver must return the axis and
surfaces. Gauge-invariant error: RMS of psi_exact(X(y)) - chi(y) over the
reference, divided by |psi_axis|; the same with KIN6D's 2D fixed-boundary GS
solution psi_2D on the same mesh is the second oracle.

- `ideal_map_derivatives_and_rotation_null`: generated gradient against central
  differences of W, reduced Hessian against differences of the reduced
  gradient, at N = 0 and N = 1 (nfp = 3); rigid rotation null to roundoff.
- `ideal_map_n0_solovev_a3_p{2,3}`, `ideal_map_n0_shaped_solovev_p{2,3}`:
  three-mesh ladders, observed order >= k + 1 - 0.4 against the exact flux and
  against the 2D solution. Measured (2026-10-13, finest pair; seed defect
  4.6e-2 / 3.6e-2):

  | case | k | boundary segments | 3D vs exact | order | 3D vs 2D | order |
  |---|---|---|---|---|---|---|
  | A3 | 2 | 32/64/128 | 2.1e-6 | 3.38 | 2.4e-6 | 3.28 |
  | A3 | 3 | 24/48/96 | 3.3e-7 | 4.04 | 3.2e-7 | 4.05 |
  | shaped | 2 | 32/64/128 | 1.3e-5 | 3.70 | 1.3e-5 | 3.70 |
  | shaped | 3 | 32/64/128 | 1.0e-6 | 4.00 | 9.9e-7 | 4.00 |

  The 2D forward solution on the same mesh is 10-60 times more accurate, so
  the 3D-vs-2D column measures the 3D discretisation error.
- `ideal_map_n3_leakage`: N = 1, nfp = 3, seed with an n = 3 interior
  perturbation and axisymmetric boundary. Map coefficients with n /= 0 are not
  gauge invariant (relabelling along B), so the test measures the phi-variation
  of the surface defect: it must fall below the axisymmetric defect and
  decrease under refinement, and the total energy must equal the N = 0 energy
  to 1e-8. Measured: seed 1.4e-2, final 1.1e-4 and 3.5e-5 (16 and 32
  segments; axisymmetric defect 2.9e-3 and 2.3e-4), energy agreement 4e-9 and
  3e-10. The variation converges at order 1.6 only, because the relabelling
  gauge (near-null) sets its pattern.

- `ideal_map_n0_seed_return_{a3_p2,a3_p3,shaped_p2}`: nested seed from the
  exact boundary, q(s), p(s) and Phi_edge, seed axis 0.19 m (A3) and 0.38 m
  (shaped) from the magnetic axis; grid-sequenced ladders; gates: return to
  1e-3 of the seed error and the k + 1 order. Each mesh stops at 1e-11 x the
  gradient of the seed X0 = id on that mesh (path independent; prolongation
  raises the start to 22-210 x that). Measured (sCluster, 4 threads, MUMPS,
  job 3225970):

  | case | k | segments | errors | orders |
  |---|---|---|---|---|
  | A3 | 2 | 32/64/128 | 2.776e-4, 3.362e-5, 5.446e-6 | 3.05, 2.63 |
  | A3 | 3 | 24/48/96 | 9.7915e-5, 6.2296e-6, 4.1551e-7 | 3.97, 3.91 |
  | shaped | 2 | 32/64/128 | 2.2243e-3, 2.0380e-4, 1.8702e-5 | 3.45, 3.45 |

  A3 P3 uses the transversal gauge weight 1e-2 (default 1e-4 for all other
  ladders) and initial damping 1e-6; Newton steps 89/84/330 (MUMPS), cap 600.
  With weight 1e-2 the end points agree between MUMPS and SuperLU to 1e-11,
  3e-8 and 3e-7 (relative, three meshes) and between damping 1e-6 and 1e-4
  to 1e-10, 1.5e-7 and 2e-6; damping 1e-4 needs 529 steps on the middle mesh
  and SuperLU hits the cap, so damping 1e-6 is the registered setting.
  With weight 1e-4 the middle mesh has two end points (6.1769e-6 and
  6.1633e-6; finest 3.956e-7 vs 3.965e-7), and damping 1e-5/1e-4 hit the cap;
  weight 1e-3 stalled. The residual spread (<= 2e-6 relative) is far below
  the order gates. The P2 gate (2.6) is marginal (2.63, 2.61 measured with
  default flags; 2.57-2.59 with -march=x86-64-v3, a known limitation of that
  flavour; the registered CTest uses default flags).
  The shaped P3 seed ladder is not registered: it converges to a
  weight-dependent end point on the first two meshes (weight 1e-4: 2.3245e-4,
  1.6077e-5; weight 1e-2: 2.3140e-4, 1.5728e-5) and the third mesh reaches
  the 600-step cap. Open. With `gauge_hessian_only` (sCluster, damping
  1e-6) the shaped P3 end points agree between weights 1e-4 and 1e-2 on the
  first mesh (2.3248620e-4, 2.3248602e-4) but need 287-307 steps there and
  the next mesh is not converged to 3e-5 between weights (1.60749e-5,
  1.60744e-5); the A3 P3 ladder with weight 1e-4 needs 82/320 steps and fails
  on the finest mesh (gradient 4e-2 after 208 steps), weight 1e-2 fails on
  the first. Damping omega less than the transversal dofs (ratios 1e-2..1e-6, P2) did
  not shorten the crawl reliably (67/52, 69, 195 steps on the finest/middle mesh). Still open: a nested solve of the omega block (W is
  exactly quadratic in omega at fixed transversal dofs, since J does not
  depend on omega) so that the flat omega valley is not driven by the
  nonlinear nu--omega coupling.

  Reproducibility: cell energies are summed serially and MUMPS factors on one
  thread by default (`factor_threads`), so a solve is identical for any
  OMP_NUM_THREADS (checked at 4 and 32 threads on sCluster); with threaded
  MUMPS the seed ladders failed intermittently. `sparse_refactor` (fortsparse
  7300414) keeps the MUMPS analysis across Newton steps: on the finest shaped
  P3 mesh the analysis drops from 12-14 s to 0.13 s and factor+solve to
  13.2 s (one thread). The SuperLU backend (fortsparse 358a2ae) keeps its
  minimum-degree column ordering across same-pattern refactorizations; the
  factor equals a fresh one, so the Newton traces are bitwise unchanged, and
  the ordering, which was most of each SuperLU factorization, is skipped:
  on the finest A3 P3 seed mesh (23269 dofs, 490 factorizations)
  factor+solve drops from 322 s to 116 s.

  Position gauge (`position_gauge`, as in N4: full (q1, q2), omega = 0): the
  pullback A3 P2 ladder is ten times more accurate (2.8e-5, 3.8e-6, 5.5e-7;
  orders 2.8, 2.8) but converges slowly (32-97 Newton steps); from a nested
  seed it does not return (1.2e-2, 8.7e-3, 5.7e-3), because redistributing
  the seed's toroidal flux within surfaces needs poloidal slides of order
  q R d(omega). A transversal phase followed by a position-gauge polish
  returns (A3 P2 2.9e-4, 3.1e-5; shaped P2 orders 3.43, 3.45) but the polish
  stalls at a relative gradient of 1e-8 on the A3 meshes; the polish is not
  needed once the seed-gradient stopping rule is used, and was removed.

  Cause of the earlier deficit (order 2.2-2.4; `seeddiag 2` separates it on
  A3 with the smooth exact map X* = P(xi(y)) of the seed labels,
  P(xi) = (sqrt(R0^2 + 2 R0 a xi1), R0 a xi2/R)): the P2 interpolant of X*
  converges at 3.1, 3.0 and quadrature degree 7/12/20 changes nothing, so
  map space, frozen data and quadrature are consistent. With b3 uniform in
  the seed disk (1/det DY), the equilibrium needs omega of 3-4.4 rad; its
  P_k error in B_phi, amplified by B_phi/B_pol, set the surface error at
  order k. Distributing the seed toroidal flux like 1/R on each surface
  (current seed) cuts omega to about 1 rad and restores k + 1. The linear
  axis-expansion terms also carried a tangential drift of 0.37-0.39 m
  (translation only: 2e-4 m), but removing them loses the P3 order (2.4,
  2.5, cone at the axis), so they stay.
  Tried without success: re-gauging after prolongation; moving normal
  gauges (centroid or nodal normals; both also degrade the pullback ladder);
  increments perpendicular to B per node (cannot change omega by radians
  without huge poloidal slides; no convergence); an anchor (harmonic-map)
  regulariser (biases, since the physical omega is O(rad)).

- `ideal_map_n0_normalised_residual_a3_p2`: the weak residual of
  equilibrium-core D3 (`ideal_map_residual`): R[xi] = -dW[xi] for every
  basis displacement (dR, dZ, R d omega), dual L2(V) norm through the
  consistent mass, times (A_ref/pi)^1/2, divided by sigma |V|^1/2 with
  sigma = <|B - F_vac/R e_phi|^2>/mu0 + <p> (1.36 here against B_phi^2 = 28).
  Seed 4.1; exact interpolant 1.4e-1, 6.0e-2, 2.9e-2 (order 1.3, 1.0 = k - 1
  for an L2 force); discrete solution 2.7e-2, 6.2e-3, 1.5e-3 (order 2.1,
  2.1 = k). The solution is stationary in the gauge space only, so the
  tangential components stay at discretisation level.

## Not yet done

Newton--Krylov with the per-n L4 preconditioner and the Hessian = L4_n
identity (equilibrium-core D5); the prescribed-q GS linearisation identity at
N = 0; reuse of the MUMPS symbolic analysis across Newton steps (needs a
refactor entry point in fortsparse); the current-constrained profile f(chi)
for forward cases; the seed order deficit above; 3D cases and free boundary;
the map -> core state adapter and the map-lane D3 gate in the core.
