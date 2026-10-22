# Non-nested equilibrium: NN2 (resistive relaxed lane)

Status: design decision, first increment (transport for a given B),
second increment (axisymmetric coupled Newton) and third increment
(Fourier-coupled Newton, N >= 1), 2026-10-13. NN2 is the
resistive relaxed lane of the three-dimensional programme
([benchmarks.md](benchmarks.md)). It fixes open question 5 of
[equilibrium-core.md](equilibrium-core.md) for this lane only. The ideal
topology-preserving lane (MRX class) and MRxMHD (SPEC/SPECTRE class) come
later and are not compared against NN2 results.

## Continuous problem

Unknowns: B = B_h + curl A on the NN1 complex ([equilibrium-nonnested.md]
(equilibrium-nonnested.md)), the pressure p, an electric potential Phi and
the loop voltage V. The harmonic part B_h = F_vac/R e_phi and the normal
trace of B on the boundary are data.

1. Pressure transport. The transported scalar is p, with uniform density,
   so p = n0 T and the heat flux is -K grad p:

       -div(K grad p) = S,   K = chi_par b b + chi_perp (I - b b),
       b = B/|B|,  p = p_b on the boundary.

   p rather than T: D3 needs p directly, and a nonuniform n would need a
   density transport model that the non-nested topology cannot borrow from
   flux labels. A nonlinear p(T) adds a constitutive law but no structure,
   so it is deferred. The source S(x) is declared input with its own
   provenance. Dirichlet p_b (default 0) on the computational boundary
   matches the Sovinec test and the HINT2 edge treatment, where p is set to
   zero on lines that leave the domain before the averaging length.
2. Force balance. J x B = grad p holds when Ampere's law is satisfied with

       J = J_perp + J_par,  J_perp = B x grad p / |B|^2,
       J_par = g_par b b.E - G grad Phi,  E = V grad(phi)/(2 pi),
       G = g_par b b + g_perp (I - b b),

   in the weak form <curl dA, B>/mu0 = <dA, J> for all dA in V1_0. The
   common D3 residual is the independent check and is evaluated, never
   replaced. Its smooth form vanishes exactly when J_perp is as above.
   g_perp acts on grad Phi only. It regularises the potential equation on
   closed surfaces and carries the O(g_perp) perpendicular current that
   the ideal limit removes; the inductive field drives parallel current
   only. g_par may be a declared profile g_par(x) = g0 w(x).
3. Current continuity. The weak Ampere system is solvable only when
   <grad chi, J> = 0 for all chi in V0. This condition is the Phi equation:

       -div(g_par b b.grad Phi + g_perp grad Phi)
           = -div(J_perp) - div(g_par b b.E).

   It is the same anisotropic operator as (1). In the limit
   g_perp/g_par -> 0 it becomes the magnetic differential equation
   B.grad(J_par/B) = -div J_perp of PIES. Its zero-loop-voltage solvability
   sets the parallel current averaged over each surviving surface to zero.
4. Fixed quantities: the boundary normal trace of B (NN1 trace; coils
   through D6 later), the harmonic F_vac, the source S, the boundary p_b,
   chi_par, chi_perp, g_par and g_perp, and the net toroidal current I_phi.
   V is the Lagrange multiplier of I_phi; the default I_phi = 0 is the
   current-free stellarator class. A q or iota profile is rejected as a
   constraint: it is undefined on islands and stochastic regions. Helicity
   is not conserved in this resistive lane.

The ideal limit is the sequence chi_perp/chi_par -> 0, g_perp/g_par -> 0
at fixed S and I_phi. It is measured as a sequence and extrapolated with
reported rates; it is never assumed. A smooth p becomes constant on
ergodic components only in that limit, so islands below the critical width
for finite chi_perp/chi_par keep a pressure gradient. [Fitzpatrick 1995]
(https://doi.org/10.1063/1.871434), Eq. (27), gives the scale width
W_c = sqrt(8) (kappa_perp/kappa_par)^(1/4) (r_s R_0/(s_s n))^(1/2)
(checked against the paper body; cylindrical large-aspect-ratio slab
analysis with constant conductivities, island width W of its Eq. (9)).
The temperature is flattened inside the separatrix for W >> W_c and not
for W << W_c. Our toroidal annulus with walls is outside those
assumptions, so the formula is a comparison for measurements, not a gate.

## Discretisation

- p and Phi lie in the Fourier-extended V0 space: P_k Lagrange in (R, Z)
  times the real modes 1, cos(n phi), sin(n phi), n <= N. A phi-dependent
  B couples all modes. FortFEM owns that assembly
  (`fortfem_meridian_fourier_scalar`), with nphi >= 2N + 1 uniform phi
  points. The isotropic forms are exact; the anisotropic tensor carries a
  declared phi-aliasing error that is part of the Fourier truncation
  ledger. The weak load uses the same rule, so manufactured consistency is
  exact.
- B stays the closed NN1 2-form; div B = 0 is an algebraic identity. b
  enters only through K and J_par at quadrature points. |B| = 0 at a
  quadrature point is a reported failure, never regularised silently.
- K and its JVP/VJP in (B_R, B_phi, B_Z) come from the generated leaf
  `anisotropic_conductivity` (gen/gen_nonnested.f90). J_perp, J_par and the
  D3 stress follow the same rule, with the stress leaf already generated.
- The Newton Jacobian on (A, p, Phi, V) has these blocks:
  - A-A: curl-curl/mu0 minus the B-derivatives of J_perp and J_par
    (generated JVPs);
  - A-p: B x grad dp/|B|^2;
  - A-Phi: -g_par b b.grad dPhi;
  - A-V: g_par b b.grad(phi)/(2 pi);
  - p-A: div(dK grad p) with dK from the conductivity JVP;
  - p-p: the transport matrix (this increment);
  - Phi-A, Phi-p, Phi-Phi: the same pattern with g;
  - V row: the current constraint.

  Derivatives are equation-level. Adjoints transpose the same factorised
  system; no solver iteration is differentiated.
- Gauge and nullspaces follow NN1: R A_phi = 0 for n > 0 and the scalar
  R A_phi at n = 0. Phi is fixed by Dirichlet data on the boundary (an
  electrode-free wall: Phi = 0).

## Perpendicular pollution

Large chi_par/chi_perp on a mesh not aligned with b pollutes the
perpendicular transport. Field-aligned meshes are rejected because the
lane exists for islands and stochastic regions. The choice is high-order
conforming P_k (k >= 3) on unaligned triangles with spectral phi.
[Sovinec et al. 2004](https://doi.org/10.1016/j.jcp.2003.10.004)
(pp. 374-375) found p >= 3 adequate up to 1e9 and bilinear elements
inadequate at 1e6. [Günter et al. 2005](https://doi.org/10.1016/j.jcp.2005.03.021)
reach anisotropy-independent error with a symmetric finite-difference
scheme on the same test. The fallback, if the ideal-limit sequence hits
the conditioning floor reported below, is an asymptotic-preserving
reformulation ([Degond, Deluzet, Negulescu 2010]
(https://doi.org/10.1137/090754200)). The Günter, Lackner, Tichmann 2007
finite-element paper (doi:10.1016/j.jcp.2007.07.016) is not yet read and
is not relied on.

## Anisotropy estimator diagnostic

The tensor patch minimisation already uses the full A^(-1) metric. FortFEM's
opt-in `correction_sweeps` adds divergence-free RT corrections on the same
vertex patches. A simultaneous sweep is damped by 1/3, since three vertex
patches overlap each triangle. Orthogonal projection and Cauchy--Schwarz give
`||e_next||^2 <= ||e||^2 - sum_a ||delta_a||^2/3` in the sampled A^(-1) metric.
Divergence is unchanged in exact arithmetic; its measured defect still enters
the bound. The correction solves use block scaling to retain the constraint
accuracy. At fixed degree and sweep count the work and storage are O(N), and
the cell and vertex passes use OpenMP. No anisotropy-uniform contraction is
proved. The sampled integrals remain approximate, and `nt_solve` retains its
existing default reconstruction.

| P2 Sovinec mesh | Ratio | Original index | Eight-sweep index | Status |
|---|---:|---:|---:|---|
| nr = 4 | 1e9 | 7.35 | 1.26 | Practical sampled bound |
| nr = 8 | 1e9 | 27.5 | 4.50 | Practical sampled bound |
| nr = 32 | 1e6 | 163 | 20.1 | Uniform-efficiency target missed |
| nr = 32 | 1e9 | 389 | 47.8 | Uniform-efficiency target missed |

`nonnested_transport_anisotropy_correction` checks containment, improvement,
index < 5 on nr = 4/8 and convergence for ratios 1 and 1e3. The larger study
is `test_nonnested_anisotropy study 32 8`; it reports containment and cost,
without asserting the unmet uniform-efficiency target. The upstream tensor
test checks exact polynomial recovery and a rotated constant-tensor ladder.

On a shared host, three Release runs at nr = 32 (12,193 DOFs), ratio 1e9,
give median reconstruction times 2.80 s at one thread and 0.781 s at four
threads (3.59x); the four-thread range is 0.768--0.791 s. The median
per-run cost is 9.44 times assembly plus the direct solve (range
8.57--10.81), excluding mesh construction and estimator input sampling.
At nr = 16 (3,105 DOFs) the four-thread median is 0.208 s: the two-point
DOF exponent is 0.97. This supports linear work at fixed sweep count but
misses the required small fraction of the solve.

The companion `nonnested_transport_anisotropy_pollution` test compares the
production tensor assembly/sparse solve with an independent directional-Gram
assembly/dense Cholesky solve, equivalent strong and weak loads, and quadrature
degree 24 instead of 8. At ratios 1e6/1e9 on nr = 4/8, solver/assembly
agreement is within 2.2e-11 and the quadrature change is below 1e-6, while the
center-pressure error exceeds 0.96. The energy error agrees with the independent
identity `||grad psi||^2 - 2(f,p_h) + a(p_h,p_h)` to 2e-8.

There is also an exact fixed-mesh explanation. If a cell polynomial q obeys
`cos(pi x) sin(pi y) q_x - sin(pi x) cos(pi y) q_y = 0` on an open triangle,
real analyticity extends the identity to the plane. At every integer x it
forces q_x to vanish for all y. A polynomial with infinitely many such roots
has q_x = 0 identically; the identity then gives q_y = 0. Thus the parallel
Gram operator is positive definite on the conforming, zero-trace polynomial
space. At fixed mesh/order `(L + (chi_par-1) G) p_h = f` therefore gives
p_h -> 0 as chi_par -> infinity, whereas the exact Sovinec solution is psi.
This proves fixed-mesh pollution for this discretization. It does not diagnose
the separately observed P4 conditioning floor or prove correctness of every
assembly path. A field-aware coarse correction or different localisation is
still required for the campaign's uniform-efficiency and cost targets.

## Verification gates

1. Manufactured 3D transport (implemented). An exactly divergence-free
   helical torus field B = grad phi x grad psi_h + F grad phi with
   psi_h = ((R-R0)^2+Z^2)/2 + eps((R-R0)^2-Z^2)cos phi, and a smooth
   non-flux-function p_e with modes 0..2. The load is weak (flux K grad p_e,
   from analytic gradients). Gates: observed L2 order >= min(k+1, 3) - 0.35
   and error <= 5 x nodal interpolant error at chi ratio 1 and 1e4.
2. Anisotropic convergence against the chi ratio (implemented). This is
   the Sovinec test on the periodic-cylinder chart:
   psi = cos(pi x)cos(pi y) on [-1/2, 1/2]^2, B = z x grad psi,
   S = 2 pi^2 psi and exact p = psi/chi_perp. It runs ratios 1, 1e3, 1e6
   and 1e9 on an affine concentric-square mesh not aligned with psi.
   Gates: backward error <= 1e-13; ratio-1 error < 1e-7; best error below
   1e-1 at 1e6 for P2 and below 1e-2 at 1e9 for P >= 3.
3. Solov'ev nested limit (implemented, n = 0). The coupled solve with
   S = -div(K grad p_GS) = -chi_perp p' lap(psi) for p_GS = p' psi must
   return the analytic Solov'ev state. B.grad p = 0 holds exactly, so the
   transport sees only chi_perp. Ohm's law constrains the case: with
   constant F the exact current is J = R p' e_phi, whose parallel part is
   J_par = p' F B/|B|^2. With Phi = 0 this is Ohmic only if
   g_par V/(2 pi R^2) = p'. A uniform g_par therefore cannot hold Solov'ev
   (the steady Ohmic current profile is not free), and the case declares
   g_par = g0 (R/R0)^2, V_exact = 2 pi R0^2 p'/g0, Phi_exact = 0 and
   I_phi = p' int R dR dZ on the discrete domain. Gates: field, pressure
   and D3 ladders, V and Phi, the gauge multiplier and offset return.
   Landreman family 1 remains the 3D nested limit; see below for what is
   missing.
4. Finite-pressure island from an offset start (implemented as a
   measurement, see the third increment). The NN1
   helical-filament vacuum field plus finite S, started from the
   axisymmetric solution, is relaxed by Newton. Measure O/X points, width,
   helical phase and the pressure jump across the island against the
   ratio. Flattening is the measured quantity: (p_O - p_X)/(p_axis - p_b)
   must decrease monotonically with chi_perp/chi_par, and the island must
   exceed the measured critical width.
5. Ideal-limit sequence (next). Use ratios 1e3..1e9 (and g likewise) at
   fixed mesh ladders, extrapolate the island observables with reported
   rates, and stop before the conditioning floor.

## Comparability

Compare only within the resistive relaxed lane; other lanes are context.

| Code | Pressure model (verified) | Admissible when |
|---|---|---|
| HINT2 ([Suzuki et al. 2006](https://doi.org/10.1088/0029-5515/46/11/L01)) | p averaged along lines over L_in; p = 0 if L_C < L_in; B relaxed at fixed p | same boundary and vessel, I_phi = 0, matched p on surviving surfaces and beta; compare in the ideal-limit extrapolation only, since L_in is not chi_perp |
| HINT ([Harafuji, Hayashi, Sato 1989](https://doi.org/10.1016/0021-9991(89)90069-7)) | algorithm not yet read from the source | none until read |
| PIES ([Reiman, Greenside 1986](https://doi.org/10.1016/0010-4655(86)90059-7)) | pressure-driven current, Ampere update, zero net current per surface | ideal limit with V = 0, I_phi = 0; same boundary and p on surviving surfaces |
| Zhang, Zhu, Hegna ([CPC 2026](https://doi.org/10.1016/j.cpc.2025.109961)) | parallel-only hyperbolic damping of p (chi_perp = 0) with artificial resistivity and flow; not a chi_par/chi_perp closure | their published case verbatim; parameter map from the paper body, not yet read |
| JOREK saturated states | full resistive MHD with chi_par, chi_perp, eta and sources | negligible flows; identical chi, eta, S, boundary and current; compare p, island width and phase |

[MRX](https://arxiv.org/abs/2510.26986) relaxes by admissible
variations and is not an NN2 reference.

## First increment (implemented)

- FortFEM `b1aca63b` (branch `meridian-fourier-diffusion`, unpushed):
  `meridian_qgrad`, `fortfem_meridian_fourier_scalar` and its test. It
  checks the isotropic blocks against the per-n Dirichlet forms to 1.6e-15,
  the loads to 3.5e-16, and exact Galerkin reproduction of P2 x trig with a
  phi-dependent tensor to 8.7e-15.
- KIN6D `kin6d_nonnested_transport`: `nt_conductivity_t` samples any
  `equilibrium_field_t`, including the NN1 free-B state, and `nt_solve`
  runs the Dirichlet solve through `nn_linear_solve`. The reported residual
  is the normwise backward error |r|/(||A||p|| + |rhs|), with up to two
  refinement steps. The generated leaf adds `anisotropic_conductivity` with
  JVP/VJP. CTest entries are `nonnested_transport_manufactured_3d` and
  `nonnested_transport_sovinec_chi_ratio_p{2,3,4}`.

Local evidence (gfortran -O3, 4 threads) gives the best |T(0,0) - 1| over
the ladder nr = 4..64:

| ratio | P2 | P3 | P4 |
|---:|---:|---:|---:|
| 1 | 9.3e-9 | 1.3e-9 | 2.1e-12 |
| 1e3 | 1.4e-5 | 3.7e-8 | 3.2e-10 |
| 1e6 | 1.8e-2 | 1.1e-5 | 5.5e-8 |
| 1e9 | 9.4e-1 | 1.6e-3 | 4.5e-5 |

P2 at 1e9 is fully polluted, as Sovinec reports for low order. The P4
ladder shows a conditioning floor: at 1e9 the error rises from 4.5e-5
(nr = 32) to 1.5e-4 (nr = 64), although the backward error stays at 8e-17.
The forward error is then bounded by cond x eps. The ideal-limit sequence
must stay above that floor or switch to the AP fallback.

Manufactured 3D (N = 2, nphi = 8): the observed orders are 2.73/3.64 (P2)
and 3.02/3.30 (P3) at ratios 1/1e4, and errors are within 2.4 x the
interpolant. The P3 rate is limited by the flux-mesh axis fan, whose
polar geometry map is not polynomial. With the r < 0.3 fan excluded, the
same solution converges at 3.98. That is a pre-existing `flux_mesh_t`
defect, recorded here and not worked around.

Claims: candidate. No inf-sup or continuous enclosure is established;
errors are reference comparisons. The phi-aliasing error of K is declared
but not bounded. Performance: assembly is O(cells x nq x nphi x (nl(2N+1))^2),
OpenMP over triangles, with no global kinetic or pairwise matrices and a
sparse direct solve. The coefficient callback of an NN1 field locates each
point in the mesh. Evaluation in the same triangle is the obvious next
optimisation, to be made only after measurement.

## Transport error estimate

`nt_solve(..., estimate_error=.true.)` retains the N = 0 tensor estimator
for a phi-independent K. For phi-dependent K and N > 0 it now reports a
practical energy majorant of the full retained, phi-sampled Fourier-Galerkin
system (`scope='phi-sampled-fourier-galerkin'`). FortFEM reconstructs
`R (K grad p - f)_pol` mode by mode in H(div), retaining the toroidal flux
and its cosine/sine derivative in the equilibrium constraint. The residual
uses the guaranteed domain-wide lower ellipticity constant chi_perp and
Parseval weights. Local divergence defects and data oscillation are
included. Neither phi/meridional quadrature, the omitted Fourier complement
nor the error of the imposed Dirichlet trace is enclosed; the claim remains
`practical-fixed-dirichlet-trace`. Derivatives of this estimator are unavailable.

The fast test `nonnested_transport_fourier_equilibrated` uses a strong
manufactured source with nonzero cosine/sine modes and phi-dependent scalar
K. P4 exact recovery leaves a 2.3e-11 bound; the P2 ladder gives energy order
1.97–1.99 and efficiency 1.21–1.23. At 83,205 DOFs and four threads,
three generic-estimator repetitions average 0.818 s against 5.64 s for assembly
and SuperLU solve (14.5%); one-thread estimator time is 2.53 s (3.09x speedup).
An earlier shared-host repeat gave 19% cost and 3.39x speedup.
These shared-host timings exclude the automatic wrapper's preliminary
phi-dependence scan and are indicative. The existing helical anisotropic
case also checks the full coupled energy: bound 0.206 versus error 0.133
at chi ratio 100. The lower-ellipticity defect norm may become loose at
large anisotropy. FFT sampling costs O(cells * nphi log nphi); reconstruction
and storage are O(cells * modes) at fixed spatial degree.

## Second increment (implemented): axisymmetric coupled Newton

`kin6d_nonnested_coupled` solves the full n = 0 system on the NN1 complex.
The unknowns are x = [A | lambda | p | Phi | V]. A = (A_pol, a = R A_phi)
is the whole V1 space, so B_phi = F_vac/R - rot A_pol is free and poloidal
current obeys Ampere. At n = 0 the Nedelec curl has the kernel grad V0_0,
so a Coulomb multiplier lambda in V0_0 gives the rows <A, grad chi> = 0 and
the column <dA, grad lambda>. The current equation is Ampere tested with
grad chi. Lambda therefore vanishes identically once that equation holds,
which is the discrete solvability of weak Ampere; the run reports it as a
check. Dirichlet data: the A trace (normal B and the toroidal flux), p_b,
Phi = 0. V is the multiplier of the I_phi row.

The leaves `coupled_current` (J above) and `anisotropic_flux` (K grad p)
and their JVP/VJP come from gen/gen_nonnested.f90. The element kernel
contracts them with FortFEM's quadrature-point bases (V2 for B, V1 for the
Ampere test, V0 values and gradients). The B columns are chained through the
FortFEM curl matrix, and the constant blocks are curl^T M2 curl/mu0, M1 grad
and its transpose. The Jacobian is exact: the unit test compares it with a
central difference of the residual (6e-11 at P2, 1e-12 at P3 for ratio
1e4). The leaf test checks JVP/VJP pairing, J x B = grad_perp p and
q.b = chi_par b.grad p. Newton uses a general sparse LU (MUMPS unsymmetric
when available, otherwise SuperLU) with backtracking on the blockwise
backward error |r|/(|J| m), where m is |x| plus the block maximum. A later
step first runs FGMRES on the current exact Jacobian, right-preconditioned
by the LU of an earlier step of the same solve, and accepts it at relative
residual 1e-10 within 20 iterations; otherwise it refactors. The Newton
histories are unchanged; the factorisations drop from one per step to
typically one or two per solve. Tangents and adjoints of a converged state
are solves with that Jacobian.

Solov'ev evidence (local, gfortran -O3, 4 threads, ratio chi_par/chi_perp =
g_par/g_perp = 1e6, seed with the axis shifted by 0.1 m and V = 0). Errors
are relative RMS over the 128 frozen points; Phi is max |Phi| in units of
V/(2 pi).

| P | nr | B | p | D3 | Phi | Newton steps |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 8 | 2.35e-3 | 1.53e-4 | 1.23e-4 | 1.5e-6 | 4 |
| 2 | 16 | 8.25e-4 | 1.06e-5 | 2.05e-5 | 3.2e-7 | 4 |
| 2 | 32 | 1.33e-4 | 2.10e-6 | 3.06e-6 | 5.6e-8 | 3 |
| 3 | 8 | 4.83e-5 | 1.34e-6 | 5.42e-6 | 2.5e-8 | 4 |
| 3 | 16 | 6.86e-6 | 8.59e-8 | 6.76e-7 | 5.6e-9 | 3 |
| 3 | 32 | 1.05e-6 | 4.51e-9 | 8.43e-8 | 3.7e-9 | 3 |

Observed orders (nr 8 to 32), B/p/Phi/D3: P2 2.07/3.10/2.36/2.67;
P3 2.76/4.11/1.38/3.00. The P3 B order is the axis-fan limit already
seen in NN1 and increment 1; Phi at P3 levels off near 4e-9. The B, flux and D3 values equal the NN1 Solov'ev ladder to the
printed digits, as they should: the coupled state is the same discrete
equilibrium. |V - V_exact|/|V_exact| is below 2e-9, and the multiplier
lambda stays below 1e-17. Newton converges quadratically from the offset
seed (P2: 5.0, 9e-5, 7e-7, 3e-11, 1e-15). From the exact nodal
interpolant with V_exact the solve returns to the offset-seed state within
1.4e-13 (P2) and 2.2e-10 (P3, in p) when both run to the rounding floor
of the backward error. That difference is cond x eps, and the transport
condition number grows like chi_par/chi_perp, so the gate is
max(1e-11, 10 ratio eps).

The errors barely depend on the ratio. At P2, ratios 1e4, 1e6 and 1e8
give the same B and D3 to the printed digits and p within 1.5%. The reason is structural. At n = 0 the discrete B_pol is
rot90(grad a_h)/R, so b_h.grad(p' a_h) = 0 at every point, and the
discrete pressure has no parallel gradient to pollute. This nested limit
therefore does not test perpendicular pollution; the Sovinec gate does.
The ratio changes mainly the conditioning: the Newton tail lengthens from
3 steps at 1e4 to 4 at 1e8.

CTest entries: `nonnested_coupled_solovev_p{2,3}` (ladder, Jacobian check,
gauge, V, Phi and D3 gates) and `nonnested_coupled_offset_return_p{2,3}`,
all at ratio 1e6. The study recipe adds the
coupled ladders at ratios 1e4, 1e6 and 1e8.

Limits. Only n = 0 is coupled: `nc_*` rejects a field with modes. Memory
peaks at 1.2 GB (P2) and 4.3 GB (P3) at nr = 32 with SuperLU, and the
local P2/P3 ladders take 45 s and 252 s. All element triplets of the
Jacobian are held before compression; chunked compression and factor
reuse in the Newton tail are the obvious savings, after measurement on
the MUMPS build. The domain is the axisymmetric chart;
the periodic-cylinder chart is rejected. Claims stay candidate: there is no
inf-sup analysis of the coupled saddle system and no enclosure.

## Third increment (implemented): Fourier-coupled Newton

`kin6d_nonnested_coupled` now solves modes 0..N (N = the field's mode
count). Every V0/V1/V2 unknown uses FortFEM's real block layout (block 0,
then cos(n phi), sin(n phi)); NN1's complex modes map by 2 B_n = c - i s.
FortFEM `a27e21f7` (`fortfem_meridian_fourier_curl`) supplies the block
curl, curl^T M2 curl and M1 grad with the identity tensor, so the constant
blocks are exact. J and K grad p are evaluated at nq x nphi points
(nphi = 4N + 1 by default, at least 2N + 1) with phi-dependent B, and the
toroidal derivatives enter through the phi components of grad p and
grad Phi. These are the existing generated leaves and JVPs: no new
symbolic terms were needed. The A-A, A-p, A-Phi and A-V Jacobian blocks
are contractions of those JVPs with the block trial and test operators,
so a phi-dependent b couples all mode pairs. The Coulomb multiplier is
used in every block, and lambda again vanishes by solvability (gate:
lambda mu0/|A| < 1e-10; observed 1e-14). The backward error now includes
the source load |<v, S>|, and the line search is scaled by |J|(|x| + |dx|)
so that a seed with p = 0 or Phi = 0 is admissible.

A shell domain (the NN1 coil annulus) has a harmonic A at n = 0: curl
grad f = 0 with f = 1 on the inner wall. `problem%floating` marks the
inner-wall V0 nodes as a floating conductor. Its mode-0 lambda and Phi
nodes are tied to one unknown each, which removes the null space and
imposes zero net current into the wall. `nn_solve`/`nn_vacuum_mode`
return the NN1 potentials (`solution`), and `field%p` holds the p blocks.
`nn_point` and the D3 residual evaluate p at phi.

Gates and evidence (local, gfortran -O3, 4 threads):

- Reduction (`nonnested_coupled_fourier_reduction_p{2,3}`). The Solov'ev
  offset seed with N = 1 stays axisymmetric: the n = 1 content is 9e-17
  relative, the mode-0 blocks equal the N = 0 solve to 3e-15 (A, p),
  2e-13 (Phi in V/(2 pi)) and 0 (V), and D3 is unchanged. The Jacobian
  matches a central difference at a non-axisymmetric state to 2e-8.
- Island flattening (`nonnested_island_flattening_p3`; P3, nr = 8, N = 1,
  90433 unknowns, 4.6 min). The setup is the NN1 coil vacuum on the annulus u in
  [0.3, 1.2], a source shell u in [0.75, 1.05] and p = 0 on both walls
  (beta ~ 1e-3). Each ratio chi_par/chi_perp = g_par/g_perp is solved from
  the vacuum seed (p = Phi = 0) in 3-4 Newton steps. The measure is G,
  |dp/dR| at the O point (on Z = 0, phi = 0, +-0.03 m) over the
  mean slope on R_O +- 0.12 m. G is 1 for an unaffected profile and 0
  when p is flat at O:

  | ratio | 1 | 1e2 | 1e4 | 1e6 |
  |---|---:|---:|---:|---:|
  | G | 0.982 | 0.973 | 0.948 | 0.143 |
  | D3 | 5.1e-4 | 4.7e-4 | 1.4e-4 | 6.7e-5 |

  The gates are: G non-increasing within 0.01, and G(1e6) < G(1)/2. The
  island width at 1e6 is 0.189 m in u at Z = 0 (edges 0.382..0.571).
  Fitzpatrick's W_c with r_s ~ 0.48 m, R_0 = 6 m, n = 1 and s ~ 2 (a line
  current's vacuum shear, not measured) is 0.34 m at 1e4 and 0.11 m at
  1e6. W/W_c therefore rises from 0.55 to 1.7 between the two ratios,
  where the flattening sets in. This is consistent with the formula but
  not a test of it: the geometry is toroidal with walls, s is estimated,
  and two ratios do not resolve a power law. The O-X pressure difference
  is reported but not gated. It changes sign (-0.38, -0.33, +0.012,
  +0.105), and p_max drops from 9.0e4 to 1.8e4 at 1e6 because source
  lines connect to the walls. At P2 the same scan loses monotonicity at
  1e6 (G stays above 0.5, chord non-monotone), as expected from Sovinec
  pollution on this mesh.
- D3 against resolution (`nonnested_island_d3_ladder_p2`; P2, N = 1,
  ratio 1e4, nr = 4/8/16, 5.3 min): 5.8e-3, 9.9e-4, 1.8e-4, order 2.5.
  G at nr = 8/16 is 0.890/0.945.
- The O and X points are found on the stellarator-symmetry line Z = 0 by
  half-turn bisection. The FEEC map is only piecewise smooth, so 2D Newton
  refinement fails below nr = 16. The one-turn closure is 2e-4 m at P3
  nr = 8 and 2e-3 m at P2 nr = 8.
- div curl A rounding is about 1e-12 at |A| ~ 20 Wb, so the island test
  gates it at 1e-10. Solov'ev keeps 1e-12.

The sCluster island study (`study_nonnested_nn2_scluster.sbatch ... island`)
runs P3 nr = 8/16 to 1e8, P2 nr = 16 and N = 2 scans, and D3 ladders.

Limits. The Jacobian holds dense element blocks of size
((n1 + 2 n0)(2N + 1))^2 for every element before compression. Memory and
time therefore grow like N^2 per element, and N >= 3 at nr >= 16 is
cluster work. S and g_par_shape are axisymmetric. The phi aliasing of K
and J is declared, not bounded.

### Landreman family 1 through the coupled solve

Family 1 (uniform iota, p = p_a - 2 psi) is posed as an exact case
(`test_nonnested_landreman`). The derivation of the data:

1. Source. B.grad p = 0, hence K grad p = chi_perp grad p and
   S = -div(K grad p) = -chi_perp lap p = -chi_perp div(J x B), because
   J x B = grad p. S(R, Z, phi) is evaluated by complex step from the closed
   form (rounding accuracy) and declared through `source_3d`.
2. Phi = 0. With g_perp > 0 the current contains -g_perp (I - b b) grad Phi.
   The exact J has no such part, so grad Phi is parallel to b and, for a
   nested state with no parallel potential, Phi is constant and equals the
   wall value 0. The solve starts from Phi = 0 and must keep it small
   (reported).
3. Ohm. J.B = g_par V B_phi/(2 pi R), so
   g_par(x) = 2 pi R (J.B)/(V B_phi). With w = R J.B/B_phi (measured in
   [3.3, 8.6] on the domain, positive, also for eps up to 0.25) take
   g_par = g0 w/w_ref (`g_par_shape_3d`) and V = 2 pi w_ref/g0. The shape
   carries the Pfirsch-Schlueter part, so no Ohmic mismatch remains.
4. I_phi is the phi mean of int J_phi dR dZ of the exact state on the
   discrete domain (32 phi points). V is the unknown and is compared with
   2 pi w_ref/g0.
5. Wall data are the L2 projection of the exact potential onto V1 in every
   Fourier block (gauge A_z = 0) and the nodal projection of p. The seed
   shifts the interior by 0.05 (1 - u^2) in R with V = Phi = 0.
6. Truncation. B_phi = a b F/R; the Fourier content of B decays like
   (eps/2)^(n/2): n = 2/4/6 carry 1.2e-2/1.9e-4/3.8e-6 at eps = 0.02.
   Errors are measured against the exact state projected on modes 0..N; the
   untruncated difference is printed as "truncation". At a large
   chi_par/chi_perp the truncation of b.grad p (not zero for truncated B)
   is amplified by the ratio, so N = 4 is needed there.

Disk wall R_c = 0.99, r_w = 0.3, eps = 0.02, p_a = 0.25, F_vac = 1,
w_ref = 5. The case needs the hooks `nc_problem_t%source_3d` and
`%g_par_shape_3d` (`nc_scalar3_i`: x(2), phi); exactly one of source and
source_3d, at most one g_par shape.

Gates (pre-declared, not loosened): B order >= k - 0.3, D3 order >= k - 0.3
with D3 falling, loop-voltage error order >= k - 0.5 on the nr ladder, and
pressure order >= k + 0.6 at ratio 1 only. CTest
`nonnested_coupled_landreman_family1_p2`: P2, ratio 1, N = 2, nr = 4/8/16.
Ratios 1e3 to 1e6, N = 4 and P3 are in the sCluster study (`landreman`
mode); a failing gate there is reported, not retuned.

Continuation in the anisotropy (`nc_continue`). The first sCluster study
(jobs 3222997/8) gave no result at ratio 1e6: from the shifted seed the
Newton iteration took damped steps (line-search factors 1/2 to 1/256)
with the backward error stuck between 1e-6 and 1e-7. This is a
globalization failure, not a linear-algebra floor. Locally, P2 N = 2 nr = 4
eventually entered the quadratic basin at step 19 and reached 2e-16 with
the same direct LU (the old print overflowed after 15 entries). N = 4 was
still at 1e-6 after 30 steps. The pressure rows carry chi_par b(A).grad p,
so an error in b enters times chi_par and the basin shrinks like
1/ratio. The discrete solution itself barely moves with the ratio.
`nc_continue` therefore solves at chi_par/chi_perp = rc^s and
g_par/g_perp = rg^s for s from 0 to 1 in steps of at most one decade.
Each stage runs Newton to the final tolerance from the previous state. A
stage that fails, or needs more than 12 steps, is retried with half the
step. Every stage converges to the discrete solution at its own s, so the
end point does not depend on the path. The test checks this on the
coarsest mesh: a second run with 2-decade steps must agree blockwise
within 1e-6 of the block maximum in A, p and V. Phi is printed only,
because it vanishes in the exact state. Local P2 nr = 4 at ratio 1e6:
N = 2 needs 8 stages and 16 Newton steps, and the path differences are
A 1.1e-9, p 9.8e-8 and V 4.8e-9. N = 4 needs 8 stages and 17 steps, with
differences A 1.2e-9, p 1.1e-7 and V 5.1e-9. The Newton-step gate now
applies per stage (at most 12; previously at most 15 for the single
solve). Ratio 1 is a single plain solve.

Pre-declared study gates with continuation (P2 and P3; ratios 1, 1e3 and
1e6; N = 4 and 6; nr = 4/8/16). The gates are unchanged: B and D3 orders
>= k - 0.3 with D3 falling, V order >= k - 0.5 and the path check.
p >= k + 0.6 applies at ratio 1 only. PASS at 1e6 means that these gates
hold at N = 4 for P2 and P3. P3 N = 6 at nr = 4/8 failed the B-order
gate (2.563). The errors there are identical to those at N = 4 on
nr = 4/8, whose 4-to-8 order is also 2.563 while the 8-to-16 order is 2.952.
The failure is therefore pre-asymptotic, and N = 6 now runs on nr = 16 too.
Local P2 N = 4 nr = 4 at 1e6 gives B 5.4e-3 (as at ratio 1), but p 4.1e-3
(5.8e-4 at ratio 1) and Phi 8.5e-2 of V. With iota = 2 every field line
closes, so Phi and p are fixed along b only up to the discrete error of
b, amplified by the ratio. The ladder shows whether B, D3 and V keep
their orders.

Local P2 results (nr = 4/8/16, 4 threads, previous revision; the gates were
then rewritten as above). The Jacobian matches a central difference to
2.7e-11. Newton converges quadratically in 4 steps.

| ratio, N | B | p | D3 | V |
|---|---|---|---|---|
| 1, 2 | 5.3e-3, 1.6e-3, 4.0e-4 (order 2.00) | 5.8e-4, 7.1e-5, 1.1e-5 (2.76) | 1.0e-3, 1.8e-4, 2.7e-5 (2.73) | 4.8e-4, 1.2e-4, 3.1e-5 (2.0) |
| 1e3, 2 | as above | 9.3e-4, 3.1e-4, 2.3e-4 (stalls) | as above | |
| 1e3, 4 (nr 4/8) | | 1.0e-3, 2.0e-4 (2.36) | | |

Phi stays at 4e-7 to 9e-6 at ratio 1.

sCluster study (job 3223744 at f08dbb1 after the full CTest job 3223743,
181/181 passed; 16 cores; nr = 4/8/16, orders from the finest pair 8/16).
Stages and Newton steps are given per nr; time is wall time and memory is
the maximum resident set of the whole ladder.

| P, ratio, N | orders B / p / D3 / V | stages (Newton) nr 4/8/16 | time | memory | gates |
|---|---|---|---:|---:|---|
| 2, 1, 4 | 2.02 / 2.76 / 2.73 / 1.96 | 1 (4) / 1 (3) / 1 (3) | 1 min 08 s | 8.9 GB | PASS |
| 2, 1, 6 | 2.02 / 2.76 / 2.73 / 1.96 | 1 (4) / 1 (4) / 1 (3) | 2 min 45 s | 17.4 GB | PASS |
| 2, 1e3, 4 | 2.02 / 2.05 / 2.73 / 1.96 | 4 (10) / 4 (6) / 4 (6) | 2 min 27 s | 9.4 GB | PASS |
| 2, 1e3, 6 | 2.02 / 2.05 / 2.73 / 1.96 | 4 (10) / 4 (7) / 4 (6) | 6 min 50 s | 18.8 GB | PASS |
| 2, 1e6, 4 | 2.07 / 4.65 / 2.75 / **1.07** | 8 (17) / 24 (193) / 8 (8) | 13 min 27 s | 9.5 GB | FAIL (V) |
| 2, 1e6, 6 | 2.09 / 5.15 / 2.76 / **0.92** | 8 (18) / 29 (253) / 8 (7) | 43 min 26 s | 19.1 GB | FAIL (V) |
| 3, 1, 4 | 2.95 / 3.89 / 3.56 / 2.96 | 1 (3) / 1 (3) / 1 (3) | 4 min 08 s | 25.7 GB | PASS |
| 3, 1, 6 | 2.95 / 3.89 / 3.56 / 2.96 | 1 (4) / 1 (3) / 1 (3) | 10 min 24 s | 52.6 GB | PASS |
| 3, 1e3, 4 | 2.95 / 4.35 / 3.56 / 2.98 | 4 (6) / 4 (6) / 4 (3) | 8 min 33 s | 25.6 GB | PASS |
| 3, 1e3, 6 | 2.95 / 4.02 / 3.56 / 2.97 | 4 (7) / 4 (5) / 4 (3) | 21 min 05 s | 52.3 GB | PASS |
| 3, 1e6, 4 | 2.95 / 4.35 / 3.56 / 2.98 | 8 (7) / 8 (6) / 8 (3) | 12 min 19 s | 25.7 GB | PASS |
| 3, 1e6, 6 | 2.95 / 5.52 / 3.56 / 2.81 | 8 (8) / 8 (6) / 8 (3) | 32 min 36 s | 52.3 GB | PASS |

Path checks (nr = 4, 2-decade steps) agree within 2e-7 in A, p and V in
every run; Phi differs by up to 6.8e-5 (P3, 1e6, N = 6), printed only. The
P3 N = 6 prediction holds: with nr = 16 the B order is 2.95 at ratios 1,
1e3 and 1e6, so the earlier 2.56 was pre-asymptotic.

Verdict at ratio 1e6: FAIL by the pre-declared rule, since P2 fails the
V-order gate (>= 1.5) at N = 4 (1.07) and at N = 6 (0.92). P3 passes every
gate at N = 4 and N = 6. B and D3 keep their orders at P2 too. The P2
failure sits at nr = 8. There Phi is 0.79 (N = 4) and 0.92 (N = 6) of V,
p is 2.0e-3 (2.0e-4 at 1e3), and the continuation needed 24 and 29 stages.
At nr = 16 the V error (3.07e-5) equals that at ratio 1 (3.13e-5), and
Phi drops to 9e-3. The V error at nr = 8 (6.4e-5) is below the ratio-1
value (1.2e-4), so the 8-to-16 ratio of 2.1 does not measure convergence.
The cause is the closed-field-line pollution predicted above. With
iota = 2 the discrete b.grad has a near-kernel of functions constant on
closed lines. The P2 error in b on the coarse mesh, amplified by 1e6,
leaves p and Phi in that near-kernel almost undetermined, and Phi feeds V
through Ohm's law. At P3 the B error at nr = 8 is 5.7e-5 instead of 1.7e-3. That is
enough at 1e6, and Phi stays below 8.1e-4 of V. Proposed next step, not a retune: add a gauge for p and Phi
along closed lines. One option is a weak constraint that sets the
line average of Phi to its wall-connected value. The other is
discretising b.grad with the same discrete B used in Ampere, so that the
exact discrete equilibrium satisfies b.grad p = 0 exactly. Rerun P2 1e6
N = 4 on nr = 8/16/32 against the same gates.

Closed-line follow-up (sCluster jobs 3224041 to 3224043 at 693f56c, which
only adds the coarsest nr as optional argument 6; 16 cores). The second
option above is already the discretisation: b in the pressure and Phi rows
is built from the same V2 coefficients as curl A in the Ampere rows. Exact
b.grad p = 0 is still out of reach, because P_k on unaligned triangles
cannot be constant on discrete field lines. The mismatch therefore enters
the current equation tested with a function constant on a closed line: the
discrete B violates the closed-line condition (oint dl/B constant on
surfaces) that family 1 satisfies exactly, and only g_perp grad Phi can
carry it, so Phi grows like g_par/g_perp. All runs P2, N = 4, gates
unchanged; Phi is relative to V.

| ratio, nr ladder | orders B / p / D3 / V | V error per nr | Phi per nr | gates |
|---|---|---|---|---|
| 1, 8/16/32 | 2.51 / 2.65 / 2.90 / 2.25 | 1.22e-4 / 3.13e-5 / 6.59e-6 | 3e-6 / 8e-7 / 4e-7 | PASS |
| 1e3, 8/16/32 | 2.51 / 1.64 / 2.90 / 2.25 | 1.21e-4 / 3.11e-5 / 6.52e-6 | 4.2e-4 / 1.0e-4 / 4.7e-5 | PASS |
| 1e6, 8/16/32 | 2.51 / 2.36 / 2.90 / 2.24 | 6.43e-5 / 3.07e-5 / 6.52e-6 | 0.79 / 9.0e-3 / 4.7e-5 | PASS |
| 1e6, 6/12/24 | 2.15 / 3.54 / 2.81 / **0.40** | 7.56e-5 / 2.40e-5 / 1.82e-5 | 0.69 / 0.26 / 1.5e-3 | FAIL (V) |
| 1e4, 4/8 | 1.75 / 2.44 / 2.53 / **1.49** | 3.42e-4 / 1.22e-4 | 1.6e-2 / 5.1e-3 | FAIL (V) |
| 1e5, 4/8 | 1.75 / 1.98 / 2.53 / **0.25** | 1.28e-4 / 1.08e-4 | 5.9e-2 / 8.0e-2 | FAIL (V) |

Wall time and memory: 4 min 14 s, 10 min 15 s and 34 min 40 s with 36 to
39 GB for the three 8/16/32 ladders; 15 min and 23 GB for 6/12/24. Newton
converged to 1e-11 in every stage, and the path checks on the coarsest
mesh agree within 3e-10 in A, p and V.

The Phi absorption is confirmed. At nr = 8 Phi grows almost linearly with
the ratio (4.2e-4, 5.1e-3, 8.0e-2, 0.79 for 1e3 to 1e6) while B and D3 do
not change, which is the g_par/g_perp amplification of a fixed closed-line
mismatch. By the rule declared for this rerun (V order >= 1.5 on 16/32 and
Phi at nr = 32 not above that at nr = 16), P2 1e6 N = 4 PASSES on
nr = 8/16/32 and FAILS on 4/8/16. The 6/12/24 ladder shows that the
pollution is not merely pre-asymptotic: Phi is not monotone in nr
(0.69 / 0.79 / 0.26 / 9.0e-3 / 1.5e-3 / 4.7e-5 at nr 6/8/12/16/24/32), and
at nr = 24 the V error is still about 1.6 times the clean nr^-2 value. It
does not decay like h^k either: at nr = 32 the 1e3 and 1e6 states agree to
the printed four digits in B, p, V, D3 and Phi, so there is no visible
ratio amplification at nr = 32, while it is large at nr = 24. The amplitude depends on
the discrete line closure for each mesh, not on h alone. The P2 verdict at
1e6 thus depends on the ladder. P3 (all ladders above) is unaffected.
Making P2 robust at 1e6 needs a design decision: a closed-line gauge for
Phi paired with a matching source so that weak Ampere stays solvable, or a
ratio-dependent mesh. Neither is a retune of the gates.
Decision (2026-10-14): NN2 comparisons at conductivity ratio >= 1e6 (NTEC)
use P3; P2 robustness at high ratio needs a closed-line Phi gauge with a
matching source and is not pursued now.

### NN2c island results (sCluster jobs 3222014, 3222015)

All 180 CTests passed. G (flattening measure) from the island study, N = 1
unless noted:

| run | G at ratio 1 / 1e2 / 1e4 / 1e6 (1e8) | island width at the top ratio (u) |
|---|---|---|
| P3 nr 8 (to 1e8) | 0.982 / 0.973 / 0.948 / 0.143 (0.148) | 0.174 (1e8) |
| P2 nr 16 (to 1e6) | 0.982 / 0.973 / 0.945 / 0.212 | 0.188 |
| P3 nr 8, N = 2 (to 1e6) | 0.982 / 0.973 / 0.976 / 0.090 | 0.189 |

D3 ladder, P2, ratio 1e4, nr = 4/8/16/32: 5.8e-3, 9.9e-4, 1.8e-4, 1.2e-4
(finest order 1.86; the 8 to 16 order is 2.4). The memory peak is 13 GB at
nr = 32. At 1e6 P3 nr = 8 N = 1 gives G 0.143 and D3 6.7e-5. Two runs
stopped with `ERROR STOP vacuum mode` (the vacuum-seed mode check) and are
not results: P3 nr = 16 to 1e8, and the P3 D3 ladder at 1e6 after nr = 8.
P3 nr = 16 is therefore open.

## Reusable island diagnostics and vacuum verification

[V] `kin6d_island_diagnostics` takes any `equilibrium_field_t`. It traces
`d(R,Z)/dphi = R (B_R,B_Z)/B_phi` with FortNum Cash--Karp, either adaptive
with declared absolute/relative tolerances or on an explicit fixed-step
ladder. This chart requires nonzero B_phi and R > 0 along the whole trace.
Sampler failures and exhausted ODE budgets propagate as unsuccessful status.
The section angle, number of full toroidal turns, centre and helicity `(m,n)`
are inputs. Phase means `m atan2(Z-Zc,R-Rc) - n phi`, modulo 2 pi, rather
than an undeclared device angle. The residue is Greene's `(2-tr DP^N)/4`.
The determinant is reported as a check of flux-area preservation at a
periodic point, without assuming canonical (R,Z) coordinates away from it.

[V] Fixed points use FortNum Newton with a centred map Jacobian at a
declared displacement. The symmetry fallback bisects the half-period Z
return on Z=0 only when the caller supplies a symmetry-line bracket. A symmetry seed that does not close is refined by complete-map Newton. Both
routes are checked against the **complete** N-turn closure. A reflected
continuous coil does not guarantee a reflected triangulation, so a small
half-turn residual alone never establishes an FEEC periodic orbit.

[V] The previous coarse-map limitation has two identifiable ingredients.
`nn_point` selects one element at
`src/mhd/equilibrium_nonnested.f90:299` and evaluates its H(div) basis at
lines 301--310. Tangential components can jump across a face. The old
fixture advanced fixed Cash--Karp stages at 512 steps/turn, while the
pinned FortNum `fortnum_multiroot.f90:650` chose a perturbation of
`epsilon^(1/3) max(1,abs(x_j))`. Stage/face crossings can therefore dominate
the differenced map. The retained `coarse` experiment compares the old
root and its Jacobian on 512/1024-step schedules before testing adaptive
tracing with a larger declared Jacobian displacement. This bounds the
claim to measured map closure, rather than asserting a differentiable
continuous FEEC flow across every face. Full field/observable derivatives
across topology changes are unavailable.

[V] The controlled-trace P3 nr=8 vacuum roots close to 2.5e-7/1.2e-6 m
(O/X). More importantly, `test_nonnested_island map 3 8 1 1e6` checks the
actual finite-pressure coupled field: O/X closures are 9.0e-7/5.5e-7 m,
with residues 0.552/-0.815. The original symmetry-only measurements on
that state close only to 1.8e-4/2.7e-4 m. Thus nr=16 is not an intrinsic
Newton threshold. On the vacuum field the legacy root itself has a tiny
fixed-schedule residual, but doubling 512 to 1024 steps changes the map by
6.4e-6 m and its small-displacement Jacobian by 0.047. A tiny root residual
on that schedule is therefore insufficient evidence of tracing accuracy.

[V] Position, residue and phase indicators compare two solves with tracing
tolerances divided by eight, maximum step halved and map displacement
halved (or a doubled fixed-step count). They are practical refinement
indicators, not certified error bars. A radius/chord label and its units are
mandatory for width reporting. Width currently means finite-time helical
libration on a declared straight section chord. Its separate indicators
are trace refinement, doubled observation time and the radial bracket
resolution. It does **not** claim an invariant separatrix in a stochastic
layer. The exact pendulum regression checks a true separatrix against this
finite-time diagnostic, N=2 points, phase, residue, determinant and error
statuses. Its square-root scaling is diagnostic evidence, not NN2 physics.

[V] NN1 retains its factor for up to three physical Newton steps when the
first floating-point update misses its residual gate. The physical residual
gate remains 1e-9 and the energy gate is checked after every update. The P3
nr=16 three-mode seed is a separate regression for the observed mode-1
residual rejection at 1.189e-9. No new sparse kernel or backend is introduced.

[V] Rung 0 uses the published synthetic NN1 192-segment filament and its
annulus. `kin6d_coil_quadrature_field` integrates the exact polygon source
with FortNum Gauss--Legendre nodes. The 8/12-node quadrature and the
independent closed straight-segment formula are compared at 64 points.
The P2/P3 ladders print `h=0.9/nr`, width, absolute error against the
quadrature trace, phase error and the individual tracing indicators. Their
fixed 32-turn definition is retained when comparing meshes. The width ladder uses 4096/8192 fixed Cash--Karp steps per turn and 16 radial
bisections. The width oracle uses the independent closed segment field after
the quadrature comparison, while periodic points use adaptive tracing. This
bounds work on long FEEC trajectories, without claiming that fixed steps
supply an enclosure. Tracing indicators still have
to resolve a reported spatial order. The 64-turn
change is a separate observation-window uncertainty, which can exceed a
spatial error and preclude a separatrix convergence claim. Tests with
substantial tracing or coupled solves are labelled `slow`.

Chris&AI

The [production scaling ledger](benchmarks.md#production-estimator-scaling)
uses an NN1 magnetic state, constant heating and chi_parallel/chi_perp = 1000
through `nt_solve(..., estimate_error=.true.)`. It verifies equality with the
separately timed estimator. The production efficiency index is unavailable
without a true-error oracle; these timings retain the practical claim.
