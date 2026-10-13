# Non-nested equilibrium: NN2 (resistive relaxed lane)

Status: design decision, first increment (transport for a given B) and
second increment (axisymmetric coupled Newton), 2026-10-13. NN2 is the
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
for finite chi_perp/chi_par keep a pressure gradient ([Fitzpatrick 1995]
(https://doi.org/10.1063/1.871434); its width scaling is not yet checked
against the paper body and is not used as a gate).

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
4. Finite-pressure island from an offset start (next). The NN1
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
| Zhang, Zhu, Hegna ([CPC 2026](https://doi.org/10.1016/j.cpc.2025.109961)) | hyperbolic parallel damping for p, dynamic resistive relaxation for B | their published case verbatim; parameter map from the paper body, not yet read |
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
backward error |r|/(|J| m), where m is |x| plus the block maximum. Tangents
and adjoints of a converged state are solves with that Jacobian.

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

### Landreman family 1: what is missing

Family 1 (uniform iota, p = p_a - 2 psi) is not reachable in this
increment. Three items are open:

1. FortFEM has the Fourier-coupled V0 diffusion but no curl-complex
   counterpart. Needed are real-Fourier V1/V2 forms with a phi-dependent
   pointwise tensor (the mixed V1-V2-V0 blocks above for modes 0..N,
   nphi >= 2N + 1) and the n > 0 gauge R A_phi = 0 inside that coupled
   system. NN1's 3D solves are mode by mode and linear. This is generic
   FEM and belongs upstream.
2. Ohmic compatibility. As for Solov'ev, the exact family 1 current must
   be Ohmic for a declared g_par(x), V and Phi. Its parallel current has
   not been derived; the case cannot be posed until it is.
3. Boundary data. The NN1 domain is an axisymmetric region in (R, Z).
   Family 1 surfaces are not axisymmetric, so the boundary is not a flux
   surface. It needs the exact normal B trace and nonconstant Dirichlet
   p_b and Phi data from the closed form, with S = -chi_perp lap(p) inside.

Chris&AI
