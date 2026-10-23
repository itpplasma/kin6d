# GLL nested-map experiment and port decision

Status: opt-in axisymmetric operator spike, not a production equilibrium solver.
KIN6D base 22d6cb5; FortFEM ac673dc6 (#110 MF1, #109 SPL1 and #117 IGAAX); FortNum
1cfe3cba. The ordinary triangle path and its dependency pins are unchanged
when `KIN6D_BUILD_GLL_SPIKE=OFF`. The option selects the stated FortFEM revision.

## Mathematical contract

The independent coordinates are a regular disk atlas `(s1,s2)` and periodic
`phi_hat`. The unknown physical map is `X=H(q)`, with
`q=(q1,q2,phi_hat+omega)`. Replace the meridional triangular interpolation by
continuous, isoparametric tensor GLL interpolation on five O-grid patches.
The O-grid changes the carrier, not the frozen-in physical model.

The reference data are the closed contravariant flux density `b` and material
pressure `p`. For an element geometry `y=G(s1,s2)`, freeze

```
c = DG^T grad_y chi
b = (-c2, c1, b3_y det DG)
J = sigma det DX > 0
B = DX b / J
W = integral [ |DX b|^2/(2 mu0 J) - p J ] ds1 ds2 dphi_hat .
```

`p` is the gamma=0 material pressure used by the existing map, rather than a
new pressure equation. The sign `sigma` is the chart orientation: cylindrical
`(R,Z,phi)` has sigma=-1. KIN6D's authoritative energy is the shared
`frozen_energy` imported by `gen/kin6d_symbolic_ideal.f90:1` and used in
`gen/gen_ideal_map.f90:49`. The reference-form pullback is already generated
there (`gen/gen_ideal_map.f90:109`). Existing freeze and energy integration
are `src/mhd/ideal_map.f90:440` and `:976` on this base.

Write the fixed gauge expansion as `q=q_ref+E U`, interpolation as `T`, and
pointwise energy as `w`. Then `g=E^T T^T w'` and
`Hv=E^T T^T w'' T E v`, with quadrature weights included. For a general chart,
`T` gathers `q` and its three reference derivatives, and `w''` includes both
`(dDX/dq)^T w_DXDX (dDX/dq)` and the chart's second derivative contribution.
Dropping the latter changes Newton's operator. The current chain rule is
`src/mhd/ideal_map.f90:886`; `:1051` assembles its global Hessian. A matrix-free
path caches quadrature leaves, not element or global Hessian blocks. Its
Euclidean coefficient pairing is the derivative of the quadrature energy.
The implicit tangent/adjoint of a converged constrained equilibrium must use
that stationary equation, not differentiate Krylov or line-search iterations.

For the cylindrical n=0 specialization, rotational invariance eliminates the
values of omega and Z from the energy. The seven inputs are
`(R,R_s1,Z_s1,omega_s1,R_s2,Z_s2,omega_s2)`. At azimuth zero,

```
DX = [ R_s1        R_s2       0
       R omega_s1  R omega_s2 R
       Z_s1        Z_s2       0 ] .
```

`gen/gen_ideal_gll.f90:29` substitutes this into the same `frozen_energy` and
derives the seven-gradient and symmetric Hessian with FortSym. Checked-in
leaves are `src/generated/kin6d_ideal_gll_leaves.f90:11`, `:33`, `:66`.
The nonlinear R factors and all their mixed Hessian terms are retained.

## Unknowns, axis and constraints

Production should preserve KIN6D's reduced unknowns: transversal displacement
nu and omega per node/mode, plus the global axis-expansion coefficients.
`ideal_map.f90:210`, `:217`, `:224` and `:263` define indexing, reduction and
expansion. The frozen transversal is `grad chi/|grad chi|` away from the axis;
it is singular at the O-point. A regular O-grid makes the axis an ordinary
interior mesh point, but does not remove this gauge issue. Preserve the
smooth axis corrections and rank filtering (`:575`, `:725`), or establish a
verified smooth quotient-space alternative before replacing them. Never
impose a polar-axis Dirichlet condition on q1 or q2.

Fixed-boundary production holds nu on the wall and preserves the frozen
fluxes and pressure law; it must not re-freeze after a map update. The rigid
omega rotation is pinned (`:427`). At n=0, the additional surface-function
nulls `omega=f(chi)` need the existing layer moments (`:499`) and their
projection (`:554`). In the nested mode, omega is eliminated at each trial
(`:1451`), so the outer action is the constrained Schur complement
`H_nunu - H_nuomega H_omegaomega^-1 H_omeganu`, not just `H_nunu`.
A matrix-free implementation needs an inner omega solve, projected moments,
and a compatible preconditioner. An ordinary scalar diffusion V-cycle does
not supply that operator.

The spike instead uses all nodal `(R,Z,omega)` coefficients, holds R and Z on
the wall, and pins one omega value. This is an unreduced physical-energy
operator check. Its identity damping regularizes remaining near-null
relabellings; it does not reproduce the nested reduced functional, its
moment constraints or its nonlinear iteration count. Wall tangential
reparametrization is also held in this experiment. These limitations are
intentional and must not be hidden in a full-port claim.

## Existing library seams and ownership

All FortFEM paths below refer to ac673dc6 unless a branch is named.

| Piece | Available public seam and evidence | Remaining work |
|---|---|---|
| GLL nodes, rules, tables | `src/discretisation/gll_line.f90:8` | Application chooses overintegration; collocation is not exact for rational map energy. |
| Regular disk/torus O-grid | `tensor_mesh.f90:242` `ogrid_tensor_mesh(dim,n,p,...,q)`; five patches, outer boundary only at `:332` | KIN6D composes the physical seed geometry and frozen profiles. |
| Gather/transpose contractions | `tensor_contractions.f90:7` `tensor_contract` | No private sum factorization in KIN6D. |
| Quadrature geometry, colouring | `tensor_operator.f90:71` `build_tensor_operator`; points/J/weights/colour arrays at `:27` | Public implementation data is used by this pinned spike; a stable coupled value/gradient callback is desirable. |
| Ready-made MF apply | `tensor_operator.f90:277` `apply_tensor_operator` | Scalar SPD diffusion plus reaction only; not the nonlinear map's coupled Hessian. |
| p-MG | `tensor_solver.f90:15`, `:49`, `:393` expose `build_tensor_multigrid`, `tensor_vcycle` | Scalar diffusion hierarchy, not block/constraint/Schur p-MG. |
| Mixed line hierarchy | Main #117: `line_hierarchy.f90:23`, `:29`; `ogrid_line_hierarchy`, `box_line_hierarchy` | Merged during GL1. Adds GLL/spline/Fourier routes and SPD coarse factor (`tensor_solver.f90:203`); no map physics. |
| Nonlinear driver | NEQ1 #115, b02d2ef1: `src/discretisation/newton_krylov.f90:23`, `:54`, `:122` | Draft; uses FortNum FGMRES. Generic driver move tracked by fortnum#110. Needs map energy/fold-aware globalization and constraint semantics. |
| Krylov | FortNum `src/linalg/fortnum_krylov.f90:69` CG, `:704` FGMRES | Procedure CG at pinned 1cfe3cba suffices; #111 object API is not required by this spike. |
| Curved unstructured mesh | UNS1 #114, 9d09504f: `tensor_mesh_gmsh.f90:14`, `island_quad_mesh.f90:21`, `tensor_mesh_tools.f90:11` | Optional future island carrier, not required for this nested O-grid. |

KIN6D owns pointwise energy and derivatives, charts, material profiles,
flux/pressure constraints, gauge choice and physical error checks, with
continuous leaves retained in `gen/`. FortFEM owns tensor interpolation,
assembly, topology and multilevel algebra. FortNum owns generic Krylov and
Newton-Krylov. NEQ1's scalar Grad-Shafranov speedup is evidence for its own
problem, not a measured speedup of this constrained map.

## Third direction and nonlinear solve

For genuinely three-dimensional states, use local periodic GLL cells in phi
by default; Fourier remains an option for low-mode/exact-symmetry cases.
The current map uses `2*ntor+1` real Fourier coefficients and a periodic
quadrature (`ideal_map.f90:349`, `:840`). A GLL phi carrier must interpolate
all three components and include their phi derivatives, preserve field-period
seams and chart orientation, and control nonlinear aliasing by quadrature
refinement. MF1's hard-coded 3D torus geometry is a fixture, not a seed for an
arbitrary equilibrium. IGAAX's `map` callback is the more general geometry
seam. The n=0 omega elimination is explicitly rejected for ntor>0.

The production step must retain a fold rejection at J<=0, a displacement and
rotation cap, damping with `D=|diag H|+epsilon*max|diag H|`, and the
roundoff-aware energy acceptance in `ideal_map.f90:1340`, `:1353`.
The spike solves `(H+mu I)s=-g` once using FortNum CG, with a unit-normalized
rhs because its tolerance at `fortnum_krylov.f90:110` is relative to
`max(||b||,1)`. It checks the actual residual afterwards. Negative curvature
at `:149` yields only breakdown, so a Steihaug policy needs an upstream
status/direction API. NEQ1's residual-norm Armijo alone is insufficient for
the map's near-null energy directions and tight seed-return ladders.

## Reproduction and evidence

Configure the ordinary project with `-DKIN6D_BUILD_GLL_SPIKE=ON`, then build
`benchmark_ideal_gll`. The option requires OpenMP and the usual native
FortSparse AMD setup. Run, for example,

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/benchmarks/gll/benchmark_ideal_gll 6 1 1000 0.01
ctest --test-dir build -j4 -R '^(ideal_gll_derivatives_p[468]|ideal_map_derivatives_and_rotation_null|ideal_map_curve_chart_landreman|kin6d_layering)$' --output-on-failure
```

Arguments are degree, cells per patch direction, timing repetitions and
positive damping. The scalar DOF count is `5(np)^2+2np+1`, tripled for the
map. Isoparametric geometry has q=p. Field quadrature uses `2p+3` Gauss points
per direction; p-MG's surrogate uses its own default rule. Generic machinery
is compiled directly from the pinned FortFEM sources. The benchmark has no
assembled map Hessian, private multigrid or private Krylov implementation.

The case is axisymmetric Landreman, eps=tau=0, delta=1/64, with the analytic
field frozen at identity. This is not the eps=0.5 three-dimensional L cold
start behind INT5's 0.445 s and VMEC++'s 0.029 s. Exact energy is independently
computed from `nested_oracle_t%averages`: `V*(<B^2>/2-<p>)`.
GLL energy refinement, triangle energy refinement and finite differences
establish operator agreement; they do not establish matched B/p/ds accuracy
of a newly solved equilibrium. Continuous bounds and implicit sensitivities
are unavailable in this experiment.

The numerical tests compare energy gradients and second directional energy
changes with fourth-order and centered finite differences, respectively, the full HVP with gradient
finite differences, and Hessian bilinear symmetry. They check masking on
both fixed rows and fixed columns, a true linear residual below 2e-8, an
unfolded Newton trial and a decrease in energy. These are independent
numerical checks, not a proof of global positive definiteness.

Raw scripts/logs belong outside the source diff. The timing table below
identifies whether the node was shared or exclusive.

## Measurements and decision

[V] sCluster pilot 3357100, shared node16, EPYC 9455, GNU 12.2, `-O3` without
`-march=native`, OpenBLAS 1 thread, three repeats, medians below. These are
provisional timings on a shared node. The exclusive node12 sweep is submitted
separately; no completed exclusive timing or VMEC++ rerun is claimed here.
All rows use mu=0.01 and CG tolerance 1e-8 on a unit rhs. MG is one FortFEM
V-cycle per R/Z block and a reaction-regularized, unpinned omega surrogate,
restricted to the physical omega pin. It is not a Hessian Galerkin hierarchy.

| n | p | map DOFs | identity energy relative error | HVP (microseconds) | CG iterations, MG / identity | one CG solve (s), MG / identity |
|---|---|---|---|---|---|---|
| 1 | 4 | 267 | 8.81e-7 | 53.6 | 425 / 402 | 0.110 / 0.0237 |
| 1 | 6 | 579 | 1.84e-9 | 109 | 708 / 623 | 0.376 / 0.0697 |
| 1 | 8 | 1011 | 1.13e-12 | 184 | 905 / 796 | 0.691 / 0.150 |
| 2 | 4 | 1011 | 9.31e-9 | 212 | 796 / 731 | 0.672 / 0.159 |
| 2 | 6 | 2235 | 1.30e-13 | 411 | 946 / 926 | 1.64 / 0.395 |
| 2 | 8 | 3939 | 6.14e-14 | 696 | 1177 / 1094 | 3.03 / 0.788 |

[V] All 18 runs satisfy the true residual check and decrease the physical
energy without folding. Scalar p-MG is slower than identity preconditioning
and does not give p-independent iterations for this operator. High-order
energy convergence is established independently: triangle P2 at boundary
counts 8/16/32 has relative errors 8.17e-4/5.18e-5/3.25e-6; P4 ends at
3.25e-10 and P6 at 1.86e-14, converging to the same exact energy. The
p=4/6/8 finite-difference, symmetry, masking and descent tests and three
existing map/chart/layer checks pass (six registered CTests). FortSym
regeneration at 8efec5bf reproduces the checked-in leaves byte for byte.

[H] A planning model for a future L solve is

```
T = Tsetup + NNewton * (Tlinearize + Ntrial * Tenergy_gradient
                      + Ncg * Titeration) .
```

[V] For n=1,p=6, measured setup including both MGs is 4.31 ms, linearization
0.226 ms, energy+gradient 0.0806 ms, and the amortized MG-CG iteration costs
0.531 ms, including about 0.418 ms of preconditioning. [H] Borrowing INT5's
12 Newton steps and assuming two energy trials per step gives 4.52 s with
the measured 708 iterations, or 0.843 s with identity preconditioning and
its measured 623 iterations. These are cost scenarios, not predictions of
a converged GLL equilibrium. Even a hypothetical seven-iteration MG solve
would cost about 0.054 s under those assumptions; the 0.029 s budget permits
about three iterations at the present n=1,p=6 costs. For p=4/p=8 the same
MG scenario is 1.32/8.30 s with the measured iteration counts.

There is no verified end-to-end time at matched L B/p/ds error. The reference
0.445 s KIN6D and 0.029 s VMEC++ times are for eps=0.5 and ntor=4, whereas
this experiment has eps=0 and no phi coupling. The production reduced gauge,
LM metric, inexact forcing and trial counts also differ. Matching one energy
integral is insufficient to select the spatial/phi resolution at the B/p/ds
budget. A defensible matched-error estimate requires a converged GLL state,
common physical samples and the existing INT5 comparison on the same node.

Decision [V on the spike, H on the port]: retain GLL p=6-8 as the carrier
candidate, but do not replace the triangle solver or promise VMEC-class time
from the scalar GS benchmark. The immediate gate is a constraint-compatible
coupled preconditioner and reduced-space solve. Kernel and setup improvements
remain useful, but cannot remove hundreds of near-null Krylov iterations.

## For FortFEM lab

- `tensor_operator.f90:277` applies scalar diffusion; `:27` exposes line and
  geometry data used by GL1. Provide a public batched multi-field value/gradient
  gather, generated pointwise residual/JVP leaf and transpose-scatter seam.
  GL1's seven-channel schedule is `benchmarks/gll/ideal_gll.f90:98`, `:173`,
  `:216`; it calls upstream `tensor_contract` rather than copying it.
- `tensor_solver.f90:49`, `:393` expose scalar diffusion p-MG. A production
  consumer needs block/constraint-aware transfer and smoothing, an omega
  Schur solve with moment constraints, and frozen-operator reuse across LM
  shifts. KIN6D supplies its gauge and physics; the reusable algebra belongs
  upstream. The measured hundreds of iterations are an observed limitation
  of this surrogate, not evidence that all p-MG formulations will fail.
- NEQ1 b02d2ef1 `newton_krylov.f90:54`, `:122` still needs landing and the
  FortNum generic-driver move (#110). Add driver hooks for energy acceptance,
  fold rejection and damping/diagonal updates before the full nested map
  consumes it. FortNum `fortnum_krylov.f90:149` also needs a distinct
  negative-curvature return if a Steihaug globalization is selected.
- PERF1: `tensor_operator.f90:174` samples coefficients in the serial setup
  loop; `:295` and `tensor_solver.f90:411` allocate apply/V-cycle workspaces.
  Measure/fix generic setup and workspace costs there. IGAAX's SPD coarse
  factor is already present at `tensor_solver.f90:203`; do not repeat the
  older coarse-LDLT diagnosis as a current-main defect.
- INT1 must integrate GEN1's emitted contractions and arbitrary pointwise
  leaves into the multi-field seam, including the value channel and energy
  gradient, rather than optimizing diffusion HVPs alone. The GEN1 emitter is
  already available through FortSym/FortGen; no private emitter or tensor
  contraction implementation belongs in KIN6D.
