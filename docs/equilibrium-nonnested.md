# Non-nested equilibrium: NN1

## Pressure decision (before implementation)

NN1 is a bounded entry to the free-B lane of the common equilibrium core.
It implements vacuum, p = 0, in three dimensions and the axisymmetric
constant-p-prime, constant-F restriction for the Solov'ev reduction gate.
The latter has p = p' (psi - psi_boundary) and B.grad p = 0 identically.
It is not a general three-dimensional finite-pressure model. Nonzero p'
with nonzero Fourier modes is rejected.

For later finite-pressure ideal relaxation, pressure must be transported
with the admissible variations or constrained weakly along field lines.
Such models admit islands only with compatible pressure and force balance;
a smooth pressure is constant on an ergodic connected component.
Parallel diffusion introduces a distinct transport model. Stepped pressure
with conserved regional helicities is MRxMHD, a different problem requiring
interfaces and traction balance. Zero beta admits vacuum and force-free
fields, but NN1 solves only the vacuum subclass, not general Beltrami fields.
These choices follow D2/D3 and R12; no arbitrary smooth 3D pressure is supplied.

## Discretisation and solver decision

Use the cylindrical chart (R, Z, phi) and FortFEM's existing meridional
complex times Fourier modes. Physical poloidal flux density R B_pol is in
RT; B_phi is a density-mapped DG field. The magnetic state is the closed
2-form B = curl_n A plus the explicitly declared harmonic F_vac/R e_phi.
The auxiliary A is internal. For n > 0 the gauge R A_phi = 0 is exact;
for n = 0 use the scalar R A_phi with Dirichlet data. Thus div_n B = 0 is
an algebraic complex identity, not a penalty or a small residual target.

The first solver minimises the quadratic magnetic energy with prescribed
boundary traces. At n = 0 a linear pressure load gives the axisymmetric
constant-p-prime energy. Its exact Newton step is the sparse FEEC solve;
an energy and equation-residual report checks the offset-start return.
The boundary trace and harmonic toroidal flux are fixed. Helicity is not
constrained: vacuum energy minimisation is not topology-preserving ideal
magneto-friction, nor MRX's full admissible-variation algorithm. An open
boundary with externally imposed currents need not conserve helicity.

Do not copy generic FEEC, Piola, sparse factorisation or ODE code into KIN6D.
The lane owns pressure restrictions, boundary/coil composition and provenance.
The common D3 Maxwell-stress residual is the physics check; the n = 0
implementation is reused directly. A three-dimensional D3 contraction
will use the same stress tensor, not a new equilibrium equation.

## Verification plan and limits

1. N = 0 Solov'ev: an independent FEEC energy solve versus the existing
   KIN6D GS solver and the analytic physical field on mesh ladders; check
   divergence, pressure compatibility, D3 and Newton return from an offset.
2. Vacuum torus with a helical filament inside an excluded coil cavity:
   compare the coil-driven FEEC field and its field-line map with direct
   straight-segment Biot-Savart, including O/X points and island width.
   The cavity makes every represented cell vacuum. Refinement separates
   spatial, Fourier and field-line integration errors.
3. Landreman analytic finite-pressure 3D remains outside NN1's pressure
   contract and is deferred to the next pressure stage.

All NN1 numerical claims are candidate/reference-comparison. No inf-sup or
continuous enclosure has been established. Numerical, coil polygon/Fourier,
model and field-line integration errors are separate; a small residual is
not a field or island-width certificate. M3/M4/M6 are not on the baseline;
NN1 uses prescribed traces and does not claim free-boundary coupling or
full migration completion. The N1 note was read from the repository object
`gs-3d-nonlinear-n1:docs/ideal-map-3d.md`, as it is absent on this branch.

Source: [MRX, arXiv:2510.26986](https://arxiv.org/abs/2510.26986), especially
sections III and IV. The NN1 vacuum minimisation has fewer constraints than
that finite-pressure algorithm and makes no matched performance claim.

## Implemented equations and interfaces

`kin6d_equilibrium_nonnested` stores V2 coefficients for modes 0:N with
B_-n = conjugate(B_n), and the harmonic F_vac/R separately. At n = 0 the
auxiliary scalar a = R A_phi satisfies

    int grad v . grad a / R = mu0 p' int R v,
    B_pol = rot90(grad a)/R,  B_phi = F_vac/R.

Thus the constant-p', constant-F GS restriction is solved independently in
FortFEM's complex; the existing KIN6D GS solver is only an oracle. Its
functional is the convex axisymmetric magnetic energy minus the prescribed
pressure load. At n > 0 the gauge a = 0 leaves a positive definite Nedelec
curl-curl operator. A boundary potential is projected with the upstream
V1 mass matrix; its trace fixes normal B. Newton solves the free increment
with zero boundary increment in one step. For coil inputs, the exact gauge
potential follows A_R = i R B_Z/n, A_Z = -i R B_R/n. Neither an interior
coil field nor a target solution is copied into the solved V2 state.

The solve uses threaded MUMPS when available, SuperLU otherwise. Derivatives
of a fixed-boundary solve are the same equation-level symmetric inverse,
with differentiated pressure and boundary loads; `nn_linear_solve` is the
shared tangent/adjoint action. Geometry derivatives, topology changes and
island classification derivatives are unavailable in NN1. The smooth stress,
pressure, gauge and tracing leaves have generated primal/JVP/VJP families.

`gen_equilibrium_core` now defines one Maxwell-stress tensor, deriving both
the original axisymmetric D3 coefficients and the new 3D leaf from it. The
new D3 evaluator integrates that tensor against zero-trace scalar P_k
components times Fourier modes. Its norm is the sum of scalar cylindrical
H1 seminorms of the three displacement components, weighted by R. It removes
F_vac/R analytically and reports sigma = mean(|B-B0|^2)/mu0 + mean(|p|).
This finite-test-space dual norm is not the continuum dual norm or a majorant.
The edge F consistency is exact for the implemented constant-F harmonic.

The core sampler now accepts (R,Z,phi), with phi = 0 for two-row queries.
The GS adapter ignores phi; the free-B adapter synthesises its modes.
Non-nested fields expose B and p, without pretending to provide global flux
labels or q. The N=0 primitive is retained solely for pressure and reduction
readback. The coordinate inverse used for mesh lookup is reentrant and
independent of magnetic surface labels. Core provenance records the free-B
lane, seed, candidate claim, harmonic data and chart.

## Frozen oracle inputs and error interpretation

Solov'ev: R0=6.2 m, a=0.62 m, F=32.86 T m, mu0=1 in benchmark units,
c=F/(12 R0^3), p'=-16c. The boundary is
R^2=R0^2+2 R0 a cos(theta), R Z=R0 a sin(theta). Three reference meshes
have 8/16/32 radial rings, with P2 or P3 geometry and FEEC order. The start
is a nonstationary interior flux perturbation, not the exact equilibrium.
Fields and flux are checked at 128 independently specified interior points.
The separate GS triangulation has 64/128/256 boundary segments. Field errors
are relative RMS poloidal errors; flux errors are relative RMS flux errors.

Vacuum: R0=6 m, harmonic F=100 T m, one 192-segment closed filament carrying
5 MA, with vertices at phi_j=2 pi j/192,
R_j=R0+0.02 cos(phi_j), Z_j=-0.02 sin(phi_j). The polygon is the declared
physical source, not an approximation claimed exact for a smooth filament.
The represented annulus has minor radii 0.3 to 1.2 m; the coil is inside the
excluded cavity. P2 FEEC/geometry, 4/8/16 radial intervals, 32/64/128 angular
vertices, modes 0:3 and 32-angle boundary Fourier quadrature are used.
No exterior FEM/BEM or free-boundary plasma motion is claimed.

Direct straight-segment Biot-Savart gives the independent field throughout
the cavity-free domain. Its potential-curl identity is independently checked.
Field error at 64 frozen points uses B-B0, not the dominant TF scale. Six
one-turn return maps are compared using FortNum Cash-Karp. O/X points are
roots of the complete two-dimensional return map, with elliptic/hyperbolic
classification from its finite-difference Jacobian; exact reflection of the
triangulation is not assumed. Island width is measured on Z=0, R<R0: the
width of initial radii whose unwrapped theta+phi range is below 2 pi over
the stated trace duration. The O-point is bracketed inside this region.
This is a finite-time libration diagnostic with a scan resolution and a
32/64-turn consistency check, not a certified infinite-time separatrix width.
Direct tracing uses the full polygon field, without Fourier truncation.

The field and map ladders combine spatial and fixed Fourier/input projection
errors. They do not independently bound every contribution. The width gates
are 4 mm between FEEC and direct tracing and 3 mm between trace refinements;
field accuracy must be better than 1% of the coil-field scale. D3 must
converge. Every scalar mode solve must decrease its quadratic functional,
meet a 1e-9 equation-residual gate and retain the exact-complex divergence
identity to rounding. No global kinetic or pairwise matrix is materialised.
Storage is linear in cells and retained modes, with sparse operators and
bounded per-cell quadrature work. The focused oracle must finish within
900 seconds at four declared OpenMP threads; no MRX speed claim is made.

FortFEM's arbitrary-point basis evaluator currently allocates scratch arrays
per query. This limits tracing performance, while the new solver loads and
D3 cell loops use cached quadrature and OpenMP. A reusable allocation-free
point-evaluation workspace belongs upstream in FortFEM. It is non-blocking
for these small readback oracles and no private replacement is introduced.

## Numerical evidence

The focused Slurm run verifies all four NN1 CTest entries. Errors below
are measured diagnostics, with no rigorous remainder or pointwise bands.
Orders are the mean over two mesh doublings, not convergence theorems.

| Oracle | Finest error | Observed order |
|---|---:|---:|
| Solov'ev P2 poloidal B / independent GS B | 1.326e-4 / 1.335e-4 | 2.07 / 2.07 |
| Solov'ev P2 flux / normalized D3 | 2.095e-6 / 3.060e-6 | 3.10 / 2.67 |
| Solov'ev P3 poloidal B / independent GS B | 1.049e-6 / 1.049e-6 | 2.76 / 2.76 |
| Solov'ev P3 flux / normalized D3 | 4.609e-9 / 8.433e-8 | 4.09 / 3.00 |
| Vacuum coil-scale B | 9.859e-4 | 2.01 |
| Vacuum one-turn map | 8.774e-4 m | 2.36 |
| Vacuum normalized D3 | 1.402e-4 | 2.70 |

The N=0 closed-complex defect is at most 7.4e-14 in the reported coefficient
normalisation. The equation-level pressure tangent error is below 9e-12.
The segment-potential curl agrees with independent Biot-Savart to 1.8e-11
relative. Stress JVP/finite-difference and JVP/VJP pairing checks pass, and
the generated 3D stress contracts to the existing axisymmetric D3 leaf.

Direct periodic-orbit R coordinates are O=5.516933 m and X=6.702895 m,
with Z=0. FEEC O/X position errors are about 1.14/0.80 mm; the complete
return-map closure defects are below 2 micrometres. Both fields retain an
elliptic O-point and hyperbolic X-point. At 32 turns, both width diagnostics
are 0.175986 m within the 0.000527 m scan resolution. At 64 turns with half
the integration step, direct/FEEC widths are 0.174756/0.174463 m, differing
by 0.293 mm. Duration/step refinement changes the widths by 1.23/1.52 mm.
The direct one-turn step-doubling difference is 5.7e-10 m; for FEEC it is
0.168 mm, below the 0.877 mm field-line-map comparison error. Tangential
FEEC jumps limit smooth-stepper accuracy; the direct-field integration rate
cannot be assigned to the FEEC trajectories.

The build/test configuration registers 149 tests, including four new NN1
entries. [The Slurm recipe](../scripts/test_nonnested_scluster.sbatch) runs
the full suite at the committed revision with the pinned FortSym build,
conda gfortran and MUMPS, recording revision and logs in the private checkout.
It follows ~/kt2.sbatch and preserves each CTest thread reservation. The
final full-suite result and job identifier are reported with the delivery.
The local check was limited to the quick symbolic-leaf test; all complete
builds, mesh ladders, long traces and full CTest runs used sCluster.

## Open items and NN2

NN1 does not implement general finite-beta 3D relaxation, helicity-preserving
admissible variations, current sheets, free-boundary coupling, geometry
sensitivities or continuous certificates. The atlas's rotating ellipse,
finite-pressure island, stochastic finite-pressure and Landreman entries
are not promoted by this vacuum fixture.

NN2 must first choose and derive pressure transport/weak field-line
constraints and their tangent/adjoint semantics, then couple the free-B
unknowns to the common nonlinear D3 equation with conserved flux/helicity
and a declared error status. It needs offset-start finite-pressure island
and pressure-flattening gates, plus the Landreman closed-form comparison.
M3/M4/M6 remain required for the general certificate, boundary and harness
contracts. FortFEM should own the point-evaluation workspace and any generic
commuting seed projection; none is privately duplicated here.
