# Analytical Benchmark Atlas

Status: catalogue and scoped executable references, 2026-10-05. Streaming,
periodic Poisson, source-balanced VPFP, bounded OU and formal polynomial Krook
cases have retained evaluators and independent tests. The [general damped-transport
case](general-transport.md) adds independent continuous state/tangent and
observable oracles; its known solutions never enter the production checker.
The [NN1 non-nested entry](equilibrium-nonnested.md) adds an independent
FEEC Solov'ev reduction and a helical-filament vacuum-island comparison
against direct Biot-Savart tracing. These are candidate/reference-comparison
results without continuous enclosures; the finite-pressure island, stochastic
and Landreman cases remain planned.
Other rows remain planned.
The audience is physicists
and implementers selecting the next independently verifiable capability.
The atlas defines its mathematics before the solver defines its tests.

## Acceptance contract

Every applicable mathematical reduction with an analytic or semi-analytic
answer becomes a permanent reference case. Build each case with the
capability that needs it, retain its executable configuration, and plot the
general model's limiting case against that reference with both interval bands.
The atlas and cross-code comparisons are equally important; agreement with
another code alone does not establish a continuous-solution enclosure.

The [catalogue](../benchmarks/atlas.json) is a planned case inventory. A row's
reference class describes the intended mathematical reference, while status
and proof level describe evidence actually produced. Proposed formulas,
literature selection, and metadata are not successful tests. Add missing
reductions as scientific cases reveal them; the seed list is not exhaustive.

| Class | Reference | Required interpretation |
|---|---|---|
| E0 | Exact closed form | Formula and data solve the named problem; evaluate rigorously |
| E1 | Exact representation needing roots, quadrature, or special functions | Include evaluation, contour, infinite-series, and truncation bounds |
| A0 | Asymptotic reference | State ordering and remainder status; a leading term alone is not an enclosure |
| M | Manufactured exact solution | Generate complete consistent forcing and boundary/initial data |
| I | Exact identity or invariant | Prove assumptions and check applicable discrete balances separately |
| X | Cross-code numerical case | Pin the external implementation and model intersection; not analytic truth |

Reference E0 is a class label, distinct from roadmap milestone E0. Keep
reference class, applicability, implementation status, run result, numerical
certificate level, and formal theorem coverage as separate fields. A mature
capability supplies its applicable E0/E1/A0/M/I/X matrix; unused classes have
a reason rather than a fictitious pass. An A0 without a proved remainder can
support a diagnostic scaling check but cannot pass a guaranteed error gate.

## Runnable reference evaluation

Closed forms use exact parameter definitions and generated rigorous
interval/ball kernels. Special functions, dispersion roots, integration,
series, and Fourier–Laplace inversion have certified error budgets. Prove
root existence/isolation for the selected branch, not merely a small residual.
Document frequency, Fourier, charge, harmonic, and collision-rate conventions.
Enclose the reference over entire plotted cells/time slabs when needed.

Use Fortran drivers and fortsym emission for all owned numerical leaves,
including analytic fields, forcing, moments, derivatives, and reference
enclosures. fortnum supplies supported verified operations. Missing generic
special-function, quadrature, root, or code-generation support is repaired
upstream when the next case needs it. No private numerical kernel or Python
postprocessor is introduced for the atlas. Lean proves the applicable
identity/reference-to-error theorem and checks certificates incrementally.

Reference and solver definitions must have independent mathematical evidence.
Sharing a symbolic DAG prevents transcription drift but does not make its
own differentiated forcing an independent oracle. Check manufactured forcing
by a separately derived identity, invariant, closed-form restriction, or
reviewed analytic derivation. Handwritten reference formulas are allowed in
the design; executable numerical formulas must be generated.

## Simple geometry before machines

Keep the geometry ladder slab, cylinder, circular torus, shaped analytic
tokamak, simple classical helical field, then real-device equilibrium. A
cylindrical domain containing the axis needs regularity and an axis treatment;
the initial toroidal domain with R>0 does not automatically support it.

For a circular toroidal domain, the field
`B_phi=B0 R0/R`, with `B_R=B_Z=0`, is an exact divergence-free, curl-free
vacuum field on R>0. It tests geometry and full-orbit frame terms but has no
poloidal rotational transform. Add poloidal fields through a stated flux or
vector potential, verifying their currents rather than calling every
hand-composed field a self-consistent kinetic equilibrium.

A simple straight periodic helical vacuum field can use cylindrical scalar
potential

$$
\Psi=B_0 z+b I_m(k r)\cos(m\theta-kz),\qquad B=\nabla\Psi.
$$

Here m is a nonnegative integer, k>0, and I_m is the modified Bessel
function. The Bessel equation gives Laplacian Psi=0, hence div B=0 and
curl B=0 on the stated regular domain. The background potential is not
periodic but its field is. This is an exact straight-helical reference, not
a toroidal stellarator equilibrium. Bend or embed helical fields using a
specified divergence-free potential construction; include induced currents
and transformation bounds. Retain field-period and broken-symmetry variants.

Solov'ev and other exact Grad–Shafranov families add circular, shaped,
limiter, and single/double-null geometry tests. They solve a reduced magnetic
equilibrium problem, not automatically stationary Vlasov–Maxwell–collisions.
[Guazzotto and Freidberg](https://doi.org/10.1017/S002237782100009X)
provide exact shaped equilibria and flux derivatives for this layer.

[Landreman's analytic families (arXiv:2609.26742)](https://arxiv.org/abs/2609.26742)
are exact 3D ideal-MHD equilibria with nested toroidal surfaces, not
perturbative in aspect ratio or symmetry breaking:

- Family 1 (uniform ι=2) is a constant linear stretch of a Solov'ev-type
  axisymmetric field, with ψ and p = p_a − 2ψ in closed form.
- Family 2 (sheared ι) is written in complex variables and avoids its branch
  cut by the choice of domain.

Both give B, J, p, ι and volume averages (|B|², β) in closed form. They test
3D force balance, div B and the surface geometry for the nested and
physical-space 3D solvers, and give the same oracle to VMEC++, DESC and
GVEC. They do not test islands or stochastic regions.

## Three-dimensional equilibrium programme (planned)

[Issan et al. (arXiv:2610.07304)](https://arxiv.org/abs/2610.07304) add a
symmetry-breaking parameter to both Landreman families, giving exact nested
benchmarks without stellarator symmetry (abstract checked only; equations are
re-derived before activation). Non-nested cases are cross-code (X) and named
after the physics problem:

| Case | Class | Model lane | Producers |
|---|---|---|---|
| landreman_3d_nested | E0 | ideal nested | KIN6D, MRX, Zhang--Zhu--Hegna, VMEC++, DESC, GVEC |
| issan_3d_nonsymmetric | E0 | ideal nested | KIN6D, MRX, VMEC++, DESC, GVEC |
| rotating_ellipse_3d | X | ideal non-nested/weak or vacuum | KIN6D, MRX, SPECTRE |
| finite_pressure_island | X | ideal non-nested/weak; resistive relaxed | KIN6D, Zhang--Zhu--Hegna, MRX/HINT where equations match |
| stochastic_finite_pressure | X | resistive relaxed | KIN6D, Zhang--Zhu--Hegna, HINT/PIES |
| mrxmhd_island_chaos | X | MRxMHD | KIN6D-MRxMHD, SPEC/SPECTRE |
| free_boundary_3d_nonnested | X | MRxMHD; ideal non-nested/weak | SPECTRE, KIN6D |

Compare codes only within one model lane (ideal nested; ideal non-nested/weak;
resistive relaxed; MRxMHD). A result from another lane is context, not a
reference.

Observables by class:

- E0: errors in B, J and p; div B; force residual; ι; V; ⟨B²⟩; β, all on an
  independent Cartesian mesh; convergence order and cost to reach accuracy.
- Islands: O- and X-point positions, island width, helical phase, ι on
  surviving surfaces, p profile across the island, energies and B·∇p. Poincaré
  plots are diagnostic only.
- Stochastic: surviving tori, island widths, field-line excursion and
  diffusion, connection lengths, classified volume fractions, p variation
  along lines, and pressure flattening.

Every island or stochastic case needs a topology-analysis specification before
activation: integration length, initial points and thresholds. The first
Zhang--Zhu--Hegna case reproduces one of its published cases verbatim (data
from the paper, pinned before activation). All entries are planned; the atlas
nulls are resolved at activation.

For large aspect ratio epsilon=a/R0 tending to zero, specify what remains
fixed: q profile, collisionality normalization, gyroradius ratios, beta,
profiles, and source/boundary data. Compare torus-to-cylinder geometry,
trapped fraction, bounce/transit motion, and neoclassical coefficients under
their own orderings. Deep-trapped bounce formulas depend on speed/pitch
definitions and q conventions; no universal coefficient is stored as truth.

## Reference families

The catalogue covers the following scientific layers. Each named member
needs its own instantiated data and validity range before activation.

| Family | Permanent targets | Scope of reference |
|---|---|---|
| Characteristics | Free streaming; constant E; uniform B; crossed E/B; harmonic potential; relativistic constant fields | Exact flow pullback with model-specific kinematics |
| Fields | Coulomb, periodic/slab Poisson, Bessel cylinder, Fourier torus, screened Poisson, current loop, wire, cylinder, solenoid, vector potential | Matching sources, gauges, domains, and singularity treatment |
| Unmagnetized linear Vlasov | Langmuir, Landau damping/growth, Case–van Kampen, two-stream, bump-on-tail, ion acoustic | Dispersion and complete initial-value response are distinct tests |
| Nonlinear phase mixing | Plasma echo | Controlled perturbative response or separately bounded nonlinear solution |
| Magnetized waves | Ion/electron Bernstein, cyclotron harmonics, upper/lower hybrid, oblique hot-plasma Harris relation, O/X/R/L and Stix cold branches | Exact linear/cold equations or ordered limits, not arbitrary nonlinear trajectories |
| Nonlinear steady Vlasov | BGK modes, Harris and force-free sheets, Bennett/cylindrical pinches, invariant distributions | Full distributions and field consistency, admissibility, boundaries |
| Simple collisions | Krook, conservative homogeneous BGK, Lorentz harmonics, Ornstein–Uhlenbeck, diffusion, slowing down | Correct fixed coefficients, reservoirs, and reduced collision assumptions |
| Landau | Maxwellian pair equilibria, invariants/entropy, Maxwellian potentials, exchange/isotropization rates, linear derivative, Lorentz/test limits | Same coefficient prescription and finite/whole momentum domain |
| Relativistic collisions | Maxwell–Jüttner equilibrium, relativistic invariants, BB potentials, c-to-infinity and low-temperature limits | Nonrelativistic correspondence at fixed physical momentum/velocity conventions |
| Classical transport | Spitzer, Spitzer–Härm, Braginskii conductivity/viscosity, Lorentz gas, weak/strong magnetization, relativistic transport | Ordered near-equilibrium transport with coefficient conventions |
| Guiding-center limits | E-cross-B, grad-B/curvature drift, mirror, trapped/passing, bounce/transit and action-angle spectra | Reduced orbit accuracy includes ordering/remainder, not exact FK equality |
| Tokamak neoclassical | Pfirsch–Schlüter, plateau, banana; heat/particle flux, bootstrap, flows, Ware pinch, trapped fraction | Geometry, collisions, electric field, and scale assumptions matched |
| Stellarator neoclassical | 1/nu, sqrt(nu), nu, superbanana plateau; Shaing–Callen, HGM, omnigenous/QA/QI and near-omnigenous formulas | Case-specific E_r, magnetic drift, resonance and boundary-layer orderings |
| Gyrokinetic | Slab/toroidal ITG, drift/kinetic-Alfvén waves, GAM, zonal residual, adiabatic/DK limit, Bessel FLR | Ordered reduced equation, initial data and time limit explicitly selected |
| Cylindrical transport | Radial classical transport, z/theta/screw pinch | Fluid, reduced kinetic, and exact VM instances identified separately |
| Manufactured | Arbitrary supported 3D fields, species mixtures, collision and field blocks, boundaries | Include every source and compatibility condition |
| Bulk MHD/hybrids | Ideal force balance, Alfvén/sound/magnetosonic waves, resistive induction, Hall limits and manufactured kinetic feedback/layers | Bulk ownership, closure, spectral or interface bounds |

Analytic eigenfrequencies alone do not enclose the initial-value solution.
Phase mixing may involve continuous spectra, transient terms, residues and
branch contributions. The [linear time-domain Landau reference](https://arxiv.org/abs/2503.22501)
and [Bernstein time-domain construction](https://arxiv.org/abs/2608.12994)
are promising E1 sources; their numerical inversions still need rigorous
evaluation here. A plasma echo is generally a controlled perturbative test,
not an exact nonlinear closed form at arbitrary amplitude.

For invariants I_j of the stated collisionless flow, f=F(I_1,...) solves its
transport equation. H and canonical toroidal momentum are useful in static
axisymmetric fields. Verify the invariants, gauge conventions, regularity,
normalizability, and boundary traces. Poisson/Maxwell self-consistency and
collision equilibrium are additional obligations. A kinetic BGK equilibrium
must also be distinguished from the BGK relaxation collision model.

Harris-sheet families and the [force-free kinetic construction](https://arxiv.org/abs/0812.1240)
are nonlinear field/distribution benchmarks. The
[cylindrical confinement constructions](https://arxiv.org/abs/2608.06780)
offer semi-explicit pinch references, whose ODEs and parameter restrictions
must be enclosed before use. Existence of a steady benchmark does not imply
its stability or a unique nearby solution without constrained stability work.

The [Rosenbluth–MacDonald–Judd formulation](https://doi.org/10.1103/PhysRev.107.1)
and [Braams–Karney potential formulas](https://arxiv.org/abs/plasm-ph/9502001)
anchor nonrelativistic and relativistic collision references. Homogeneous
BGK has a fixed Maxwellian target only when its conserved moments and
constant relaxation frequency make that target fixed. Lorentz rates must
state whether the angular Laplacian includes a factor of one half.
Slowing-down profiles are exact for a stated reduced drag/source problem,
not general full Landau with diffusion and all species response.

Use the [near-omnigenous bootstrap implementation rabe](https://github.com/itpplasma/rabe)
only within its documented Lorentz, monoenergetic, and near-omnigenous
assumptions. The local source inspected at revision
3282d659cb44108a61d537c11a5d1dd8a5526c92 is cross-code evidence, not a
certified oracle or a universal superbanana solver. Pin its matched case;
independently evaluate the underlying formulas for an A0 reference.

The [BSL6D slab-ITG example](https://bsl6d.pages.mpcdf.de/bsl6d/examples/slab_itg.html)
provides full-kinetic/reduced dispersion comparisons under its background and
Boussinesq assumptions. Cyclone Base Case is a standard numerical X case,
not an analytic E0 case. No benchmark result is inferred from its popularity.

## Limiting-case lattice

Store each limit as a directed edge between named problems, with parameter
path, unchanged data, scaling of the observable, admissible norm, and
remainder or diagnostic status. Compare complete interval bands at finite
parameters. Limits need not commute, and singular limits can change norms,
domains, constraints, or solution branches.

| Limit | Intended comparison | Required qualification |
|---|---|---|
| E or B amplitude to zero | Magnetic-only or unmagnetized FK | DK/GK coordinates become singular at B=0 |
| nu to zero/infinity | Collisionless or local-equilibrium/fluid regime | Time/source scaling and boundary/initial layers |
| T to zero | Cold response | Distributional limit; generally not strong convergence in a thermal L2 norm |
| Ion mass to infinity | Immobile response or Lorentz scattering | Fixed ion response and collision-reservoir limit are separate claims |
| c to infinity | BB to Landau, consistent kinematics | Parameter/energy normalization held consistently |
| k_perp rho or omega/Omega small | DK/GK ordered response | Neither single ratio alone specifies every reduction |
| lambda_D/L to zero | Quasineutral limit | Well-prepared data and possible unresolved layers/fast fields |
| Perturbation amplitude to zero | Linear response | Fixed reference, forcing, norm, and field/collision derivative |
| Nonzero toroidal amplitudes to zero | Axisymmetry | Integer mode labels do not tend continuously to zero; N=0 is an exact configuration restriction |
| beta to zero | Vacuum-dominated response | External currents, drive, and self-field scaling remain specified |
| a/R0 to zero | Large-aspect-ratio cylindrical/tokamak asymptotes | Profile, pitch, q, and collisionality path |

These are metamorphic behavioral tests only after the target relationship is
proved or analytically justified. A resolution sweep is not a remainder
proof. Do not use zero toroidal modes, cold fluid response, a simple helical
field, or an invariant identity as proof of the final nonlinear machine case.

## Catalogue and implementation sequence

Version the JSON catalogue now; implement its smallest Fortran reader only
when R0 needs executable metadata. Common defaults and family definitions
are explicitly inherited by each case, with case fields taking precedence.
Each activated instance freezes equations, species/pair/field models,
geometry, parameters/validity, initial/boundary data, formula and source
version/DOI, observables, generated reference implementation, certificate
method and proof level, expected order, resources, and plot recipe.
Null fields are unresolved obligations and forbid promotion to runnable.

Implement R0 plus characteristic/Poisson cases first. Add simple harmonic
collisions with B0–B3, circular/helical field checks with P0/P1/G0, exact
nonlinear equilibria with field coupling, Coulomb references with C0–C6,
neoclassical asymptotes with device geometry, and reduced/transient responses
with their consumers. Expand the catalogue before implementing each new
capability; do not build every special-function evaluator at once.

Keep exact mathematical expectations independent of discretization order.
Record the numerical scheme's proposed order only when selected and prove
or diagnose it in the relevant norm. Production promotion requires the
applicable reference matrix, the continuous certificate, and the physical
device gates. Plot generation remains on demand; no output figures are
committed. [Demonstrations](demonstrations.md) defines the band contract.

Chris&AI
