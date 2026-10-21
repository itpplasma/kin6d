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
cases remain planned. The nested Landreman/Issan field evaluators and boundary
exports are implemented as [scoped exact references](nested-oracles.md); solver
comparisons remain pending.
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

Both give closed-form B, J, p and nested surface maps. Family 1 also gives
closed-form ι and volume averages; family 2 uses an ODE rotation number and
quadratures for these integrated observables. They test
3D force balance, div B and the surface geometry for the nested and
physical-space 3D solvers, and give the same oracle to VMEC++, DESC and
GVEC. They do not test islands or stochastic regions.

## Three-dimensional equilibrium programme (planned)

### Nested-reference input admission, 2026-10-21

[V] The integer-transform Landreman and nonsymmetric Issan cases now have
matched pressure, iota and toroidal-flux inputs, emitted from the retained
FortSym definitions. An independent surface integral checks toroidal flux to
4.07e-11 relative. The affected native checks pass on sCluster with
gfortran 12.2. The sheared families still require pressure/flux conversion.

[V] The Fourier boundary ladder below uses eps=0.5, delta=1/64, tau=0 or 0.4,
64 angular samples and eight threads. Errors are checked on a staggered grid.
The minor scale is 0.125 m, so the requested 1% radius is 0.00125 m. Each
listed resolution passes. Full raw timing and memory
records are in [the input ledger](../benchmarks/nested/boundary-inputs.tsv).

| Case | Modes m,n | Sampled boundary error (m) | Export time (s) | Peak RSS (KiB) |
|---|---:|---:|---:|---:|
| Landreman integer | 8 | 4.165e-6 | 0.07 | 15,164 |
| Landreman integer | 12 | 2.141e-8 | 0.08 | 15,812 |
| Landreman integer | 16 | 1.117e-10 | 0.11 | 15,616 |
| Issan integer | 8 | 1.271e-6 | 0.07 | 15,232 |
| Issan integer | 12 | 3.318e-9 | 0.10 | 15,900 |
| Issan integer | 16 | 8.776e-12 | 0.24 | 15,508 |

[V] Temporary renders of the input ladder:
[Landreman](https://box.sloppy.at/f0cf6.png) and
[Issan](https://box.sloppy.at/7c605.png). The retained
[plot recipe](../benchmarks/nested/plot_inputs.py) regenerates both figures.

[V] These are input-preparation measurements, including point intervals and
I/O, with one timing per point. They establish neither solved equilibria nor
a speed ranking. Boundary errors are sampled indicators, not continuous
bounds. The [comparison adapter](../benchmarks/nested/README.md) retains
external states and uses a Fortran verifier for errors in B, pressure and
surface labels and for the strong force defect. State bounds and estimator
effectivity remain unavailable. A solved KIN6D epsilon>0 initializer and D8
case admission are still required before a production time-to-accuracy claim.

[V] The first reference runs exposed a coordinate defect in the original
integer-family boundary adapter: its material poloidal angle winds twice with
toroidal angle, while iota=2 uses the conventional untwisted angle. The adapter
now uses material angle 2 phi-theta. An axisymmetric-limit test fails with the
old mapping and passes with the correction. The physical field and pressure
definitions are unchanged. GVEC additionally needs both theta and zeta
reversed for these clockwise boundaries, consistently transforming iota and
flux. VMEC's internal theta reversal changes the sign of its stored iota;
compare it after the declared coordinate conversion, not by taking its magnitude.

### Catalogue

[Issan et al. (arXiv:2610.07304)](https://arxiv.org/abs/2610.07304) add a
symmetry-breaking parameter to both Landreman families, giving exact nested
benchmarks without stellarator symmetry. [V] The equations are retained in
FortSym, with exact factored proofs and independent numerical checks described
in [the oracle contract](nested-oracles.md). Both families and both symmetry
choices have executable inputs. Integrated/rotation/boundary bounds and
matched external solver comparisons remain unavailable. Non-nested cases are cross-code (X) and named
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

Model match of the resistive comparators (checked against the papers). The
Zhang--Zhu--Hegna code (NTEC) damps p by a parallel-only hyperbolic equation
with artificial resistivity and flow; it has no chi_par/chi_perp closure.
HINT computes p by field-line tracing with a connection-length cutoff. The
closest published match to NN2's model (anisotropic p transport, momentum,
resistive induction) is
[Hudson and Nakajima 2010](https://doi.org/10.1063/1.3431090), which proposes
it for HINT and has no published results. The transport operator alone is
tested against Hudson and Breslau
([PRL 100, 095001](https://doi.org/10.1103/PhysRevLett.100.095001)) in
`test_nonnested_hudson_breslau`: exact integrable limit gated; the chaotic
field (k = 2.1e-3) is an unconverged ladder.

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
from the paper, pinned before activation). Entries without an implementation note remain planned; the atlas
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

## Production estimator scaling

Verified measurements, 2026-10-15: Release `-O3`, GNU Fortran 16, OpenMP MUMPS,
Ryzen 9 5950X (16 physical cores, 32 logical CPUs), shared host. The
[raw ledger](../benchmarks/error-control/production-consumers.tsv) retains both
samples at every size/thread point, the baseline and optimized revisions,
and a three-repeat paired NN2 comparison. These are modest toroidal
production-path cases, not a large-device scaling limit. No production
true-error oracle is available, so production efficiency indices are unavailable.
The existing independent efficiency/containment CTests pass; claim levels stay
practical/approximate.

The table uses medians at four threads. `E/S` divides estimator time by solve
time. The exponent is `log(E_large/E_small)/log(N_large/N_small)`; it is a
measurement, not a complexity proof. The final free-boundary mesh has 5,481 or
5,489 DOFs at the largest size because the refit rebuilds its triangulation.

| Consumer | Three DOF counts | E/S, small / medium / large | Large E (s) | Measured DOF exponent |
|---|---|---|---:|---:|
| NN1, P2 mapped Solov'ev, 16/32/64 rings | 3,105 / 12,193 / 47,777 | 1.41 / 1.39 / 1.07 | 0.573 | 0.90 |
| NN2, NN1 field, heating 1, chi ratio 1,000, 8/16/24 rings | 865 / 3,105 / 7,201 | 1.24 / 1.39 / 1.24 | 0.535 | 1.01 |
| Ampere n=2, m=3 ring sheet, setup + first solve, 4/8/12 rings | 592 / 2,064 / 4,624 | 0.028 / 0.045 / 0.119 | 0.070 | 1.30 |
| Same Ampere estimator / reused sheet solve | same | 3.88 / 2.98 / 5.38 | 0.070 | 1.30 |
| D3, solved GS field, P2 displacement, 48/96/144 boundary nodes | 2,378 / 9,386 / 21,010 | 3.55 / 3.64 / 3.40 | 0.342 | 1.07 |
| Coupled/refitted gs_free, 48/72/96 boundary nodes | 1,377 / 3,177 / about 5,485 | 0.52 / 0.70 / 0.85 | 2.217 | 0.99 |

| Largest case: estimator speedup relative to one thread | 1 | 2 | 4 | 8 | 16 |
|---|---:|---:|---:|---:|---:|
| NN1 | 1 | 1.81 | 3.48 | 5.52 | 9.00 |
| NN2 | 1 | 1.85 | 2.51 | 4.40 | 6.16 |
| Ampere sheet | 1 | 1.66 | 1.83 | 1.06 | 1.77 |
| D3 | 1 | 2.31 | 3.10 | 6.13 | 8.10 |
| gs_free | 1 | 1.73 | 2.68 | 5.19 | 7.71 |

NN1/NN2 exclude mesh/FEEC construction. NN2 exercises
`nt_solve(..., estimate_error=.true.)` and checks its bound against the timed
`nt_axisymmetric_error` call. Ampere includes BEM assembly/factorization in the
fresh denominator and separately times a reused `ampere_n_solve`, as in
`gs_response_3d`. D3 estimator time is the paired difference of calls with and
without `error_estimate`; its denominator is assembly/sampling plus the Riesz
solve. `gs_free` uses its existing `estimator_s` and `wall_s`: the latter
includes mesh setup, Newton, refits and readback. Setup is not silently
amortized into the sheet-solve comparison.

Contention is material: at the largest size, optimized max/min repeat ratios
reach 1.06 (NN1), 1.15 (NN2), 1.44 (Ampere), 1.58 (D3) and 1.25 (gs_free).
Consult the individual samples before comparing thread counts. Baseline gs_free at 16 threads ranged from 0.86 to 7.47 s at
the largest size. The table does not establish hardware-independent speedups.

`perf record` on NN2 identified triangle Lagrange evaluation (54.7% of sampled
cycles), with allocation and a dense Bernstein-to-nodal transform per query.
The FortFEM campaign branch removes that allocation/transform, reuses retained
generated P1/P2 jets and generated cardinal-product leaves, and preserves node
ordering. In the paired four-thread, 3,105-DOF check, median estimator cost
fell from 0.554 to 0.227 s (2.44x); bounds changed by less than 4e-12 relative.
KIN6D additionally parallelizes meridian/D3 geometry hulls, GS profile
conversion and D3 stress-load conversion. Across the full driver sweeps,
optimized bounds differ from baseline by less than 9e-10 relative. For
`gs_free` the 48/72-node bounds agree within 8e-11; the 96-node bound spans
0.051% across the two final triangulations, already present on the baseline.
No identical-mesh claim is made for that adaptive trajectory.

The small-fraction cost target remains unmet for NN1, NN2, D3 and reused sheet
solves. FortFEM `equilibrated_hcurl_meridian_n.f90:239` rejects patch recovery
for sheet/free-trace data, and `solve_global` at line 387 ends with
`sparse_solve_once` at line 427. FortSparse's default solver selects serial
SuperLU (`fortsparse_solver.F90:49,592`): the sheet path retains a global serial
factorization and has no O(N) claim. The math-owning workers must close that
path. `ideal_map%factor_threads=1` is absent from these consumers' call paths
and is unchanged; `gs_direct` uses the OpenMP maximum thread count instead.

Reproduce with the opt-in `benchmark_estimator_consumers` target and
`benchmarks/run_estimator_consumers.sh BUILD OUT`; the coupled application
recipe is `benchmarks/run_estimator_gs_free.sh BUILD OUT`. Both retain two
measurements at 1/2/4/8/16 threads and cap individual processes at 60 s. The
optimized build uses `FETCHCONTENT_SOURCE_DIR_FORTFEM` pointing to FortFEM
`eb/scale` at `a0f6fd49ba11b22dc019c518a94937dd54cd75a5`;
KIN6D's dependency pin is unchanged. Source revisions are recorded in the
campaign handoff. No plots or large generated artifacts are tracked.
