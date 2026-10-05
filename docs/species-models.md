# Species models and hybrid coupling

Status: required long-term feature envelope. Only implement a reduced model
when a concrete scientific case and its coupling/certificate obligations are
defined. The first runnable configuration is a steady 1X1V laboratory;
the ultimate parent reference remains full-orbit kinetics.

## Model selection axes

A separate [bulk-plasma model](model-family.md#bulk-plasma-and-mhd) can
aggregate thermal species into multi-fluid/MHD state and couple to selected
kinetic populations. It is not an additional per-species label. Population
ownership and field/exchange accounting must prevent double counting.

Select the species model independently for each population. Separately select
pair collision physics, discrete structure treatment, field model, geometry,
stationary/transient consumer, and certificate scope. Independence of these
choices in configuration does not imply that every combination is admissible.
Reject unsupported or inconsistent combinations before solving.

| Species model | State | Retained response | Required reduction or closure |
|---|---|---|---|
| Full kinetic (`FK`) | `f(x,p)` | Full orbits and gyroangle within the chosen parent equations | Collision/interaction closure only |
| Gyrokinetic (`GK`) | `F(R,p_parallel,mu)` or specified equivalent | Ordered gyrocenter motion, trapping, and retained FLR effects | Low-frequency/magnetization ordering, gyrocenter transform, field action |
| Drift kinetic (`DK`) | Guiding-center distribution | Parallel kinetics and ordered drifts/trapping | Selected long-wavelength guiding-center order and field coupling |
| Parallel kinetic (`PK`) | `f(x,p_parallel)` or field-line equivalent | Parallel streaming, current, and resonances supported by its equation | Perpendicular distribution/energy closure and geometry assumptions |
| Bounce/orbit averaged (`BA`) | Distribution in chosen orbit invariants | Ordered slow orbit/collision evolution | Bounce/transit averaging, orbit topology, and separatrix treatment |
| Fluid/moments | Declared moment fields | Selected density, flow, pressure, heat-flow responses | Moment set, constitutive closure, and entropy/energy assumptions |
| Gyrofluid/Landau-fluid | Declared low moments and nonlocal/kinetic closures | Selected FLR, phase-mixing and thermodynamic responses | Derived closure, ordering, field coupling and free-energy ledger |
| Gyromoments | Hermite–Laguerre coefficients of a named GK model | Progressively resolved parallel/perpendicular velocity dependence | Retained index set, basis/measure, closure and controlled moment tails |
| Adiabatic | Algebraic density/response | Specified instantaneous response | Boltzmann or other response law, normalization, zonal/mean treatment |
| Fixed background | Prescribed distribution or charge/current | No dynamical response | External reservoir definition and input uncertainty |

These models form a collection of reductions and closures, not a universal
ordered accuracy ladder. GK and DK depend on retained order; parallel closure
and bounce averaging discard different information. Full-orbit energetic
particles can coexist with cheaper thermal populations only when the combined
equations are derived consistently.

Examples of intended supported configurations are all-FK reference;
DK electrons with FK ions; DK electrons, GK D/T, and FK energetic alphas;
parallel-kinetic electrons with FK ions; and kinetic electrons with prescribed
ions. Gyrofluid/Landau-fluid and gyromoments are explicit priority branches,
with [closure and representation contracts](moment-models.md). Relativistic
FK kinematics do not automatically make a DK/GK/fluid closure relativistic;
each model states its Hamiltonian and invariant definitions.

## Three production tiers

1. Certified full-kinetic reference for homogeneous and local problems,
   annuli/full surfaces, global axisymmetric steady states, and affordable
   nonaxisymmetric steady responses.
2. Certified reduced/hybrid production for larger domains, ion-scale
   turbulence, and long-time evolution within an explicitly ordered model.
3. Certified comparison and eventually reduction-error bridges between
   those models on common observables and physical cases.

The [global full-f GENE-X gyrokinetic code](https://genex.ipp.mpg.de/about)
is a precedent for device-scale tokamak/stellarator turbulence. The
[ORB5 electromagnetic gyrokinetic code](https://arxiv.org/abs/1905.01906)
documents multiple gyrokinetic ions, drift-kinetic electrons, and alternative
electron responses. These demonstrate useful model choices, not a theorem
that arbitrary FK/GK/DK combinations are already available or certified here.

## Reduced equations and field coupling

For DK electrons require, among other ordering assumptions, small electron
orbit scale relative to field/distribution variation and low characteristic
frequencies compared with electron cyclotron motion. `k_perp rho_e << 1`
and `omega/Omega_e << 1` are useful diagnostics, not a complete reduction
theorem. Ion full-orbit equations may retain ion gyromotion, but a hybrid's
admissible fields/frequencies remain constrained by its reduced electron model.

DK/GK distributions use transformed coordinates and Jacobians. Their field
response includes the appropriate push-forward, polarization, magnetization,
and boundary terms at the chosen approximation order. Simply retaining local
`m=0` harmonics of an FK state does not implement these equations.
[Guiding-center variational principles](https://arxiv.org/abs/1602.05030)
provide a reference for deriving the corresponding Maxwell response.

Use a common field infrastructure with contributions from each species;
choose one mutually consistent combined field action. Full Maxwell, Poisson,
or a reduced low-frequency/quasineutral field system are separately named
models. Do not mix field equations from incompatible orderings and call the
result exact Maxwell. Static polarization charge and magnetization remain
relevant even when the time derivative of polarization current vanishes.

A PK equation needs an explicit perpendicular closure. A one-dimensional
parallel velocity state does not generically retain magnetic-moment trapping;
mirror forces or averaged trapping need the corresponding closure/orbit
information. Adiabatic electrons require a prescribed temperature, reference
density, gauge/normalization, and treatment of surface-constant or zonal
components. Their charge law alone does not define an electromagnetic current
response or a consistent dissipative collision partner.

## Mixed collision pairs

For representations `A,B`, a useful schematic construction is

$$
C_{ab}^{A\leftarrow B}=P_a\,Q_{ab}[L_a u_a,L_b u_b],
$$

where `L_s` lifts to a named parent particle phase space and `P_s` returns
the reduced action. A lift is a model-dependent reconstruction; information
removed by a reduction is not recovered exactly. Include transform order,
phase-space measure, particle/gyrocenter location shifts, and field dependence.
In nonuniform fields, collisions can induce reduced spatial transport;
gyrophase averaging a momentum update at the same grid point can miss it.

Derive matching contributions for both partners. Prove the appropriate
particle and paired momentum/energy identities using the parent/reduced
Hamiltonian and field coupling. A low-dimensional conservation correction
cannot repair an inconsistent physical lift, omitted polarization, or lost
asymptotic terms. It also cannot establish positivity or entropy automatically.

For fluid partners supply compatible moment exchange closures. Adiabatic and
fixed partners have no evolved distribution onto which a pair update can be
projected; prescribe their exchange/reservoir interpretation separately.
Use full Landau/BB primitives where mathematically applicable, rather than
inventing a bespoke formula for every pair or treating every partner as an
array of the same shape.

The species/pair descriptor records supported combinations, transform and
measure definitions, ordering assumptions, field compatibility, invariant
tests, and certificate theorem. A configuration with a missing pair bridge
fails explicitly. Introduce the smallest representation seam with the first
implemented hybrid, not a hierarchy of speculative classes now.

## Electron kinetics with fixed ions

This is an early useful branch. Ions can supply only a prescribed charge
density/current, or a physical distribution `F_i(x,p)` for electron–ion
collisions. Those are distinct inputs: an immobile neutralizer alone does
not specify a finite-temperature collision background.

Evolve electron kinetics with `C_ee[f_e,f_e]+sum_i C_ei[f_e,F_i]+S_e`.
Choose fields independently: prescribed `E,B`, self-consistent Poisson, or
a supported Maxwell field problem. With electron charge `-e`,

$$
\rho=-e\int f_e\,d^3p+\rho_i^{\rm fixed},\qquad
j=-e\int v_e f_e\,d^3p+j_i^{\rm fixed}.
$$

For reduced electrons replace these particle moments by the chosen model's
consistent charge/current coupling. Prescribed ion flows must obey the
stationary or time-dependent charge-continuity and boundary constraints of
the field problem. Any imposed reservoir source is part of the model.

Record electron–ion momentum and energy exchange as transfer to the external
ion reservoir. Do not force it to vanish with closed-system projection.
Electron self-collisions retain their own invariants. When ion response is
later enabled, require the paired exchange and field/source balances of the
combined model. A finite-mass prescribed ion distribution is not silently
replaced by the infinite-mass Lorentz limit.

Preserve fixed-ion configurations and overlay their certified observables
with the corresponding responsive-ion model on common inputs. These cases
support electron transport/current and energetic-electron studies in real
fields before the cost of evolving every species is justified. Runaway
interpretation still needs the declared hard-collision/radiation/atomic terms
for the regime, rather than BB alone being labeled complete runaway physics.

## Certified reductions

Each model first receives an enclosure of its own continuous solution.
For the same observable, `I_hybrid-I_full` then encloses the difference of
the certified solutions. Branch correspondence, consistent source/initial
data, particle-versus-gyrocenter pullbacks, and observable normalization
must be stated. This difference is not automatically a proved asymptotic
remainder or a causal isolation of one physical mechanism.

A reduction-error theorem needs verified ordering and regularity, a controlled
transform/lift defect, initial/boundary layers, collision and field reductions,
and stability over the requested domain/time. It may then combine a model
remainder with the numerical certificate in a common observable/norm. Keep
input and artificial-domain uncertainty separate. Long-time error can grow
even when an instantaneous ordering parameter is small.

Eventually choose the cheapest admissible species/pair model whose proved
numerical plus reduction error meets the requested tolerance. Without the
remainder theorem, report a model comparison rather than automatic fidelity
selection as a guarantee. Switching a model during a trajectory requires
conservative state transfer and a certified switching/projection defect.

## Memory and scale savings

An angular `L=12` FK expansion has `169` coefficients; a gyroangle-independent
pitch expansion has `13` in a suitable local frame. This illustrates an
electron angular-state saving only, conditional on comparable radial/order
choices and the transformed model. It is not a GK implementation by itself
or a full-code factor-of-thirteen speed claim.

Removed cyclotron dynamics and relaxed ordered spatial/temporal requirements
can matter more than state count. Certification of the reduced model bounds
its numerical error; certification against the parent requires the separate
reduction theorem. A reduced electron model can leave ion waves, fast
parallel streaming, or electromagnetic stiffness requiring temporal treatment.

## Populations beyond the initial kinetic spine

[Coupling roles](coupling.md) are independent of the phase-space model:
a FK test population may be passive while DK is active. Neutral/charge-state
populations and their closures are specified in [edge physics](edge-physics.md).
MHD aggregates retain explicit population ownership; they are not a per-species
transport flag. [Space scope](space-plasmas.md) reserves CGL/nonlocal fluid
closures, arbitrary distribution shapes and minor populations when needed.

Mass is positive and charge may be signed or zero. Pair-plasma or nuclear
product cases need compatible collision/reaction definitions; general species
parameters do not imply all interactions are implemented. Relativistic fluid
reductions and alternative velocity bases are later derived capabilities.

Chris&AI
