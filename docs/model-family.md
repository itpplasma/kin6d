# One mathematical engine, selectable plasma models

Status: canonical long-range model map. The implemented code is still only the
first 1X1V slice.

## Principle

KIN6D changes state space, ordering, closure, representation or coupling while
retaining common geometry, fields, residuals, diagnostics, derivative contracts,
verification and certification machinery.

The main rule is:

\[
\boxed{\text{do not create a new numerical universe when a model reduction is enough}.}
\]

Names such as NEO-2-like, GENE-X-like, MARS-like or SOLPS-like identify
comparison points in model space. They never define executable architectures.

## Model axes

| Axis | Intended choices |
|---|---|
| Per-species physics | fixed, adiabatic, fluid, gyrofluid, gyromoment, PK, bounce/orbit averaged, DK, GK, FK |
| Bulk plasma | none, ideal/resistive/Hall/XMHD, multi-fluid, kinetic-MHD |
| Distribution form | full-f, delta-f, control variate |
| Carrier | Eulerian, semi-Lagrangian, particles, blobs, hybrid/remapped |
| Particle-field coupling | PIC, PIF, particle-FE, general particle-in-basis |
| Sampling | deterministic, QMC/randomized QMC, Monte Carlo |
| Collisions | none through Landau/Beliaev--Budker, deterministic or justified stochastic realization |
| Fields | prescribed, Poisson, quasineutral EM, Darwin, full Maxwell |
| Coupling | test/passive, fixed background, linear response, fully active |
| Geometry | slab through arbitrary mapped 3D/open/wall-bounded domain |
| Problem type | equilibrium, steady transport, forced response, eigenproblem, transient |
| Evidence | analytic, structural, cross-code, interval certificate, statistical confidence |

The Cartesian product is not automatically valid. Each supported combination
declares the equations, measures, boundary/source assumptions and certificate
scope.

## Full kinetics as parent, not religion

All-species FK is the high-fidelity charged-particle reference. It is not the
default cost paid for every problem.

A practical production calculation may use, for example,

\[
e:{\rm DK},\qquad D,T:{\rm GK},\qquad \alpha:{\rm FK},
\]

or

\[
\text{thermal bulk}:{\rm MHD},\qquad
\text{energetic particles}:{\rm DK/FK}.
\]

The correct question is not "which model is best?" but "which model gives the
requested observable within the requested total error at minimum cost?"

Before rigorous reduction bounds exist, this remains a comparison programme,
not an automatic guarantee.

## Bulk MHD and multi-fluid models

One-fluid MHD combines populations into bulk variables. It is not the same as
assigning each species a fluid closure. The model records exactly which
populations contribute to bulk density, inertia, pressure, current and energy.

The intended bulk branch includes:

- static ideal equilibrium;
- ideal linear response and eigenmodes;
- resistive response;
- Hall/XMHD;
- multi-fluid;
- kinetic-MHD with active kinetic species;
- spatially hybrid MHD/kinetic layers;
- nonlinear extended-MHD applications when a real target demands them.

Ideal/resistive/MHD field equations are low-frequency/quasineutral reductions
and must not be presented as full Maxwell merely because they use the same
field infrastructure.

Kinetic feedback uses one derived current- or pressure-coupling formulation.
A kinetic population must not be counted both inside the MHD pressure and
again as an explicit kinetic stress.

## Equilibrium hierarchy

External equilibria are first-class inputs before every solver is reproduced
internally. Adapters may support axisymmetric G-EQDSK/EFIT/CLISTE/CHEASE and
3D VMEC/VMEC++/DESC/GVEC/SPEC/SPECTRE candidates, with explicit units, signs,
topology and provenance. These are interchange representations or named-model
candidates, not the ontology of equilibrium.

Internal progression is:

\[
{\rm imported\ candidate}
\rightarrow
{\rm axisymmetric\ Grad\!-\!Shafranov}
\rightarrow
{\rm free\ boundary}
\rightarrow
{\rm nested\ 3D\ specialization\ and\ general\ physical\!\!\!-space\ 3D}
\rightarrow
{\rm resistive/XMHD/relaxed\ stationary\ states}
\rightarrow
{\rm flowing/anisotropic/kinetic\ equilibrium}.
\]

Nested flux surfaces are exploited where they exist but are not required by the
general physical-space model. HINT/PIES-like nonnested full-pressure methods and
SPEC/SPECTRE relaxed-region models are comparators with different physical
constraints, not alternate file readers for one identical problem.

An imported equilibrium is not certified merely because the producer reports a
small force residual. KIN6D may certify an external candidate against its stated
continuous equilibrium model and, once the native model exists, re-solve it.
Distinct stationary solutions remain distinct equilibria. TIAGO, not this model
layer, owns reconstruction/posterior weights over several candidates.

A kinetic equilibrium is broader than a Maxwellian or ideal-MHD state. The
important target is

\[
{\rm transport}+{\rm collisions}+{\rm sources}+{\rm losses}+{\rm fields}=0,
\]

including strongly non-Maxwellian steady states.

## Gyrofluid, Landau-fluid and gyromoments

Gyrofluid is a first-class bridge, not a side model. Where possible, low-order
gyrofluid closures are truncations of a Hermite--Laguerre or equivalent
gyrokinetic moment hierarchy.

The progression is

\[
{\rm gyrofluid}
\rightarrow
{\rm Landau\ fluid}
\rightarrow
{\rm gyromoments}(N)
\rightarrow
{\rm GK}.
\]

The closure, ordering and free-energy structure matter as much as the number of
moments. Increasing \(N\) without controlling recurrence, streaming and
truncated tails is not a convergence theorem.

Scott/GEM and Kendl TIFF/T3FF are important structural/reference landmarks;
their exact assumptions must be matched before comparison.

## DK, GK and FK

DK/GK distributions use their own transformed phase-space measures, locations,
Hamiltonians, polarization/magnetization and field push-forwards. They are not
implemented by deleting dimensions from an FK array.

A mixed FK/GK/DK system must derive the common field action and pairwise
collision/exchange rules. Missing reduced physics is reported as model error,
not hidden in a numerical tolerance.

The same DK/GK/FK model may use continuum, delta-f particle or other
representations as described in [representations.md](representations.md).

## Coupling roles

Representation and coupling are independent. A population may:

- evolve;
- source fields;
- feel induced response.

This separates:

- passive diagnostics/test particles;
- prescribed reservoirs;
- active species;
- fixed bulk + induced plasma response;
- fully evolving coupled systems.

Likewise fields are decomposed into background, external drive and induced
response where useful.

## Spatial hybrid physics

A long-range KIN6D capability is model fidelity that varies by region, for
example

\[
{\rm MHD\ bulk}
\rightarrow
{\rm GK\ turbulent\ region}
\rightarrow
{\rm FK\ resonant/sheath\ region}.
\]

This requires derived overlap/interface equations, conservative transfer and a
bounded model/interface defect. A spatial mask between unrelated residuals is
not a valid hybrid.

Adaptive model fidelity is separate from numerical \(h/p/L/n\) refinement.

## Field hierarchy

Fields are selected independently:

\[
{\rm prescribed}
\rightarrow
{\rm electrostatic}
\rightarrow
{\rm quasineutral\ EM}
\rightarrow
{\rm Darwin}
\rightarrow
{\rm full\ Maxwell}.
\]

The field reduction is a physical model choice, not merely a preconditioner.
Conversely, reduced field equations may be used as a preconditioner while the
target residual remains full Maxwell.

## Edge and open-system models

The same model axes extend from core through SOL/divertor to wall. Cross-field
transport may be prescribed, closure-based or resolved by gyrofluid/GK/FK
turbulence. Neutrals may be fluid, deterministic kinetic or Monte Carlo.
Atomic/reaction, radiation, sheath and wall models are independent but coupled
through explicit particle/charge/momentum/energy exchanges.

Open-system physics does not generally preserve isolated-plasma invariants.
The correct ledger includes reservoir, radiation and wall terms.

## Space-plasma models

Solar-System applications reuse the same hierarchy with open boundaries,
gravity/rotation/radiation pressure, CGL/ten-moment/Landau-fluid closures,
pickup ions, charge exchange, photoionization and surfaces.

Useful mixed configurations include FK ions + fluid electrons, MHD solar wind
+ kinetic pickup ions + kinetic neutral H, and local FK reconnection embedded in
a larger reduced model.

## Stability as a common consumer

Stability is orthogonal to the physical model. A selected MHD, gyrofluid,
gyromoment, DK, GK or FK state may expose a linearized generalized eigenproblem,
forced response, nonmodal growth, continuation/bifurcation and nonlinear
stability analysis when its equations support them.

The linearized action is the same \(R_Uv\) used by Newton and tangents; the
adjoint action \(R_U^\dagger v\) enables resolvent/optimal-growth analysis.
See [stability.md](stability.md).

MHD targets include kink, interchange/ballooning, tearing, RWM and Alfvén
continua/eigenmodes, with NTM and other finite-amplitude branches treated as
nonlinear models rather than labels on a linear eigenvalue. Reduced/kinetic
models similarly cover their supported drift/ITG/TEM/ETG/KBM/MTM and
space-plasma instabilities.

## Macro transport and multiscale modelling

Transport dimensionality and closure fidelity are independent of the
microphysics model. A 1D/1.5D ASTRA-like profile evolution, 2D/3D global
transport and direct full-dimensional kinetic evolution all belong in one
hierarchy.

The closure-provider contract may return fluxes, Onsager matrices, local
coefficients, nonlocal kernels, linear mode data, electromagnetic
cross-spectra or direct nonlinear microphysics fluxes. KAMEL/QL-Balance is an
important internal precedent for response -> quasilinear coefficient ->
profile evolution.

Global orbit transport supplies an intermediate level between local radial
coefficients and global turbulence. It may use GORILLA-like characteristics to
construct 2D/3D transition/nonlocal kernels and may consume turbulence
cross-spectra rather than only scalar diffusion.

See [transport-multiscale.md](transport-multiscale.md).

## Controlled model comparisons

To attribute a difference to physics, freeze as much as possible:

- geometry and physical inputs;
- normalization and source/boundary data;
- observable definition;
- numerical error target.

Then change one model assumption or one inseparable group of terms. If two
certified models give observable intervals \(I_A\) and \(I_B\), their difference
lies in \(I_B-I_A\). That is not automatically a theorem that either brackets
the real plasma.

For different state spaces, compare through explicit lift/projection maps and
physical observables.

## Adaptive fidelity objective

Eventually KIN6D should estimate the observable effect of each unresolved
defect using an adjoint:

\[
\delta Q_k\approx-\lambda^\dagger d_k.
\]

Candidate defects include spatial resolution, velocity/harmonic tails,
collision fidelity, GK/FK reduction, neutral treatment and atomic fidelity.

Refinement may then be ranked by

\[
|\lambda^\dagger d_k|/\text{additional cost}.
\]

This becomes a guarantee only after the corresponding stability and remainder
bounds exist.

## Admission discipline

A new model is implemented only when a concrete case needs it. It must have:

- a complete equation/closure/coupling statement;
- a limiting or structural oracle;
- an external reference where available;
- differentiability semantics;
- certification status;
- performance target.

No model becomes production merely because the generic architecture could
represent it.

Chris&AI
