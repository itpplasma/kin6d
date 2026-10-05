# Moment hierarchies and deterministic structure

Status: required high-priority model branch; all implementations are planned.
Direct steady kinetic science remains the primary development spine.
Gyrofluid and small turbulence cases are early demonstrations using shared
operators, not prerequisites for steady neoclassical or full-kinetic work.

## Parent models and retained moments

Add gyrofluid, Landau-fluid and gyromoment descriptions explicitly to the
[model family](model-family.md). Distinguish a physical closure from its
velocity representation. A Hermite–Laguerre expansion of a named GK model is
a spectral representation; truncating it with a closure changes its finite
equations. A converged hierarchy targets that GK model, not arbitrary FK.

Define the dimensionless parallel velocity and magnetic-moment variables,
Maxwellian/basis weights, Hermite/Laguerre convention, species scales and
phase-space measure before writing coefficients `N_s^{pj}`. Full-f states
need not be close to Maxwellian, but convergence/positivity with sharp tails
is not automatic. If the reference density, flow, temperature or B changes,
include all basis/measure derivatives and frame terms in the moment equations,
JVP/VJP and transfers. Count moments independently of physical-space DOFs.

[Frei–Jorge–Ricci](https://doi.org/10.1017/S0022377820000100) derives a
Hermite–Laguerre gyrokinetic moment hierarchy with variational electromagnetic
fields for plasma-periphery conditions. The
[flux-tube moment study](https://doi.org/10.1017/S0022377823000715) compares
linear gyromoment results with GENE. These give independent model/formula
references; reproducing a finite linear case is not a nonlinear global-GK
convergence or reduction theorem.

## Closures and limits

The selected hierarchy declares retained indices, closure of each omitted
moment, collision projection, gyroaveraging/polarization and field coupling.
Support low-moment gyrofluid, nonlocal/rational Landau-fluid closure, and
explicit higher moments as progressively tested choices. Low moment count
alone does not reproduce every established gyrofluid model: its closure and
ordering must match. Require independently checked projection/reconstruction
and free-energy/entropy transfers when changing the retained set.

[Hammett–Perkins](https://doi.org/10.1103/PhysRevLett.64.3019) motivates
fluid closures matching kinetic Landau damping; [higher-order closures](https://doi.org/10.1103/PhysRevLett.121.135101)
provide further linear-response targets. Such damping represents unresolved
phase mixing; it is not automatically Coulomb entropy production. Nonlocal
parallel closures need a geometry/boundary realization, not a periodic Fourier
sign multiplier silently applied to an open SOL field line.

Isothermal, adiabatic, Braginskii and Landau-fluid closures retain different
assumptions and are not a universal monotone accuracy ladder. Ordinary
pressure moments also support CGL, six/ten-moment anisotropic descriptions
and heat-flux closures, with their own reductions toward FK. A ten-moment
pressure tensor requires an explicit third-moment closure. Increasing a
Hermite cutoff without resolving streaming, closure reflections or recurrence
does not establish kinetic convergence.

Gyrofluid/gyromoment tests vary retained parallel/perpendicular moments and
compare common physical observables at fixed model parameters. Certifying a
finite hierarchy proves that hierarchy only. A GK claim additionally needs
a moment-tail/closure remainder; an FK claim needs the GK reduction bridge.
Nonnegative reconstructed distributions and invariant balances are separate
properties, including for multispecies collision projection.

## Independent gyrofluid landmarks

[Scott's GEM formulation](https://arxiv.org/abs/physics/0501124) derives
coupled gyrofluid terms from fluctuation free-energy consistency. Its GEM
gyrofluid name is distinct from the GEM reconnection challenge in the atlas.
Retain its actual electromagnetic, thermodynamic and FLR assumptions rather
than introducing a single named-code mode.

[Kendl's TIFF paper](https://doi.org/10.1016/j.cpc.2023.108953) supplies a
full-f/full-k gyrofluid reference. The author's
[publication ledger](https://www.uibk.ac.at/de/ionen-angewandte-physik/compsys/team/kendl/)
identifies [T3FF](https://doi.org/10.1016/j.cpc.2026.110137) as a 2026 toroidal
3D **isothermal** gyrofluid drift-Alfvén code. It is not a thermal ITG oracle
merely because it is electromagnetic. Freeze actual releases/equations before
comparison. Full-k denotes the selected FLR treatment; it does not remove
low-frequency, magnetization or other gyrofluid/GK orderings.

[Hamiltonian gyrofluid space-plasma work](https://doi.org/10.1016/j.fpp.2024.100055)
provides additional wave/reconnection references within its specified regime.
Together these guide cheap-to-expensive configurations without importing
JAX/Python or legacy named-code solver branches into kin6d.

## Structure before transport algorithm

Deterministic continuum is the production default for charged distributions,
collision/reaction quadrature, native neutrals, radiation and loss integration.
PIC/Monte Carlo remains an independent oracle or justified optional specialist
backend. EIRENE interoperability is practical but carries statistical or
conditional coupling status; it does not supply deterministic neutral bounds.
Determinism of representation is distinct from reproducible floating reductions
and does not guarantee differentiability across topology/events.

Structure preservation governs the mathematical method. Characteristics,
conservative semi-Lagrangian DG, characteristic Galerkin, variational maps,
exact subflows and splitting are choices for operators where they are useful.
The certified steady residual is a defined continuous/weak problem with a
reconstruction; an interpolation/time-marching algorithm does not define
stationary physics. Characteristics may accelerate transport/preconditioning.

For a closed metriplectic formulation `dot U=J grad H+M grad S`, require
skew J, Poisson/Jacobi structure where claimed, symmetric nonnegative M, and
the degeneracies `J grad S=0`, `M grad H=0`. These yield conservation of H
and nonnegative entropy production. Antisymmetry/positivity alone does not
establish all these consequences. Derive subsystem energy/free-energy and
invariant forms in their natural pairings; different reduced models need
not share one quadratic energy.

Sources, reactions, boundaries and reservoirs have declared transfers.
Krook/OU with fixed targets are open-system closures: OU can exchange energy
with a heat bath and its ordinary differential entropy need not increase.
Test the appropriate relative entropy or complete reservoir ledger instead
of imposing an isolated energy/entropy theorem. Ordinary symplectic stepping
does not generally conserve the exact Hamiltonian; report temporal and
roundoff defects separately from semi-discrete invariants.

## Teaching and turbulence branch

The [student laboratory](student-laboratory.md) starts with a stationary
1X1V configuration and retains transient characteristics/phase mixing.
Its early turbulence door adds HM, HW, interchange/blob, minimal ITG,
then drift-Alfvén/gyrofluid examples. These are reusable configurations of
shared Laplace/bracket/parallel/closure operators. Every selected model first
passes linear dispersion/growth and a free-energy/invariant/source ledger.
Numerical movies alone do not establish model or certificate validity.

Implement a minimal two-dimensional bracket from one FortSym definition,
with primal/JVP/VJP, independent analytical antisymmetry and integration
oracles, alias control and boundary conditions. Continuous bracket identities
do not prove a finite discretization's Jacobi identity or full-time energy
conservation. A sparse/FFT or tensor action uses the same defined projection
and pairing. Small Fourier slab examples need no general 3D/FCI framework.

## Steady scientific progression

The main sequence is a certified steady 1X1V laboratory, 0X3V collision box,
1X3V classical transport, then real-device DK neoclassical coefficients with
NEO-2/SFINCS and applicable rabe asymptotes. Gyromoments and stationary GK
follow, then multispecies nonlinear Landau/BB and global FK/Maxwell. Reactions,
deterministic neutrals and walls extend the same steady formulation.

Keep early HM/gyrofluid dynamics as a secondary structure, derivative and
performance demonstration; it does not delay the publication-oriented steady
DK target or force a time-step-centric core. Reserve compatible moment
interfaces now; implement each concrete closure when its case arrives.
Later comparisons use matching geometry/observable definitions and the
[adjoint defect contract](adjoints-structure.md), with explicit lifting and
tail/reduction proofs before advertising controlled model error.

Chris&AI
