# One mathematical engine, selectable plasma models

Status: long-range model map and coupling obligations. The first executable
product remains the small certified steady reference. No MHD, reduced kinetic,
hybrid, equilibrium importer, or transient implementation exists in kin6d.

## Chameleon philosophy

Use the cheapest mathematically supported model that demonstrably answers the
scientific question at the requested tolerance. Full kinetics is the parent
reference, not a requirement to pay its cost for every calculation. Keep one
engine for physical actions, constraints, steady/transient consumers,
reconstruction, observables, and certification. Different mathematical models
have different state spaces and actions within that engine.

| Axis | Intended choices |
|---|---|
| Species | Fixed, adiabatic, ordinary/gyro/Landau fluid, gyromoments, PK, orbit averaged, DK, GK, FK |
| Bulk plasma | None, multi-fluid, ideal/resistive/extended MHD, kinetic–MHD |
| Collision pairs | None, Krook/BGK, Lorentz/LBD, model linearized, Landau, BB, later research closures |
| Field coupling | Prescribed, electrostatic, full Maxwell, explicitly reduced field action |
| Geometry/domain | Slab, cylinder, flux tube, annulus/surface, global torus, supported general 3D domain |
| Formulation | Equilibrium, linear forced response, spectral stability, nonlinear evolution |
| Time consumer | Steady or transient, using the same selected physical actions |
| Execution | CPU reference, optional generated CUDA kernels, distributed execution when needed |
| Evidence | Analytic references, cross-code comparisons, numerical certificates, formal proof coverage |

These are configuration axes, not a claim that every Cartesian product is
valid or implemented. Derive and register supported combinations with their
ordering, conserved quantities, field action, boundary data, and certificate
theorem. Reject missing bridges explicitly. Cross-code validation and
certification are complementary evidence, not mutually exclusive modes.

NEO-2-like means a specified steady neoclassical DK/linearized-FP case;
GENE-X-like means a specified global full-f GK transient case. These names
describe model intersections for comparisons, never executable architectures
such as `neo2_mode`, `genex_mode`, `mars_mode`, or `jorek_mode`. Thin application consumers may eventually select mathematical configurations
from shared libraries. Existing codes
remain independent validators, not presumed interchangeable implementations.

The strongest invariant is shared mathematical/numerical infrastructure:
change state space, derived reduction/closure or physical brackets instead of
creating another numerical universe. The [structure policy](moment-models.md)
combines a reversible Hamiltonian core, thermodynamically constrained metric
dissipation and open-system transfers, with explicit proof hypotheses.
The [student gallery](student-laboratory.md) and later real-device states
use the same residual, derivative, structural and certificate contracts.

## Controlled physics comparisons

Freeze physical inputs, equilibrium provenance, domains, source/boundary
data, normalization, and physical observable definitions. Change one physical
assumption at a time when possible: species reduction, field approximation,
collision action, Ohm law, or kinetic feedback. Retain the same geometry and
field infrastructure where the model allows it, with convergence-independent
continuous bounds for both results.

Different models may require different phase-space measures, meshes, or
polarization/magnetization. Match the physical observable through proved
pullbacks rather than forcing identical arrays. If one change requires
several coupled terms, report that complete change; a six-percent difference
cannot automatically be attributed to one named mechanism.

For certified observables Q_A in I_A and Q_B in I_B, the difference is in
I_B-I_A. Branch correspondence and a model-reduction remainder are additional
claims. Adaptive choice of model or spatial fidelity becomes guaranteed only
when numerical plus reduction/coupling errors can be bounded. An uncertainty
estimate or visually small difference alone does not justify promotion.

## Bulk plasma and MHD

One-fluid MHD aggregates populations into bulk density, velocity, pressure
and fields. It is not another value of a single species' `fluid` flag.
Record which populations belong to each bulk, their equation of state,
transport/Ohm closures, and which populations are evolved independently.
Every charged population is accounted for once in charge, current, pressure,
inertia, energy, and collision exchange.

The intended branch covers static ideal equilibrium, linear ideal forced
response, dynamic ideal modes, linear resistive response, kinetic–MHD,
extended/two-fluid MHD, and eventually nonlinear extended MHD. These are
distinct models, not a universal increasing-accuracy chain. Ideal and
resistive field actions are quasineutral/low-frequency reductions; they do
not silently retain all full-Maxwell waves or arbitrary relativistic physics.

For static isotropic no-flow ideal equilibrium, the selected model is
`J cross B=grad p`, `div B=0`, with matching Ampere law, boundary/exterior
data, and flux/current constraints. Flows, anisotropic pressure, electric
forces, and external/source momentum add their own terms. Static kinetic
moment balance is not generically that simple no-flow ideal equation.

Linear forced response differentiates the constrained model about an enclosed
equilibrium. It must include equilibrium uncertainty and background residual.
Dynamic eigenmodes add inertia, spectral parameter and eigenvector
normalization. Resistive response additionally needs its selected Ohm law,
induction, thermal/source closure, and wall/exterior conditions. A resistivity
switch alone does not implement tearing layers or resistive-wall modes.

Extended Ohm laws can include Hall, electron pressure, inertia and transport.
Derive the complete chosen law and energy balance from the specified bulk or
multi-fluid equations. A schematic list of terms is not an executable closure.
Nonlinear MHD/XMHD, disruption and ELM regimes are long-range cases, triggered
by a concrete scientific need rather than a general plugin framework.

## Hybrid coupling

A species hybrid can use thermal e+i in an MHD bulk and energetic alphas in
DK/GK/FK. Choose a consistent current-coupling or pressure-coupling formulation
and account for energy, momentum, charge, equilibrium subtraction, and
kinetic feedback. Adding kinetic pressure to a bulk that already includes
that population double-counts its response. A frozen-field kinetic diagnostic
and a self-consistent kinetic–MHD solve carry different claims.

A spatial hybrid can use fluid/MHD bulk physics and kinetic resonant layers.
Its interface or matched overlap must preserve the required fluxes, charge,
energy and asymptotic matching, with a bounded interface/model defect.
Layer location, width, rational-surface uncertainty, and singular limits
enter the certificate. Applying a spatial mask to unrelated residuals is
not a consistent coupled model. Species and spatial hybridity are separate
axes; implement neither without a derived first case.

Imported MHD equilibrium can initialize kinetic reconstruction or numerical
continuation toward DK/GK/FK. Such a lift needs a distribution/pressure closure,
constraints, and source compatibility; many distributions share its moments.
Continuation across changing state spaces needs defined transfer maps and
branch assumptions. It is not an automatic proof that every MHD equilibrium
is the limit of a stationary collisional kinetic solution.

## Equilibrium candidates

Imported equilibria are supported before internally solving every equilibrium
family. Add the next needed adapter with units, COCOS/sign conventions,
topology, input hashes, geometry accuracy and boundary assumptions:

- Axisymmetric G-EQDSK/EFIT/CLISTE/CHEASE-style inputs.
- VMEC wout/VMEC++, DESC or GVEC nested-surface candidates.
- SPEC/MRxMHD candidates with explicitly represented interfaces, pressure
  jumps, islands or non-nested regions and their own weak equations.

An imported equilibrium is a numerical input or untrusted candidate. It is
not certified merely because an external code reports a small force residual.
An internal axisymmetric Grad–Shafranov consumer is the first proposed
equilibrium solver: fixed boundary before free boundary, coils/passive
structures, X-points and eventually inverse/evolving problems, one difficulty
at a time. Exact Solov'ev/shaped cases and named external solvers provide
independent oracles. Python-based external reference tools do not create a
required Python layer in kin6d.

For 3D equilibrium, first evaluate and certify the external candidate against
its precisely stated continuous model. VMEC nested surfaces and SPEC relaxed
interfaces require different equations and spaces. Gauge, flux, relabeling,
nullspaces, resonance/regularity and stability/complement estimates must be
handled before a useful local enclosure exists. This certification is an
analytical research milestone, not a generic residual-to-radius conversion.
Write an internal 3D iteration only when the external-candidate route leaves
a demonstrated scientific or performance need.
[VMEC++](https://proximafusion.github.io/vmecpp/) is an external fixed/free-boundary
ideal-MHD candidate source. [SPEC](https://github.com/PrincetonUniversity/SPEC)
provides a different relaxed model and topology scope.

## Independent references and certification

The local ARES README at revision
533866a74d2d30909a6a62be2913fb041c796462 describes linear ideal forced
MARS-F response and particular Shaing/frozen-field drift-kinetic NTV paths.
It excludes eigenvalue/RWM, resistive-wall, time-evolution, and energetic
branches in that scope. Keep it an independent, versioned oracle; do not
claim all MARS-K feedback is already implemented there or migrate its code
before a shared-library need and equivalence evidence exist.

[GPEC's documented suite](https://github.com/PrincetonUniversity/GPEC/blob/develop/CLAUDE.md)
includes ideal, resistive and kinetic response components.
[JOREK](https://github.com/iterorganization/JOREK) supplies nonlinear extended
fluid/hybrid precedent. Each comparison freezes the actual equations,
closure, wall, frequency/rotation, equilibrium, and diagnostic definitions.
Their breadth does not establish the kin6d implementation or certificates.

Static equilibria, forced response, spectral stability and trajectories need
different theorems. Ideal force operators may have neutral or unstable modes;
kinetic/resistive operators can be non-Hermitian or non-normal. Near resonance,
certification needs constraints, parameter isolation and rigorous resolvent
or stability bounds, not a generic small residual. Eigenvalue enclosures do
not by themselves prove bounded forced response or a long-time trajectory.

The architecture keeps only mathematical seams needed by the next case:
active state, selected closures, evolution/linearized actions and constraints,
physical observables, and independent continuous certificate obligations.
Every owned numerical leaf remains fortsym-generated; CPU/CUDA choice does
not change the model. Preserve the first runnable demo and all earlier
limiting cases as these configurations appear.

## Expanded scientific map

[Geometry](geometry.md) makes Cartesian physical vectors fundamental and
retains cylindrical/Fourier as a useful initial fusion chart. Figure-8 is a
planned mandatory geometry gate; FCI is a separate map service. [Coupling](coupling.md)
distinguishes the fixed equilibrium, external drive and induced response.

[Edge physics](edge-physics.md), [Solar System scope](space-plasmas.md) and
[extensions](extensions.md) reserve neutral/reaction/wall, wave/heating,
field-reduction, fusion-source, circuit and diagnostic actions. The combined
state only contains selected active blocks. A shared library may acquire thin
application consumers as needed; one engine does not require one executable.
Unsupported combinations fail explicitly. Every feature earns production
status with a scientific use case, derived coupling, an independent oracle
and an actual certification/performance target.

Chris&AI
