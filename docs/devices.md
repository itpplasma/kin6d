# Real-device applicability

Status: acceptance policy and planning estimates, 2026-10-05. No device
benchmark, equilibrium file, or hardware performance result is implemented.
Full kinetics is the reference model; reduced models make additional domains
and evolution problems practical. [Species models](species-models.md) defines
the permissible reductions and their mathematical obligations.

## Scientific domains

| Domain | Intended full-kinetic use | Reduced/hybrid use |
|---|---|---|
| Homogeneous or 1X3P | Collision reference and certification | Limiting-case validator |
| Local physical patch | Operator and reduction comparisons | Matched bridge problem |
| Annulus or full magnetic surface | Orbit, transport, and response reference | Practical surface calculations |
| Global axisymmetric core | Primary steady full-orbit target | Cheaper transport and continuation |
| Stellarator field period | Three-dimensional steady reference when affordable | Full-surface/global production |
| Sparse nonaxisymmetric response | Full-orbit reference with certified Fourier complement | Reduced response comparison |
| Whole-device turbulence | Selected small reference problems | Gyrokinetic production objective |
| Wall and Debye sheath | Local resolved reference | Explicit edge/sheath model with its own scope |

This table is a development target, not a feasibility theorem. Flux tubes can
be useful reduced turbulence benchmarks; they are not the flagship geometry
for full-orbit steady kinetics. Gyrokinetics changes the equations and scale
requirements in addition to removing a velocity angle. GENE-X demonstrates
the practical role of global full-f gyrokinetics in tokamak and stellarator
core/edge/SOL calculations; it does not supply kin6d's certificates.
[GENE-X description](https://genex.ipp.mpg.de/about)

## Device envelope

Use published machine scales to select cases, then freeze a particular
equilibrium and profiles. Catalogue maxima do not define one simultaneous
discharge or its species temperatures and densities.

| Device | Published scale | Planned benchmark |
|---|---|---|
| TCV | Major radius 0.88 m; maximum toroidal field 1.54 T in the cited presentation | Low-field, strongly shaped tokamak; cheap reduced case before FK promotion |
| ASDEX Upgrade | Major radius 1.6 m; magnetic field 3.2 T; plasma density up to 2e20 particles/m3; temperature about 100 million degrees | Axisymmetric core, then a selected RMP response |
| SPARC | Major radius 1.85 m; on-axis field 12.2 T in the 2020 design | Compact high-field burning-plasma design case |
| W7-X | Major radius 5.5 m; magnetic field 3 T; plasma density 3e20 particles/m3; temperature 60–130 million degrees | One field period when symmetry permits, then a symmetry-breaking case |
| LHD | Major radius 3.9 m; field 3 T in the published machine table | Helical-coil geometry and transport |
| ITER | Major radius 6.2 m; toroidal field up to 5.3 T in the cited design | Reactor-like axisymmetric core |
| ARC | Major radius 3.3 m; on-axis field 9.2 T in the 2015 conceptual study | Compact reactor tokamak; distinct from future ARC revisions |
| DEMO-like | A named, dated design point must be selected | Later reactor comparison; no universal DEMO parameter set assumed |
| Proxima Stellaris | QI, four field periods; major radius 12.7 m; average field 9 T in the 2025 concept | Reactor-scale quasi-isodynamic stellarator |
| Thea Helios | QA, two field periods; major radius 8 m; aspect ratio 4.5; on-axis field 6 T in revision v3 | Reactor-scale quasi-axisymmetric stellarator |
| Additional compact stellarators | Public, versioned equilibria selected when available | Fill size/topology gaps without inventing unavailable inputs |

The AUG and W7-X scales come from their official technical descriptions.
These density entries are catalogue plasma densities, not automatically the
electron density of a chosen case. The ITER values identify a published
design envelope, not an available experimental discharge.
[AUG](https://www.ipp.mpg.de/16208/einfuehrung),
[W7-X](https://www.ipp.mpg.de/16931/einfuehrung),
[ITER design presentation](https://www.iter.org/sites/default/files/media/2025-07/i-1_loarte.pdf)

The other entries identify published experimental descriptions or specific
conceptual designs. They establish a planning span of about 0.88–12.7 m and
1.54–12.2 T; field definitions differ and are retained per case. This is not
a rectangular box of simultaneously attainable plasma parameters. QA is a
property of the magnetic-field strength in suitable coordinates, not exact
axisymmetry of the three-dimensional geometry.
[TCV](https://crppwww.epfl.ch/~weisen/TCV_PRESENTATION.HTML),
[SPARC overview](https://doi.org/10.1017/S0022377820001257),
[LHD](https://www-lhd.nifs.ac.jp/pub/LHD_Project.html),
[ARC v2](https://arxiv.org/abs/1409.3540v2),
[Stellaris concept](https://publikationen.bibliothek.kit.edu/1000179851/172386752),
[Helios v3](https://arxiv.org/html/2512.08027v3)

Maintain this envelope as a versioned set of physical missions, not merely
optional validation names. The progression spans TCV/AUG/LHD/W7-X,
SPARC/ITER/ARC/DEMO, and QA/QI reactor stellarators. A later design revision
adds a new case rather than overwriting the archived design point.

Before a device benchmark becomes a gate, its manifest records:

- Input provenance, access/redistribution conditions, exact file hashes,
  equilibrium convention, units, and coordinate/frame transformations.
- Physical masses, charges, species profiles, sources and reservoirs;
  dimensionless parameters computed from those inputs and stated definitions.
- Domain, symmetry, fields, species and pair models, boundaries, exterior
  treatment, and whether inputs are fixed values or enclosed uncertainties.
- Requested observables and certified tolerances, independent limiting-case
  mappings, runtime/memory budget, and expected plots with interval bands.
- Certificate theorem, norm, assumptions, and complete numerical/input error
  scope. Missing experimental uncertainty remains explicitly outside the claim.

Select concrete inputs as soon as the first electron-in-prescribed-field
case needs them. Keep the reference solver small before that trigger; do not
build speculative importers or download a catalogue of equilibria now.

## Production gate

For every selected machine keep three independently reported verification
levels: geometry smoke test, small physical-parameter kinetic case, and
production-resolution certified case. Geometry checks include independently
known field values, derivatives, periodic maps, measures, boundary normals,
and divergence/field representations where applicable. They do not certify
a plasma solution. Small cases can initially use adiabatic electrons and DK
ions where that field/closure combination is mathematically supported.

CI runs affordable geometry and small-case gates from pinned public inputs;
scheduled or explicitly launched resource jobs run larger production gates.
Do not make workstation CI depend on a full-device solve, private downloads,
or automatic access to a protected host. Report unavailable inputs and
unsupported reductions explicitly. Promote the same frozen equilibrium
through species and collision models, retaining earlier plots and intervals.

A component reference milestone can pass its analytic continuous-model gate
without solving a whole machine. Production promotion additionally requires
a relevant AUG-like and W7-X-like benchmark at physical mass ratios and
physical dimensionless parameters, with certified requested observables and
meaningful plots. Local patches or surfaces are acceptable only when they
match the capability's declared scientific scope. A claim of global operation
requires a global benchmark.

Record each algorithm's scope and the two device cases before promotion.
Artificial masses or enlarged gyroradii may diagnose bugs but cannot replace
these gates. A feature useful only in one geometry remains scoped and
experimental until its advertised cross-device gate exists. ITER/DEMO points
extend the envelope later; they do not delay the first useful electron case.

The first practical branch evolves or solves kinetic electrons with prescribed
ions, first in prescribed fields and then with electron-self-consistent
fields. Fixed ions are an external reservoir, not a closed multispecies
conservation test. Promoting the ions to fluid or kinetic response preserves
that earlier case and plots the certified response difference.

## Geometry and field periods

The physical geometry contract is generic three-dimensional field and domain
evaluation: `B_vac(x,t)`, `E_ext(x,t)`, derivatives required by the selected
model, and wall/domain geometry with normals and boundary labels. Certified
evaluations additionally provide enclosures and interpolation/tail bounds.
Tokamak equilibrium formats are import adapters, not fundamental physics
types. Initially implement the smallest analytic evaluator and domain needed
by P0/P1; add a data adapter only for the next frozen machine case.

Symmetry is declared and checked case metadata. Axisymmetric tokamaks may
use n=0; LHD uses its selected helical/toroidal periodicity, and QA/QI designs
use their native field periods. Breaking symmetry activates a larger domain
and spectrum. Orbit and field operators never infer symmetry from a machine
name, QA/QI label, or the existence of nominally repeating coils.

W7-X has fivefold device symmetry. A one-period calculation requires the
selected field, boundary, sources, reservoirs, and initial data to share
that symmetry. In full-torus Fourier notation it retains multiples of five;
on a period coordinate, differentiation includes the corresponding factor.
All transforms, momentum-frame rotations, and normalization must respect
the period map. [IPP symmetry description](https://www.ipp.mpg.de/ippcms/eng/presse/archiv/05_05_pi)

A certificate in a symmetric subspace establishes its stated local solution
claim there. It does not exclude symmetry-breaking branches or certify their
stability. Nonperiodic perturbations require the larger domain and modes.
Never assume W7-X needs only a few toroidal harmonics because tokamak RMP
examples do; its equilibrium spectrum must itself be resolved and bounded.

The initial cylindrical Fourier representation handles three-dimensional
fields on simple domains. A realistic toroidally varying plasma boundary
may need mapped elements or another boundary representation. Choose that
seam when the frozen W7-X case requires it, deriving transformed transport,
field, measure, and certificate terms together. This is not an axisymmetric
solver rewrite or a promise that a rectangular mesh fits every boundary.

## Debye scales and numerical stiffness

An asymptotic-preserving formulation or preconditioner can remove a numerical
stability restriction without replacing full kinetic equations. Published
semi-Lagrangian Vlasov–Maxwell work demonstrates this opportunity in the
quasineutral limit. Its scope does not establish a bound for our collision
system or permit unresolved sheaths and waves.
[Liu et al., JCP 2025](https://www.sciencedirect.com/science/article/pii/S0021999125001238)

For a full-model claim the certificate still controls the continuous spatial
and temporal complements. Small Debye scales with negligible solution
amplitude need not impose uniform refinement; an actual thin boundary layer
does. A quasineutral, Darwin, or gyrokinetic field model is instead an explicit
physics choice with its own equations and model-difference status. Never
identify implicit stability with fast-wave accuracy.

## Memory envelope

For equal species grids, real-equivalent FP64 kinetic state bytes are
`8 Ns N_RZ Np Nangle (1+2N)`, where `N` is the highest positive toroidal mode
for a contiguous spectrum. Spatial counts are actual scalar degrees of
freedom, including element order. Independent complex storage may cost more.

For four species, 48 radial coefficients, and 128 x 256 spatial degrees:

| Positive modes through N | 39 angular coefficients | Initial full L=8: 81 angular coefficients |
|---|---:|---:|
| 0 | 1.83 GiB | 3.80 GiB |
| 1 | 5.48 GiB | 11.39 GiB |
| 2 | 9.14 GiB | 18.98 GiB |
| 4 | 16.45 GiB | 34.17 GiB |
| 8 | 31.08 GiB | 64.55 GiB |

The 39-coefficient column assumes L=8 with a restricted angular order and
is a later optimization scenario. A locally B-aligned truncation needs frame
connection terms and a certified angular complement; it is not the initial
full-harmonic implementation. These are storage calculations, not evidence
of adequate physical resolution or of workstation fit.

At 256 x 512 degrees, N=4 costs 65.81 GiB with 39 coefficients and 136.69 GiB
with 81. N=16 costs 241.31 or 501.19 GiB respectively. An ideal eight-way
distribution gives 30.16 or 62.65 GiB per worker for that last state alone.
Eight devices with 96 GiB each are a planning scenario, not verified hardware
or a performance claim. Fields, collision potentials, padded transforms,
Krylov vectors, preconditioners, communication buffers, and the certifier add
material memory. Distribution does not make these costs disappear.

Select a scientific tolerance and hardware budget together. Verified
residual and stability estimates decide the required spatial, radial,
angular, toroidal, and temporal refinement. Resolution alone cannot certify
that a machine case is affordable. Profile the complete solve/certify/plot
workflow before introducing compression, GPU kernels, or distributed storage.

Chris&AI
