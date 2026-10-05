# Solar System plasma scope

Status: long-range application map. No space-plasma solver is implemented.

## Physical scope and reuse

The intended envelope includes solar atmosphere/corona, solar wind and
heliosphere, planetary magnetospheres, ionospheres/exospheres, comet/moon
plumes and dusty/embedded-body environments. It uses classical or
special-relativistic electromagnetic kinetics, fluid/MHD reductions,
partially ionized/reacting matter and specified neutral/atomic physics. General relativity,
pair creation, nuclear networks and comprehensive astrophysical radiation
transport are outside the planned scope.

Physical Cartesian vectors and general domain maps support these applications
without a second plasma architecture. Cartesian boxes, spherical shells and
body-conforming meshes are chart choices. [Geometry](geometry.md),
[coupling roles](coupling.md), [model family](model-family.md) and
[reaction/surface physics](edge-physics.md) retain the same contracts.
Global calculations still require their own affordable discretization,
consistent closures and certificate theorems; geometric generality alone does
not make global full kinetics practical.

| Application | Possible supported configuration, once derived |
|---|---|
| Reconnection/current sheet | FK electrons and ions, Maxwell |
| Collisionless shock/magnetosheath | FK ions, FK or fluid electrons |
| Solar-wind turbulence | GK/FK ions, DK/FK electrons, kinetic minority alphas |
| Global magnetosphere | MHD bulk, kinetic ring-current/energetic populations |
| Heliosphere | MHD/multi-fluid wind, kinetic pickup ions and neutral hydrogen |
| Radiation belts | DK/orbit-averaged electrons with declared wave diffusion |
| Ionosphere/exosphere/plume | Partially ionized fluid/kinetic populations, neutral winds, reactions and photon drivers |
| Dusty/body interaction | Charged grains/finite bodies with sheath, charging, emission and external forces |

Species model, active/test/reservoir role and field self-consistency remain
independent choices subject to compatibility. Collisionless `C=0` is an
important configuration; phase mixing, sharp beams and instability can still
make its computation and certification difficult.

## Forces and frames

The parent characteristics permit a specified external force,

$$\dot{\mathbf x}=\mathbf v(\mathbf p),\qquad
\dot{\mathbf p}=q(\mathbf E+\mathbf v\times\mathbf B)+\mathbf F_{\rm ext}.$$

Gravity, radiation pressure and rotation enter only with a stated frame,
force law and energy/source bookkeeping. Nonrelativistic rotating-frame
Coriolis/centrifugal formulas are not automatically valid relativistic
momentum equations. Derive any such extension in its intended regime.

For a distribution defined against physical phase-space volume, arbitrary
forces require the conservative phase-space equation

$$\partial_t f+\nabla_x\cdot(\dot{\mathbf x}f)
 +\nabla_p\cdot(\dot{\mathbf p}f)=C+S.$$

The simpler advective form is valid only when the relevant phase-space flow
has zero divergence (or the correct invariant measure is used). Radiation
reaction or velocity-dependent drag must not silently use a Hamiltonian
volume-preservation assumption. Generate force, Jacobian and enclosure
kernels through fortsym.

## Boundaries and closure choices

Support mathematically specified inflow, characteristic outflow, absorbing,
reflecting, periodic, conducting and material-boundary conditions. A
nonreflecting approximation has its own defect. A distant finite boundary
does not prove enclosure of an infinite-domain solution without a separate
truncation/decay estimate. Mapped unbounded domains also require rigorous
tails and regularity.

Space fluid closures extend through isotropic MHD, anisotropic/CGL,
Landau-fluid/nonlocal heat flux and kinetic species. CGL double-adiabatic
invariants require their collisionless magnetized ordering and neglected
heat-flux assumptions; they are not exact invariants of arbitrary FK states.
Coronal Spitzer conduction, optically thin radiation and prescribed heating
are distinct constitutive/source choices. Photospheric/chromospheric or other
optically thick work may move up the common radiation hierarchy through
opacity, diffusion and multigroup/moment transport; highest-fidelity non-LTE
spectral transfer remains a provider boundary. No Maxwellian reference is
structurally mandatory for the FK state.

Initializers and imported states can contain beams, loss cones, rings,
bi-Maxwellians, kappa distributions and pickup shells. Positivity, normalization,
regularity and momentum-tail integrability are explicit input obligations.
Sharp anisotropy may require many harmonics. Add rotated/adaptive angular or
alternative velocity representations only when measured accuracy/performance
requires them, with a reconstruction/certificate path; there is no claim that
spherical harmonics are optimal for every distribution.

## Wave, neutral and surface interactions

A quasilinear/wave operator is separate from Coulomb collisions. Specify its
coordinates/actions, measure, diffusion tensor, boundaries, wave spectrum and
ordering. A symmetric positive-semidefinite diffusion tensor alone does not
prove conservation of the full plasma/wave system. An imposed spectrum is an
external energy/momentum reservoir unless its dynamics are coupled. Partition
resolved fluctuations and modelled wave scattering to avoid double counting.

Kinetic charge exchange supplies velocity redistribution, energetic neutral
atoms and pickup-ion sources. ENA detector observables require defined
line-of-sight, attenuation, instrument response and their integration bounds.
Photoionization adds photon-driven reaction rates from a specified incident
spectral photon flux, with shielding/geometry and energy accounting. A scalar
unattenuated rate is a closure, not a photon-transport solution.

Planetary bodies and spacecraft use the same surface operator for absorption,
reflection, charging, photo/secondary emission, sputtering and desorption.
Capacitance/current conventions and the external electrical reservoir must
be defined. Dust/grain charging, finite embedded bodies and radiation pressure are common
particle/surface/force capabilities admitted when a real case needs them; if charge is `q=Z e`, its evolution uses `e dot Z=I`, with
consistent current units. They do not justify a speculative particle framework.

## Spatial hybrids and evidence

Long-range targets include MHD globally with kinetic bow-shock/reconnection
patches, and fluid solar wind with kinetic minorities/neutrals. Interfaces
must preserve the appropriate population ownership, field constraints and
flux balances. Model switching needs a controlled reduction and interface
error estimate. Kinetic patch existence does not certify a whole global model.

The permanent atlas should cover Parker wind and spiral in their specified
reduced equations; Harris/force-free sheets; linear whistler, Alfvén/KAW,
Weibel, mirror/firehose and ion-cyclotron problems; and pickup-ion motion.
GEM reconnection, Orszag–Tang, Brio–Wu and collisionless-shock comparisons are
standard numerical families, not automatically analytic truth or certified
weak solutions. Shock/instability cases need a norm, solution concept and
appropriate existence/selection theorem before assigning proof status.

Planned cross-code references include Vlasiator, Gkeyll, Smilei, WarpX,
PIConGPU and Athena++/PLUTO within their actual model overlap. Activate a
reference only after fixing equations, boundary conditions, units and inputs.
Each promoted capability generates candidate/reference bands and an earlier
limit overlay. Numerical, data/model and statistical uncertainties retain
separate labels. The cheapest supported model is selected only with evidence
that it answers the requested observable to the stated tolerance.

Chris&AI

See [reacting and partially ionized plasmas](reacting-partially-ionized.md),
[radiation](radiation.md), and the [capability catalogue](capability-catalogue.md)
for shared cross-domain primitives.
