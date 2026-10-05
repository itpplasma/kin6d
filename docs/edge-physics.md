# Edge, reactions, neutrals and surfaces

Status: long-range specification. No edge capability is implemented.

## Independent choices

Edge/SOL/divertor calculations use the same species, bulk, field and time
operators as core calculations. Add independent choices for cross-field
transport, neutral transport, reactions, radiation, sheath treatment and
material surfaces. There is no separate `edge_mode` physics implementation.
Geometry is physical Cartesian space with optional charts and field-line
maps, as specified in [geometry](geometry.md). Nested surfaces are optional;
open lines, islands, stochastic regions, separatrices and X-points are valid
physical domains. Field-line maps must terminate correctly at material walls.

| Axis | Intended choices |
|---|---|
| Plasma | Fluid/Braginskii, PK, DK, GK, FK; compatible bulk models |
| Cross-field transport | Prescribed flux; diffusive/convective closure; resolved fluctuations |
| Neutrals/gas | Fixed, reaction-diffusion, multispecies/compressible fluid, deterministic kinetic, DSMC/Monte Carlo |
| Reactions | Effective/EEDF rates, collisional-radiative/internal-state networks, distribution-aware rates, differential kernels |
| Radiation | Optically thin; opacity/escape; grey/multigroup/moment transport; angular/provider path |
| Sheath | Analytic/kinetic boundary closure; locally resolved FK/Poisson layer |
| Surface | Conductor/dielectric; surface charge; absorption/reflection; recycling/pumping; sputtering/emission; reactive coverages |

A prescribed turbulent transport closure and resolved turbulence must have a
defined scale partition; adding both indiscriminately double-counts transport.
Their closure error is separate from numerical error. A fluid exhaust model
and a kinetic sheath patch must exchange matching particle, charge, momentum,
energy and electromagnetic boundary fluxes. Spatial model adaptation requires
a derived interface and a reduction bound, not merely a changed array layout.

## Neutral transport and practical references

Neutral atoms and molecules are ordinary populations with zero charge and
appropriate mass/internal states. The same hierarchy also supports reacting
heavy gas through reaction-diffusion and multispecies compressible flow before
a kinetic/DSMC description is needed. A deterministic kinetic neutral equation
uses streaming plus elastic, charge-exchange and reaction operators. It shares
physical geometry, reconstruction, sources and output with charged species;
it does not contain a Lorentz force.

Use EIRENE initially as an optional coupling backend or independent oracle
when a real exhaust benchmark needs it. It is not a mandatory dependency.
A native deterministic Boltzmann neutral backend is a later certificate goal.
The [public EIRENE repository](https://github.com/iterorganization/EIRENE)
provides the implementation and database/manual entry points.

Monte Carlo sampling produces statistical uncertainty. A confidence interval
is never labelled a deterministic enclosure. A coupled plasma result cannot
claim a deterministic certificate for the full neutral/plasma problem unless
all neutral contributions used by the theorem have deterministic bounds.
Possible partial results include certification conditional on an explicitly
prescribed neutral source, with the source's statistical uncertainty reported
separately. A fixed random seed improves reproducibility, not mathematical
certainty.

## Species and reaction graph

Explicit charge states `X^0,...,X^Z` form the parent representation. Each can
have its own transport model. Coronal equilibrium, bundled charge states and
other population reductions are declared closures with reconstruction rules;
there is no requirement to evolve 75 full tungsten distributions.

A reaction has reactants/products, stoichiometry, internal-state changes,
rate or differential cross-section data, product distributions and reservoir
exchanges. With signed stoichiometric coefficients `nu_ar`, number-density
sources obey

$$S_a=\sum_r\nu_{ar}R_r.$$

Check charge and nuclei conservation algebraically from every reaction record.
Momentum and energy require the velocity-resolved products and internal
energies as well. Stoichiometry alone cannot prove those balances.
Charge exchange may leave each population's particle count unchanged while
exchanging its velocity distribution; a zero net stoichiometric source must
not suppress its kinetic collision action.

Hydrogen isotope benchmarks begin with atoms and ions. Introduce molecules,
molecular ions, negative ions and vibrational/electronic states when required
by the chosen dissociation, MAR/MAD or detachment case. Reaction definitions
must distinguish isotopes and specify units, thresholds and binding energies.
Impurities enter through the same species/reaction/transport path.

## Atomic-data fidelity

| Level | Data and action | Acceptance requirement |
|---|---|---|
| A0 | Effective rates `n_e n_X K(T_e,n_e)` | Declared electron distribution assumption and table domain |
| A1 | State-resolved collisional-radiative population equations | Level coverage, transition balance and stiffness bounds |
| A2 | Cross sections integrated against distributions | Correct measures, thresholds and rigorously enclosed quadrature |
| A3 | Differential reaction kernels | Product distributions and pair/reservoir momentum/energy exchange |

OPEN-ADAS [ADF11](https://open.adas.ac.uk/adf11) supplies derived effective
ionization, recombination, charge-exchange and radiated-power data classes.
[ADF04](https://open.adas.ac.uk/adf04) supplies resolved level/rate collections
whose included channels vary by file. Inspect the actual dataset comments.
An effective Maxwellian rate is not an arbitrary-distribution cross section;
it cannot generally be inverted uniquely to recover one. Differential product
kernels need additional data when the chosen source does not contain them.

Support adapters for ADAS/OPEN-ADAS, AMJUEL, HYDHEL, H2VIBR and user data as
needed, subject to their distribution terms. Normalize records without
silently changing their physical assumptions. Record dataset version/hash,
source, units, validity range, isotope, interpolation/fit definition and
uncertainty information for every run. Generate owned interpolation,
quadrature and reaction kernels through fortsym; repair missing generic
capabilities upstream. Keep bulk database assets outside the solver repository.

For stationary targets, a distribution-aware rate can have the form

$$R=n_X\int f_e(\mathbf p)\sigma(E(\mathbf p))|\mathbf v(\mathbf p)|\,d^3p.$$

For moving targets use the appropriate relative velocity and two-population
integral; a relativistic process needs its derived flux/cross-section
convention. Non-Maxwellian electrons do not automatically justify evaluating
a Maxwellian rate at their moment temperature.

## Sheath and surface operators

Separate the plasma-to-sheath coupling from the material response. Cheap
Bohm/sheath and transmission closures, intermediate kinetic closures and
locally resolved FK Debye layers are distinct models. Global core resolution
need not resolve a local Debye sheath, but omitted sheath/interface effects
need explicit bounds before certifying a larger-domain observable.

A surface operator maps incoming traces and material state to outgoing
traces, reservoir exchanges and material evolution. Use oriented physical
normals: incidence has `v·n<0` for the outward plasma-domain normal, with
incident flux weighted by `|v·n|`. Momentum distributions use the same physical
measure as the plasma. The surface state can initially be a material label
and fixed temperature; inventory, composition and temperature evolve later.

Introduce absorption, specular/diffuse reflection, neutralization, recycling,
pumping, sputtering, secondary/thermionic/field/photoelectron emission,
dielectric surface charge, reactive coverages, deposition/etching and
desorption only when a benchmark requires them. A sputtering yield depends on incident energy
and angle and must be paired with a normalized emitted-particle distribution.
Emitted atoms become ordinary neutral species that can ionize and return to
the wall. The reaction/transport code has no special tungsten-only branch.

Recycling bookkeeping counts nuclei: one returned molecule contains two
hydrogen nuclei. Track returned, implanted, retained and pumped inventory,
along with external replenishment. Charge deposited in a wall must appear in
a declared current balance or external circuit. An evolving wall potential
uses a specified capacitance and total currents, or consistently uses their
per-area counterparts.

## Energy and uncertainty accounting

Include binding/excitation energies, photon energy, wall deposition,
reflected/sputtered energy and external heating/pumping/open-boundary fluxes.
An optically thin loss is energy exported to a named radiation reservoir.
The closed-system balance, when its assumptions hold, contains plasma,
neutrals, internal states, radiation and material energy. Driven/open cases
report the complete ledger rather than forcing every subsystem to conserve
its energy independently.

A certificate encloses the solution of the stated reaction/table/surface
model. Roundoff and interpolation/quadrature error can be bounded separately
from physical data uncertainty. Missing experimental uncertainty is not zero
uncertainty. A theorem robust over an input-data interval must establish its
premises uniformly over that interval; a numerical band for one nominal table
does not cover all physically plausible rates or surface yields.

## Benchmark progression

Start with an exactly balanced homogeneous reaction or neutral free-streaming
case. Then add a slab recycling boundary, charge exchange, a sheath benchmark,
a reduced transport/neutral coupling case, and finally a physical-device
exhaust comparison. Each increment retains earlier limit curves and their
whole enclosures in generated plots; statistical error bars use a different
legend and cannot substitute for guaranteed bands.

The broader weakly/partially ionized formulation is specified in
[reacting-partially-ionized.md](reacting-partially-ionized.md), and radiation
fidelity in [radiation.md](radiation.md).

SOLPS-ITER, EMC3-EIRENE, SOLEDGE3X-EIRENE, GBS, GRILLIX and BOUT++/Hermes
are planned references for their overlapping fluid/exhaust models. XGC,
GENE-X, Gkeyll, PICLS and COGENT are planned kinetic/GK comparisons. This list
does not assert identical equations or existing certified solutions.
[GENE-X](https://genex.ipp.mpg.de/about) and
[XGC](https://xgc.pppl.gov/html/general_info.html) already provide important
kinetic edge precedents. TCV-X21 is a planned experimental comparison whose
measurement uncertainty remains distinct from numerical certification.

Chris&AI
