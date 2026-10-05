# Species models and hybrid coupling

Status: canonical species-model policy. Only the first 1X1V laboratory is
implemented.

## Independent species model selection

Each population selects its physical model independently, subject to compatible
field and collision/reaction coupling.

| Model | State | Retained physics |
|---|---|---|
| FK | \(f(x,p)\) | full-orbit kinetic dynamics |
| GK | gyrocenter distribution | ordered gyrocenter dynamics and retained FLR |
| DK | guiding-center distribution | parallel kinetics, trapping and drifts |
| PK | reduced parallel kinetic state | selected parallel response |
| bounce/orbit averaged | invariant/orbit distribution | slow averaged kinetics |
| gyromoment | retained Hermite--Laguerre-like moments | systematic kinetic moment content |
| gyrofluid/Landau-fluid | low moments plus closure | selected FLR/phase-mixing response |
| fluid | declared moments | chosen moment closure |
| adiabatic | algebraic response | instantaneous closure |
| fixed | prescribed state | external reservoir/background |

These labels are not array-shape shortcuts. Every reduction declares its phase
space, measure, Hamiltonian, fields, closure and observable mapping.

A separate bulk MHD/multi-fluid model may aggregate explicitly owned
populations. Bulk ownership is not another per-species label.

## Representation is separate

For any kinetic/reduced distribution state, also select the representation
described in [representations.md](representations.md):

- full-f, delta-f or control variate;
- Eulerian, semi-Lagrangian, particles or blobs;
- PIC/PIF/particle-in-basis if particles source fields;
- deterministic/QMC/MC sampling where applicable.

Thus "GK particles" and "GK continuum" are two representations of a named GK
model, not two unrelated physical models.

## Intended mixed configurations

Examples include:

\[
e:{\rm DK},\quad D/T:{\rm GK},\quad \alpha:{\rm FK},
\]

\[
e:{\rm adiabatic},\quad i:{\rm GK},
\]

\[
\text{thermal bulk}:{\rm MHD},\quad
\text{energetic species}:{\rm DK/FK},
\]

and

\[
\text{thermal continuum}+\text{particle energetic tail}.
\]

Each requires one complete field and exchange accounting. Unsupported pairings
fail before solving.

## Coupling roles

Species model and coupling role are independent. A population may be:

- active and field-sourcing;
- test/passive but evolved;
- fixed reservoir;
- part of a bulk model.

The role determines whether it evolves, sources fields and feels induced
response. Test-particle status does not imply DK/FK/PIC and vice versa.

## Reduced field coupling

DK/GK states require consistent transformed charge/current, polarization and
magnetization at the chosen ordering. A reduced species cannot simply insert
raw moments into an incompatible full-Maxwell action.

Conversely, full Maxwell is not required for every kinetic species. A
quasineutral/Darwin field model may be the named target model, while a reduced
field solve may also serve only as a preconditioner for a fuller target
residual. Keep those cases distinct.

## Mixed collision pairs

For representations/models \(A,B\), a useful schematic is

\[
C_{ab}^{A\leftarrow B}
=
P_a\,Q_{ab}[L_a u_a,L_b u_b],
\]

where \(L\) lifts to a named parent phase space and \(P\) projects the collision
action back. This notation is only a design pattern: each actual lift/projection
must be derived and checked.

Reduced-space collisions must account for phase-space measure, location shifts,
polarization/magnetization where relevant and pairwise momentum/energy
exchange. A low-dimensional conservation correction cannot repair an
inconsistent physical lift.

Fluid, adiabatic and fixed partners require explicitly defined exchange
closures/reservoir behavior rather than pretending every population owns a
distribution of the same shape.

## Fixed-ion electron branch

A useful early case evolves electron kinetics against prescribed ions.

Distinguish:

1. a fixed neutralizing charge/current background;
2. a prescribed physical ion distribution used in electron-ion collisions.

They are not interchangeable.

For active electrons,

\[
\rho_e=-e\int f_e\,d^3p,\qquad
j_e=-e\int v_e f_e\,d^3p,
\]

with the selected field model adding the fixed-ion contribution consistently.
Momentum/energy transferred to fixed ions is a reservoir exchange, not an
internal invariant that must vanish.

## Gyrofluid and gyromoments

Gyrofluid/Landau-fluid and gyromoment species are explicit priority branches.
A low-order state should, where possible, arise from the same parent
gyrokinetic moment hierarchy that supports later higher-order calculations.

Changing retained moments is a model change until a controlled tail/closure
remainder exists. The same species can be compared at several moment counts,
then against GK.

## Bulk MHD and kinetic-MHD

A bulk state aggregates selected populations into mass, momentum, pressure and
fields. Record population membership explicitly.

An active kinetic species may couple through a derived current- or
pressure-coupling formulation. Never add kinetic pressure to a bulk pressure
that already contains that population.

Spatial MHD-to-kinetic coupling requires a matched/interface formulation with
conservative transfer and a bounded interface/model defect.

## Relativistic scope

Relativistic physical momentum is the FK reference coordinate. This does not
automatically make every DK/GK/fluid reduction relativistic. Each reduced model
must state its own relativistic ordering and Hamiltonian before such a claim.

Mass is positive and charge may be signed or zero. Pair plasmas and unusual
species therefore need no hard-coded taxonomy, only compatible interactions.

## Reduction evidence

Each model first receives a numerical certificate for its own equations.
Comparing two certified models gives a certified difference of their predicted
observables, not a theorem that either is closer to physical reality.

A rigorous reduction bridge additionally needs:

- ordering bounds;
- regularity;
- lift/projection control;
- initial/boundary layer control where relevant;
- field and collision reduction error;
- stability over the domain/time of interest.

Only then may numerical and model errors be combined into a guaranteed total
observable bound.

## Adaptive species fidelity

Long term, fidelity may vary by population and region. The target is the
cheapest supported model whose total proven observable error is below the
request.

Before that theorem infrastructure exists, adaptive model selection is an
experimental heuristic and must be labeled accordingly.

State transfer during a model switch requires conservative projection and a
measured/bounded switching defect. It is not merely changing an enum.

Chris&AI
