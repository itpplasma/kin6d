# Background, drive and induced response

Status: configuration semantics. Bulk/background evolution and induced
perturbation feedback are independent choices, subject to model consistency.

Define fields as B=B0+delta B_ext+delta B_plasma and similarly for E.
The background contains its already-declared plasma and external components;
the new drive and induced response must not recount those sources. Each
component records whether it is prescribed, solved, enclosed or time-varying.
External electric/magnetic drive must satisfy the chosen Maxwell/boundary
relations, including induction when applicable.

| Background | Drive | Induced response | Problem |
|---|---|---|---|
| Fixed | Prescribed | None | Test-field experiment |
| Fixed | Prescribed | Self-consistent | Forced linear or stated perturbative response |
| Fixed | None | Self-consistent | Eigenmode/stability with frequency and normalization |
| Evolving | Optional external | Self-consistent | Nonlinear trajectory |
| Solved simultaneously | Optional external | Self-consistent | Coupled steady equilibrium/state |

For an RMP case the orbit experiences delta B_ext+delta B_plasma when the
selected model includes both. Freezing B0 and bulk profiles does not remove
the self-consistent perturbation solve. Linearization gives the constrained
operator L[U0] delta U=S_ext, with background residual/uncertainty accounted
for. It is not the same problem as solving the full nonlinear equilibrium.

## Population roles

Separate representation (FK/GK/DK/fluid/bulk) from population coupling role.
A solved FK minority can be passive; a reduced DK population can back-react.
Describe whether it is solved/evolved, contributes charge/current/pressure
or exchange feedback, and feels the induced response. Fixed sources can
contribute prescribed charge/current without an evolved distribution.

Conceptual flags `evolve`, `sources_fields`, `feels_response` express these
questions, but flags alone are not a physical model. An algebraic adiabatic
population can respond and source fields without an evolution equation;
collision/reaction momentum and energy feedback need separate declarations.
A frozen distribution that 'feels' fields requires an explicit prescribed
functional response or is just a diagnostic interpretation. It does not
silently acquire dynamical feedback.

| Population | Typical semantics | Required balance |
|---|---|---|
| Active kinetic | Solved, sources response, feels response | Coupled charge/energy/momentum |
| Test kinetic | Solved, no back-reaction, feels selected fields | One-way energy/exchange with external reservoir |
| Fixed reservoir | Not evolved, prescribed sources/collision background | Declared continuity and reservoir transfer |
| Algebraic closure | Instantaneous selected response | Consistent field/closure action |

For thermal MHD plus energetic kinetic species, declare which populations
are counted in the bulk and whether the kinetic response is a diagnostic or
feeds back. Species collisions, atomic reactions and surfaces follow the same
source/exchange accounting. A passive diagnostic must not secretly source
the fields or get a closed-system total-energy certificate.

Each run records background, drive, induced response, population roles, field
equations, constraints and observable pullbacks. Comparisons can hold the
background fixed while changing response feedback, then separately allow
bulk relaxation. Preserve test-field and frozen-background limits in the
atlas and their plots. The certificate names the actual partially or fully
self-consistent mathematical problem; 'self-consistent' alone is insufficient.

Chris&AI
