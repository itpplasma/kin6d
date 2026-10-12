# Equilibria, imported candidates and TIAGO inference

Status: canonical architecture policy, 2026-10-05. Most capabilities described
here are planned; the implemented baseline remains the 1X1V laboratory.

## One stationary-state concept

KIN6D defines an equilibrium or stationary state through the selected physical
equations,

\[
R(U;\theta)=0.
\]

"Equilibrium" is not synonymous with EQDSK, VMEC, nested flux surfaces or any
other file/coordinate representation. Open systems with sources, collisions,
walls or transport may have stationary states that are not thermodynamic or
ideal-MHD equilibria, but they use the same residual/implicit-derivative
architecture.

The map from inputs to stationary state need not be unique. If two valid
stationary solutions are found, they are two equilibria. Keep them separate.
Each result records its initial state, external data, continuation/evolution
history and solution identifier when those facts matter. KIN6D does not need a
general probabilistic branch object.

## Imported equilibria are candidates

Imported data serve four roles: immediate usable field/state data, an initial
guess for a native stationary solve, an independent cross-code comparator, and
a source whose residual/interpolation/convention uncertainty can be checked.

Import must never silently promote the source code's coordinate assumptions into
KIN6D's physical model. Normalize physical units, orientation, flux/toroidal
angle conventions, field periods and vector-component conventions at the
adapter boundary. Prefer libneo when it already owns a mature reader/converter.

The libneo boundary is intentionally replaceable, but its generic implementation
is improved on demand. When a KIN6D use case exposes a reusable defect or
missing capability in parsing, convention normalization, interpolation, field
evaluation, spatial derivatives or representation-level error diagnostics, fix
and test it in libneo rather than creating a KIN6D fork. If an old evaluator is
too entangled to make differentiable/certifiable economically, use libneo only
to normalize the source data and construct a native KIN6D representation for
the differentiable/certified computation.

## Axisymmetric MHD

KIN6D should own a native Grad--Shafranov path for fixed and free boundary
problems. For an imported EQDSK/gEQDSK: read and normalize it, reconstruct the
physical field/profile data, evaluate divergence/Ampere/Grad--Shafranov and
boundary residuals, recover the represented inputs, re-solve natively, and
compare/certify imported and native observables.

A free-boundary equilibrium cannot in general be reconstructed from the EQDSK
alone. External coil geometry/currents or an equivalent vacuum/boundary
specification are additional physical inputs.

## Three-dimensional MHD without a flux-surface ontology

Nested surfaces remain useful when they exist. A VMEC/DESC/GVEC-like state may
be the cheapest representation of a smooth nested equilibrium and can be used
directly, as an initializer, or as a preconditioner. It is not the general
KIN6D equilibrium definition.

The general physical-space formulation must remain valid through islands,
stochastic/ergodic regions, separatrices and wall contact. Reference families
with different assumptions include VMEC/DESC/GVEC for nested equilibria,
HINT/PIES-like nonnested full-pressure methods, and SPEC/SPECTRE for
multi-region relaxed MHD.

Dropping nested surfaces does not drop ideal-MHD force balance. For smooth
static ideal MHD,

\[
\nabla p=J\times B,\qquad B\cdot\nabla p=0.
\]

An ergodic magnetic component therefore cannot support an arbitrary smooth
pressure gradient. Pressure flattening, ideal barriers/current sheets, relaxed
regions or nonideal evolution are physical/model choices.

Smooth non-symmetric equilibria with nested surfaces do exist. Grad doubted
that smooth toroidal equilibria with non-constant pressure exist without a
continuous symmetry. Two 2026 preprints give counterexamples:

- [Landreman, arXiv:2609.26742](https://arxiv.org/abs/2609.26742) gives
  closed-form analytic families. They have two field periods, are
  stellarator-symmetric, and have exact nested surfaces, with ∇p vanishing
  only on the axis. One family has uniform ι=2, the other sheared ι. This
  refutes the formulation in Cardona et al. 2025, Conj. 2.5.
- [Gómez-Serrano, Liehr and Taylor, arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
  give smooth solutions with cyclic symmetry, certified in Lean 4.

Existence is therefore not an obstruction in principle. Whether a given
boundary and set of profiles admits a smooth nested solution, and what happens
on rational surfaces with sheared ι, stay case-specific questions. The
analytic families serve as exact 3D benchmarks ([benchmarks](benchmarks.md)).

## Topology and history

History dependence is represented by the selected equations and initial/
continuation data, not by a universal topology flag. Ideal evolution preserves
its frozen-in constraints; MRxMHD-like models retain a finite declared invariant
set; resistive MHD/Hall/XMHD can reconnect according to their equations.

Continuation is both a solver and scientific provenance. Increasing and
decreasing a control parameter may reach different stationary states; keep both
if they are distinct.

## Kinetic and open-system stationary states

The same contract extends from MHD through DK/GK/FK. Collisionless kinetic
equilibria are highly nonunique. With collisions, sources, transport, neutrals
and walls, the scientifically relevant object is often a balanced stationary
solution rather than a collisionless equilibrium.

## TIAGO owns reconstruction and inference

KIN6D provides the differentiable forward model: stationary/evolution solve,
physical state and observables, JVP/VJP and implicit tangent/adjoint actions,
numerical certification, and model-reduction/error information.

TIAGO owns measurements and diagnostic likelihoods, correlated noise and
calibration models, nuisance parameters and priors, MAP/least-squares
reconstruction, posterior/Laplace covariance, and comparison of candidate
equilibria. TIAGO's existing VMEC++ exact-adjoint reconstruction remains an
independent provider/reference.

## Several equilibria and local Laplace UQ

If the same data are compatible with equilibria A and B, run two local inverse
problems:

\[
p(\theta\mid d,A)\approx N(\hat\theta_A,C_A),\qquad
p(\theta\mid d,B)\approx N(\hat\theta_B,C_B).
\]

If useful, TIAGO may attach evidence/prior-derived weights and report a finite
mixture of these local Laplace approximations. The physical states remain
separate; do not average A and B into a state that may satisfy neither set of
equations. More expensive general sampling is deferred until a concrete
within-solution posterior requires it.

## Orthogonal uncertainty ledger

Keep separate: measurement/calibration uncertainty; posterior/reconstruction
uncertainty conditional on a model/equilibrium; diagnostic-model uncertainty;
KIN6D numerical error; KIN6D physical model/reduction error; stochastic
sampling uncertainty; and the existence of several distinct stationary
solutions. Do not combine deterministic intervals, posterior standard
deviations and model-error bounds by root-sum-square without a justified joint
probability model.

Chris&AI
