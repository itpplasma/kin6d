# Collision models and discrete structure

Status: required scientific feature set; all implementations are planned.
This extends the [design](../DESIGN.md) and the [certification contract](certification.md).
Every project-owned numerical leaf follows the [fortsym policy](tooling.md#kernel-generation-policy).

## Independent choices

Select the physical collision model, its nonlinear/linearized action, the
discrete conservation/entropy treatment, and the certificate scope separately.
Each choice changes a stated mathematical or numerical contract. A conservative
discrete update is not a continuous-solution certificate; an accurate enclosure
for a simple model is not a bound on omitted physical effects.

The gyrokinetic GENE-X code is a useful precedent: its
[collision module](https://phoenix-public.pages.mpcdf.de/genex/module/op_coll_m.html)
documents a common collision base, BGK, LBD, Lorentz, and FPL types, with
grid/spectral implementations. Borrow the narrow operator-hierarchy idea.
Its gyrokinetic equations and discretizations have a different scope from
the full-orbit relativistic target here. The complete feature combination's
novelty remains a prior-art question.

## Operator family table

`N_s` means each species' particle number, `P` total kinetic momentum, and
`K` total kinetic energy of an isolated collision system. Entropy statements
refer to positive admissible distributions and the exact specified closure.
Sources, prescribed reservoirs, momentum boundaries, and fields can change
the balances. Application cost is qualitative; certificate difficulty also
depends on the domain, norm, and stability regime.

[Structural contracts](adjoints-structure.md) declare weighted adjoints,
collision sign, invariant complement, entropy and Onsager–Casimir scope for
each construction. Removing a nullspace does not alone prove a uniform
Coulomb whole-velocity spectral gap. Fixed-target Krook's derivative can be
self-adjoint/dissipative without physical pair conservation; nonlinear models
need their actual linearization or bracket, not blanket adjoint flags.

| Model or action | Purpose and physics | Conditional structure | Omitted physics or restriction | Work and certification hurdle |
|---|---|---|---|---|
| `none` | Collisionless reference | Zero collisional moments | All collisions | No collision work; kinetic/field certification can still be hard |
| `krook` | Fixed-target relaxation | Target-dependent balances | Coulomb rates and reaction | Local; source/reservoir theorem |
| `bgk` | Moment-matched relaxation | `N_s,P,K`, entropy for the chosen conservative construction | Detailed drag/diffusion | Moments and local action; multispecies target/entropy proof |
| `lorentz` | Angular scattering in a prescribed heavy background | Test-species number and speed-dependent energy | Slowing down and field-particle reaction; momentum reservoir | Sparse angular action; frequency/domain bounds |
| `lbd` | Lenard–Bernstein/Dougherty drag and isotropic diffusion | Pair conservation and entropy for a specified valid closure | Coulomb high-energy rates and anisotropic diffusion | Local FP; moment coefficients and discrete correction bounds |
| `hirshman_sigmar` / `sugama` | Reduced linearized moment/friction model | Conservation and adjoint/entropy properties under the selected formulation | Exact nonlinear Coulomb kernel | Moment response; model-specific linear stability |
| `landau_test` | Nonrelativistic differential/test-particle action | Number; exchange with fixed reservoir | Field-particle reaction and nonlinear correction | Fixed-background coefficients; continuous transport/diffusion bounds |
| `landau_field_only` | Integral/field-particle diagnostic | Meaningful in its paired decomposition | Not a standalone physical evolution model | Perturbation potentials; diagnostic identity bounds |
| `landau_linear` | Full first-order Coulomb action | Linearized invariants and entropy dissipation around an appropriate equilibrium | Quadratic perturbation term | Potential response; constrained linear stability |
| `landau_nonlinear` | Full nonrelativistic small-angle closure | `N_s,P,K`, H-theorem under its assumptions | Relativity, detailed dressing, hard close encounters | Nonlinear potentials; derivatives and infinite complements |
| `bb_test` / `bb_field_only` | Relativistic test/field components | Reservoir or diagnostic structure as above | Other decomposition components | Relativistic coefficient/potential bounds |
| `bb_linear` | First-order Beliaev–Budker action | Relativistic linearized invariants under stated equilibrium assumptions | Quadratic perturbation term | Relativistic constrained stability |
| `bb_nonlinear` | Relativistic small-angle production closure | Relativistic `N_s,P,K`, H-theorem | Detailed dielectric dressing and hard collisions | Nonlinear relativistic potentials and complements |
| `balescu_lenard` | Dynamically screened weak correlations | Invariants/H-theorem for the specified admissible model | Strong close collisions; not automatically relativistic | Dielectric/resonance work and stability bounds |
| Matched screened/hard | Defined overlap-subtracted closure | Must be derived and proved | Higher correlations still closed | Screening, scattering, overlap, and closure analysis |
| Mean-force closure | Effective equilibrium correlations | Model-specific equilibrium and balances | Nonequilibrium higher correlations approximated | Pair data and applicability estimates |
| Truncated BBGKY | Explicit low-order correlations in a small reference problem | Depends on truncation and boundary model | Higher hierarchy terms | High-dimensional reference and remainder theorem |

Krook means `C=-nu(f-F_target)` with a prescribed target. Matching a target's
normalization does not conserve number for arbitrary `f`; the relevant moment
must agree or a reservoir supplies the difference. BGK recomputes a target
from moments. Relativistic BGK needs a compatible relativistic target and
invariant construction; a nonrelativistic Maxwellian target is not silently
reused in that branch.

For Lorentz scattering choose and document the normalization, for example
`C=(nu_D(p)/2) Delta_Omega f`. Then its degree-`l` eigenvalue is
`-nu_D(p) l(l+1)/2`. A bounded manufactured frequency provides a clean early
oracle. Physical low-momentum behavior and a `v^-3` asymptote need their own
domain treatment. Angular diffusion conserves radial energy but exchanges
momentum with the prescribed background.

The [improved Sugama model](https://arxiv.org/abs/1906.07427) is a reference
for reduced linearized friction-flow relations and conservation. Choose its
precise equilibrium, species, and collisionality assumptions before adopting
an adjointness or entropy claim. Keep reduced models as later named choices,
rather than prerequisites for reaching nonlinear Landau.

## Shared pair primitive

For Landau or Beliaev–Budker with fixed masses, charges, pair coefficients,
and Coulomb-logarithm prescription, define a bilinear primitive
`Q_ab(g_a,k_b)`: action on species `a`, with potentials generated by species
`b`. Separate species arguments and their physical scales in its contract.

For `f_a=F_a+h_a` and `f_b=F_b+h_b`, bilinearity gives exactly

$$
Q_{ab}(f_a,f_b)=Q_{ab}(F_a,F_b)+Q_{ab}(h_a,F_b)
+Q_{ab}(F_a,h_b)+Q_{ab}(h_a,h_b).
$$

Define test-particle `Q_ab(h_a,F_b)`, field-particle `Q_ab(F_a,h_b)`, linear
action as their sum, and nonlinear correction `Q_ab(h_a,h_b)`. Sum over
partners for the species action. The full nonlinear operator uses
`sum_b Q_ab(f_a,f_b)`; it does not pass one species array as both arguments
for a cross-species pair.

The background term vanishes only for a reference that is actually a
collisional equilibrium. Unequal-temperature/drift species generally exchange
momentum and energy, even when each is individually Maxwellian. Around a
nonequilibrium reference retain `Q_ab(F_a,F_b)`. Spatially varying equilibrium
backgrounds can also have nonzero orbit/field/source residuals; the perturbation
form must retain those terms.

The test action is differential in its first argument because the background
potentials are fixed. The field action depends nonlocally on the perturbation
through its potentials, even if those are obtained by differential radial
solves. This is the differential/integral distinction to expose for NEO-2-like
comparisons in an agreed overlapping regime. `field_only` is diagnostic and
cannot be selected as an accepted physical collision closure.

If a Coulomb logarithm or other pair coefficient depends on `f`, the complete
operator is no longer bilinear as written: include coefficient variations in
the Jacobian and certification. Balescu–Lenard's self-consistent dielectric
and moment-dependent LBD/BGK coefficients likewise do not inherit this
decomposition without additional terms.

For fixed coefficients, verify

$$
DC(F)h=\sum_b[Q_{ab}(h_a,F_b)+Q_{ab}(F_a,h_b)].
$$

Compare bilinear composition, a supported fortad JVP path through fortsym
emission, and a finite-difference diagnostic. Agreement detects wiring or
derivative errors; it does not independently prove that the common formula
is the intended continuous operator. Retain exact harmonic/equilibrium and
pair-exchange oracles as well. An entropy theorem for linearized perturbations
is distinct from a full nonlinear H-theorem for arbitrary signed `h`.

## Discrete conservation treatment

Ulbl, Michels, and Jenko's
[GENE-X LBD work](https://onlinelibrary.wiley.com/doi/full/10.1002/ctpp.202100180)
uses a number-conserving finite-volume form and modifies mixing quantities
to enforce discrete momentum and energy conservation to arithmetic precision.
Its construction motivates a separate discrete-structure choice; its
gyro-averaged formula is not copied into full-orbit momentum coordinates.

Prefer a flux/metric-bracket discretization with the right nullspace by
construction. For a permitted correction, measure the raw invariant defects
`delta I_k=<psi_k,C_h>` with the same discrete moment functionals used by the
solver. Choose correction functions `chi_j`, form
`B_kj=<psi_k,chi_j>`, and solve

$$
B\lambda=-\delta I,\qquad
\widetilde C_h=C_h+\sum_j\lambda_j\chi_j.
$$

Define the norm/metric and admissible subspace before calling a correction
minimal. Prove the small constrained system is solvable with sufficient
conditioning; unequal species grids require paired correction and scaling.
Species particle invariants, total vector momentum `p`, and energy use the
physical quadrature measure. Energy is `p²/(2m_s)` for nonrelativistic models
and `(gamma_s-1)m_s c²` for relativistic models.

With a complete harmonic basis, number and energy occupy degree zero and
momentum degree one. A correction may be supported there if radial and species
constraints permit it. This is a small solve, not an excuse to set each
species' exchange moments to zero: only the paired total must vanish.

A test-particle or Lorentz reservoir model intentionally exchanges momentum
or energy. Enforcing closed-system invariants on it changes the physics.
Metadata must prevent inappropriate conservation modes. Conservative
projection does not automatically preserve positivity, entropy dissipation,
adjointness, or consistency. Prove each claimed property for the actual
corrected operator and treat a failed property as a separate capability gap.

The certifier bounds the correction relative to the intended continuous
operator, including its derivative where needed. A conservative corrected
discrete root is not itself a root of the continuous equation. Include the
projection/correction defect in the residual and stability budget.
Floating moment cancellation to machine precision is a diagnostic; rigorous
arithmetic encloses the remaining balance error.

## Narrow interface and metadata

Expose model selection and batch `apply`; add `jvp` and certificate-bound
actions when the next capability needs them. Select the concrete implementation
outside hot loops. A collision kernel knows pair data, basis/grid information,
and local distribution inputs, but not Newton, MPI decomposition, temporal
driver, plotting, or file formats. Avoid independent nonlinear/linear/test/
field implementations of the same pair physics.

At the first executable collision stage, introduce a small Fortran descriptor
as the single source for executable capability reporting and generated operator
documentation. The table above is a design specification until then. Record:

- Model identity/version, pair coefficient prescription, reference state,
  action selection, kinematics, and physical applicability.
- Expected invariants and entropy/adjoint claims with their assumptions and
  theorem identities; reservoir exchanges and diagnostic-only restrictions.
- Discrete scheme/correction selection, its proven properties, and measured
  balance defects separately from continuum identities.
- Generated-kernel provenance, CPU/CUDA support, and independent test status.
- Certificate support for a named domain/norm and, for the actual run, the
  verified proof level, radius/tube, observables, and unmet obligations.

A backend-level `certified=true` flag is insufficient. A supported certifier
can fail for particular data, resolution, reference, or stability margin.
Keep physics fidelity, discrete structure, and run-specific proof strength as
separate fields; do not print a check mark for an unimplemented or failed case.

## Mixed species and reservoirs

Species representation and pair physics are separate choices with a
compatibility gate. FK–DK/GK pairs require the selected coordinate transforms,
phase-space measures, spatial/velocity effects of the reduction, and field
action. A lift/parent-pair/projection construction is a useful starting
expression, not a proof that arbitrary projections preserve the coupled
physics. Pair balances use the model's invariant functionals and any field
contributions. Conservative correction cannot repair an inconsistent model.
[Hybrid coupling](species-models.md#mixed-collision-pairs)

Fixed finite-temperature ions provide a prescribed collision distribution;
a charge-only neutralizing background does not. Electron self-collisions
retain their closed invariants, while electron–reservoir momentum and energy
transfer is reported explicitly. Do not project that transfer to zero or
call it a failed all-species conservation test. Promoting the ions to evolved
response retains the reservoir limit as a plotted comparison.

## Certified model differences

For matching inputs and observable definitions, compute independent certified
intervals `I_i` and report differences `I_j-I_i`. This supports comparisons of
Lorentz, TP, full linearized, nonlinear Landau, and BB predictions, and later
screening/matching choices. A branch correspondence theorem is needed to
claim a tracked solution family rather than a difference of two solutions.

Changing only TP versus TP+FP isolates that model choice if all other
definitions are held fixed. Nonlinear versus linearized comparison must
retain the same reference, forcing, and perturbation conventions. Landau
versus BB can change kinematics as well as collisions; name both changes.
Ordinary nonrelativistic Balescu–Lenard versus BB is not automatically a pure
screening comparison. These are certified model differences, not bounds on
real-plasma error without an additional model-reduction theorem.

Chris&AI
