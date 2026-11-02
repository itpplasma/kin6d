# Relaxation comparators and what kin6d adopts

Status: design note (2026-10-28). Nothing on this page is implemented yet
unless a module or test is named.

This page relates the published three-dimensional relaxation codes without
nested surfaces to kin6d's equilibrium lanes and lists the algorithms kin6d
adopts from them. Sources and the reading state are in
[references](references.md#three-dimensional-equilibrium-comparators);
benchmark use follows the [source policy](literature-policy.md#benchmark-sources).

## Relation to kin6d

**Frozen-in map and Eulerian relaxation.** kin6d's nested 3D lane minimizes
the energy over a map from a reference domain ([ideal map](ideal-map-3d.md)):
field lines move with the map, so the topology is preserved exactly. MRX
evolves B on a fixed mesh by δB = curl(v × B), the linearized form of the same
frozen-in transport; topology changes only by grid-scale reconnection. Both
discretize the same second variation (the δ²W energy principle), so the map
Hessian and MRX's Newton operator can be cross-checked on one nested case. The
map tangles as islands grow, which is why non-nested states use the
[FEEC lane](equilibrium-nonnested.md).

**Pressure closure.** The codes differ mainly in how pressure is determined:

| Closure | Codes | Effect in islands |
|---|---|---|
| adiabatic transport with the displacement | kin6d map lane, SIESTA | pressure follows compression |
| Lagrange multiplier for divergence-free flow | MRX | final profile depends on the initial state |
| parallel diffusion towards B·∇p = 0 | NTEC, HINT | flattening at a diffusion-set rate |
| stepped pressure, Taylor state per region | SPEC, SPECTRE | MRxMHD |

These are different models, not numerical variants. kin6d treats the closure
as an explicit model choice of each equilibrium, and comparisons hold only at
the same closure or in limits where the closures agree (vacuum, low β with
nested surfaces, small islands).

**Metriplectic relaxation.** kin6d's [metriplectic relaxation](metriplectic-relaxation.md)
dissipates energy and conserves Casimirs such as helicity by construction.
MRX obtains the same structure with an H⁻¹-smoothed gradient flow, a rank-one
correction that restores discrete helicity, and a resistive phase guided by a
prescribed current. MRX's result that seeded and unseeded states relax to one
equilibrium after the resistive phase is a case kin6d's metriplectic and
resistive lanes should reproduce.

**Islands.** kin6d represents prescribed islands through its island closure
provider ([equilibria with islands](island-equilibrium.md)). MRX seeds resonant
perturbations at rational surfaces and chooses their amplitudes by minimizing
a quadratic energy model, a* = −G⁻¹∇E. That criterion is a principled source
of initial states and closure parameters.

**NTEC.** C¹ Bell triangles times Fourier modes are not structure-preserving
(div B needs a cleaning term), and the hyperbolic parallel diffusion relaxes
towards B·∇p = 0 rather than imposing it. kin6d keeps FEEC and solves the
constraint directly. NTEC's published SIESTA and HINT comparisons are useful
reference data for the resistive lane.

**SIESTA, SPEC/SPECTRE, PIES, HINT.** SIESTA assembles the Hessian by finite
differences and factors it block-tridiagonally, which is close to kin6d's
current factored non-nested Newton. SPEC and SPECTRE correspond to the planned
native MRxMHD lane; their Taylor states are kin6d's
[Beltrami closures](equilibrium-beltrami.md). PIES and HINT are Picard and
explicit relaxation methods and serve only as published references.

## Algorithms to adopt

Each is reimplemented in kin6d's own stack (FortFEM operators, FortSym-generated
kernels), in order of expected payoff:

1. **Matrix-free Newton–Krylov for the non-nested FEEC lane.** Newton-MR
   (Liu and Roosta, ref. 100 of the MRX paper) with MINRES on the
   second variation, applied through sum-factorized FortFEM actions, with
   preconditioners built from metric lumping (Kronecker products and sums of
   1-D matrices, a dense axis block) and geometric multigrid on the div-div
   stiffness. This replaces the assembled Hessian and FortSparse factorization
   in the [harness](equilibrium-harness.md) as the default path.
2. **Parallel-only Levenberg–Marquardt term.** Damping on the v ∥ B component
   only, which is the near-null mode of the second variation, in place of a
   diagonal damping of all unknowns.
3. **Harmonic deflation and the potential-velocity force.** Compute the
   harmonic forms per geometry and project them out of every solve; evaluate
   the force through curl (curl curl)⁻¹ curl(J × B).
4. **Energy-criterion island seeding** as an initializer for island cases.
5. **Residual floor as a diagnostic.** Current sheets at rational surfaces set
   a force-residual floor that does not decrease under h or p refinement;
   acceptance gates report it instead of chasing it with smooth spaces.

Not adopted: explicit damped iteration, Bell elements, hyperbolic pressure
relaxation as a production method, and accelerator-specific frameworks.

## Benchmarks

- Analytical and limiting cases: [Landreman and Issan nested references](nested-oracles.md),
  Cerfon–Freidberg and Solov'ev, MRX's manufactured vacuum field, Taylor
  states, the pendulum island-width estimate, vacuum and low-β limits.
- Open source, run at a pinned revision: MRX (`ncsx_relax_seed_relax` in the
  [atlas](benchmarks.md)), SIESTA, SPECTRE, VMEC++, DESC, GVEC.
- Published only: NTEC and HINT figures and tables, recorded with their model
  assumptions and caveats.

## Implementation status

- Item 2: `parallel_damping` of the ideal map Newton (non-nested solve),
  [ideal map](ideal-map-3d.md), tests `ideal_map_parallel_mass` and
  `ideal_map_parallel_damping_p2`. The NN2 coupled Newton has A, not a
  displacement, as unknown: v || B maps to dA = v x B = 0, so that system
  has no parallel near-null mode to damp.
