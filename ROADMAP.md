# MILCA roadmap — active work first

Status: 2026-11-12. **This is the only active priority list.** Issues not listed here are
backlog, not instructions to implement. Physics specifications and long-range
scope remain in [DESIGN.md](DESIGN.md), [the capability catalogue](docs/capability-catalogue.md),
[the model graph](docs/model-graph.md), and the
[previous detailed milestone plan](https://github.com/itpplasma/milca/blob/94caec95e474bd23640278c5e03b50204943f369/ROADMAP.md).

## Mission and division of work

MILCA is the group's fast-moving, independently verifiable **physics and numerics
laboratory**, not a requirement that everyone use one executable. Established
specialized codes (GORILLA, SIMPLE, NEO-2, MARS, BEAMS3D, VMEC++, etc.) are
independent products, competitors, consumers and scientific comparators.
Improve them through focused, independently evidenced upstream PRs; their
maintainers decide adoption and review. Do not transfer MILCA implementations
into comparator codes just to obtain agreement.

The scientific target is **time and cost to a named observable at matched
physical model, accuracy and uncertainty**, not amount of code, number of
features or raw runtime alone. Keep the group's human effort on physics.

## Current work (ordered by priority)

| Priority | Work package and next gate | Tracked work |
|---|---|---|
| **P0** | **Restore trustworthy gates:** finish the adaptive NTV mesh error/size-limit defects (refinement and band joins are capped since #372 and #382; start meshes are not yet); investigate open nightly and equilibrium timeout failures. Do not accept an accuracy or resource claim when its gate was bypassed. | [#363](https://github.com/itpplasma/milca/issues/363), [#362](https://github.com/itpplasma/milca/issues/362), [#259](https://github.com/itpplasma/milca/issues/259), [#145](https://github.com/itpplasma/milca/issues/145) |
| **P1a** | **NTV/neoclassical physics:** converge AUG 30835 and ITER TC24 observables with explicit model/convention match, resolution, error and cost; focus on unresolved resonant-ion discrepancies. Finish only the currently exercised collision/source work, not blanket feature parity: the energy-coupled kernel with finite orbit width on the hp orbit-class mesh (stacked branches `ntv/a01` to `ntv/a05`) and the torque density in real space. | [#233](https://github.com/itpplasma/milca/issues/233), [#197](https://github.com/itpplasma/milca/issues/197), [#198](https://github.com/itpplasma/milca/issues/198), [#369](https://github.com/itpplasma/milca/issues/369) |
| **P1b** | **3D equilibrium and exterior:** the coil-driven BEM exterior agrees with independent Biot–Savart and virtual-casing references (#350, #352); next is the finite-pressure coil-driven case with a displaced seed. Preserve NN/island and Beltrami convergence gates. Avoid advancing derived applications on an unvalidated field. | [#385](https://github.com/itpplasma/milca/pull/385), [#391](https://github.com/itpplasma/milca/pull/391) |
| **P1c** | **Drift-wave turbulence for the resonant-layer response:** spectra, correlation scales and transport coefficients at one rational surface ([drift-turbulence.md](docs/drift-turbulence.md)). Stages L and HW with the saturated-flux benchmark and the homogeneous drift-Alfven step of stage DA are in; next is stage DA beyond the homogeneous slab. | [#371](https://github.com/itpplasma/milca/pull/371), [#379](https://github.com/itpplasma/milca/pull/379), [#383](https://github.com/itpplasma/milca/pull/383) |

These are **existing active campaigns**, not a proposal to restart them. Each
campaign closes a specific observable-level verification gate before adding
adjacent scope. Document unresolved discrepancies explicitly (e.g. the AUG
resonant-ion per-surface differences in
[drift-kinetic lessons](docs/lessons/drift-kinetic-transport.md)).

## Next work (start small, after a P0 gate is stable)

1. **P1d — First independent Monte Carlo/continuum check:**
   [#366](https://github.com/itpplasma/milca/issues/366). Use a common
   Lorentz or Ornstein–Uhlenbeck collision model, an analytic limit, a MILCA
   continuum implementation and an *independent* GORILLA_APPLETS stochastic
   implementation. Converge mesh, timestep and sample count; compare the same
   weak observable with deterministic and statistical errors kept separate.
   Do not assume test-particle trajectories implement a full linearized
   multispecies Landau operator.
2. **P2 — Evidence-to-fix automation:** make a failed analytic, cross-method or
   cross-code case produce a minimal reproducer, discrepancy classification,
   evidence record, targeted test and optionally a proposed upstream PR.
   Reuse the [benchmark atlas](docs/benchmarks.md), fo and Fortyard's
   [testing-at-scale design](https://github.com/lazy-fortran/fortyard/blob/main/doc/testing-at-scale.md);
   avoid a new universal verification framework.
3. **P2 — Expand comparisons only where physics needs them:** NTV regime limits
   [#232](https://github.com/itpplasma/milca/issues/232), established transport
   codes [#227](https://github.com/itpplasma/milca/issues/227)–[#235](https://github.com/itpplasma/milca/issues/235),
   collision hierarchy [#188](https://github.com/itpplasma/milca/issues/188),
   and native-field orbit coupling [#331](https://github.com/itpplasma/milca/issues/331).
   Promote each independently validated slice on scientific demand.

## Verification and integration rules

1. **Different evidence, not duplicate implementations:** challenge a result by
   analytic/asymptotic limits, conservation/adjoint/reciprocity identities,
   independent symbolic derivation, discretization, Monte Carlo versus
   continuum, established codes and ultimately observations as applicable.
   Shared equations, generated kernels, inputs and historical derivations are
   *correlated evidence*. Record exactly which physical models intersect.
2. **A test must be able to fail:** use failing-before/passing-after
   reproducers, convergence and resolution sweeps, mutation/counterexample
   searches, boundary and sign checks. A golden record or two solvers using
   the same symbolic expression is not independent truth.
3. **Benchmark observable-first:** pin model, conventions, solver revision,
   input, geometry, tolerance, random seed, hardware and actual result status.
   Report continuous error bounds, empirical discretization estimates,
   statistical uncertainty and model mismatch separately. A disagreement is
   allowed to remain *unresolved*.
4. **AI investigates continuously:** generate tests, minimize failures,
   check paper formulas and errata, compare independent CAS systems (including
   FortSym/SymPy/Mathics), diagnose regressions and propose upstream repairs.
   No agent may alter production code *and* silently redefine its authoritative
   benchmark or mathematical acceptance condition.
5. **Review scientific risk, not lines or AI authorship:** automatically
   integrate low-risk, well-tested mechanical changes; audit samples. Demand
   independently checkable evidence for numerical changes. Escalate new physics,
   unresolved high-impact result changes, reference/oracle changes and
   publication claims. External maintainers always choose their own process.
6. **Continuous verification, bounded cost:** affected tests immediately,
   cross-method/downstream tests after integration, expensive reference and
   device tests in a resource-limited background queue, and complete
   evidence for published claims. A passing cached result is valid only for its
   correctly identified full input/toolchain closure. Track stale and failed
   tests; repair the test harness rather than silently relaxing the gate.

## Admission and delivery

Implement the **smallest complete scientific slice** with a named question,
an explicit mathematical model, an independent oracle or clearly scoped
limiting case, convergence/error status, and time/memory at matched accuracy.
Keep derivations, reference inputs, reproducible commands and provenance;
report missing certificates honestly.

Generic numerics belong in FortNum/FortSparse/FortFEM/FortSym/FortIO and
interchange in libneo as defined in [AGENTS.md](AGENTS.md). Maintain independent
specialized reference implementations; do not require established codes to
depend on MILCA.

## Future, not active

All other roadmap families (full NEO-2 parity, comprehensive 3D kinetic
response, gyrokinetics, edge/neutrals, reactor integration, generalized orbit
framework, turbulence beyond the drift-turbulence branch above, radiation, microscopic/quantum extensions and broad
multifidelity automation) remain **future capability choices**, even where
pieces already run. Their requirements and milestone detail remain in the
[previous staged plan](https://github.com/itpplasma/milca/blob/94caec95e474bd23640278c5e03b50204943f369/ROADMAP.md) and domain
specifications. Start one only when a live physics question and available
verification evidence justify moving it into the current table.

Chris&AI
