# Software-engineering and clean-architecture policy

Status: canonical repository policy, 2026-10-05.

This document applies ordinary software-engineering discipline to an unusually
broad research programme. The breadth of KIN6D is a reason to defer detail, not
a reason to build speculative infrastructure.

Clean Architecture here means **dependency direction, separation of policy from
details, and delaying irreversible choices until a real use case needs them**.
It does not mean reproducing a ceremonial set of concentric folders, abstract
base classes or plugin registries.

## KISS — simplest thing that serves the current slice

For the active vertical slice, prefer the simplest complete design that
preserves the scientific contract.

- A concrete procedure/type is preferred to a framework when one consumer is
  enough.
- Do not introduce a registry, factory hierarchy, service locator, callback
  system or plugin mechanism merely because future variants are imaginable.
- Keep data flow explicit. Avoid hidden control flow and configuration magic.
- Prefer small ordinary derived types and procedures to deep inheritance.
- A sophisticated algorithm is justified by the mathematical/scientific need
  or measured performance, not architectural fashion.

"Simplest" does not mean scientifically incomplete. Model assumptions,
conservation laws, derivatives and error status remain explicit.

## YAGNI — catalogue is not implementation

The capability and numerical-method catalogues record admissible future design
space. They do not justify code.

A feature, abstraction, configuration option, compatibility hook or extension
point enters production only when the **current prioritized vertical slice**
needs it.

In particular:

- do not add placeholder state fields or namelist options for planned features;
- do not implement a Cartesian product of the documented model axes;
- do not add an interface solely because a future catalogue entry might use it;
- do not preserve an accidental young API for hypothetical future users;
- delete obsolete scaffolding rather than maintaining speculative fallback
  paths.

Architecture should avoid obvious dead ends, but implementation decisions are
taken at the latest responsible moment.

## DRY — one authoritative source of knowledge

Duplicated *knowledge* is forbidden even when the text/code looks different.

There is one authoritative production definition for:

- each physical/model equation or operator;
- sign/unit/normalization conventions;
- configuration semantics and defaults;
- state ownership and conservation/exchange ledgers;
- imported-format interpretation;
- derivative semantics and error/certificate status.

Project-owned primal/JVP/VJP/enclosure/performance specializations should be
generated from one retained mathematical source wherever practical.

Do not copy a generic implementation from FortFEM/FortNum/FortSym/FortAD/
FortIO/libneo into KIN6D to modify it privately; improve the owning dependency
when the capability is generic.

### Deliberate DRY exceptions

Duplication is allowed only when independence or execution requires it:

1. **Independent verification/oracles.** An analytic/reference implementation
   may intentionally not share production code, because sharing would destroy
   independence.
2. **Performance/reference implementations.** A simple reference action may
   coexist with an optimized CPU/GPU action. Prefer generated specialization
   from one mathematical source; if handwritten duplication is unavoidable,
   cross-test it and document why it exists.
3. **Semantically different physics.** Do not create a wrong abstraction merely
   because two models currently contain similar algebra.

The goal is one source of truth, not maximum textual deduplication.

## SRP and separation of concerns

A component should have one coherent reason to change.

Keep separate:

- physical laws/model composition;
- discretization and representation;
- nonlinear/linear/eigen solver orchestration;
- hardware scheduling and optimization;
- configuration resolution;
- external-format adapters;
- scientific I/O/provenance;
- diagnostics/inference;
- visualization/UI.

A module that parses user input, defines plasma physics, runs PETSc and writes a
result file is already too broad.

Mutable scientific state is explicit and owned. Avoid module-global/singleton
state in reusable library paths; pass state/context explicitly so multiple
instances, threads and tests can coexist.

## Dependency inversion and the clean dependency rule

Scientific policy must not depend on replaceable outer details.

Conceptually the dependency direction is:

    physics/model policies
        <- application/use-case composition
            <- capability ports/adapters
                <- PETSc/libneo/FortIO/CUDA/MPI/CLI/web/site details

This is a direction rule, not a required directory layout.

KIN6D physics must not import PETSc, libneo, FortIO, CUDA, web or scheduler
internals. Outer adapters may depend on KIN6D-owned contracts.

Not every replaceable detail requires an abstract interface on day one. A seam
is extracted when a real second implementation/consumer, testing boundary or
high-cost external dependency demonstrates the need. Until then, a small direct
adapter is preferable.

## Interface segregation

Prefer small capability-specific contracts to god objects with many optional
methods.

Examples:

- field values;
- spatial derivatives/JVPs;
- rigorous enclosure;
- writable result sink;
- equilibrium tangent/adjoint;

are separate capabilities unless a concrete consumer proves they belong
together.

A backend that cannot satisfy a semantic capability does not return a dummy or
silently degrade; the capability is absent and the unsupported combination is
rejected early.

## Open/closed without speculative plugin machinery

The common core should normally gain a new model/representation/provider by
composition with existing contracts rather than by editing named-code branches
throughout the solver.

However, do not build a generic plugin architecture before concrete variants
exist. First implement the vertical slice; extract the stable common seam after
real implementations reveal it.

## Substitutability and semantic contracts

Two implementations of one capability are substitutable only if they share the
same documented physical semantics: units, orientation, normalization, domain,
pairing, error status and failure behavior.

Use contract tests across alternative backends. If semantics differ, represent
a different capability/model instead of hiding the difference behind a flag.

## Composition over inheritance

Favor composition of orthogonal physics/numerics capabilities and shallow
interfaces. Deep class hierarchies are strongly discouraged. Shared behavior
comes from reusable operators and explicit composition, not inheritance trees.

## Make invalid states hard to represent

Unsupported combinations are rejected during configuration resolution or
construction, not discovered deep inside a nonlinear solve.

State ownership prevents double counting of species, current, pressure and
energy. Units/signs/model assumptions belong in typed metadata/contracts rather
than comments that callers may ignore.

## Testing and verification

Tests protect behavior and scientific contracts, not implementation structure.

- A bug fix starts with a reproducer at the owning layer.
- Prefer analytic/manufactured/structural oracles over self-comparison.
- Alternative implementations share contract tests but retain independent
  oracles.
- Generated-code equality is provenance, not independent correctness evidence.
- Keep fast tests small; large scientific/performance cases run only where
  their evidence is needed.
- Do not mock away the mathematical property the test is intended to verify.

## Performance without premature optimization

Correct/simple/reference code comes first unless performance is itself the
current vertical-slice requirement.

Optimize measured bottlenecks at their owning layer. Performance-specific
specialization must not fork the scientific definition. Keep a simple reference
action when it adds verification value.

Do not generalize an optimization before the profile demonstrates reuse.

## Error handling and observability

Fail early at boundaries with enough context to diagnose the violated contract.
Do not silently repair scientifically meaningful invalid input unless the
normalization/correction is explicitly part of the adapter contract and
recorded in provenance.

Solvers expose convergence/failure status and relevant residuals; adapters
expose transformations; certificates state precisely what is and is not proven.

## API evolution

KIN6D is young. Internal APIs may change aggressively when that reduces
complexity or reveals a better abstraction.

Backward compatibility becomes a constraint only for a real external consumer
or explicitly stabilized public contract. Do not accumulate deprecated paths
for hypothetical users.

For shared mature dependencies such as libneo, preserve existing consumers with
bounded shims/migration when practical while fixing generic behavior upstream.

## Boy-scout rule

Leave every touched production path modestly better than it was, but only in
the bounded area needed by the current task. Fix nearby clear defects, tests,
naming, explicit state or ownership problems when doing so reduces future risk.

Do not turn the boy-scout rule into an invitation for unrelated cleanup or
repository-wide refactoring.

## Review checklist

Before merging a nontrivial change, ask:

1. What current vertical-slice requirement justifies this change?
2. Is this the simplest complete solution?
3. Did we create an abstraction before two real uses demonstrated the seam?
4. Is any scientific/configuration knowledge duplicated?
5. Does each changed component have one coherent reason to change?
6. Do dependency arrows point from outer details toward inner policy?
7. Could a smaller capability interface replace a god object/flag?
8. Is mutable state explicit and reentrant?
9. Is the new behavior independently testable?
10. Is optimization based on evidence, and does the reference definition remain
    intact?
11. Did we accidentally promise or implement catalogue scope not needed now?
12. Can obsolete scaffolding be deleted instead?

Chris&AI
