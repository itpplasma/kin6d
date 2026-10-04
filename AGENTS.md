# kin6d agent instructions

## Load before work

Read README.md, DESIGN.md, and ROADMAP.md for the project map, current
milestone, and acceptance gates. Read docs/tooling.md before implementation
and docs/certification.md before proposing a certified numerical result.
The repository currently contains design documents only. Do not claim solver,
certifier, formal-proof, CUDA, or build evidence that has not been produced.

## Required policy

- Scientific and technical communication is English. Sign messages drafted
  or sent on Chris's behalf to other people with `Chris&AI` on its own final
  line unless explicitly instructed otherwise.
- Scientific/application code, generator drivers, and the certifier are
  Fortran. No required Python layer. Lean is the trusted proof layer;
  optional native CUDA lives behind ISO_C_BINDING and a thin C ABI.
- Every project-owned numerical kernel, including reference, derivative,
  enclosure, and CUDA leaves, must be emitted by fortsym. Retain generators,
  exact definitions, assumptions, provenance, and independent oracles.
  Never edit generated numerical code or hide formulas in handwritten glue.
- Generic dependency defects and missing features are repaired in their
  owning repositories. The user authorizes the controller to push locally
  verified fixes to owned upstream main promptly. Coordinate with active
  owners, retain a focused reproducer, pin the fix, and recheck kin6d.
  No permanent consumer workaround or copied generic numerical kernel.
- Follow the Fortran skill for implementation. Preserve CMake/CTest
  presets/profiles and use resident fo Gremlin feedback with a verified
  driver. Missing support is an upstream issue, not a reason to omit a gate.
- Make the smallest complete change. Defer abstractions until required by
  the next capability, existing duplication, or measured cost/accuracy.
- The first executable commit and every later increment must retain a
  runnable end-to-end default demo with meaningful Fortran postprocessing.
  Read docs/demonstrations.md. Plot candidates with full certified intervals
  and overlay general models' limiting cases with earlier reference cases.
  Pointwise bands require pointwise/cellwise bounds, not a generic L2 radius.
  Never commit generated plot images, PDFs, SVGs, frames, or animations.
- Tests require independent behavioral oracles. Regeneration equality,
  patch matching, and a small algebraic residual are not continuous tests
  or continuous-PDE certificates.
- Collision physics, discrete conservation/entropy treatment, and actual
  certificate strength are separate contracts. Read docs/collisions.md;
  compose fixed-coefficient Landau/BB actions from shared pair primitives.
  Never enforce closed-system invariants on an incompatible reservoir model
  or advertise a diagnostic field-only action as a physical closure.
- Read docs/species-models.md and docs/devices.md before adding a reduction
  or machine case. Species representation, collision pairs, and field models
  require explicit compatibility. Generic 3D geometry has no fundamental
  tokamak type. Fixed species are reservoirs; fields can remain self-consistent.
  Keep numerical certificates, model differences, and proved reduction
  remainders distinct. Production needs scoped physical-parameter device gates.
- Preserve unrelated edits and stage explicit paths. One controller owns
  integration and promotion; workers return evidence. Freeze review inputs
  by exact base commit and patch digest when evidence is uncommitted.
- Never touch `faepop*` or `faepcr*` hosts without an explicit user request
  naming the host. The faepmac1/faepmac2 exception remains as instructed.
- Keep every response and tool call below about 300 lines / 8k tokens.
  Write long files incrementally and report summaries. Repeat this limit
  in any authorized subagent prompt; do not delegate without authorization.
