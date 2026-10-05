# Configuration, UX, I/O and reproducibility

Status: canonical interface policy, 2026-10-05. The current executable slice
does not yet implement this full interface.

## User goal

The normal workflow should state the physical problem/model, desired observable
and requested accuracy/certificate target before exposing low-level numerical
knobs. KIN6D should choose or adapt ordinary numerical settings when it can
measure their error. Experts retain an advanced numerical view.

Certification does not remove model judgment. Model selection remains explicit
and may be assisted by documentation or AI, but AI-generated configurations are
only proposals passed through the same KIN6D validator as hand-written input.

## Canonical configuration

Fortran owns typed configuration semantics. Fortran namelist is the canonical
human-editable syntax. Every run preserves both the exact original input and a
fully resolved input in which defaults and automatic choices are materialized.

A machine-readable schema is generated from the same Fortran definitions for
documentation, AI proposals, Python convenience and browser forms. Do not
maintain independent scientific schemas in Fortran, Python and TypeScript.

The eventual reference CLI should provide small operations such as validate,
resolve, run, inspect, extract-input, doctor and schema. Scheduler/resource
settings stay outside the plasma-physics namelist.

## Browser and Python convenience

There is no mandatory desktop application and no mandatory Python runtime.
After the schema/result contract is stable, an optional small Python server and
TypeScript frontend may provide generated forms plus a synchronized raw
namelist editor, validation, scheduler submission and result views.

An SSH tunnel to a user-local server is sufficient for personal cluster use;
sites may integrate the same service with an existing portal. KIN6D does not
own SSH credentials. An optional Python package may provide notebooks and
parameter scans but does not own scientific semantics.

## Scientific result format

Use NetCDF-4 semantics stored in the HDF5 subset supported by FortIO. Do not use
arbitrary HDF5 features that make ordinary NetCDF tooling unable to read the
file. A result should contain groups for input, provenance, geometry/state,
observables/diagnostics and certificates.

A valid KIN6D result always contains enough information to recover its effective
KIN6D input. Embed original and resolved namelists. Record code/dependency
revisions, compiler/version/flags, strict-vs-fast build identity and material
hardware/runtime information.

Large external files need not be embedded verbatim. Record their logical/source
name, detected format/version, units/coordinate/sign conventions,
transformations applied, content hash and access/redistribution limitations.

## Interchange

Use libneo as the shared adapter/convention layer when it already supports a
format, including EQDSK/gEQDSK, VMEC wout, Boozer/boozmn/chartmaps,
SPECTRE/SPEC data, JOREK data and coil formats. KIN6D consumes a normalized
physical state/field/domain contract, not a named-code mode.

## Compiler and build portability

Production support targets GNU gfortran, Intel LLVM ifx, NVIDIA nvfortran and
LLVM Flang. CMake must select or probe compiler-specific flags; GNU-only flags
are not generic project policy. Keep fast candidate, conservative/reference and
rigorous-certifier build contracts distinct.

A compiler is supported only when relevant tests actually compile and execute
with it. Releases should eventually support offline/pinned cluster builds
without mandatory outbound dependency downloads.

## Result UX

Lead with the requested observable and its uncertainty ledger: numerical
certificate, model-reduction status, physical-input uncertainty and stochastic
uncertainty. Advanced residuals, spectra, meshes and solver diagnostics remain
available. TIAGO reconstruction/posterior uncertainty is displayed as an
orthogonal layer rather than folded silently into a KIN6D deterministic
certificate.

Chris&AI
