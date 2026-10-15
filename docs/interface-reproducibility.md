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

## Implemented contract (campaign F5)

Status 2026-10-05. Every application in `app/` uses `kin6d_io`
(`src/apps/result_io.f90`, linked against FortIO pinned in `CMakeLists.txt`):

- `--config FILE` names an optional namelist file; all other arguments
  (`--snapshot [frame]`, benchmark sizes, paths) are unchanged; a command-line
  override is applied before the resolved text is captured. Each
  application owns one group with its defaults in one place. gfortran rejects
  a group named like the program, so a demo's group drops `_demo`
  (`&debye_sphere`) and a benchmark keeps its name (`&benchmark_fp`) with
  the program renamed `benchmark_fp_main`. The two periodic-pair demos also read the shared `&periodic_pair`
  group (`src/apps/periodic_pair_config.f90`: length, mass, coupling, time_max,
  frames). `classical_pair_demo` rejects a mass or coupling other than the
  declared electron pair rather than silently simulating another system.
- Unknown entries and malformed values fail with the group name; a group in
  the file that the application does not read is an error. A missing group
  keeps its defaults.
- The result is `output/<application>.nc` (NetCDF-4/HDF5). It holds named
  double variables (scalar to rank 3) with `units`, `long_name` and, where the
  application has one, `kin6d_status` (candidate, certified-discrete,
  sampling interval, measurement, ...). Global attributes hold the literal
  file text, the resolved namelist after automatic choices, the command line,
  the build-time git commit and tree state and a runtime recheck of both, the
  build type/flags, `compiler_version()`/`compiler_options()`, the
  FortSym/FortNum/fortplot/FortIO revisions and the UTC creation time.
- The resolved text is itself a valid input and resolves to itself.
- `result_io_contract` (CTest) reads the file back through FortIO's HDF5
  reader and, when they are installed, through `ncdump`/`h5dump`. It also
  checks the rejection cases. The gallery renderer records each `.nc` with
  its sha256.

Integrating a new application takes a few lines: `call run%start(name)`, read
its group between `run%open_config`/`run%close_config`, write the group to
`run%open_resolved`/`run%close_resolved`, then `result%create(run)`,
`result%write(...)` and `result%close()`. Finally add the target to the
`kin6d_io` application list in `CMakeLists.txt`.

Planned and not yet implemented: groups for input/provenance/state/
certificates inside the file, a generated schema, the validate/resolve/
inspect CLI, hashes of external inputs, hardware/runtime records, and
derivative or ledger variables beyond the per-variable status string.

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
