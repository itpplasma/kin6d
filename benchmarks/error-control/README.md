# Error-control measurement receipts

These fixed-input CPU measurements support the policy in
[error-control-strategy.md](../../docs/error-control-strategy.md). They compare
claim strength, width and cost without changing the mathematical inputs.
The CSVs retain every trial; JSON and the atlas retain medians, ranges,
machine/compiler, source revision and timing boundaries. Rendering is excluded.

- `baseline.csv` / `baseline.json`: original 21-path sweep.
- `measured.csv` / `measured.json` / `table.md`: 23 paths, including standalone
  Poisson and streaming comparison, with embedded Q1/Q7 stages instrumented.
- `modes.csv` / `modes.json`: five tiers on identical frozen Q0 candidates,
  125 and 343 modes, five trials each; audit/workhorse order alternates.
- `reference.csv` / `reference-summary.txt`: binary128 diagonal-exponential
  oracle at three step sizes and seven fractions, all D/S/A radii and efficiency
  indices. Worst D midpoint index is 5.34; all-fraction maximum is 49.3 near
  the right endpoint. Both containment and budget close; the ~10 efficiency
  target fails there. Unknown true error on large fixtures is explicitly
  unavailable, not replaced by the strict radius.
- `inputs/` and `crosschecks/`: existing benchmark app recipes and numeric
  receipts. They exercise additional times/temperatures/goals, so their widths
  need not match the dedicated drivers. Their exercised app code is unchanged
  from `63ff062`. The transport goal and derivative targets fail on its
  coarse affine-source fixture and remain reported as failures.

## Reproduction

FortNum additions remain unpushed on local branch `error-strategy`, commit
`a431cf220efebc583c5fa39d032437e36089cfe3`, based on `e97a493`.
The isolated checkout is `/tmp/fortnum-error-strategy`; the original FortNum
main checkout was not switched. A clean configure needs this local source
(or a new checkout of that local branch). Fetching the unpushed pin remotely
cannot work.

```
FC=gfortran cmake --preset cpu \
  -DKIN6D_FORTSYM_BUILD=/Users/ert/code/lazy-fortran/fortsym/build \
  -DFETCHCONTENT_SOURCE_DIR_FORTNUM=/tmp/fortnum-error-strategy
cmake --build build -j 10
ctest --preset cpu -j 8
build/benchmark_error_continuum > /tmp/continuum.csv
build/benchmark_error_quantum > /tmp/quantum.csv
build/benchmark_error_microscopic > /tmp/microscopic.csv
build/benchmark_error_modes > /tmp/modes.csv
```

Run benchmarks sequentially without competing builds/tests. The executables
are registered build targets, not CTest performance gates. Tiny solve timings
are averaged where supported and still have clock-resolution noise; compare
absolute cost before interpreting enormous ratios. Paired-difference rows
subtract separately measured candidate work. Embedded RK E indicators are
already paid for by the candidate: their query overhead is not a complete
standalone RK defect cost.

Q1 divides initialization into eigensolve/candidate, residual/orthogonality
and reconstruction diagnostics. Q7 optional timers separate low-rule work
(including the low outer inversion-node table), inner low-rule work at high
nodes, and outward tail evaluation. The final inversion query is assigned to
estimation; it includes the candidate observable and indicator propagation.
Probe overhead is retained, not subtracted. The uninstrumented original Q7
whole-table time is in the baseline; profiling is not the production default.
Strict tail time does not imply the full Q7 response is certified.

Q0 D/S radii enclose a finite supplied Hamiltonian, not the spatial continuum
or assembly/model uncertainty. A reports epsilon_b=5% and separately checked
omission allowance; its acceptance uses 1.05 times the displayed radius.
P uses 24 samples, delta=.05, seed 27182, and an ideal iid coordinate law with
a Threefry realization. It bounds one residual norm only; its current widths
103/368 are unusable for the stated scale. E quantities cannot accept a
certificate. H1 intervals retain batch-means/mixing assumptions. The 1% width
comparison alone does not upgrade any of these statuses.

The existing FP benchmark evaluates the same checker again for acceptance;
this doubles its checker cost. This is reported separately rather than hidden
in the dedicated driver's single-check total. Reuse requires immutable data
or a checked identity key because the FP state has public mutable fields.
Fermi, transport and streaming receipts likewise preserve their own setup,
reuse, derivative and band boundaries.

## Verification and development tools

Release and runtime-checked Debug use native CMake/CTest. The independent
oracles cover generated allowance value/JVP/VJP, binary128 D/S containment,
A's true-error ratio and failed epsilon gate, plus sampled norm coverage.
FortNum's CMake gates cover cancellation, interval vertices, underflow and
nonfinite rejection, and 2000 streams each on constant, heterogeneous and
sparse vectors. Coverage experiments supplement the concentration theorem.

The installed `/opt/homebrew/bin/fo` rejects `fo gremlin` as an unknown command
(SHA256 `f80316b9922aa4d216fcaeaf8159507c7b442ce244b679d37403c38415be2dec`).
The recorded `/var/tmp` resident build is absent on this Mac. Next action:
install/rebuild the current CMake-capable Gremlin delivery and attach the
native project preset; this run used native build/CTest feedback throughout.
FortNum's attempted `fpm test --target test_verified_matvec --profile release`
hits the existing C smoke test's duplicate `_main` discovery/link defect.
Its reproducer and next action are retained upstream; native CMake passed.
Neither limitation changes the scientific certificate claim.

Chris&AI
