# Steady Poisson and FP laboratory

The CPU laboratory now includes a self-consistent source-balanced 1X1V
Vlasov–Poisson–Fokker–Planck reference, a bounded homogeneous OU/FP steady
state, and a separately formalized prescribed-field Krook reference. Each
case declares its own boundaries and proof level. The original Gaussian
free-streaming demonstration remains available.

## Build, check and plot

```sh
./scripts/demo.sh --snapshot
# output/final-snapshot.png, output/fp-steady.png,
# output/poisson-fp-steady.png, output/steady-observables.png,
# output/steady-bounds.csv, output/steady-certificate.txt
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

GNU Fortran with binary128 support and CMake/Ninja are required.
Lean and mathlib are optional and pinned in [the proof package](../proof/README.md).
Initialize that package once with `cd proof && lake update && lake exe cache get`.
The cache download is optional when building the pinned mathlib sources.
At the reference-infrastructure checkpoint, `cpu` and `debug` ran three
independent numerical gates. The [general transport increment](general-transport.md)
adds a fourth gate for source-driven continuous state and tangent bounds.
The explicit `formal` preset additionally checks the Lean library, rational
fixtures, corrupted-certificate rejection and a Fortran-produced certificate:
`cmake --preset formal`, `cmake --build --preset formal`, `ctest --preset formal`.
A numerical-only build does not supply formal evidence.

## Source-balanced self-consistent 1X1V case

Use dimensionless `q=m=epsilon0=1`, spatial period `2*pi`, a positive integer
mode `k`, and finite velocity interval `[-1,1]`. A fixed charge background
has density `-1`. The active distribution and induced electrostatic field are

```
B(v) = 3/4*(1-v²),  integral[-1,1] B(v) dv = 1,
f(x,v) = B(v)*(1+a*cos(k*x)),
rho(x) = a*cos(k*x),
phi(x) = a/k²*cos(k*x), E(x) = a/k*sin(k*x).
```

The direct solver takes the cosine moment of the prescribed source, `drive`.
For Krook rate `nu>=0`, FP rate `beta>=0`, and heat-bath temperature
`theta>0`, it solves

```
(nu+3*beta*theta)*a = drive,
v*d_x f + E*d_v f + nu*(f-B) - beta*d_v(v*f+theta*d_v f) = S.
```

The denominator must be positive and `|a|<1`. The source is generated from
the continuous profile and includes the nonlinear second spatial harmonic
from `E*d_v f`. The field is computed from the candidate density through
Poisson, with exact zero-mean charge compatibility and zero-mean potential
gauge. Incompatible net charge is rejected, including a nonzero charge too
small to change a floating-point sum near unity.

Both velocity endpoints have Dirichlet `f=0`. With FP diffusion these are
open particle boundaries. Integrating the FP action gives
`-3*beta*theta*(1+a*cos(k*x))`; the generated source supplies the matching
particle transfer. The case has a heat bath and optional Krook reservoir.
Closed-system particle, momentum and energy conservation are not its claims.

This solver resolves the specified finite Fourier/polynomial source family
directly. The explicit exact profile proves existence for that source and
boundary data. The interval checker independently encloses its amplitude
and Fourier field coefficients, then converts their discrepancies to
uniform pointwise bounds over the entire spatial/velocity domain using
`0<=B<=3/4` and `|sin|,|cos|<=1`. It also encloses the spatial mean electric
energy `(a²)/(4*k²)`. These are reviewed continuous bounds; general nonlinear
branch uniqueness and a formal coupled error theorem remain open.

The tests retain independent rational distribution, collision and source
fixtures, including a nonzero second harmonic, analytical derivatives,
Poisson signs/scales, endpoint traces, pure FP and zero-FP limits, and
corrupted-state/certificate rejection. A small sampled residual is a
diagnostic alongside these oracles, rather than the certificate.

## Bounded homogeneous OU/FP case

On `v∈[-2,2]`, solve the stationary heat-bath equation
`d_v J=0`, `J=theta*d_v f+v*f`, with `J(-2)=J(2)=0` and prescribed positive
mass. The candidate comes from a centered zero-flux recurrence followed by
mass normalization. It requires positive recurrence coefficients and rejects
a grid too coarse for that condition. Work and storage grow linearly with
the number of velocity nodes; two vectors hold the velocity and candidate.

Independently, zero boundary flux implies `J=0` everywhere, so
`f*=mass*exp(-v²/(2*theta))/Z`. The exact finite-domain Gaussian integral
`Z` fixes the unique equilibrium. The certifier encloses `Z` with interval
midpoint quadrature and the whole-cell bound
`|G''|<=V²/theta²+1/theta`, where `G=exp(-v²/(2*theta))` and `V` is the
largest velocity magnitude in that cell. Its integral remainder is
`h³*sup|G''|/24`; the linear-interpolation remainder is
`h²*mass*sup|G''|/(8*Z_lower)`.

The sum of the largest nodal discrepancy and largest interpolation remainder
encloses the reconstructed candidate at every velocity. Reference
normalization, arithmetic and interpolation uncertainties are included.
Independent high-precision erf fixtures check the normalized equilibrium and
normalizer; parameter, derivative, symmetry and invalid-data cases are also
retained. The plotted mass interval integrates the linear candidate exactly
with interval arithmetic and widens it by `4*radius`, the domain length times
the uniform error bound.

This is an auxiliary finite-domain 0X1V FP certificate at the reviewed
continuous level. It supplies neither whole-line tails nor a Lean FP theorem.

## Formal continuous reference

Lean proves a continuous maximum principle for periodic differentiable
solutions of `v*d_x f+nu*f=s`, with arbitrary real `v` and positive `nu`.
The theorem does not assume a discrete inverse or continuous stability
constant. At a periodic global extremum the spatial derivative vanishes;
applying the argument to both signs bounds the error by the residual divided
by `nu`.

The concrete profile is `3/4*(1-v²)*(1+a*cos(x))`, prescribed `E=0`,
`nu=1`, and `v∈[-1,1]`. Lean discharges the source identity, regularity,
periodicity, existence, uniqueness, nonnegativity and continuous error bound.
Exact rational certificate fixtures include two approximate amplitudes and
independently justified rejection of corrupted residuals and radii.

`formal_krook_certificate` computes an amplitude with a generated direct-solve
leaf and writes its actual binary64 coefficient as an exact dyadic rational.
The writer checks the required IEEE representation before interpreting its
bits. Lean checks the emitted numerical premises and applies the continuous
theorem to the precisely declared polynomial reconstruction. The default
source amplitude `1/2` gives an exact candidate and zero error radius.

The formal certificate applies to that polynomial prescribed-field family.
The Gaussian sampled streaming state, FP interval arithmetic implementation
and self-consistent kinetic error bounds have separate proof status. Lean
also checks the unit Poisson Fourier equation, derivative signs, charge
neutrality, gauge and field-energy integral as component identities.

## Figure and provenance contract

The FP plot shows both 32- and 64-cell candidates with full pointwise bands
and rigorous Gaussian cell ranges. The coupled plot shows the distribution
slice and induced field with whole-cell ranges expanded by their checked
pointwise radii. It retains the zero-FP and zero-drive cases. Cell envelopes
cover trigonometric variation between plotted nodes; joining tiny endpoint
error bars would not provide that coverage.

All numerical, reference, enclosure and derivative leaves are emitted by
the retained Fortran drivers in `gen/`. Normal builds use committed emissions.
The generating FortSym revision and local verification evidence are recorded
below. The Fortran postprocessor writes generated figures under ignored
`output/`; no figures or movies are tracked.

The demo retains FP node values and their exact binary64 interpretation,
coupled samples, certificate parameters and exact binary bounds under
`output/`. Decimal columns are readable diagnostics; the hexadecimal fields
preserve interval endpoints without inward decimal rounding. The helper also
records the base commit, worktree patch digest and source-content manifest.

## Verification and upstream repairs

The original integrated release and runtime-checked debug builds passed all
four registered tests, including the computed Lean certificate. With the
general transport increment, `cpu` and `debug` each pass four numerical tests;
`formal` additionally retains the proof gate. At 32 and 64
velocity cells, the bounded FP sup radii are `5.585456e-3` and `1.387730e-3`.
The demo's mixed coupled case has distribution/potential/field radii about
`2.50e-16`, `6.66e-16`, and `5.55e-16`. These are measured run values rather
than universal thresholds.

A clean release build fetched the public FortNum and FortPlot pins without
local source overrides. Build, all registered tests, all three native
generators, both demos, the certifier and the streaming movie passed with
Python absent from PATH. Final plots were inspected in color and grayscale.
The figures retain readable limiting-case styles, complete interval bands,
separate subplot legends and the requested observable axis limits.

An independent mathematical review passed against base
`7141d9b163834c6a4c10d241e4b55edc173695ad` and frozen integration patch
`555dc059794d47c879b839b05326e4380c48e5ea6479789158f690565e9f4a39`.
It checked source moments, boundary transfers, nonlinear fixtures, coefficient
bounds, observable energy and band interpolation. Subsequent changes added
benchmark reporting and plotting dependency repairs; these mathematical
definitions and reconstruction contracts did not change.

The optional `./build/benchmark_fp [cells]` audit defaults to 4096 cells, five
trials, 100 solves and ten certifications per trial. It checks the achieved
bound against an independent finite-domain erf fixture and prints every raw
timer sample. On the local AMD Ryzen 9 5950X host with GNU Fortran 16.2.1,
the clean release median solve and certificate CPU times were about `4.86e-5` and `2.07e-2`
seconds. The achieved sup radius was `3.37e-7`; the two state vectors occupied
65552 bytes. The measured process peak RSS was 3460 KiB. Raw samples, resource
measurements and the executable digest are retained under ignored `output/`.
This is a first fixed-input baseline, without a speedup claim.

`./scripts/benchmark.sh` runs 16, 32, 64, 256 and 4096 cells by default; explicit
cell counts may be supplied. It records a declared unit observable scale and
requested 0.01 radius, achieved continuum sup radius, state bytes, raw solver,
certificate and acceptance CPU samples, executable/source identity, compiler,
hardware and process resource usage. Acceptance currently recomputes the
strict certificate and is therefore measured separately and included in the
reported total. Oracle and corrupted-bound validation run outside timers.
The 16-cell case may miss the target; that is reported rather than suppressed.
This known-Gaussian certificate remains a reference benchmark alongside the
[arbitrary represented-source transport checker](general-transport.md).

FortSym revision `102c1cbe392e1ba70d169fd75aeda0f2ffe8ab88` is the required
regeneration pin. It adds the native emitted-kernel audit while retaining the
published integer-division repair `ef981ff`
and the literal polynomial derivative repair `e171c32`. Independent compiled
oracles reproduce both original failures and check their repairs. Repeated
derivatives of raw unsimplified composed constant exponents remain an explicit
upstream limitation with a retained manual reproducer; the laboratory uses
literal exact exponents.

The upstream-first practical certification increment is now implemented as
benchmark infrastructure: FortNum `b91091e26b613af0b4c0ef2a7d2e1eaa412f6a69`
provides strict/aggregate sum and dot experiments, exact cancellation/subnormal
fixtures, a reviewed aggregate-rounding argument and strict exponential timing.
Release/debug validation and the existing rounding/interval gates pass.
Aggregate reductions are benchmark-local candidates, with no production switch.
On a pinned CPU4 Ryzen 5950X run, 4096-element positive sum/dot medians were
0.954/1.210 ns per element versus 8.991/12.519 for strict intervals; cancellation
can widen the aggregate radius. Strict interval exponential evaluation on
[-2,0] measured about 1545 ns per element. These are primitive time/width
tradeoffs, not end-to-end certification speedups.

FortSym's shared-definition audit emits primal/JVP/VJP and reference enclosure
leaves, checks independent real128 values, directional differences and adjoint
duality, and records 66 raw timing samples separately from its short CTest gate.
Native focused gates pass 4/4 and pinned Fo focused tests pass 2/2; generator
execution completes. The interval runtime is a reference emitter fixture, not
a newly validated general arithmetic backend. Broad resident checks still have
known asset-inventory gaps; FortNum's resident cold build remains pending, so
no full resident verification is claimed. The general damped-transport residual
certificate is now implemented with compact Bernstein interval evaluation.
Aggregate production arithmetic and forced-FP certificates remain open gates.

FortPlot pin `e589583a790e1b2dec662c10cca57e03739a0f3d` repairs selected
subplot limits, local legends, legend row spacing and bottom label clearance.
Independent requested-value and raster-pixel oracles, native focused tests,
Fo focused tests, `make test-ci` and `make verify-artifacts` pass. The original
demo was rerendered against the published fix.

Fo revision `48fc48dde87da78320e43c6dbcaaa08b46cc2758` repairs frozen
capture of flattened transitive FPM dependencies. Offline compiled fixtures
retain the captured answer 42 after changing the live source. Focused
inventory, manifest and execution-view gates pass with immutable task driver
SHA256 `317142a367bab97ed19b58bfa6d6d8924887098db23da4839d1150d695be6ed7`.
The original FortSym resident recheck reports `BUILD_PASS` and a passing
required `test_fortsym_kernel_emit` gate. Its frozen acquired-source chain now
includes FortPlot, FortArray and FortIO without refetching. External system
libraries remain outside Fo's modeled inventory, so this is source-chain
repair evidence rather than a claim of complete dependency capture.

Regenerate all numerical leaves with:

```sh
cmake --preset cpu -DKIN6D_FORTSYM_BUILD=/path/to/pinned/fortsym/build
cmake --build --preset cpu --target gen_laboratory gen_poisson gen_fp
./build/gen_laboratory
./build/gen_poisson
./build/gen_fp
cmake --build --preset cpu
ctest --preset cpu
```

CMake checks the FortSym source revision before enabling regeneration.
Normal builds do not require FortSym. The installed Fo driver remains
`6e54f805af13272e859ebf6c7ec48d205a22d7fdac7c84d435fdffc922bc76ee`.
The shared kin6d Gremlin development lane still reports the documented
CMake-only input-capture blocker; native CMake/CTest supplies the evidence.

Chris&AI
