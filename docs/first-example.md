# First 1X1V laboratory

The first executable is a deterministic prescribed-field kinetic laboratory.
It solves a source-balanced stationary problem directly and animates free
streaming with an exact characteristic oracle. This is the first implemented
slice of the [student laboratory](student-laboratory.md); self-consistent
Poisson, the full collision hierarchy and Lean verification remain open.

## Build and run

Use GNU Fortran with `real128`, CMake 3.22+, Ninja and FFmpeg for MP4 encoding.
The first configure downloads the pinned public FortPlot and FortNum sources.
Subsequent builds are local. Python and a CAS are unnecessary for normal builds.

```sh
./scripts/demo.sh
# output/free-streaming.mp4, output/final-snapshot.png, output/streaming.csv
./scripts/demo.sh --snapshot
# Only the final frame and its snapshot-bound record; no encoder required.
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

CMake presets and registered CTest tests are the build/test authorities.
The selected, digest-verified Fo driver also supports `fo test` and
`fo exec kin6d_demo`. The shell helper uses the native CMake/CTest route
explicitly while resident Gremlin's CMake input capture remains defective.

## Mathematical problem

Use dimensionless thermal velocity, streaming time and spatial coordinate:
$x\in\mathbb R/(2\pi\mathbb Z)$, $v\in[-4,4]$, $E=0$. The movie solves

\[
\partial_t f+v\partial_x f=0,\qquad
f(x,v,0)=e^{-v^2/2}(1+A\cos x),\quad A=1/2.
\]

The exact solution is $f^*(x,v,t)=e^{-v^2/2}[1+A\cos(x-vt)]$.
Velocity does not advect, so the finite velocity boundary needs no artificial
inflow condition. There is no claim about the omitted infinite-velocity tails.

Store the three Fourier coefficients $m(v),a(v),b(v)$ and use piecewise-linear
velocity reconstruction. The exact streaming subflow rotates $(a,b)$ at each
velocity node. The implementation reuses phase factors for a fixed timestep;
changing the timestep recomputes them. The initialized velocity grid is immutable
for evolution: reinitialize the state to change grids. This tiny configuration
needs no FortFEM spatial mesh; its single Fourier harmonic is exact in space.
Future general spatial discretizations remain owned by FortFEM.

The direct stationary example solves

\[
 v\partial_x f+\nu(f-M)=S,\quad M=e^{-v^2/2},\quad
 S=\nu A M\cos x-v A M\sin x.
\]

For each velocity node, solve the two-by-two system for $(a,b)$. Its exact
solution is $a=AM,b=0$, with $\nu=3/4$. This is a manufactured open-system
steady state, not a self-consistent kinetic equilibrium. It is obtained without
time marching and keeps the direct-residual scientific spine visible.

## Continuous snapshot bound

The candidate denotes the continuous Fourier/piecewise-linear reconstruction,
not just the sampled heatmap. FortNum encloses the exact node coefficients;
the largest sum of their absolute discrepancies is $\eta_{\rm node}$.
For each cell of width $h$, let $V=\max |v|$ over the cell, and let
$M_+$ enclose its largest Maxwellian value. Differentiating the exact flow gives

\[
 |\partial_v^2 f^*|\le
 M_+\{(1+A)(V^2+1)+2A|t|V+A t^2\}=B_K.
\]

The linear-interpolation remainder is bounded by $h^2B_K/8$. Therefore

\[
 \sup_{x,v}|f_h-f^*|\le
 \eta_{\rm node}+\max_K h_K^2B_K/8.
\]

The radius and band endpoints use outward interval arithmetic, including
Taylor/remainder transcendental enclosures from FortNum; they do not assume
that a system libm result is exact. Generated Fortran is compiled without fast
math or floating-point contraction. The plotted lower-profile band encloses
the continuous slice between velocity nodes. The heatmap itself is a sampled
view with fixed color normalization, not a plotted interval at every pixel.

This is a reviewed analytical bound evaluated with interval arithmetic.
It is not yet a Lean-checked theorem. It bounds each recorded snapshot at its
recorded floating-point time, not every instant between movie frames.

## Symbolic source and derivatives

`gen/gen_laboratory.f90` defines the numerical expressions. FortSym emits the
primal leaves and JVP/VJP families for initialization, interpolation,
reconstruction, streaming factors/application, steady sources/solve/residual, strict `real128` characteristic reference
and interval leaves into `src/generated/kin6d_leaves.f90`. Handwritten Fortran
owns orchestration, state, tests and plotting. Normal builds use committed
outputs; regeneration is explicit:

```sh
cmake --preset cpu -DKIN6D_FORTSYM_BUILD=/path/to/fortsym/build
cmake --build --preset cpu --target gen_laboratory
./build/gen_laboratory
cmake --build --preset cpu
ctest --preset cpu
```

Use FortSym revision `90b68132dbecb83b982826dc8b193f3a627d94e1` for regeneration.
It adds generic optional exponential/trigonometric enclosure emission, with
compiled behavioral tests. Generator consistency does not replace physics tests.

## Evidence and performance

The registered test checks characteristic values against `real128`, quadratic
velocity refinement, the manufactured steady coefficients, norm preservation,
JVP/VJP pairing and an independent centered parameter perturbation. It rejects
zero/NaN bounds, invalid amplitudes, NaN time and corrupted candidates, and
checks changed-timestep/reverse-flow cache behavior. The steady tangent also
satisfies the implicit residual equation and a closed rational oracle.

The animation contains 121 frames at 24 fps, 960 by 720 pixels; physical time
runs from zero to 7.5 in exactly representable increments of 1/16. On the
2026-10-05 local GNU Fortran run, the direct steady residual was
$1.11\times10^{-16}$ and the final 256-cell continuous radius was
$3.76868\times10^{-3}$. These numbers are run evidence, not portable thresholds.

Run the separate performance audit with `fo exec benchmark_streaming` or
`./build/benchmark_streaming`. It alternates reference/cached order over five
trials, includes state copies, checks output agreement and reports median,
minimum and maximum CPU time. For 131073 nodes and 200 applications, local
medians were 0.246049 s uncached and 0.042891 s cached: a measured 5.74-fold
improvement. One-time factor generation is outside these amortized timings.
Production evolution allocates no per-step scratch and stores no trajectory,
phase-space operator, global Jacobian or derivative tape. Rendering uses one
bounded image tile; FortPlot stages PNG frames on disk for FFmpeg.

## Upstream repairs and workflow limits

FortPlot required selected-subplot mesh/fill routing, explicit color-limit
handling, legend geometry and subplot-aware animation capture. A new pixel
oracle exercises capture without pre-rendering. Its focused tests, essential
CI suite and required `make verify-artifacts` rendering gate pass. Changes live
upstream; kin6d pins the repaired revision in CMake.

The session used Fo driver SHA256
`6e54f805af13272e859ebf6c7ec48d205a22d7fdac7c84d435fdffc922bc76ee`.
Public delegated test/exec commands pass. Resident capture reports
`cannot parse fpm.toml for generation inputs` on this CMake-only project.
The upstream CMake plan owns that unresolved input-inventory/generation gap;
no fake FPM manifest was added. Native CMake/CTest results are reported as such.

FPM 0.13.0 also reused a stale default-profile plotting executable after a
source edit and a build in a second compiler profile. Preserving that profile
and rebuilding the missing default profile removed the stale binary; the same
rendering gate then passed. This is recorded as an unresolved FPM cache defect,
not evidence of a remaining FortPlot rendering failure.

Published upstream revisions:

- FortSym: `90b68132dbecb83b982826dc8b193f3a627d94e1`.
- FortPlot: `3abf42c9a3c1c7026e1eda5e4aea8d607f2b8f1e`.
- Fo reproducer/report: `d31db5c` (records the open defect; no backend-fix claim).

The completed movie took 10.88 s wall time in the measured local run, including
rendering and encoding; the command-tree peak RSS was 443676 KiB. This includes
FFmpeg's own buffers and threads and is not the kinetic state size. The solver
holds eight vectors of 257 doubles (16448 bytes) for the fine-grid state,
including velocity, factors and work arrays. Encoded-frame inspection checks
first, middle and last frames; the final figure was also reviewed in grayscale.

A clean build fetching the exact public dependency pins also passed the
registered characteristic/steady test, without local source overrides. The
`free_streaming` atlas row now records this implemented scope; other atlas
rows remain planned.

The external literature sync also completed: 75 curated DOI records and 31
verified PDF attachments, with offline empty-search/stale-key/access-error
oracles. Imported literature remains marked unreviewed. Bibliographic assets
and their receipts are kept outside the solver repository.

The literature path also needed two published Helpy repairs: `6ebb1fb`
(typed Zotero write payloads) and `8bfe74e` (dry-run write safety). Their focused
behavioral tests and installed CLI smoke checks passed; no credentials or
library exports are stored in kin6d.

Chris&AI
