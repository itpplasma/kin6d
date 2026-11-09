# MILCA

<a href="https://cloud.tugraz.at/index.php/s/cMDgzfagwJnEyiq/download?path=%2F&amp;files=milca.png"><img src="https://cloud.tugraz.at/index.php/s/cMDgzfagwJnEyiq/download?path=%2F&amp;files=milca_small.png" alt="MILCA" width="420"></a>

MILCA is a Fortran research code for plasma equilibrium, stability and kinetic
models. Current work centres on axisymmetric and three-dimensional
magnetohydrodynamic (MHD) equilibria, with small kinetic, classical-particle
and quantum reference problems alongside them.

<a href="https://fortran-lang.org/"><img src="https://raw.githubusercontent.com/fortran-lang/assets/cc6fafc67918ea2daef46b965dcec67ed19ef585/powered_by_fortran/wide/powered_by_Fortran_wide.svg" alt="Powered by Fortran" width="240"></a>

This is research software under active development: interfaces, inputs and
outputs change frequently. Each documented case states its verification and
error-bound status. Numerical kernels and their derivatives are generated with
[FortSym](https://github.com/lazy-fortran/fortsym); generated sources are
included, so ordinary builds do not require FortSym. Supported error bounds
use interval arithmetic from [FortNum](https://github.com/lazy-fortran/fortnum).

## Build

You need GNU Fortran with `real128` support (version 13 or newer recommended),
a C compiler, CMake 3.22 or newer, Ninja, BLAS/LAPACK and SuiteSparse AMD.
The first configure needs network access to fetch pinned dependencies.

```sh
cmake --preset cpu
cmake --build --preset cpu
ctest --preset cpu -R '^characteristic_and_steady_oracles$'
```

Select tests by name or model class; do not run the full suite. For example:

```sh
ctest --test-dir build -j 32 -L equilibrium -LE slow
```

OpenMP is enabled by default. Build options, per-class development presets and
optional CUDA, FFmpeg and Lean support are described in
[the tooling guide](docs/tooling.md) and [CMakePresets.json](CMakePresets.json).

## Run an example

Solve a fixed-boundary Grad–Shafranov equilibrium for a circular plasma:

```sh
cmake --build --preset cpu --target grad_shafranov_demo
./build/grad_shafranov_demo --config examples/grad_shafranov_circle.in
```

Applications read Fortran namelists and write NetCDF-4 results with resolved
configuration and provenance. For the first kinetic laboratory, including
free streaming with a continuous error band:

```sh
./build/milca_demo --snapshot
```

## Models and documentation

- [Documentation index](docs/README.md) and [figure gallery](docs/chapter-gallery.md)
- [Grad–Shafranov equilibrium](docs/grad-shafranov.md) and
  [free boundary](docs/gs-free-boundary.md)
- [3D ideal equilibria](docs/ideal-map-3d.md),
  [3D free boundary](docs/ideal-free-boundary.md) and
  [equilibria with islands](docs/island-equilibrium.md)
- [Linear MHD stability](docs/stability.md)
- [First kinetic laboratory](docs/first-example.md) and
  [steady Poisson–Fokker–Planck laboratory](docs/steady-laboratory.md)
- [Scope](docs/scope.md), [design principles](DESIGN.md) and [roadmap](ROADMAP.md)

The numerical libraries are [FortFEM](https://github.com/lazy-fortran/fortfem),
[FortSparse](https://github.com/lazy-fortran/fortsparse),
[FortIO](https://github.com/lazy-fortran/fortio) and
[FortPlot](https://github.com/lazy-fortran/fortplot), alongside FortSym and
FortNum. Coil readers and field utilities come from
[libneo](https://github.com/itpplasma/libneo).

MIT license; see [LICENSE](LICENSE).
