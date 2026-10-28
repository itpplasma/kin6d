# kin6d

kin6d is a Fortran research code for plasma equilibrium, stability and kinetic
models, written so that each numerical result comes with derivatives and, where
the theory allows, a checked error bound. Most of the current work is on
axisymmetric and three-dimensional MHD equilibria, including free-boundary and
magnetic-island configurations, alongside a set of small kinetic, microscopic
and quantum reference problems.

Status: research code under active development; interfaces, input files and
outputs change frequently. Each documented case states how far its result is
verified.

The long-term aim is a single code base in which fluid, kinetic and microscopic
plasma models share operators, derivatives and error estimates, so that the
cheapest adequate model can be chosen for a given question. The scope and
design principles behind that aim are in [docs/scope.md](docs/scope.md) and
[DESIGN.md](DESIGN.md).

## What it does today

- Axisymmetric Grad–Shafranov equilibria on curved triangular finite elements:
  fixed boundary (circular, shaped or supplied), prescribed-q inverse problems,
  and free boundary by FEM–BEM coupling with external coils and a coil-current
  inverse. Checked against analytic Cerfon–Freidberg solutions and a FreeGS
  reference case.
- Three-dimensional ideal equilibria computed as a frozen-in energy map, with
  fixed or free boundary (boundary-element vacuum, coil circuit currents as
  parameters), and evaluators for the exact nested equilibria of Landreman and
  Issan et al. as references.
- Non-nested equilibria with islands on a finite-element exterior calculus
  (FEEC) discretization: vacuum and Solov'ev cases, mapped 3D domains, Taylor
  (Beltrami) closures, metriplectic relaxation, and an island closure provider
  with tangent and adjoint.
- Linear ideal-MHD stability and response: displacement-energy Hessian,
  fixed-wall kink spectra, and per-toroidal-mode response.
- A small 1D–1V kinetic laboratory: free streaming, Krook and Fokker–Planck
  steady states, a self-consistent Vlasov–Poisson–Fokker–Planck case, damped
  transport with an a-posteriori certificate, and particle versus continuum
  comparisons with interval-arithmetic error bands.
- Microscopic classical and quantum reference cases: Coulomb encounters,
  one-component plasma oscillation, Debye shielding and correlation functions,
  finite-basis quantum pair and jellium models, Hartree/Hartree–Fock
  comparisons, ideal Fermi gas, and a proton–proton fusion teaching case.

Pointwise numerical kernels are generated with
[FortSym](https://github.com/lazy-fortran/fortsym) as primal, JVP and VJP
routines; the generated Fortran is committed in `src/generated/`, so FortSym
is needed only to regenerate it. Implicit solves are differentiated at the
equation level. Error bounds use outward-rounded interval arithmetic from
FortNum; a few continuous reference bounds are also checked in Lean
([proof/](proof/README.md)).

## Building

Requirements:

- GNU Fortran with `real128` support (version 13 or newer recommended) and a C
  compiler
- CMake 3.22 or newer and Ninja
- BLAS/LAPACK and SuiteSparse AMD
- network access on first configure: CMake fetches pinned revisions of
  FortNum, FortFEM, FortSparse, FortPlot, FortIO and two coil modules of libneo
- optional: FFmpeg for MP4 animations; Lean/Lake for the `formal` preset;
  CUDA for experimental GPU solves

```sh
cmake --preset cpu
cmake --build --preset cpu
ctest --preset cpu -R '^characteristic_and_steady_oracles$'
```

The full test suite is large; select tests by name (`-R`) or by area label
and tier, e.g. `ctest --test-dir build -L equilibrium -LE slow`. The `dev-fast`
preset adds per-area test presets; see [CMakePresets.json](CMakePresets.json)
and [docs/tooling.md](docs/tooling.md). OpenMP is on by default
(`-DKIN6D_ENABLE_OPENMP=OFF` gives a serial build).

## Examples

Fixed-boundary Grad–Shafranov equilibrium of a large-aspect-ratio circular
plasma. Applications read a Fortran namelist given by `--config` and write
`output/<application>.nc` (NetCDF-4) with the resolved configuration and
provenance:

```sh
cmake --build --preset cpu --target grad_shafranov_demo
./build/grad_shafranov_demo --config examples/grad_shafranov_circle.in
```

First kinetic laboratory (free streaming against the exact characteristic
solution, with a continuous error band):

```sh
./build/kin6d_demo --snapshot   # writes output/final-snapshot.png
```

`scripts/demo.sh` builds and runs the kinetic laboratory demos together.
Further runnable cases, namelists and expected outputs are described in the
pages listed below.

## Documentation

- [docs/README.md](docs/README.md): index of all documentation pages
- [First 1D–1V laboratory](docs/first-example.md) and
  [steady Poisson and Fokker–Planck laboratory](docs/steady-laboratory.md)
- [Fixed-boundary Grad–Shafranov](docs/grad-shafranov.md),
  [free-boundary Grad–Shafranov](docs/gs-free-boundary.md)
- [3D ideal map equilibrium](docs/ideal-map-3d.md),
  [3D free boundary](docs/ideal-free-boundary.md)
- [Non-nested equilibrium](docs/equilibrium-nonnested.md) and
  [equilibria with islands](docs/island-equilibrium.md)
- [Animation and figure gallery](docs/chapter-gallery.md)
- [Scope and design principles](docs/scope.md), [DESIGN.md](DESIGN.md),
  [ROADMAP.md](ROADMAP.md)

## Related projects

- [FortSym](https://github.com/lazy-fortran/fortsym): symbolic generator of
  the kernels in `src/generated/`
- [FortNum](https://github.com/lazy-fortran/fortnum),
  [FortFEM](https://github.com/lazy-fortran/fortfem),
  [FortSparse](https://github.com/lazy-fortran/fortsparse),
  [FortPlot](https://github.com/lazy-fortran/fortplot),
  [FortIO](https://github.com/lazy-fortran/fortio): numerical, finite-element,
  sparse-solver, plotting and NetCDF/HDF5 libraries used by kin6d
- [libneo](https://github.com/itpplasma/libneo): MAKEGRID coil reader and
  Biot–Savart routines

## License

MIT; see [LICENSE](LICENSE).
