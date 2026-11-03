# GNU Fortran host-associated automatic extents

`gfortran_host_extent.f90` isolates the compiler defect behind the spectral
matcher's automatic-array failure. Compile it separately from the test suite:

```sh
gfortran -O3 -g -std=f2018 -pedantic -fcheck=all -fbacktrace \
  -finit-real=snan -finit-integer=-99999999 -fdump-tree-original \
  gfortran_host_extent.f90 -o host_extent
./host_extent
```

The required extents are `3 3`. GNU Fortran 16.2.0 on macOS arm64 instead
reports `3 0` and rejects `candidates(1)`. The original tree dump evaluates
`size(oracle)` using the compiler's upper-bound temporary before assigning
that temporary from `n`. The same failure occurs at `-O0`; the larger matcher
exposes different results as optimization changes the uninitialized extent.
Removing the internal function, or declaring `candidates(n)`, restores both
extents to three at both optimization levels.

The matcher keeps allocatable scratch arrays and allocates `oracle` before
inquiring its size. Its analytic folding and coarse-frequency tests pass with
identical output at `-O0` and `-O3`, with runtime checks and with
`-fsanitize=address,undefined`. This workaround does not change branch
isolation or multiplicity rules.
