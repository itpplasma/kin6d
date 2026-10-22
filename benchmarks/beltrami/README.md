# Synthetic Taylor/SPEC comparison

These public inputs describe a single fixed-boundary plasma volume with
constant pressure zero, mu=1.1, toroidal flux 1, R0=6 and minor radius 1.
`axis.sp` is axisymmetric. `helical.sp` has nfp=2 and

    R = 6 + cos(theta) + 0.1 cos(theta - 2 phi),
    Z =     sin(theta) - 0.1 sin(theta - 2 phi).

The FEEC map is the corresponding affine deformation of each reference
section. The magnetic 2-form is pushed forward with its contravariant Piola
map. No nested magnetic surfaces are assumed in the solve.

From the KIN6D worktree, with `out` on a storage filesystem:

```sh
out=/mnt/storage/kin6d-v2beltrami
cmake -S ~/code/SPEC -B "$out/spec-build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS=-std=legacy
cmake --build "$out/spec-build" --target xspec -j 2
python3 benchmarks/beltrami/run_spec.py \
  "$out/spec-build/build/bin/xspec" "$out/spec-ladder"
cmake --preset cpu -B "$out/build" -DKIN6D_BUILD_EXTRAS=OFF \
  -DKIN6D_BELTRAMI_SPEC_SAMPLES="$out/spec-ladder"
cmake --build "$out/build" --target test_nonnested_beltrami_spec \
  test_nonnested_beltrami test_nonnested_beltrami_pairing -j 4
ctest --test-dir "$out/build" -j 4 --output-on-failure \
  -R '^nonnested_beltrami_'
python3 benchmarks/beltrami/analyze.py "$out/spec-ladder"
```

SPEC was tested at b23a2408fc3f1dfbc14896e16d73a57f7d1d9e88. GCC 16 needs
`-std=legacy` for SPEC's existing FORMAT descriptors. No SPEC source changes
are needed. The run script records the executable SHA256, retains native
HDF5 files and logs, and uses one OpenMP thread for each SPEC run. The optional
Python tools only run the external code, extract its native grid and reduce
sample arrays. All FEEC solves, field evaluations and residuals are Fortran.

The SPEC ladder uses `(Mpol,Lrad)=(4,8),(6,10),(8,12)`, with Ntor=0 for the
axisymmetric case and Ntor=Mpol/2 for the helical case. The diagnostic grid is fixed at 192 poloidal and 96 toroidal points (one
toroidal point for the axisymmetric case). SPEC Ndiscrete is 12/8/6 on
the three resolutions (SPEC `src/preset.f90:1292`). The sampled angular indices and
radial indices 3/7/11 (r=3/16, 7/16, 11/16) are fixed across resolutions
and avoid FEEC facets. The native comparison checks this before assembly.
There are 24 axisymmetric and 192 helical sample points.
`extract_spec_grid.py` multiplies the native grid's B^phi by R to obtain the
physical B_phi; see SPEC `src/sphdf5.f90:662-666` at the revision above.

The FEEC P2 mesh ladder has 2/4/8 radial rings, and the helical mode ladder
retains n=0:2, 0:4 and 0:6 at four rings. Reported component errors use the
common full-vector RMS scale so nearly zero components remain meaningful.
The magnitude error is separate. The sampled element-curl residual uses
finite differences wholly inside each selected element. It does not include
facet jumps. SPEC's native `Lcheck=1` residual is retained in each run log;
its curvilinear component norms differ from the FEEC element-curl norm.

`analyze.py` reports the reference, mesh and mode refinement differences
separately. They are observed uncertainty diagnostics, not rigorous error
bounds or additive certified bars. They do not establish arbitrary-geometry
P3 convergence. The P2/P3 manufactured and Lundquist order tests instead use
an exact affine square section, as specified in
[the closure contract](../../docs/equilibrium-beltrami.md).
