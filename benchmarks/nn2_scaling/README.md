# NN2 scaling measurements, 2026-10-22

[V] These are derived summaries of the existing synthetic helical-coil island
fixture in [KIN6D](https://github.com/itpplasma/kin6d). No machine equilibrium
or non-public coil data is included. Methods, interpretation and failed gates
are in [the scaling report](../../docs/equilibrium-nonnested-scaling.md).

## Provenance

The development base is ff665b0 on `islands/v3scale`. The dependency pins are:

- FortFEM: 7e3e2d281887081c7280b5c4fafc050be7c633b1, adding only the identity forms to the original pin.
- FortNum: 1cfe3cba327e9a32743ed10292e75a1311a6ced4.
- FortSparse: 4f716f68152a1db44e4b7c58521b52535d6585fa.

Local builds use GNU Fortran 16.2.1, RelWithDebInfo (`-O2 -g`), OpenMP,
`-fno-fast-math -ffp-contract=off`, and an AMD Ryzen 9 5950X. Each run sets
OPENBLAS_NUM_THREADS=1 and at most eight solver threads, unsets DISPLAY and
limits virtual memory to 25165824 KiB. Cluster builds use the guarded native
CMake preset, the installed conda Fortran compiler and system C compiler.
The cluster measurements below ran on the same node, with eight solver threads.
Local and cluster timings are not interchangeable.

| Summary | Raw directory under `/mnt/storage/kin6d-v3scale` | Job / reservation |
|---|---|---|
| local-dense.tsv | local-dense | local, 8 threads |
| local-threads.tsv | local-threads | local, 1/2/4/8 threads |
| local-krylov.tsv | local-krylov | local, 8 threads |
| scluster-dense.tsv | scluster-dense | 3354209, 16 cores, 188.8 GiB guard |
| scluster-krylov.tsv | scluster-krylov | 3354257, 8 cores, 94.4 GiB guard |
| scluster-dense-highchi.tsv | scluster-dense-supplement | 3354343, 16 cores, 188.8 GiB guard |

Job 3354209 ends FAILED at the pre-existing P3/nr=16 vacuum seed rejection,
after nine successful coupled solves. Job 3354343 attempts the remaining
ratio-100 cells and every ratio-1e6 cell, preserving per-case exits. Job
3354257 completes the full Krylov scan and four high-anisotropy equivalence
attempts; the individual failures remain failures even though the batch job
completes. Their raw logs are in `scluster-equivalence`.
FortFEM's 69/69 focused Debug checks pass in job 3354199.

Remote originals are under `/home/ert/runs/kin6d-v3scale-20261022` on sCluster.
Immutable source exports there are `source-mf14` for the first dense sweep
and `source-mf13` for the final-tolerance sweep. The former's tighter
matrix-free tolerance does not affect its dense-path measurements. SHA256
digests of their complete source manifests are respectively:

```text
dd7b3f0f3655dad01b80dcd54ff9a29ab0d32c60eec68780f24901a5814e77cd
f9acc3f01d5b45bf6ef720260b33f586040a1d8fecdd58f9d63bb4e18e731bc7
```

The local timing binaries, before later trace-helper and formatting edits,
have SHA256 digests:

```text
a07ee4bdaec7c1d146fceef977a3c69a2d5a62e8d57ae1bc811b92021593682c  test_nonnested_island
070ed6c3e5ff20f6ca8fb8b9d96ee56f2f96cb083124316936de0871b1a6bbdc  test_nonnested_coupled
```

Raw logs, source manifests and the gprofng experiment are checksummed in
`/mnt/storage/kin6d-v3scale/measurements.sha256`; its SHA256 is
`66e2d3fc1303d3bbade8f0bf056c1d2aed3e740d45d202587b24bb0fddabd904`. Campaign outputs stay outside Git.

## Reading the TSVs

[V] `seconds` sums same-named wall-time samples across the solve, including
the terminal Jacobian rebuild. Parent and child phases overlap.
`max_end_rss_KiB` is the largest RSS sampled at that phase's end, while
`max_process_hwm_KiB` is cumulative process HWM, not an isolated phase peak.
`process-total` uses `/usr/bin/time -v`. Zero RSS there means no phase sample,
not zero process memory. Failed setup has solve_status `unavailable`.

The Krylov arrays include the terminal Newton slot, so a trailing zero is
not a zero-iteration linear solve. Early failures print only attempted slots.
Chi scans in `*-krylov.tsv` hold g_parallel/g_perp=100. Dense and thread
scans use chi=g=100; `chi1e6-*` dense cases use chi=g=1e6.

## Reproduce

Run from the source checkout, choosing fresh build and result paths under
`/mnt/storage`. The following limits local work to eight threads and 24 GiB:

```sh
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
unset DISPLAY
ulimit -v 25165824
cmake --preset dev-fast -B /mnt/storage/kin6d-v3-repro/build \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build /mnt/storage/kin6d-v3-repro/build -j8 --target \
  test_nonnested_island test_nonnested_coupled test_nonnested_leaves
scripts/study_nonnested_scaling.sh /mnt/storage/kin6d-v3-repro/build \
  /mnt/storage/kin6d-v3-repro/local local
scripts/summarize_nonnested_scaling.py /mnt/storage/kin6d-v3-repro/local
ctest --test-dir /mnt/storage/kin6d-v3-repro/build -j8 --output-on-failure \
  -R '^(nonnested_matrixfree_.*|nonnested_symbolic_derivatives_and_d3_identity|nonnested_coupled_fourier_reduction_p[23]|kin6d_layering)$'
```

The focused CTest command is expected to retain the documented P3/N=1
failure. Use the driver's `threads` scope for strong scaling and `krylov`
scope for all N and chi ratios. For large dense cases, export the source
to sCluster and submit `scripts/study_nonnested_scaling_scluster.sbatch`
with arguments SOURCE RUN MEM_GUARD; do not run those cases locally.

CPU sampling, separate from wall-time benchmarking:

```sh
gprofng collect app -p 1 -o /mnt/storage/kin6d-v3-repro/dense.er \
  /mnt/storage/kin6d-v3-repro/build/test_nonnested_island profile 2 8 2 1e2 dense
gprofng display text -functions /mnt/storage/kin6d-v3-repro/dense.er
gprofng display text -lines /mnt/storage/kin6d-v3-repro/dense.er
```
