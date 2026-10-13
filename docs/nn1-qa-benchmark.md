# NN1 mapped 3D: precise-QA vacuum benchmark

Status: candidate. Program `test/test_nonnested_qa.f90`, data and checks in
`benchmarks/nn1_qa/`, Slurm script `scripts/test_nonnested_qa_scluster.sbatch`.
Solver: [nn1-mapped-3d.md](nn1-mapped-3d.md).

## Sources

- Boundary: `tests/test_files/input.LandremanPaul2021_QA` (VMEC input, nfp 2, 61
  `RBC/ZBS` modes, m <= 5, |n| <= 5) of simsopt, https://github.com/hiddenSymmetries/simsopt,
  commit 81436b4f7c65 (master, 2026-10-13; release v1.11.1, MIT). Surface of
  Landreman and Paul, PRL 128, 035001 (2022), arXiv:2108.03711.
- Coils: `examples/1_Simple/inputs/biot_savart_opt.json`, same commit: 16 coils
  (4 base Fourier curves of order 5, nfp 2, stellarator symmetry), 1e5 A each,
  written by `examples/2_Intermediate/stage_two_optimization.py` for this
  boundary. A stage-2 fit, not the coils of the PRL. `qa_data.py` converts both
  files to `qa_boundary.dat` and `qa_coils.dat` (192-point polygons).

## Oracle

Biot-Savart of the polygonal filaments, exact per segment: `B` and the vector
potential `A` (checked: `|curl A - B|/|B| ~ 1e-10`, net flux `F = -0.32 T m`
`= mu0 I/(2 pi)` equal at two circles to 12 digits, Fortran and numpy agree to
15 digits at three points). The prescribed Neumann data are the tangential traces of
`A`; the interior field is compared with `B`. The vacuum residual of the published
surface with these coils is `<|B.n|>/<|B|> = 4.0e-4` (max 1.7e-3), `<|B|> = 0.299 T`,
minimum coil-surface distance 0.304 m. The surface is therefore not exactly a coil
vacuum surface, which does not matter here: the coil field is an exact vacuum field
inside any surface that excludes the coils.

## Wall map

The published boundary is represented exactly (no truncation) as the image of a
disk (centre R = 1, radius 0.2). With `z = ((R-1) + iZ)/0.2` and
`w = z^m exp(-i nfp n phi)`,

    R = sum (RBC - b (1-|z|^2) m ZBS) Re w,   Z = sum (ZBS - b (1-|z|^2) m RBC) Im w.

For `b = 0` this is the harmonic polynomial extension of the VMEC parametrisation,
which folds at phi = 0, theta = pi (min det = -0.0085, inboard dent). The term
vanishes on the boundary and rotates the interior normal derivative; `b = 0.1` gives
det in [0.0071, 0.0962] (`qa_check.py`). The arclength reparametrisation also removes
the fold but needs more than 16 toroidal harmonics; it was rejected. Map degree in
phi: 10.

## Rule (declared before the first solve)

E = ||B_h - B||_L2 / ||B - F/R e_phi||_L2 over the mapped domain.

1. div B <= 1e-12 and FGMRES residual <= 1e-9 (hard stops).
2. Mode ladder at nr = 4: E(N = 12) <= 0.5 E(4) and |E(12) - E(10)| <= 0.2 E(10).
3. h ladder at N = 10, nr = 2, 4, 8, 16: observed order >= 1.2 for every pair and
   E(nr = 16) <= 0.05. The target 2 is reported, not gated.

## Results

See the results section appended after the Slurm job.
