# Synthetic nfp = 3 free boundary: KIN6D against VMEC++

A public synthetic case for the 3D free boundary on the curve chart
(decision 19, step a to b). No device data are involved.

- `make_coils.py` writes 24 planar circular modular coils (8 per period,
  radius 0.30 m) centred on the helical axis R = 1 + 0.15 cos 3phi,
  Z = 0.15 sin 3phi, one circuit, MAKEGRID format, |B| about 1 T.
- `run_vmecpp.py` runs VMEC++ free boundary at beta = 0 and zero net current
  on the makegrid response table of the same coils file, and writes `wout.nc`,
  `iota.dat` and `timing.json`.
- `fb3d_synthetic` (app/fb3d_synthetic.f90) reads the same coils file through
  libneo, builds the helical curve chart, the coupling surface traces, the
  field-period R2 exterior and solves the free-boundary residual. The frozen
  iota profile is VMEC++'s `iota.dat` and the toroidal flux is phiedge.
- `compare.py` measures the distance between the two boundaries and |B| at the
  foot points. `vacuum_reference.py` is the independent beta = 0 oracle: the
  coil-field flux surface that encloses phiedge, from Biot-Savart field-line
  tracing, compared on the phi = 0 symmetry plane with both codes.

```sh
python make_coils.py coils.fb3d_synth
python run_vmecpp.py coils.fb3d_synth vmec --nzeta 96 --grid 121 121 96 --mpol 8 --ntor 6 --ns 16 32 64
fb3d_synthetic coils.fb3d_synth vmec/iota.dat 0.0314159 k 12 16 2 12 18 384 1 200
python compare.py k vmec/wout.nc
python vacuum_reference.py 0.0314159 --kin6d k --wout vmec/wout.nc
```

The results and their limits are recorded in the pull request, not here.
