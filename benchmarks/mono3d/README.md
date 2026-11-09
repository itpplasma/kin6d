# Monoenergetic 3D coefficients on published fields

Recipes behind the MONKES and rabe comparisons of `mono_coefficients_3d` in
[drift-kinetics.md](../../docs/drift-kinetics.md) (K1). The field files stay
outside the repository:

- W7-X standard configuration (EIM): `boozmn_wout_w7x_std_v_0.0.txt.nc` and
  the MONKES database of `examples/W7X-EIM/boozer_xform_input` in MONKES
  (Escoto et al., Nucl. Fusion 64, 076030, 2024), commit add6c3f;
- Landreman-Paul (2022) reactor-scale QH equilibrium
  `wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc` and rabe's
  golden output `test/golden/expected/rabe.nc` in rabe, commit 3282d65.

Steps:

```sh
python boozer_surfaces.py WOUT qh.json 0.25510204,0.5 48 24  # booz_xform 0.1.0
python surface_grid.py qh.json 1 41 51 qh_41x51.txt          # b on theta x zeta
python trapped_fraction.py qh_41x51.txt                       # f_t of that surface
gfortran -O3 -fopenmp -I build/modules mono_scan.f90 -o mono_scan \
    build/libmilca_continuum.a build/libmilca_core.a build/libmilca_fft.a \
    <fortnum libraries> -lopenblas
./mono_scan qh_41x51.txt 240 13.635429024906703 -0.0735110806949427 1e-2 1e-3
```

`mono_scan GRID N_LEG R D31_LC NU*...` prints D31 at nu = nu*/(pi R) and its
ratio to D31_LC = f_t iota/(3 (iota + 4)), the Landreman-Catto value in the
normalisation of `dk_mono`; R is rabe's major radius. For MONKES, MONKES
writes the modes it used (after its own truncation) to
`monkes_Magnetic_configuration.dat`; the grid is built from those modes with
B = sum B_mn cos(m theta + n N zeta) and its iota, B_theta, B_zeta, and the
coefficients compare as D11_M = 2 a^2 D11, D31_M = -2 a D31, D33_M = 2 D33,
a = G/(iota B00 psi').

Results (faepmac1, 8 threads):

| case | grid | Legendre | wall | memory |
|---|---|---|---|---|
| QH s = 0.5, nu* = 1e-2 to 1e-4 | 33 x 41 | 240 | 294 s | 7.1 GB |
| QH s = 0.5, nu* = 1e-3 | 41 x 51 | 240 | 314 s | 14 GB |
| W7-X s = 0.2, nu = 1e-3 to 1/m | 25 x 25 | 160 | 142 s (workstation) | 1.0 GB |
