# Circular CHEASE fixture: stored and reconstructed data

This synthetic circular tokamak (aspect ratio 10, nominal q-axis 1.55) entered
public [KIN6D](https://github.com/itpplasma/kin6d) in
71bccfa289d128c79b351f88f546c4e55f913350. That revision contains only the four
profile arrays and a 65 by 65 poloidal-flux grid (902 lines). Its title also
extends into the G-EQDSK integer fields. It contains no q, LCFS or limiter
records. No private source file was used to complete it.

The current file preserves all 901 original numeric records byte for byte.
Only the title/layout is corrected, followed by a reproducible trailer:

- q on uniform normalized flux s = 0..1: the axis value is
  F/(R_axis sqrt(det Hessian(psi))), evaluated with `gs_axis_safety_factor`.
  For s > 0 it is minus the geometric signed q from `surface_build`, with 512
  poloidal samples. This is positive source COCOS 1 q (q-axis 1.5499435744,
  q-edge 3.1006248991). These are reconstructed field integrals, not retained
  original CHEASE q values and not an independent cross-code q oracle.
- LCFS: 128 counter-clockwise rays about the spline-recovered magnetic axis,
  intersected with psi = psi_edge using `flux_ray_point`. No limiter is known;
  its count is zero. The contour is reconstructed from the stored grid.

`complete_shield_torus_fixture.f90` regenerates the header and trailer through
libneo's writer, preserving the original numeric prefix. With the optional
adapter enabled:

```sh
cmake --build --preset cpu --target complete_shield_torus_fixture
build/complete_shield_torus_fixture test/fixtures/shield_torus_a10.eqdsk completed.eqdsk
```

The source declares COCOS 1; libneo estimates its GS prefactor as 0.9994255092
and classifies it as 1 without profile repairs. Import standardizes to 3,
reversing flux, derivative profiles and q. The legacy-reader equivalence gate
compares physical fields, F and p after this convention change.

The stored edge pressure is -3.4873077e-20 Pa. The native re-solve explicitly
chooses zero edge pressure and records the 3.4873077e-20 Pa offset; derivative
profiles and original imported data are preserved. The native comparison uses
three P2 meshes and reports axis, flux and q differences alongside 65/33-grid
interpolation estimates. Exceeding those estimates is printed as a finding,
not hidden by changing tolerances. These estimates are not rigorous bounds.

The final empty record represents the zero-length limiter array emitted by
libneo. It is required by its formatted reader and intentionally retained.
