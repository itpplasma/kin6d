# Circular CHEASE fixture: stored and reconstructed data

This synthetic circular tokamak (aspect ratio 10, nominal q-axis 1.55) entered
public [MILCA](https://github.com/itpplasma/milca) in
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

The completed fixture is interpreted as COCOS 1; the original source version
and convention were not retained. Its CHEASE title alone does not establish
them. Libneo estimates its GS prefactor as 0.9994255092
and classifies it as 1 without profile repairs. Import standardizes to 3,
reversing flux, derivative profiles and q. The legacy-reader equivalence gate
compares physical fields, F and p after this convention change.

The stored edge pressure is -3.4873077e-20 Pa. The native re-solve explicitly
chooses zero edge pressure and records the 3.4873077e-20 Pa offset; derivative
profiles and original imported data are preserved. The native comparison uses
three P2 meshes and reports axis, flux and q differences alongside 65/33-grid
interpolation estimates. The smooth fixture explicitly selects
`curved_lcfs=.true.` to retain the periodic cubic LCFS curve. Interior mesh
refinement alone does not control the boundary or derivative-profile
reconstruction errors. Exceeding a grid estimate remains a printed finding;
the interior q regression also requires 1% agreement. These comparisons do
not establish an original CHEASE q oracle or a rigorous bound.

Conditioning of the prescribed-profile re-solve. With p'(psi) and FF'(psi)
prescribed, this fixture lies close to a fold of the native problem:
scaling both profiles by 1-1e-4 changes the native axis flux by 0.90%
(amplification 90 on the level-2 mesh, CTest `geqdsk_chease_conditioning`),
and scaling them by 1+1e-4 gave no converged solution in a one-off run.
The source itself is consistent: a fourth-order finite-difference GS residual
of the stored grid with the stored profiles (normalized flux below 0.9) has
relative rms 4.4e-5, and finite differences of F^2/2 and p reproduce FF' and
p' to about 1e-4. A one-off scan also showed that resampling the LCFS with
128 or 512 rays from the stored grid leaves the mismatch unchanged. The
native/source axis-flux mismatch (0.39% at level 2, 0.49% at level 3) and
the 0.6% interior q difference are therefore consistent with that data
consistency times the measured amplification. The conditioning test gates
the level-2 axis-flux mismatch against amplification times 1e-4. This is a
sensitivity explanation, not an enclosure; a well-posed qualification needs
a current- or q-constrained native re-solve.

The final empty record represents the zero-length limiter array emitted by
libneo. It is required by its formatted reader and intentionally retained.
