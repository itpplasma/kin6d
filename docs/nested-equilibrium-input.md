# Fixed-boundary nested equilibrium from a boundary and profiles

`nested_equilibrium --config FILE` solves the fixed-boundary nested 3D map
equilibrium (the D8 boundary/profile lane, [ideal-seed-3d](ideal-seed-3d.md))
for a boundary, a toroidal flux and prescribed profiles, and writes the solved
state's global quantities to `output/nested_equilibrium.nc`. It is the native
entry for cross-code comparison with VMEC-class codes; converting a VMEC
`&INDATA` deck into this namelist is the adapter's job, not KIN6D's.

## Input: namelist `nested_input`

| Name | Meaning | Default |
|---|---|---|
| `nfp` | field periods | 1 |
| `mpol`, `ntor` | largest poloidal and toroidal boundary mode | 1, 0 |
| `rbc(n,m)`, `rbs(n,m)`, `zbc(n,m)`, `zbs(n,m)` | boundary R, Z = sum C cos(m theta - n nfp phi) + S sin(m theta - n nfp phi), indexed (n, m) as in VMEC | 0 |
| `phiedge` | total toroidal flux (Wb), positive | 1 |
| `iota(k)`, `pressure(k)` | coefficients of s^k, s the normalized toroidal flux; pressure in Pa | 0 |
| `boundary_nodes`, `degree` | disc mesh of the final level | 16, 3 |
| `ladder_nodes(:)`, `ladder_degree(:)` | coarser levels solved first; each corrects the next seed (`ideal_map_prolong_correction`) | none |
| `angle_mpol` | poloidal modes of the boundary field-line angle; -1 pins the boundary angle | -1 |
| `tolerance`, `max_iterations` | free reduced-gradient reduction and Newton limit per level | 1e-8, 200 |
| `prescribe_current` | match the enclosed toroidal current `current(k)` (A, coefficients of s^k) instead of keeping `iota` | false |
| `max_current_iterations`, `current_tolerance` | iota updates and the largest current error relative to the largest prescribed current | 30, 1e-5 |

With a prescribed current, `iota` is the starting profile. For a fixed map the
current enclosed by each surface is affine in iota there, so iota is updated on
12 surfaces with secant slopes (the first from a step of 0.05), refitted as a
degree-8 polynomial and re-solved from the previous map.

The boundary must run counterclockwise in (R, Z) for increasing theta. A
boundary given in an arbitrary poloidal angle (every VMEC input) needs
`angle_mpol >= 0`: with the boundary angle pinned, the field lines on the wall
are tied to that angle and the solve converges to a different, constrained
state (on the circular test case the axis is then 8-20 mm off).

## Output

`volume_p`, `aminor_p` (sqrt of the mean cross-section area / pi), `rmajor_p`
(volume / (2 pi mean area)), `aspect`, `raxis_cc` (toroidal mean of the
O-point's major radius), `iotaf_axis`, `iotaf_edge`, `pressure_axis` (Pa),
`wb` (int B^2/2 dV / (2 pi)^2 with mu0 = 1, VMEC's normalisation), `betatotal`,
`force_residual` (free reduced-gradient norm relative to the final level's
uncorrected seed), `tolerance`, `iterations`, `ctor` (enclosed toroidal
current at the edge, (1/mu0) of B.dl along the discrete wall at phi = 0; inner
surfaces use the image of the reference circle), with a prescribed current also
`current_error` and `current_iterations`, and the global attribute `converged`. Volume integrals use collapsed Gauss-Legendre on each triangle and
the trapezoidal rule in phi; the mean area is int dV/R / (2 pi).

These are candidate values: no error bound is reported.

## Checks

`equilibrium_nested_circular-tokamak` solves
`examples/nested-equilibrium/circular-tokamak.nml` (R0 = 3 m, a = 1 m,
iota = 0.4 + 0.2 s, p = 0) and requires the exact volume 2 pi^2 R0 a^2 to 1e-4,
the axis within 1.5 mm of 3.0366 m and `wb` within 2e-5 of 0.0750949; PARVMEC
(ns = 101, mpol = 8) gives 3.03634 m and 0.0750949, GVEC gives 3.03683 m.
`equilibrium_nested_circular-tokamak-current` prescribes I(s) = 340876 A * s on
the same boundary, starting from iota = 0.3; PARVMEC with CURTOR = -340876 A and
AC = 1 (its current sign convention is opposite) gives raxis 3.04266 m and wb
0.0755050, KIN6D 3.04281 m and 0.0755056 after 4 iota updates.

