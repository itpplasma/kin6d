# Energy-criterion island seeding

Status: implemented, opt-in (2026-10-28). Item 4 of the
[relaxation comparators](relaxation-comparators.md#algorithms-to-adopt).
Module `src/mhd/island_seed.f90` (`kin6d_island_seed`), generated leaves
`src/generated/kin6d_island_seed_leaves.f90` from `gen/gen_island_seed.f90`,
test `test/test_island_seed.f90`. Nothing calls it by default.

## What it does

`island_seed(field, options, report, potentials)` acts on an NN1 state
([free-B lane](equilibrium-nonnested.md)): an axisymmetric background in mode
0 plus any n > 0 modes, for example a converged nested state with shielding
currents. For every coprime (m, n) with 1 <= n <= `n_max` (default: all modes
of the field) and 1 <= |m| <= `m_max` it finds the circles where the background
rotational transform equals n/m. For each one it builds two seed potentials

    dA = alpha((r - r_mn)/w) B0/|B0| {cos, sin}(m theta - n phi),
    alpha(xi) = (1 - xi^2)^2,  w = min(1/|m|, r_mn - r_inner, r_outer - r_mn),

with r in units of the radial scale. Each potential is projected into V1 with
zero boundary trace and moved to the R A_phi = 0 gauge by subtracting
grad_n(a/(i n)). This leaves its curl unchanged exactly (the measured gauge
defect is 0). The magnetic energy W = int |B|^2/(2 mu0) is exactly quadratic
in the seed amplitudes a:

    W(B + sum a_i curl dA_i) = W(B) + a.g + a.G a/2,
    g_i = <J, dA_i> = (curl dA_i, B)/mu0,   G_ij = (curl dA_i, curl dA_j)/mu0,

and the seed is a* = -G^{-1} g (dense Cholesky; G couples seeds only within
the same n). The state's V2 modes and, if given, the caller's gauge
potentials are updated in place, so `nn_solve` or `relax_run` can continue
from the seeded state. `island_resonances` returns the resonances without
seeding.

For each chain the report gives r_mn, iota, d iota/dr, the support, the
resonant radial-field amplitude dB^r before and after seeding, the O-point
angle at phi = 0, and the pendulum width

    w_mn = 4 sqrt(R dB^r / (|m| |d iota/dr| B_phi)),

where R is the metric radius (R0 of the axis, or the cylinder radius) and
B_phi is averaged on the circle. This is the separatrix width of
H = m iota' x^2/2 - (R dB^r/B_phi) sin(chi + phase), chi = m theta - n phi,
a large-aspect-ratio, single-harmonic estimate.

## Numerical choices

- iota = 1/<dphi/dtheta> and dB^r are computed on geometric circles about
  `field%centre`. The background must therefore have nearly circular nested
  surfaces; the largest background B_r/B_theta on each resonant circle is
  reported as `circularity`.
- The mean over a circle uses Gauss-Legendre rules on arcs split where the
  circle crosses element edges. Represented fields are only piecewise smooth
  along a circle, and a uniform rule over their jumps converges at first
  order: with uniform sampling, dB^r of the seeded state still changed by 4%
  between 512 and 8192 samples. With the split rule it agrees to all printed
  digits for 256, 512 and 1024 scan samples.
- d iota/dr is a least-squares slope over 17 points in a window of at least
  two ring intervals. The represented iota oscillates on the element scale:
  a central difference with a 0.02 step was 8-12% off on 8 rings.
- Resonances inside half the first ring of a disc are not seeded, because
  the axis cells do not resolve iota.
- A purely axisymmetric state gives g = 0 for every n > 0 seed by Fourier
  orthogonality, so it is never seeded. Seeding acts only on the resonant
  n > 0 content of the state.

The pointwise kernels (`circle_rotation`, `seed_potential`, `pendulum_width`)
are FortSym leaves with JVP/VJP. They were regenerated with the pinned FortSym
8efec5b, built here against FLINT 3.6.0 and the system MPFR 4.2.1. The same
build reproduces `kin6d_nonnested_leaves.f90` and
`kin6d_metriplectic_leaves.f90` byte for byte. Sensitivities of the seed with
respect to the state (through G^{-1}) are not exposed yet.

## Verification case

The case is a periodic cylinder with R0 = 3, radius 1, B_z = 1, P2 FEEC and
P2 geometry, modes 0..1. The prescribed nested background has
iota(r) = 0.4 + (0.1/0.5625) r^2, so the only resonance with m <= 4 is
q = 2 at r = 3/4, which is a mesh ring. The n = 1 boundary trace has
B = grad psi x e_z with psi = (eps/2) (x - i y)^2 and eps = 1e-4. The relaxed
state is the vacuum mode (`nn_solve`), which has a q = 2 island. The
independent "converged nested state" is the shielded minimiser: the least
energy with the same trace and zero resonant B^r on r = 3/4, obtained with
one Lagrange multiplier. Its shielding sheet lies on element edges.

`island_seed_seed` (16 rings; 24 s on a 4-core VM, ladder 26 s, none 4 s):

| Check | Result |
|---|---|
| resonances found | one, (m, n) = (2, 1), r = 0.7519 (exact 0.75), iota' = 0.2670 (exact 0.2667) |
| G symmetry defect / min eigenvalue ratio | 6.2e-18 / 1.0 |
| internal model vs recomputed energy change | -3.670526035363349e-7 / -3.670526035363444e-7 |
| quadratic model vs independent `relax_state` energy at 0, 1, 2 x a* | defect 3.7e-15 absolute, 1.2e-16 relative to W |
| energy ordering | W(shielded) 30.053959143 > W(seeded) 30.053958776 > W(vacuum) 30.053958612 |
| dB^r shielded / seeded / vacuum | 2.1e-6 / 1.044e-4 / 1.491e-4 (seeded = 0.70 of vacuum) |
| pendulum vs traced width, eps | 0.09688 vs 0.10051 (ratio 1.037) |
| pendulum vs traced width, eps/4 | 0.04844 vs 0.04962 (ratio 1.024) |
| width(eps/4)/width(eps) | 0.500000000000000 |
| Newton hand-off (`nn_solve`) | returns the vacuum to 1.9e-11, energy decreases |
| ideal metriplectic hand-off, 2 steps | status 0, energy non-increasing |

The shielded state's residual dB^r of 2.1e-6 is the offset between the
shielding circle r = 3/4 and the represented resonance r = 0.7519. The traced
width is the radial extent, on the phi = 0 O-point ray, of field lines whose
helical angle librates over 160 transits (RK4, 32 steps per transit), with
16 bisections per edge. On 8 rings the field-line tracing alone differs by
6% from the pendulum width of the smooth vacuum field. On 16 rings that
difference is 0.4%, so the width comparison uses 16 rings. The ratio
approaches 1 as eps decreases, as expected of an estimate that is
asymptotic in the island width. The ideal hand-off also relaxes the
background: the prescribed iota profile carries a current that is not
force-free.

`island_seed_ladder`, the predicted width under refinement (the test gate
uses 8/16/32 rings; the 4-ring level is from a separate run):

| rings | r_mn | iota' | G_11 | dB^r seeded | seeded/vacuum | width |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.76223 | 0.2754 | 554.1 | 1.2217e-4 | 0.812 | 0.103181 |
| 8 | 0.75565 | 0.2682 | 511.5 | 1.0918e-4 | 0.730 | 0.098846 |
| 16 | 0.75186 | 0.2670 | 503.3 | 1.0440e-4 | 0.700 | 0.096882 |
| 32 | 0.75049 | 0.2668 | 500.7 | 1.0361e-4 | 0.696 | 0.096533 |

The successive width differences are 4.3e-3, 2.0e-3 and 3.5e-4, giving
observed orders of 1.1 and 2.5. The first version shielded on the
represented resonance circle, which cuts through elements. Its gradient g did
not converge (1.632e-2, 1.611e-2, 1.580e-2 on 8/16/32 rings) because the
shielding sheet was not resolved. Aligning the sheet with a ring fixed this:
g is then 1.939e-2, 1.922e-2, 1.921e-2.

`island_seed_none` covers two cases. With iota in [0.55, 0.75], no n/m with
m <= 4 lies in range: no seed, and the state and potentials are bitwise
unchanged. With a resonant but axisymmetric state, g = 0 and a* = 0 exactly.
`island_seed_leaves` checks the leaf identities, the JVP against central
differences (2e-10 or better) and the JVP/VJP pairing (4e-16 or better).
`island_seed_exact_regeneration` checks the generated leaves against
the generator when `KIN6D_FORTSYM_BUILD` is set.

## Limits

- Model: the magnetic energy only, as in the NN1 vacuum restriction. No
  pressure term, and no helicity or flux constraint enters the quadratic
  model (the seed has zero boundary trace, so fluxes through the boundary are
  unchanged).
- The seed space has two functions per chain. The seeded resonant field is a
  fraction of the relaxed one (0.70 here), which is consistent with
  W(seeded) >= W(vacuum).
- Geometric circles and theta, not straight-field-line coordinates. Shaped
  or strongly shifted backgrounds need flux-surface coordinates first.
- Overlapping chains (Chirikov) are not detected. Widths are reported, but
  nothing limits the amplitudes.
- No certificate: all numbers are candidate diagnostics.
