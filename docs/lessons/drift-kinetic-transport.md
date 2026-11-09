# Lessons: drift-kinetic transport

Physics area: neoclassical toroidal viscosity (NTV), the torque that a
non-axisymmetric perturbation delta B of a tokamak field exerts on each
species, with its heat flux and offset rotation; the axisymmetric transport
matrix (bootstrap current, parallel flows, radial electric field E_r from the
force balance); and finite orbit width on the canonical Maxwellian. The
equations and switches are in [drift-kinetics.md](../drift-kinetics.md) and
[drift-kinetics-extensions.md](../drift-kinetics-extensions.md).

References: [NEO-2](https://github.com/itpplasma/NEO-2) in its quasilinear
tokamak version NEO-2-QL for the collisional NTV and axisymmetric transport,
and [NEO-RT](https://github.com/itpplasma/NEO-RT) for collisionless resonant
transport with thin orbits; the POTATO orbit code in the NEO-RT repository for
exact guiding-centre orbits. All three are MIT-licensed. Devices: the ITER
scenario TC24 and ASDEX Upgrade (AUG) discharge 30835, whose Boozer files and
profiles stay outside the repository.

Notation: s is the normalised toroidal flux, theta and phi the Boozer angles,
b = B/B_0 with B_0 the (0, 0) harmonic, xi = v_par/v the pitch, nu the
deflection frequency, q = 1/iota, G and I the covariant Boozer field
components, Omega_tE the E x B rotation frequency. The *kernel* is MILCA's
monoenergetic NTV solve at one speed: the Legendre-Fourier kernel
(`milca_dk_ntv`), its spectral-element counterpart on the orbit-class mesh
(`milca_dk_hp`, "hp"), and the finite-width kernel, which runs the hp solver
on exact orbits (`milca_ntv_orbit_geometry`). The *thin* kernel neglects the
orbit width. f_10 is the axisymmetric (n = 0) solution that the perturbed
operator acts on. Test names in backticks are ctest names.

Outcome: on ITER TC24, torque, heat flux, flows, bootstrap current and E_r
agree with NEO-2-QL within 0.2 to 5 % once the causes below were removed. On
AUG 30835 the electron torque agrees within 2 to 6 % and the ion torque summed
over the profile within 4 %, while the ion torque on single flux surfaces near
energy resonances still differs by up to 70 %. Of the 30 causes below, six
were wrong or missing terms in MILCA, seven were conventions or switches of
the reference, five were resolution, three were checks against the wrong
limit, three were input data and six were implementation defects.

## 1. Wrong formula

- **Mirror source differentiated along the wrong direction** (#264). NEO-2-QL
  solves two quasilinear source terms besides the toroidal-derivative drive:
  the rest of the perturbed radial drift (the *drift source*) and the
  perturbed mirror and collision operators acting on f_10 (the *f_10
  source*). The mirror part of the f_10 source took i m, the theta derivative
  of delta B, where the equation needs the derivative of delta B/B_0 along the
  field line at fixed field-line label, i (m - n q). With it, the two sources
  removed 49 to 86 % of the torque on a ladder of circular surfaces, although
  NEO-2-QL, which contains both sources, agreed with the source-free kernel
  to 2 %. Found by an independent Galerkin transcription of the source
  equation with its own n = 0 solve, and by the 1/nu limit, where the two
  sources must add up to a parallel derivative. The earlier identity test of
  #262 had only checked the implementation against itself.

## 2. Missing terms

- **Change of the drift factor in the drift source** (#264): the term
  -2 (delta B/B) d_theta ln B was missing; dropping it from the oracle changes
  the torque by 3e-3 relative.
- **Surface measure in the torque of the hp kernel** (#322). In Boozer
  angles b.grad theta = k b, and a flux-surface average weights with
  1/(b b.grad theta) = 1/(k b^2). The hp kernel already carried 1/(k b) in its
  phase-space measure but lacked the further 1/b in the torque weight: 2 %
  against the Legendre kernel. It was hidden inside the AUG finite-width
  comparison until constant and Boozer streaming were compared separately.
  `dk_hp_ntv_reference` now checks Boozer streaming to 2e-9 relative.
- **Boundary term of a moving quadrature split** (#237). The collision
  operator between species a and b splits its speed quadrature at
  v = v_th,b, which moves with the temperature of b. The temperature
  derivative of that block missed the boundary term of the moving split and
  was wrong by O(1) between species. Found against centred differences.
- **Current and shear term of the precession** (#318). On a Boozer file with
  magnetic shear, the exact orbit precession equals the thin bounce average
  only with NEO-2-QL's current and shear drift term; without it the
  precession is off by up to a factor 2. `boozer_orbit_shear` checks it.
- **f_10 counterpart of the canonical drive** (open). The finite-width kernel
  drives with the canonical Maxwellian. To first order in the orbit width this
  drive differs from the local one by exactly the sum of the two sources of
  section 1, so the finite-width kernel reproduces the thin kernel with both
  sources off, which was checked (section 5). The collisional f_10 source has
  no finite-width counterpart yet.

## 3. Different model or convention in the reference

- **Drift weight**: NEO-2-QL and NEO-RT weight drift and drive with
  (1 + xi^2)/b, grad-B plus low-beta curvature drift. MILCA's first NTV test
  case used the grad-B part (1 - xi^2)/b alone (#238). Both are switchable.
- **Streaming**: Boozer angles have b.grad theta = k b, not a constant (#238).
- **Collision model switch**: in NEO-2-QL `isw_lorentz = 1` selects the
  monoenergetic mode, not pitch-angle scattering at every energy. A reference
  run chosen by the flag's name therefore solved another model. NEO-2 pull
  request 302 (itpplasma/NEO-2) replaced the three flags by one
  `collision_model` switch that rejects mixtures.
- **Normalisation**: NEO-2 normalises the flow coefficients of every species
  with the electron temperature (#308). Its D_3j, the parallel-flow responses
  to the gradients A_j, follow from MILCA's flow coefficients
  d(a, b, 3, j) as D_3j = -2 kappa q G d(a, b, 3, j), with kappa = iota B_0/
  (G + iota I) (#309).
- **Toroidal direction**: NEO-2 counts Omega_tE and the toroidal rotation
  V_phi along the direction in which the covariant toroidal field component is
  -G of the Boozer file (#310), on both devices.
- **Perturbation convention**: the perturbation files hold
  delta B = Re sum_m (bmnc - i bmns) e^{i(m theta + n phi)}, the kernel
  Re sum_m c_m e^{i(m theta - n phi)}, so c_{-m} = (bmnc + i bmns)/B_0 and the
  rotation enters the kernel as -Omega_tE (#271).
- **Left-handed files**: q < 0 occurs and must keep its sign (#271).

## 4. Resolution

- **Energy resonances** (#273): in strong rotation the resonance of E x B and
  magnetic precession enters and leaves the trapped region inside the ion
  energy range and puts kinks into the energy integrand. The AUG deuterium
  torque changed by up to 57 % between 48 and 96 Gauss-Laguerre nodes;
  adaptive panels converge it to 2e-3 relative.
- **Kink detection** (no pull request): the first speed grid aligned to these
  kinks for the energy-coupled collision operator found no kinks, because the
  detection threshold was too strict, and silently reproduced the uniform
  grid. A refinement that changes nothing has to be checked for having acted.
- **hp rows across the precession reversal** (open): the trapped rows of the
  hp mesh do not follow the levels of the conserved pitch and cut across the
  resonance where the precession reverses.
- **Low collisionality** (#211, #266): the cost sits in theta, not in
  Legendre degree, and the finite-nu corrections go like sqrt(nu) or nu, so
  limits are extrapolated from several nu.
- **Reference resolution** (#212, #214): NEO-2 reached MILCA's monoenergetic
  coefficients only at four times its default resolution, and its pitch grid
  was not converged at inverse aspect ratio 0.01.

## 5. Wrong limit or reference in a check

- **Paired sources** (#323): the finite-width kernel at small orbit width was
  compared with the thin kernel with only the drift source on and gave 0.78 of
  its torque on AUG at s = 0.5. With both sources off the two agree to 2e-7
  relative; switching both on changes the thin torque by only 0.2 %.
- **Scaling to the thin-orbit limit** (#323): scaling the mass at fixed
  temperature keeps the poloidal Mach number v_E B/(v B_pol) fixed, so the
  E x B energy change across the orbit never vanishes and the two kernels
  differed by 5 to 10 % at any mass. The thin limit scales temperature and
  E_r together. The relative deviation then falls linearly with the orbit
  width over the minor radius: 6e-4 at 1e-2 and 3e-5 at 1e-3, for both signs
  of the rotation.
- **Sign of the rotation**: on the up-down symmetric test equilibrium the
  thin torque barely depends on the sign of Omega_tE (4.049e-11 and
  4.061e-11), so agreement for both signs there does not fix the sign of
  E x B in the orbits. On AUG 30835 at s = 0.5 the thin torque changes by
  0.76 % between the two signs; the finite-width kernel matches it to 2e-4
  with the device's sign convention and misses it by 0.76 % with the opposite
  one.

## 6. Inconsistent input data

- **AUG Boozer file** (no pull request): |B| computed from the file's
  geometry and the tabulated |B| differ by 1.4 to 9.5 %, and G by up to 9 %, so
  dB/ds cannot be validated on this file.
- **ITER TC24** (#240): <B^2> from the geometric Jacobian is 8e-4 below the
  flux-Jacobian average of the same file.
- **Cancelling outputs** (#309): at s = 0.975 the ion flows cancel to 1 % of
  their terms, so relative differences there do not measure agreement.

## 7. Implementation and toolchain defects

- Allocation from an expression does not keep its lower bounds in gfortran 16
  (#260).
- The Jacobian-vector product of the collision assembly in the zero
  direction returned the primal matrix (#255).
- A Boozer regularity factor used the signed poloidal mode number and broke
  the symmetry of equivalent harmonics (#256).
- An internal function passed its own name as an actual argument; GCC 12 then
  built a trampoline and an executable stack (#323).
- Test checks that ran in parallel wrote one shared file and raced (#323).
- Root finders by bisection, secant and finite differences were replaced by
  Newton steps on generated residuals with exact tangents (#323), as the
  differentiability principle of [DESIGN.md](../../DESIGN.md) requires.

## Rules from this campaign

1. Every term gets an oracle written from its defining equation and a limit
   in which its effect is known; a test of the implementation against itself
   does not count.
2. Terms that belong together are switched and checked together.
3. Thin-orbit limits scale energy and E_r together, never the mass alone.
4. A reference's switches, weights, normalisation and directions are read
   from its source, not inferred from flag names.
5. MILCA and the reference both carry a resolution ladder; energy integrands
   with resonances use adaptive panels; a refinement is checked for having
   acted.
6. Device input is checked for internal consistency before a difference is
   attributed to the physics.
7. Sign checks use a case whose result changes with the sign.
