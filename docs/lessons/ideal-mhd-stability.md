# Lessons: global ideal-MHD stability

Physics area: global ideal-MHD stability of fixed- and free-boundary
equilibria (roadmap item P1e), specified in
[stability.md](../stability.md#global-3d-ideal-mhd-stability-p1e). delta W is
the ideal-MHD potential energy, xi the plasma displacement, and the hybrid
space the radial pair of a continuous normal component and discontinuous
tangential components of one degree less.

References so far: exact uniform-cylinder spectra (fast and slow roots from
the zeros of the Bessel derivative J_m'), the independent test-side cylinder
oracle of [s0-spectrum.md](../s0-spectrum.md) and its analytic m = 1 limit
omega^2 -> -3.

Outcome: in the periodic cylinder the hybrid space reproduces the exact
uniform spectrum to 1e-7 without spectral pollution, and the fixed-wall m = 1
kink of the independent oracle to 1.3e-10. On a Solov'ev torus the PEST and
Hamada charts give the same internal-kink eigenvalue to 2.6e-5.

## Causes

### Implementation and toolchain defects

- **Axis scaling of Fourier coefficients.** Dividing harmonic m of R and Z
  by rho^m before the radial fit amplified rounding noise by up to rho^-32
  and destroyed the chart beyond m = 10; in the geometric angle about an
  elliptic axis the odd harmonics vanish only like rho, so the chart left
  its flux surfaces by 10 %. The divisor rho^k with k = min(m, 2 - mod(m, 2))
  keeps the parity and is valid in every regular angle (`ideal_flux_tokamak_field`).
- **Banded numbering.** Numbering the unknowns harmonic by harmonic gave a
  bandwidth of the matrix size, and the factorization of a 1712-row problem
  did not finish in minutes; radial-major numbering makes it banded.

- **Spectral pollution from physical components.** With the physical
  components (xi_r, xi_theta, xi_z) as unknowns, the uniform cylinder
  showed 22 eigenvalues between the Alfven point and the first fast mode,
  where the exact spectrum is empty. The divergence (r xi_r)'/r then is not
  in the space of the tangential terms. Taking the contravariant densities
  sqrt(g) xi^i = (r xi_r, xi_theta, r xi_z) as unknowns makes
  r div xi = X' + i(m T + kz Z) exact in the derivative space and removed
  every spurious eigenvalue. Found by the gap count (observable O5) of
  `ideal_cylinder_uniform`, which now guards it.
- **Lagrangian form outside H1.** The Hessian of the frozen-in energy uses
  the full gradient of xi, including radial derivatives of the
  discontinuous tangential densities, and gave omega^2 = -2.3e4 instead of
  -0.195 in the hybrid space. The Eulerian form needs only div xi and
  tangential derivatives and is the one assembled there; the Lagrangian form
  remains an independent check for H1 fields (`ideal_dw_pinch_equivalence_and_derivatives`).

- **Node-wise slip on a curved wall.** In the position gauge of the map
  lane a fixed wall keeps the displacement component along the wall. Letting
  each wall node slide along its own tangent gave 26 to 40 eigenvalues below
  -10 on the Solov'ev internal-kink case (oracle: one mode at -1.3e-3), more
  on finer meshes. Two causes. A point that slides along a curved wall stays
  on it only with a second-order normal shift, on which the wall force does
  work; without this curvature term the stiffness of the slip is wrong at
  order one. With the term, 2 to 8 spurious modes near -1e-2 remained
  whenever the radial spacing at the wall was finer than the poloidal one:
  between the nodes the interpolated slip is not tangent, the discrete wall
  bulges, and nothing restores a plasma that expands. The slip is now
  tau(s) X_s on each wall edge X(s) with linear tau, which is tangent to the
  discrete wall everywhere; all spurious modes are gone and the eigenvalue
  approaches the oracle from above. Clamping the wall also removes them but
  leaves the eigenvalue 5 to 17 % above the oracle on meshes where the
  tangent slip is 1 to 7 % above it.
- **Toroidal displacement as third unknown.** Adding the toroidal
  displacement with the frozen-in Hessian on all three components for the
  adiabatic term gave an eigenvalue 7 % below the gamma = 0 value in the
  limit gamma -> 0. A parallel displacement beta B, on which the frozen-in
  Hessian does not act, returns the gamma = 0 value to 5e-10. [H] The cause
  of the first result was not established.

### Resolution of MILCA

- **Unstructured mesh for a near-marginal internal kink.** On quality
  triangulations of the Solov'ev case the fixed-boundary kink converges at
  second order only (7.7 % off the oracle at 41 000 cubic nodes; isotropic
  packing at q = 1 gives 1 % at 79 000). A flux-aligned ring mesh, whose
  radial spacing is independent of the poloidal one, gives 0.46 % at 11 000
  nodes and 5e-5 at 120 000. The radial spacing away from q = 1 mattered
  more than the packing at q = 1: halving the largest spacing reduced the
  error by 3.6, halving the smallest by 1.4.

- **One-sided harmonic window.** With iota < 0 the resonant harmonics of
  n = 1 are positive m. Windows m = -2..16 (PEST) and m = -8..2 (Hamada)
  showed a stable spectrum, m = -2..8 (Hamada) an internal kink, and only the
  wide window m = -10..28 brought both charts to the same unstable
  omega^2 = -1.3033e-3. Truncation is Ritz-monotone: a missing sideband can
  turn an unstable equilibrium stable. Found by comparing two charts; the
  invariance gate `ideal_flux_tokamak_invariance` uses a window with
  sidebands on both sides.
- **Flux-grid spacing.** The current of an EQDSK-type input comes from the
  second derivatives of the interpolated flux. On a 161 x 201 grid the force
  balance of the Solov'ev chart was off by 6 %; it converges as h^3.5 and is
  2e-4 on 641 x 801.

### Inconsistent input data

- **Under-resolved equilibrium current.** The circular-torus states of
  [vrmhd-torus.md](../vrmhd-torus.md), sampled from cubic Grad-Shafranov
  elements of size 0.05 onto a flux grid, gave unstable n = 1 eigenvalues in
  the low-beta control, whose number grew under refinement and differed
  between PEST and Hamada charts. The chart's current, from second
  derivatives of the sampled flux, violated force balance by up to 1.3 at
  mid-radius and 147 near the axis relative to grad p. On elements of size
  0.0125 the residual is 0.013 at mid-radius and the control is stable. A
  stability result now reports the force residual of its chart; an unstable
  count that grows under refinement and differs between charts is a sign of
  the input, not of physics.

### Wrong limit or reference in a check

- **Symbolic proofs blind to the current term.** The exact identities of
  `gen_ideal_dw` (uniform medium, plane waves) do not contain the current
  and pressure-gradient terms; flipping the sign of the current term left
  them all proved. Integrating both energy forms on a force-balanced screw
  pinch exposes it with an 8 % mismatch, and the test fails without force
  balance (1.5 %).

## Rules

- A slip condition on a curved wall is imposed on the whole discrete wall,
  not at its nodes, and carries the curvature term of the wall force.

- Discretize the displacement through its contravariant densities, so that
  the divergence is exact in the tangential space.
- Every observable with an exact gap in its spectrum gets a gap count.
- A symbolic identity covers only the terms it contains; check the remaining
  terms with an integral identity on an equilibrium.
