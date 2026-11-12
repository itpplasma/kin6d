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
kink of the independent oracle to 1.3e-10.

## Causes

### Implementation and toolchain defects

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

### Wrong limit or reference in a check

- **Symbolic proofs blind to the current term.** The exact identities of
  `gen_ideal_dw` (uniform medium, plane waves) do not contain the current
  and pressure-gradient terms; flipping the sign of the current term left
  them all proved. Integrating both energy forms on a force-balanced screw
  pinch exposes it with an 8 % mismatch, and the test fails without force
  balance (1.5 %).

## Rules

- Discretize the displacement through its contravariant densities, so that
  the divergence is exact in the tangential space.
- Every observable with an exact gap in its spectrum gets a gap count.
- A symbolic identity covers only the terms it contains; check the remaining
  terms with an integral identity on an equilibrium.
