# Frozen fitted profiles in the non-nested lane

[V] `kin6d_nonnested_profile` implements the I5a seed, fit, freeze, solve
sequence of [the island design](island-equilibrium.md#b-prescribed-profiles-frozen-fitted-sources-then-projection).
It supplies data to the existing NN2 A+D equations. It introduces neither a
profile penalty nor a projected equilibrium solve. Numerical results remain
candidates, and continuous error bounds and geometry derivatives are unavailable.

The caller implements `nn_profile_t%evaluate` with owned reference-label and
profile data, returning pressure, its physical gradient, the parallel-current
target lambda = J.B/|B|^2 and the seed potential gradient. Points have order
(R,Z,phi), while vectors have physical cylindrical order (R,phi,Z). The caller
records the reference label, units, normalisations, seed, boundary traces,
transport coefficients and current constraint in the case provenance. The label
is reference data and need not label surfaces of the final field.

[V] `nn_fit_source` copies the seed and profile into an `nn_fitted_source_t`
whose data are private. The implementations supplied by callers must own their
data, including anything reached through pointer components, and support
concurrent read-only sampling. Later edits to caller-owned value data cannot
change the frozen fit. Null fields, failed sampling, nonfinite targets and
invalid coefficients return failure rather than a field-magnitude floor.

On that fixed seed it defines

    q_seed = K(B_seed) grad p_target,
    S = -div q_seed,
    d_seed = lambda_target
             - g_par (B_seed,phi V_seed/(2 pi R)
                      - B_seed.grad Phi_seed)/|B_seed|^2.

The pressure source is represented weakly as `<grad v,q_seed>`. This avoids
requiring second derivatives of a supplied profile and retains the toroidal
metric in the weak load. Dirichlet pressure data are still the caller's fixed
boundary entries. A Neumann heat-flux interpretation is not implemented here.
The signed drive adds `d_seed B` to NN2's current. The live B, pressure,
potential and loop voltage still solve the original coupled equations with
that additional explicitly declared parallel drive. In particular, the fit
does not replace the pressure-driven current or waive current continuity.

A typical current-driven call is:

```fortran
allocate (problem%fitted)
call nn_fit_source(seed, target, problem%chi_par, problem%chi_perp, &
    problem%g_par, problem%g_perp, voltage_seed, problem%fitted, status)
! Check status, supply the fixed trace entries of x and problem%current.
call nc_solve(field, problem, x, report, status)
```

[V] Scalar source hooks and conductivity-shape hooks are rejected when a fit
is present. This increment keeps uniform positive base conductivities and
fits the signed current drive instead of trying to infer a potentially
negative conductivity. The frozen functions are evaluated on their copied
seed during assembly, with no fitting on the live state. Changing live
transport coefficients is a different declared model at the same frozen source.
Equation-level state derivatives hold the fit fixed. Generated primal/JVP/VJP
leaves define the drive fit, current correction and point diagnostics. The
source assembly and existing coupled element loops use OpenMP.

`nn_state_profile_t%init` adapts an existing uniform-conductivity A+D solution,
including its potential and voltage. It applies FortFEM's Fourier gradient
operator once and samples the resulting V1 fields. It therefore retains all
Fourier derivatives and orientation signs when fitting an existing solution.
It does not infer current from an unqualified curl of an H(div) field.

`fit%defects` reports physical weighted RMS pressure/target and lambda/target
mismatches, the seed current-law fit error and the parallel-pressure gradient.
The caller supplies the solved gradients, positive weights and voltage. Its
conductivities must match those of the frozen fit. Evaluate `nc_residual` on
the seed and inspect the Ampere and continuity rows separately: an exact local
fit can still be incompatible with force balance or charge conservation.
`nc_report_t` retains the solved block residuals and candidate status. A small
fit defect alone is never an equilibrium test.

[V] `kin6d_nonnested_projection` implements I5b only for **declared** invariant
components and positive quadrature weights. Its orthogonal projection is the
weighted component mean of the supplied pressure, corresponding to
Pi_B[p_target(s_ref)], rather than p_target(Pi_B[s_ref]). Closed-line weights
must represent dl/|B|, invariant-surface weights their invariant measure and
explicitly mixed volumes their volume measure. Internal island surfaces stay
separate unless the caller explicitly declares whole-island mixing. Open-line
components are rejected. The operator does not detect topology, prove
ergodicity or certify that a finite trace resolved a component. Its fixed-
component tangent is the same projection and its adjoint uses the declared
weighted pairing. Derivatives across component changes are unavailable.

[V] The registered tests retain the requested numerical thresholds:

- `nonnested_fitted_manufactured_p2/p3`: a smooth non-flux-function profile,
  exact Fourier content and an independent strong cylindrical Laplacian check
  of the fitted weak load, with L2 orders at least k+0.6. This is transport
  verification with force-incompatible manufactured data.
- `nonnested_fitted_nested_p2/p3`: prescribed Solov'ev pressure/current on
  straight-sided domains, offset return, volume L2 pressure order at least
  k+0.6, poloidal-field order at least k-0.4 and pressure error below 1e-3.
  Exact nodal traces belong to the declared polygon, isolating profile order
  from curved-boundary approximation. A blocking seed-sampling defect in `flux_mesh_locate` is repaired by
  searching adjacent strips before projecting a point to a cell. The
  `nonnested_fitted_sampling_p2/p3` tests compare physical sampling against
  direct element quadrature, so this repair is checked by its real consumer.
  The physical profile is not copied
  into the solved interior.
- `nonnested_fitted_reproduce_p2/p3` and `reproduce3d_p2/p3`: fit a solved A+D
  state, freeze and return from an offset, with relative state error below
  1e-8. Frozen-source Jacobians are compared with directional differences.
- `nonnested_profile_projection_g1b`: integrable surfaces, closed rational
  lines, helical-island internal surfaces and explicit volume mixing, with
  idempotence, component means and weighted orthogonality at 1e-9. Trace
  quadrature ladders compare rational means with a closed form and island
  means with independent complete elliptic integrals at 1e-6.
- `nonnested_profile_forced_island_g4`: the known Hamiltonian
  H=x^2/2+epsilon cosh(x) cos(alpha), an independent outer boundary-value
  solution, time-integrated separatrix classification and phase checks.
  Four amplitudes test the square-root width law in canonical coordinates
  and toroidal geometries at aspect ratios 6 and 60, with explicit physical
  radial conversion. The same field admits an invariant pressure that the
  frozen-source transport solve recovers within 1% on a mesh ladder.

[H] The G4 field is a prescribed divergence-free Hamiltonian oracle with a
regular outer perturbation. It is not a nonlinear saturated tearing solution,
a free-boundary solution or evidence that arbitrary pressure/current targets
admit a static MHD equilibrium. Full coupled island return and D3/ideal-limit
qualification remain separate programme gates. The component diagnostic has
no general field-line topology classifier.

[V] I1 integration boundary: `islands/i1closure` was read while its branch
is at 5fda1b3654c718fd290135fc7085545aa0475fa1. Its declaration names
`B_frozen_fit` and returns `ISLAND_DELEGATE_B`. This lane changes only NN2's
source/current hooks, with no edits to I1's declaration, reader or report
modules. During integration, I1's validator must accept that delegation only
when a ready `problem%fitted` is present, retain rejection of `B_projection`
as a solve, and use the weak-source/current-drive path in its balance report.
Its voltage-driven path can reuse the same frozen current drive. This lane
uses NN2's current-driven interface pending that integration.
