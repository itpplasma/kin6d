# Frozen fitted profiles in the non-nested lane

`kin6d_nonnested_profile` implements the I5a seed, fit, freeze, solve
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

`nn_fit_source` copies the seed and profile into an `nn_fitted_source_t`
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
problem%closure%pressure = 'B_frozen_fit'
allocate (problem%fitted)
call nn_fit_source(seed, target, problem%chi_par, problem%chi_perp, &
    problem%g_par, problem%g_perp, voltage_seed, problem%fitted, status)
! Check status, supply the fixed trace entries of x and problem%current.
call nc_solve(field, problem, x, report, status)
```

Scalar source hooks and conductivity-shape hooks are rejected when a fit
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

## I1 declaration and report integration

Declare `pressure_closure='B_frozen_fit'` and supply a ready `problem%fitted`.
The case reader remains transactional. Missing or uninitialised fits stay
explicitly delegated, conflicting scalar source or conductivity-shape hooks
are rejected, and `B_projection` remains a diagnostic delegation rather than a
solve. Readiness does not prove valid sampling; point-data failures propagate.

The frozen parallel current drive enters the assembled residual and I3
JVP/VJP actions with the fit held fixed. Net toroidal current keeps its phi-mean
meridian measure, while boundary charge uses the full-torus factor 2*pi.

For the weak fitted pressure source the raw pressure reaction is the live-minus-
seed flux reaction. The report exposes `weak_heat_reaction` and its absolute
compatibility defect `weak_heat_balance`; they are not physical heat. Physical
`heating`, `outward_heat` and `heat_balance` are NaN with an explicit unavailable
`heat_status`, because no qualified seed normal-flux trace is supplied.
Scalar-source A+D retains its existing physical balance semantics. Neither
readiness nor these balance diagnostics establishes nonlinear saturated
islands, a general topology classifier, or a continuous equilibrium bound.

The public fitted-profile regressions cover manufactured pressure convergence,
nested-limit pressure and field convergence, existing A+D state reproduction,
and curved-cell sampling at P2/P3. Additional I1 regressions exercise readiness,
transactional rejection, the relative weak heat reaction, signed toroidal
current, residual finite differences and assembled/JVP/VJP consistency at
prescribed nonstationary states. Such states are not promoted equilibria.

For a pure harmonic field at zero pressure, the retained D3 normalization has
sigma = mean(|B-B0|^2)/mu0 + mean(|p|) = 0. Its absolute finite-test-space load
and the other report quantities remain available, while `force%normalized` is
NaN and `force_status` explicitly says that normalization is unavailable. The
standalone D3 evaluator returns `NN_RESIDUAL_ZERO_SCALE`; the island report
recognizes only that named availability condition and propagates other errors.
This does not replace the declared perturbation scale with total magnetic
energy, introduce a denominator floor, or promote a prescribed current state
to an equilibrium.
