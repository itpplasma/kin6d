# D8 equilibrium verification harness

`kin6d_equilibrium_harness` owns the M6 case registry, seed descriptions,
resolution validation and report contract. `run_equilibrium_case` invokes a
lane driver and returns `equilibrium_state_t` plus `equilibrium_run_report_t`.
The driver is an explicit procedure argument: the equilibrium layer does not
depend upward on the MHD map implementation. Its native context retains the
mesh, map/frozen labels or free-boundary system needed to reproduce the result.
The drivers export through `map_equilibrium` and `gs_equilibrium`; the psi
export uses the solved boundary, plasma-cell mask and resolved profile law.
Missing field/chart exports or inconsistent lane provenance fail the run.

| Lane | Registry cases | Independent oracle |
|---|---|---|
| map | `solovev`, `shaped` | exact flux and a separate 2D GS solve |
| map | `seed`, `seedshaped` | exact flux from a displaced nested-axis seed, and 2D GS |
| psi | `limiter-jump`, `limiter-smooth` | manufactured compact flux bump plus exact coil flux |
| psi | `targets-jump`, `targets-smooth` | the same oracle with target-mean LCFS |
| psi | `xpoint-jump`, `separatrix-smooth` | manufactured two-bump FB-cut/separatrix |

The existing `ideal_map_n0_*` and `gs_free_boundary_t3_manufactured_newton`
CTest ladders run through this interface without changing their Newton
settings or analytic refinement gates. Grid sequencing still uses the seed
gradient of the current mesh as its stopping scale. Map Newton steps retain the
canonical symmetric Hessian and assembly-roundoff check, with the masked matrix
cached per pattern so that a trial only rewrites the diagonal. The damped Hessian
uses FortSparse's native supernodal pivoted LU action, as the free-boundary KKT
systems do. A nonpositive Cholesky pivot is not a rejection of the Newton trial:
raising the damping before the energy line search changes the nonlinear path
and can change the attained oracle error at the requested gradient tolerance.
No inertia is consumed by this map driver.
`verify_factor` of the map checks every factorisation and the cached input.
The factor follows the caller's OpenMP budget by default (`factor_threads = 0`);
a positive `factor_threads` still supplies an explicit cap.
An offset seed must be
corrected relative to the independent exact flux; stationarity alone cannot
pass the retained return gates.

`equilibrium_d8_psi_offset_{reference,rp,rm,zp,zm,el}` runs the P2 limiter
ladder at 24, 48 and 96 boundary nodes. The starts displace the manufactured
plasma bump by +/-15 mm in R or Z, or elongate it by 2%, while holding the
forcing, limiter, coils and profile law fixed. The finest result must have
less than 3% relative maximum nodal error and remove at least 90% of the
nonzero seed error; the original convergence-order gate also applies.
The nested seed uses FortNum's FFT for the supplied boundary and retains
coefficients down to roundoff. `ideal_seed_boundary_label` checks its label
on independent exact LCFS points between the Fourier sampling nodes, against
the unchanged `1e-12` physical-state export precondition.
The legacy T3 start retains its larger displacement and amplitude change,
and its quadratic-Newton-phase gate. A near-solution offset need not traverse
that same quadratic phase.

Each report records the case/lane/start/resolution, convergence status,
iterations, independent oracle error, seed error and an a posteriori residual.
The map residual is the retained sigma-normalised D3 L2 force measure;
sigma excludes the curl-free toroidal background. The manufactured psi
residual is the discrete forced-equation residual divided by the manufactured
load norm. These are different norms and must not be compared as equal
quantities. In particular, the manufactured body load is not unforced ideal
MHD, so its residual is not labelled as D3 force balance.

All these results remain **candidate** states. State and derivative error
bounds are explicitly **unavailable**: sampled oracle errors and algebraic
residuals are not deterministic certificates. The harness rejects missing,
negative or nonfinite evidence even if the driver reports success. The
registry contract test checks rejected combinations, returned native state,
NaN/missing evidence and the absence of automatic certificate promotion.

Focused execution (native CMake/CTest through Fo):

```sh
FO_CMAKE_BUILD_DIR=build-dev-fast fo test --all \
    equilibrium_d8_registry_contract 'equilibrium_d8_psi_offset_*' \
    gs_free_boundary_t3_manufactured_newton 'ideal_map_n0_*'
```

The remaining D8 cross-lane transfers, common D7 oracle comparisons, free-B
reductions and D5 identities are separate gates; this M6 migration does not
claim they have all been implemented by the harness.
