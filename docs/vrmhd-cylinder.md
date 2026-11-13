# Helically symmetric straight-cylinder flux-pumping case

This case is the 2D part of RM-B8 and RM-B10 in
[the MHD design](nonlinear-resistive-mhd.md) (milestone RM2): a periodic
straight cylinder driven by a fixed loop voltage, restricted to the helicity
of the m/n = 1/1 mode. It is the same physical setup as a set of JOREK
cylinder runs and of an independent single-helicity reduced-MHD model, both
kept outside this repository; this page states the setup so that their
numbers can be compared.

On [the model ladder](nonlinear-resistive-mhd.md#model-ladder-on-the-resistive-time-scale)
this case is rung L1, the cheap mechanism model: it isolates the saturation
of the 1/1 mode, its electromotive force, `eta(T)` and anisotropic conduction
at the cost of a 2D problem, and it contains no toroidal effect. Its planned
diagnostic is the loop-voltage ledger RM-B13 of the ladder.

## Geometry, units and symmetry

Minor radius `a = 1`, axial period `2 pi R0` with `R0 = 10`, axial field
`Bz(a) = 1`, `mu0 = mi = 1` and mass density one, so that time is in Alfven
times `a sqrt(mu0 rho)/B0`. Fields depend on `r` and `theta - z/R0`, the screw
symmetry `xi = (-k y, k x, 1)` with `k = 1/R0` of `gs_symmetry_t`. The section
uses the closed coordinate coframe of [the section contract](vrmhd-section.md),
which the Fourier operators of RM3 also need. The wall at `r = a` is an ideal
no-slip conductor with fixed temperature.

## Profiles

`vrmhd_cylinder_profiles` reimplements the recipe of the JOREK case generator
(uniform perpendicular conductivity `K`, heat source
`S0 (1 - tanh((psi_N - 0.3)/0.1))/2`, Ohmic current `j ~ T^(3/2)`, `psi_N`
of the same current, iterated to a fixed point on a uniform radial grid with
trapezoidal integrals). With `q0 = 0.9` on axis the edge value is
`q_a = 3.352` and `q = 1` lies at `r = 0.264`, `rho_p = 0.355`. The
temperature is scaled to the poloidal beta `beta_p` and enters as
`Te = Ti = T/2`, with `T = p` for unit density.

The full-MHD state needs radial force balance, which the reduced models do not.
With `p = T`, the axial field follows from
`d(p + Bz^2/2)/dr = -B_theta J_z` and `Bz(a) = 1`; the resulting poloidal
current `J_theta = -dBz/dr` is kept by the prescribed source
`Js = J_theta e_theta`. An optional central current source
`S_j(psi_N) e_z` adds to `Js`, at fixed total current, as in the JOREK
current-source scan. The loop voltage drives `E_z = eta (J_z - S_j)`, which
is uniform for `eta = eta0 (Te/T0)^(-3/2)` with `eta0 = 1/S` at the axis.

Here `q = r Bz/(R0 B_theta)`; the reduced models use `Bz = 1`. The two differ
by `|Bz - 1|`, about `1e-3` at `beta_p = 1`.

## Closures and sources

| quantity | value |
|---|---|
| resistivity | `eta0 (Te/T0)^(-3/2)`, `eta0 = 1/S = 1e-6`, `T0 = Te(0)` |
| viscosity | constant, `Pm/S` with `Pm = 10` |
| perpendicular conduction | `chi = 2.4e-5` in `d_t T = chi lap T`, hence `1.5 chi` in the internal-energy form |
| parallel conduction | constant, `chi_par/chi_perp = 1e6` |
| Ohmic heating | `ohmic_fraction = 0` by default, as in the JOREK runs |
| heat sources | `Qe = Qi = -div(n chi grad Ts)` of the initial profiles |

The heat sources make the initial temperature stationary under perpendicular
conduction; a helically symmetric-free state has no parallel gradient. With
`ohmic_fraction = 0` the Ohmic power leaves the system as an explicit ledger
entry ([weak form](vrmhd-weak-form.md#mixed-current-and-energy-exchange)).

## Mesh

`vrmhd_cylinder_mesh` triangulates the disk in rings at given radii. The
number of points per ring starts at a given count and doubles where the arc
length exceeds a given multiple of the radial step, so radial refinement at the
`q = 1` layer does not force the same azimuthal resolution. The wall is a
polygon through the ring points.

## Stationary solve

`vrmhd_steady_solve` (module `milca_vrmhd_steady`) finds driven stationary
states of the section, `R(U) - V h.W = 0` with the Coulomb multiplier of
[the stationary wrapper](vrmhd-section.md#boundary-and-stationary-gauge), by
pseudo-transient continuation: each step solves
`(R_U + M/dtau) dU = -R` with a supernodal direct factor whose ordering and
matching are analysed once and reused, since the Jacobian pattern is fixed.
The pseudo step grows by a constant factor unless the residual rises by more
than a factor of ten, and is halved otherwise. Beyond `dtau = 1e8` (a hundred
resistive times) a step counts as Newton, and the solve is converged when a
Newton step changes no field by more than `1e-8` of its scale (the Alfven
speed for the velocity): the Ampere rows reach their roundoff floor long
before the others, so a residual-norm test alone does not terminate.

At low beta the density is weakly determined: the pressure must absorb the
discretisation error of `J x B`, whose scale is `beta` times smaller than the
magnetic pressure. Pure Newton steps then change the density by tens of
percent on coarse meshes and leave the region of the linear model; the
pseudo-transient path does not. A small particle diffusivity
(`density_diffusion = 1e-6`, JOREK's `D_perp`) is part of the case.

`vrmhd_steady_spectrum` returns eigenvalues of `M x_t = -R_U x` about a
stationary state nearest a real shift `sigma` by Krylov-Schur iteration on
`(R_U + sigma M)^(-1) M`. Fixed coefficients, optionally frozen fields and the
gauge multiplier carry no mass. Each `m >= 1` eigenvalue of the section is
double, because the stationary state is invariant under rotation of the
section; the pair is returned as a nearly real conjugate pair.

## Measurements

Mixed family with thermodynamic degree two, rings clustered at the `q = 1`
radius (strength 3, uniform fraction 0.3), 24 rings (1264 triangles, 90817
unknowns) or 48 rings (4976 triangles, 356657 unknowns).

| quantity | 24 rings | 48 rings | reference |
|---|---|---|---|
| m = 1 growth rate, frozen temperatures, `beta_p = 0.01` [1/tau_A] | 7.1706e-4 | 7.1723e-4 | 7.2187e-4 (visco-resistive reduced MHD, Chebyshev eigenvalue solver) |

The remaining 0.6 % is converged in the mesh; the full-MHD state has
`Bz = 1 + O(B_theta^2)` from radial force balance, while the reduced model
has `Bz = 1`, and the growth rate changes by about 1 % per 0.1 % change of
`q` (computed with the same reduced model).

With `chi_par/chi_perp = 1e6` the stationary state on 24 rings differs from
the interpolated profiles by 0.9 % in the electron temperature and 1.3 % in
`q(0.1)`; with `chi_par/chi_perp = 1` the difference is below `1e-4`. The
parallel conduction along the discrete field, whose radial component is a
discretisation error, adds a spurious perpendicular conductivity of order
`chi_par (B_r/B)^2`; at `chi_par/chi_perp = 1e6` a relative error of `1e-4`
in `B_r` already matches `chi_perp`. This is the RM-B3 question of the design
on an unaligned mesh.

## Gates

| test | check |
|---|---|
| `vrmhd_cylinder_profiles` | `dT`, `q_a`, `q(0.3)` against the reference generator on the same 2001-point grid to `1e-12`; force balance of the interpolated profiles |
| `vrmhd_cylinder_stationary` | stationary solve from the interpolated profiles with isotropic conduction; temperature change below `1e-3` of T(0), velocity below `1e-5` (Alfven units), `beta_p = 0.01` (slow) |
| `vrmhd_cylinder_linear` | m = 1 growth rate with frozen temperatures within 5 % of the reduced-MHD eigenvalue `7.218665e-4` of the same profiles (slow) |

## Probes

`test_vrmhd_cylinder` takes `key=value` options (`rings`, `degree`, `cluster`,
`dtau`, `beta_p`, `chi_ratio`, `frozen`, `seed`, `t_end`, `dt_out`, `dt_max`,
`state`, `check`). The
`evolve` mode solves the stationary state, seeds BDF2 with its unstable
`m = 1` eigenvector scaled to a peak velocity `seed`, and prints `q` at
`r = 0.1` (annulus 0.09-0.11, from the mean of `B_theta/r`) and the magnetic
and kinetic energy of the `m = 1` part (Fourier projection on 400 annuli).
`dt_max` caps the adaptive step separately from the output interval `dt_out`
(the step also ends on every output time). `check=2` prints every
linearization and factorization. With `state=FILE` the final section state is
written to `FILE`. The `pumped` mode evolves as `evolve` (or reads `FILE` if it
exists), then solves the stationary helical state from it by the same
pseudo-transient Newton iteration and prints the six eigenvalues of its
linearisation nearest zero: a stable pumped state has all of them in
`Re < 0`.
These trajectories are candidates; the comparison with the reduced model and
with JOREK is not a registered gate.

### Cost of the evolution

The BDF2 Newton systems are solved by FGMRES on the exact Jacobian-vector
product, preconditioned by a lagged supernodal factor of the assembled
Jacobian. The factor is replaced when the BDF shift changes by more than
`shift_tolerance`, after `lag` linearizations, or when the previous solve
needed more than `krylov_refactor` (8) products. Perturbed or grown pivots are
accepted (`verify_perturbed`, `pivot_growth_limit`), since FGMRES iterates on
the original matrix; otherwise fortsparse replaces such factors by its serial
scalar LU, which takes minutes at 1e5 unknowns. The residual and its
Jacobian-vector product evaluate the pointwise form once per test slot and
contract it with the basis.

On 24 rings (90817 unknowns), 16 OpenMP threads on one sCluster node, the
eta(T) case with `chi_ratio=1e6` reaches t = 3500 tau_A in 441 s, against
35741 s before these changes, where the lagged factor went stale and FGMRES
hit its 200-product limit. A solve then needs about 1.5 FGMRES products per
Newton iteration, and the run covers about 2e4 tau_A per wall hour through
the first saturation of the m = 1 mode.
