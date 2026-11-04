# Wall-slab island benchmark

The wall, equilibrium and initial-induction checks admit a setup for nonlinear
island calculations. They do not establish the Rutherford or plasmoid gate
RM-B5 in [the MHD design](nonlinear-resistive-mhd.md).

## Equilibrium and walls

All quantities below are nondimensional. The Cartesian domain is
`-1 <= x <= 1`, `0 <= y <= 2 pi`, with translation symmetry in `z` and
periodicity in `y`. The permeability and ion mass are one. The equilibrium is

```text
v = 0,             A = (0, 2 x, -x^2/2),
B = (0, x, 2),     J = (0, 0, 1),
Te = Ti = 1,       rho = 1-x^2/4,       p = 2-x^2/2.
```

Here `v` is velocity, `A` the vector potential, `B=curl A` the magnetic field,
`J=curl B` the current, `rho` mass density, and `Te,Ti` electron and ion
thermal energies. Pressure satisfies `p=rho(Te+Ti)`. Thus `grad p=J cross B`.
The prescribed current source is `Js=J`; both the resistive electric field
and the source-accounted Ohmic heat vanish at equilibrium. Momentum, density
and thermal sources are zero. FortSym checks the force balance and generates
the equilibrium and seed evaluations.

The walls at `x=+/-1` have zero velocity and fixed tangential perturbations of
both `A` and `J`. The latter is consistent with Ohm's law at a stationary
conducting wall with fixed tangential `A`. Density and temperatures retain their
natural equations. Independent checks interpolate the radial coordinate to
show that opposite walls are distinct, compare periodic tangential traces,
and verify zero tangential wall perturbations through physical point evaluation.

The thermodynamic polynomial degree is two. The equal family uses velocity
and second-kind Nedelec magnetic degree two; the mixed family uses velocity
and magnetic degree three. Both families represent this equilibrium exactly.
A full mass-matrix residual checks the potential projection. The stationary
residual, consistent initial derivative and BDF2 drift are checked separately.

## Seed, induction and width

For seed amplitude `epsilon`, add

```text
delta A_z = epsilon (1-x^2)^3 cos(y).
```

The perturbation and its first two radial derivatives vanish at the walls.
Initially `v=0` and both temperatures are one, so Ohm's law independently gives

```text
d/dt [A_z(0,0)-A_z(0,pi)] = -14 eta epsilon,
```

where `eta` is electrical resistivity; with unit permeability it is also the
magnetic diffusivity. The initial-induction check compares this strong-form
value with the consistent derivative of the assembled differential-algebraic
system, using a 2% relative error budget. The [native refinement data](
data/vrmhd-island-induction.csv) use `eta=epsilon=0.01` and eight or sixteen
cells in each direction, compiled with GNU Fortran 16.2 and OpenBLAS in Release
mode. The mixed family's relative errors are 1.01% and 0.067%; the equal
family's are 7.37% and 2.61%. The [plot source](plot_vrmhd_island_induction.py)
reproduces the comparison; these two levels do not establish an asymptotic
convergence theorem.

For a reflection-symmetric single island, the separatrix value is
`A_z(0,pi)`. The full cut width `W` is the distance between its two crossings
on `y=0`, found by bisection of the physical finite-element potential. A
separate quadratic-flux example has exactly `W=4 sqrt(epsilon)`; adding a
constant to the flux leaves the measured width unchanged. This diagnostic
assumes the stated reflection symmetry and locations of the elliptic and
hyperbolic critical points. An evolved field that loses these properties
requires critical-point tracking before the cut can be interpreted as an
island width.

For the unperturbed constant-shear slab, the ideal outer equation is
`psi''-psi=0`, with conducting-wall perturbation `psi(+/-1)=0`. Its normalized
logarithmic derivative jump is `Delta'=-2 coth(1)`. The reference Rutherford
slope is `1.22 eta Delta'`, using the full-width convention of
[the independent island references](linear-mhd-oracles.md#rutherford-growth-and-saturation).
Comparison requires a resolved resistive layer much narrower than the island,
an island much narrower than the equilibrium scale, and `abs(Delta' W)<<1`.
Initial resistive relaxation does not by itself establish that regime.

## Native checks and candidate trajectories

Build `test_vrmhd_island` and `gen_vrmhd_island`, then run
`ctest --test-dir build -R '^vrmhd_island_' --output-on-failure`.
The existing periodic BDF2 checks remain separate regression tests.

The executable also accepts
`probe FAMILY NY NX ETA DT STEPS [EPSILON [scalar|supernodal]]`, where
`FAMILY` is `equal` or `mixed`, `NY,NX` are cell counts, `DT` the constant time step, and `STEPS` the
number of steps. The default amplitude is 0.01. The viscosity coefficient is
set equal to `ETA`; because density varies, the kinematic viscosity and
magnetic Prandtl number vary across the slab. Output gives the cut width,
critical-point flux difference and maximum velocity coefficient at every step,
followed by the physical mass and energy balances described below.
The default `supernodal` and the `scalar` control select the
[native direct factors](vrmhd-section.md#time-evolution). The lagged factor
is only a right preconditioner: Newton and Krylov iterations still use the
full nonlinear residual and tangent. A nonlinear retry or an internal
subdivision fails this constant-step probe rather than admitting a trajectory
on a different time grid. The `vrmhd_island_decimal_steps` regression checks
eight steps of length `0.1`, including decimal endpoints that previously
caused the BDF driver to add a microscopic internal step. It checks time-grid
continuation, without admitting the coarse spatial mesh for island rates.
Output also gives iteration counts and the last
factor's family, fallback status, retained owned array storage, and available
pivot and last-solve backward-error statistics. Retained array storage excludes
process memory and the peak workspace of sparse graph ordering.

These trajectories are candidates until spatial and temporal refinement,
critical-point checks, physical ledgers and the required regime checks admit
a rate comparison. No Rutherford agreement is claimed by the setup checks.

The [direct-factor comparison](data/vrmhd-island-preconditioner.csv) uses
the mixed family with four periodic and eight radial cells, `eta=0.001`,
`epsilon=0.005625`, `DT=0.25` and twenty steps. Both factors accept exactly
the requested time grid. The full cut widths differ by at most `1.3e-14`
relative and the flux differences by `5.5e-15`; a separate physical
quadrature norm checks all state fields against the scalar trajectory.
With GNU Fortran 16.2, Release optimization and OpenBLAS, two executions
take 84.7–86.9 s with the scalar factor and 58.8–58.9 s with the supernodal
factor. Both use `OMP_NUM_THREADS=4`; this OpenMP OpenBLAS build also reports
four threads despite `OPENBLAS_NUM_THREADS=1`. These are whole-run timings
under that configuration, not an isolated factor-kernel scaling result.
The [timing data](data/vrmhd-island-preconditioner-time.csv) and
[plot source](plot_vrmhd_island_preconditioner.py) retain the comparison.
The unchanged energy, decay, wave and steady-state gates check independent
analytical references; these timing trajectories do not admit island rates.

## Physical balances of candidate trajectories

The mass diagnostic integrates density. The energy diagnostic integrates
kinetic energy `rho |v|^2/2`, magnetic energy `|B|^2/(2 mu0)` and the two
species' internal energy `3 rho (Te+Ti)/(2 mi)`. Their instantaneous rates
are calculated by differentiating these physical expressions, independently
of the residual leaves. Analytic polynomial fields check every integral and
rate. A separate seeded test with electron thermal energy `1+0.1 cos(y)`
checks the source-accounted balances and exercises a nonzero wall power for
both field families. This test is distinct from the uniform-temperature
initial-induction oracle.

A fixed nonzero tangential current has an affine constraint power, `P_wall`,
as defined in [the current-wall contract](
vrmhd-weak-form.md#conducting-wall-current-and-its-power). Its lift has the
prescribed boundary coefficients and zero interior coefficients. The
benchmark evaluates this power from its separate rate, ideal and resistive
pairings, and checks it against the constrained induction reaction. It is
not inferred from an energy difference, and it is not identified with
continuum Poynting flux. The instantaneous balances for this setup are

```text
dM/dt = 0,
dE_total/dt = P_wall,
d(E_kinetic+E_magnetic)/dt = C - Q_ohmic - D_viscous + P_wall.
```

Here `C` is the integral of pressure times velocity divergence, `Q_ohmic`
is the signed, source-accounted heating `eta J.(J-Js)`, and `D_viscous`
is the integral of `2 nu S:S`, with `S` the symmetric traceless velocity
strain. The current source remains the equilibrium current. No density,
heat, voltage or transport port is added.

Each output row includes mass, total energy, instantaneous wall power,
the total energy change minus the trapezoidal wall-power integral, the
mechanical power, and the mechanical energy change minus its trapezoidal
power integral. Time integration and power quadrature both contribute to
these integrated defects, so they need a time-step refinement study. A
nonzero wall power must be retained when assessing conservation and its
spatial limit must be checked before a continuum benchmark claim. All
reported quantities use the nondimensional units stated above.

The [retained time-step study](data/vrmhd-island-time.csv) uses the mixed
family with four cells in each direction, `eta=epsilon=0.01` and final
time 0.5. Reproduce each row with `probe mixed 4 4 0.01 DT STEPS 0.01`,
where `DT,STEPS` are its time step and step count. These are candidate
trajectories for the time-error check, not spatially admitted island rates.
The coarser steps are outside the asymptotic time regime. On the final two
halvings, the total-energy defect decreases with observed orders 2.11 and
2.19, while the mechanical-energy defect gives 2.12 and 2.20. The largest
absolute mass drift is below 7e-14. These measurements use GNU Fortran 16.2
and OpenBLAS in Release mode. The [instantaneous accounting data](
data/vrmhd-island-ledgers.csv) and [plot source](plot_vrmhd_island_ledgers.py)
retain the independent wall-power comparison.
