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
`probe FAMILY NY NX ETA DT STEPS [EPSILON]`, where `FAMILY` is `equal` or
`mixed`, `NY,NX` are cell counts, `DT` the constant time step, and `STEPS` the
number of steps. The default amplitude is 0.01. The viscosity coefficient is
set equal to `ETA`; because density varies, the kinematic viscosity and
magnetic Prandtl number vary across the slab. Output gives the cut width,
critical-point flux difference and maximum velocity coefficient at every step.
The existing lagged sparse factor is only a right preconditioner: Newton and
Krylov iterations still use the full nonlinear residual and tangent.

These trajectories are candidates until spatial and temporal refinement,
critical-point checks, physical ledgers and the required regime checks admit
a rate comparison. No Rutherford agreement is claimed by the setup checks.
