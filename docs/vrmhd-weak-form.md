# Mixed-current visco-resistive MHD weak forms

The Cartesian leaves define the volume integrands of
`M(U) U_dot + R(U; theta) = 0` for one singly charged ion species and
quasineutral electrons. Electron mass is neglected in the mass density, so
number density is `n = rho/mi` and pressure is `p = n(Te + Ti)`.
The state `U` consists of density, velocity, vector potential, both
temperatures and the auxiliary current. `M` is its mass operator, `R` its
spatial residual, and `theta` denotes the closure parameters and prescribed
source data.

All quantities use SI units. Temperatures `Te`, `Ti` and reference temperature
`T0` are energies per particle in joules, not kelvin or electronvolts;
conversion is an input responsibility. Density `rho` is in kg/m³, ion mass
`mi` in kg, velocity `v` in m/s, vector potential `A` in T m, magnetic field
`B = curl A` in T, current `J` and prescribed current `Js` in A/m², and
permeability `mu0` in H/m. Resistivity `eta` is in ohm m and dynamic viscosity
`nu` (despite its symbol) is in Pa s. Heat diffusivities `chi` and particle
diffusivities `D` are in m²/s; heat flux therefore includes number density.

The fields are `rho, v, Te, Ti` in H1 (square-integrable values and first
weak derivatives), and `A, J` in H(curl) (square-integrable values and curls).
Finite-element exterior calculus (FEEC) uses the discrete de Rham complex
`V0 -grad-> V1 -curl-> V2 -div-> V3`, with scalars in `V0`, velocity in
`[V0]^3`, and `A, J` in `V1`. The strong discrete curl places `B` in `V2`;
the complex identity `div curl = 0` then enforces discrete `div B = 0`.
Here `V0` is in H1, `V1` in H(curl), `V2` in H(div) (square-integrable
values and divergence), and `V3` in L2 (square-integrable values).
Scalar tests are `r, qe, qi`, the velocity test is `w`, and the vector-potential
and auxiliary-current tests are `W, K`. The leaves receive only values,
gradients and curls at one quadrature point. `grad_v_ij = partial_j v_i`;
matrix arguments are ordered by column. No second trial derivatives occur.

## Mass and residual

We use conservative momentum `rho v` and species internal energies
`3 rho Ts/(2 mi)` while retaining primitive fields as coordinates. Their
chain-rule mass terms retain density coupling explicitly and avoid dividing
momentum or energy residuals by density.

The mass integrands, where a dot denotes a time derivative, are

```
rho_dot r
(rho v_dot + rho_dot v) . w
A_dot . W
3/(2 mi) (rho Ts_dot + rho_dot Ts) qs,  s = e, i.
```

The current equation is algebraic and has a zero mass row. With `:` denoting
sum over both matrix indices, the volume residual integrands are

```
Rrho = -rho v . grad r + (D grad rho) . grad r - Srho r
Rv   = -rho (v tensor v) : grad w - p div w
       + tau : grad w - (J x B) . w
RA   = -(v x B) . W + eta (J - Js) . W
RJ   = mu0 J . K - B . curl K
RTs  = -3/2 n Ts v . grad qs + n Ts div v qs
       + n (Ks grad Ts) . grad qs - Hs qs
He   = eta J . (J - Js) + Qei + Qe
Hi   = tau : grad v - Qei + Qi
Qei  = n rate (Ti - Te).
```

`Srho` is a mass source in kg/(m³ s); `Qe, Qi` are heat sources in W/m³.
The nonnegative `rate` is an energy-equilibration rate in 1/s, with its
normalization defined by `Qei` above. These are source inputs, not extra
species. Density diffusion and mass injection are prescribed model ports;
a future total-energy ledger must include their kinetic-energy exchanges.
No closed-system total-energy claim is made for these ports here.

Volume integration uses `dx`. These forms use periodic boundaries; they also
apply to a stationary ideal conducting wall with vanishing tangential electric
field, compatible zero-tangential tests for `A`, and no flow, diffusive or
thermal flux through the wall. Nonzero boundary data require their surface
terms and work ledger. Gauge fixing and symmetry reduction are separate work.

## Transport and closures

Define `b = B/sqrt(B.B + epsilon_B²)`, where positive `epsilon_B` in tesla
is an explicit model parameter, not a numerical clipping operation. For either
species, and separately for particle diffusion,

```
Ks = chi_parallel(Ts) b tensor b + chi_perpendicular (I - b tensor b)
D  = D_parallel b tensor b + D_perpendicular (I - b tensor b)
eta(Te) = eta0 (Te/T0)^(-3/2)
nu(Ti)  = nu0 (Ti/T0)^(-3/2)
chi_parallel(Ts) = chi_parallel0_s (Ts/T0)^(5/2).
```

`I` is the Cartesian identity. Each perpendicular coefficient is constant.
For the viscous stress, define the symmetric traceless velocity strain `S`
and a regularized traceless magnetic tensor `P`:

```
S   = (grad v + transpose(grad v))/2 - I div v/3
P   = b tensor b - I (b.b)/3
tau = 2 nu S + 3 nu_parallel (P:S) P.
```

The optional dynamic parallel viscosity `nu_parallel` is in Pa s and defaults
to zero. This definition stays traceless and gives nonnegative viscous heating
even when regularization makes `b` shorter than one. The isotropic stress has
zero bulk viscosity. The caller supplies positive `rho, mi, mu0, Te, Ti, T0,
epsilon_B`, nonnegative transport/viscosity references and rate, and finite
remaining inputs. Leaves do not clamp invalid states.

`vrmhd_parameters_t` holds `mi, mu0, t0, eta0, nu0`, electron/ion arrays
`chi_parallel0, chi_perpendicular` (electron first), `d_parallel,
d_perpendicular, nu_parallel, equilibration_rate, epsilon_b`. The caller must
set `mi` and `mu0`; dissipative coefficients and rate default to zero,
`T0` to one joule and `epsilon_B` to 1e-12 T.

## Mixed current and energy exchange

The auxiliary current is defined by weak Ampère's law `RJ = 0`, rather than
by a second pointwise derivative of `A`. With `W = J` in induction and
`K = A_dot` in Ampère, the discrete magnetic-energy derivative is

```
d/dt integral |B|²/(2 mu0) dx
  = -integral (J x B).v dx - integral eta J.(J - Js) dx.
```

Thus Lorentz work cancels the ideal induction work, and the resistive loss
equals electron heating using the same discrete `J`. At `Js = 0` the heating
is exactly `eta |J|²`. For nonzero `Js` we include `-eta J.Js` in electron
heating, treating it as net nonideal electrical work; it can be negative.
A prescribed drive's separate power deposition, if desired, belongs in the
heat-source inputs. Weak Ampère's law has no boundary contribution for the
periodic or compatible ideal-wall tests just specified.

The viscous force contribution is `-tau:grad w`, so with `w = v` it cancels
the ion heating `tau:grad v`. Both use the same stress, and its dissipation
is `2 nu S:S + 3 nu_parallel (P:S)^2`, which is nonnegative for nonnegative
coefficients. Pressure work cancels the two thermal compression terms.
For constant thermal tests the volume advection and conduction terms vanish,
while electron-ion equipartition cancels between species.

For a closed system without density diffusion or sources, these balances give
conservation of kinetic plus internal plus magnetic energy when the scalar
space contains constants and the advection assembly preserves the kinetic
continuity product. In the present conservative form that product is the
density equation tested with `|v|^2/2`: subtracting it from momentum tested
with `v` gives the kinetic-energy time derivative, and the advective integrands
cancel pointwise. Constants alone do not ensure that `|v|^2/2` lies in a
particular finite-element space; an assembled solver must preserve this
product through compatible spaces or a conservative split. The tests check
these local identities; full semidiscrete total-energy conservation remains
an assembly gate.

## Leaf interface and verification scope

The generator is [gen_vrmhd_leaves.f90](../gen/gen_vrmhd_leaves.f90); committed
[Fortran leaves](../src/generated/kin6d_vrmhd_leaves.f90) need no symbolic
runtime. Each routine is named `vrmhd_<term>` and has `_jvp` (Jacobian-vector
product) and `_vjp` (vector-Jacobian product) companions. Inputs, including
parameters and test data, precede the output `value`; tangent inputs follow
all primal inputs, and the adjoint weight precedes input duals. These are
Euclidean scalar leaf duals. A finite-element adjoint must supply quadrature
weights, assembly and the physical Riesz maps separately.

`mass_induction` uses generic vector argument names `v_dot_1:3` for `A_dot`;
`ampere` uses `w_1:3` for `K` and `curl_w_1:3` for `curl K`. `pressure_work`
returns `-p div w`. `transport` returns the tensor contraction for particle
diffusion; `conduction` additionally includes `rho/mi`. The Lorentz, viscosity,
ideal induction, resistive induction, Ohmic heating, viscous heating,
equipartition and source leaves return the signed right-hand-side terms above.
`RA` negates the sum of both induction leaves, and `Rv` negates the viscosity
and Lorentz leaves. Advection and compression leaves have residual signs;
the electron equation subtracts equipartition and the ion equation adds it.
Scalar `source` is used
for mass and both external heat sources. Temperature leaves use `temperature`
for either `Te` or `Ti`; closure leaves use the corresponding reference.

The `rm1` tests compare every primal with an independent numerical expression,
all directional and coordinate derivatives with centered differences, and
all adjoints with the Euclidean pairing. They verify a quadratic finite-
difference error decrease, proper-rotation covariance, conduction transverse
to `B`, stress dissipation, pressure work, and ideal, resistive and viscous
energy exchange.
For weak Ampère's law, an analytic trigonometric `A` and independently
differentiated current/test curls are integrated on a periodic cube: the
uniform eight-point rule in each direction is exact for these degree-two
trigonometric products, up to floating-point rounding. Regeneration is checked
against the pinned FortSym revision. These are pointwise and weak-form gates;
no assembled nonlinear solve or rigorous continuous certificate is claimed.

## References

The compatible magnetic spaces and weak-current energy pairing follow the
structure-preserving FEEC approach of
[Hu, Ma and Xu, Numerische Mathematik 135, 371–396 (2017)](https://doi.org/10.1007/s00211-016-0803-4)
and [Hu, Lee and Xu, Journal of Computational Physics 436, 110284 (2021)](https://doi.org/10.1016/j.jcp.2021.110284).
[Gawlik and Gay-Balmaz, Journal of Computational Physics 450, 110847 (2022)](https://doi.org/10.1016/j.jcp.2021.110847)
treat energy and helicity preservation for incompressible MHD. Their
compressible, barotropic variational formulation is instead in
[Journal of Plasma Physics 87, 835870501 (2021)](https://doi.org/10.1017/S0022377821000842).
These references motivate the discrete structure; they do not certify this
two-temperature extension or its future spatial and temporal assembly.
