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
terms and work ledger.

With outward unit normal `N`, the boundary residual terms are

```
Rrho_boundary = (rho v - D grad rho).N r
Rv_boundary   = [rho v (v.N) + p N - tau N].w
RJ_boundary   = -(N cross B).K
RTs_boundary  = [3/2 n Ts v - n Ks grad Ts].N qs.
```

Induction has no integration-by-parts boundary term. A no-slip ideal wall
sets `v=0`, fixes `N cross A`, and fixes scalar Dirichlet values, so their
tests vanish at the boundary. Auxiliary-current tests remain unrestricted:
the Ampere surface term is required when the magnetic tangential trace is
nonzero. Manufactured nonzero Dirichlet data use the same trace lifting.
These prescribed ports need their work ledger before a closed-system
energy statement applies.

For a stationary solve, the mixed Coulomb gauge adds
`integral grad(lambda).W dx` to induction and
`integral A.grad(q_lambda) dx = 0`. Multiplier `lambda` and test
`q_lambda` belong to H1 with zero boundary trace. A compatible continuous
physical solution has `lambda=0`. The time-dependent problem retains Weyl
gauge; its A mass row removes the stationary gradient kernel.

Two driven-state conventions are specified. The current-source convention
prescribes `Js` and zero loop voltage, as in JOREK comparisons. The
loop-voltage convention prescribes voltage `V` and an invariant, closed,
non-exact one-form `h` with loop integral one. Rotation uses
`h=grad(phi)/(2 pi)`; translation and screw symmetry use `h=grad(z)/L_z`,
with axial period `L_z`. In Weyl gauge set `A=Atilde-t V h`, giving

```
B = curl Atilde
RA = [Atilde_dot - V h - v cross B + eta (J - Js)].W.
```

Thus `Js=J` gives stationary Ohmic induction at `V=0`; with `Js=0` it
requires `eta J=V h`. A full stationary state must also balance electron
heating: the latter convention retains `eta |J|^2` when no thermal sink is
prescribed. Voltage drive supplies magnetic power
`integral V h.J dx` without changing the specified Ohmic heating term.
The screw field `xi/|xi|^2` has nonzero curl and cannot replace `h`.
Discrete harmonic coefficients are edge integrals of the multivalued
coordinate divided by its period, and their discrete curl must vanish to
roundoff. A projection with a small nonzero curl is insufficient. The
section wrapper's implementation and remaining harmonic-drive construction
are stated in [the section contract](vrmhd-section.md).

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

## Flux-limited parallel heat transport

`flux_limited(s)` selects a harmonic-mean limiter for species `s`, with
`g = b.grad(Ts)`, positive `gradient_floor(s)` in J/m, positive
`alpha = flux_limit(s)`, and `v_th = sqrt(2 Ts / mass_s)`. The electron mass
is `electron_mass`; the ion mass is `mi`. With `chi_Sp` the Spitzer value,

```
chi_eff = chi_Sp / (1 + chi_Sp sqrt(g^2 + gradient_floor^2)/(alpha Ts v_th)).
```

Thus the parallel part of the conductive flux has magnitude at most
`alpha n Ts v_th`, including the regularized `|b| <= 1`. Density cancels in
the ratio of classical to free-streaming flux. Perpendicular heat transport
retains its separate coefficient. The positive gradient floor makes the
coefficient differentiable at zero parallel gradient; it is a declared model
parameter, not a clipping operation. The harmonic-mean interpolation follows
[Miceli et al., A&A 458, 213 (2006), Section 2.2](https://arxiv.org/abs/astro-ph/0608063),
which attributes it to Dalton and Balbus; the cap here is specified per species
as `alpha n Ts v_th`, rather than their hydrodynamic sound-speed cap.
The `rm2` closure gate checks the harmonic mean, flux bound and derivatives.

The invariant [section volume assembly](vrmhd-section.md) calls these same
primal and JVP leaves and retains the conservative mass coupling.

## Mixed current and energy exchange

Weak Ampère recovery makes `J_h` a discrete curl-curl of `A_h`. Its curl
is an additional recovered discrete derivative, for which this formulation gives
no convergence guarantee. Here `h` denotes mesh spacing. The measured `curl J_h`
L2 error converges at approximately order 0.4–0.5 both over all cells and away
from the boundary, so this slow rate is not a boundary-layer effect; no term
of the formulation contains `curl J`. Induction, Lorentz force and Ohmic
heating use only L2 pairings of `J`. The solved-current gate therefore uses
L2, with the same design order as `curl A` in L2 (one for the lowest-order
spaces). The curl error remains a diagnostic, reported over the whole mesh
and after excluding the two cell layers next to the boundary.

Stationary stability concerns the primal fields
`p=(rho,v,A,Te,Ti,lambda)`, with auxiliary current eliminated by the exact
algebraic recovery `mu0 M_J J = C^T B` including prescribed boundary terms.
Here `M_J` is the L2 current mass matrix and `C^T B` denotes the weak curl
load. If `D=R_U` is the stationary Jacobian, the reduced Jacobian is
`S=D_pp-D_pJ (mu0 M_J)^(-1) D_Jp`. Measuring independent oscillatory current
perturbations through the dual H(curl) induction norm instead introduces an
artificial factor of `h`. For constant resistivity the eliminated resistive
block is `eta/mu0 (curl A,curl W)`, coercive in H(curl) modulo gradients,
which the Coulomb multiplier constrains. The numerical gate uses H1 Riesz
matrices for the scalar, velocity and multiplier fields and H(curl) for `A`
on both sides of `S`; it retains a positive floor and rejects more than 25%
loss of the minimum singular value per refinement. It also reports the
maximum singular value, the minimum mode's field composition and the
condition number of `M_J`. This is finite-mesh evidence, not a proof of a
uniform continuous inverse bound.
The mass gate checks the condition number of `D_J^(-1/2) M_J D_J^(-1/2)`,
where `D_J=diag(M_J)`, because section edge and orbit nodal bases have different
mesh scaling; the raw condition number is diagnostic only.

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
Closest in purpose are the structure-preserving compressible visco-resistive
MHD scheme of [Fambri and Sonnendrücker (arXiv:2407.14141)](https://arxiv.org/abs/2407.14141)
and the FEEC fields of the hybrid code Struphy (Holderied, Possanner and Wang,
J. Comput. Phys. 433, 110143 (2021)).
These references motivate the discrete structure; they do not certify this
two-temperature extension or its future spatial and temporal assembly.
