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
`V0 -grad-> V1 -curl-> V2 -div-> V3`, with `A, J` in `V1`.
Density, temperatures and velocity use a separate
continuous Lagrange space `V_s` and `[V_s]^3`, respectively; the scalar
space `V0` of the magnetic complex remains the domain of its exact gradient.
The degree convention, quadrature and failed spectral qualification of the
proposed pair are stated in [the section contract](vrmhd-section.md).
The strong discrete curl places `B` in `V2`;
the complex identity `div curl = 0` then enforces discrete `div B = 0`.
Here `V0` is in H1, `V1` in H(curl), `V2` in H(div) (square-integrable
values and divergence), and `V3` in L2 (square-integrable values).
Scalar tests are `r, qe, qi`, the velocity test is `w`, and the vector-potential
and auxiliary-current tests are `W, K`, both in H0(curl), with zero
tangential trace at a conducting wall. The leaves receive only values,
gradients and curls at one quadrature point. `grad_v_ij = partial_j v_i`;
matrix arguments are ordered by column. No second trial derivatives occur.

## Mass and residual

Momentum uses the skew-symmetric kinetic-energy form, while species internal
energies remain `3 rho Ts/(2 mi)`. Primitive fields are coordinates. The density
rate in momentum is the same discrete derivative used in continuity; its
coupling to momentum is half the conservative chain-rule coupling.

The mass integrands, where a dot denotes a time derivative, are

```
rho_dot r
(rho v_dot + rho_dot v/2) . w
A_dot . W
3/(2 mi) (rho Ts_dot + rho_dot Ts) qs,  s = e, i.
```

The current equation is algebraic and has a zero mass row. With `:` denoting
sum over both matrix indices, the volume residual integrands are

```
Rrho = -rho v . grad r + (D grad rho) . grad r - Srho r
Rv   = rho/2 [(grad v)v . w - (grad w)v . v] - p div w
       + tau : grad w - (J x B) . w
       + Srho v.w/2 - (D grad rho).grad(v.w)/2
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
species. Density diffusion and mass injection are prescribed model ports. Write
`F = D grad rho` and `Q = Srho + div F`, so continuity is
`rho_dot + div(rho v) = Q`. The default convention is particles born at rest:
the conservative momentum equation has no accompanying momentum injection.
The skew left side equals the conservative left side minus `v Q/2` whenever
continuity holds, so its right side includes `-v Q/2`. Integration by parts
of that explicit diffusion coupling gives the last two volume terms in `Rv`.
The kinetic source power is `-Srho |v|^2/2 + F.grad(|v|^2)/2`, with boundary
power `-F.N |v|^2/2`. These exchanges are explicit ledger entries; they do not
change the channel flux, stress or electromotive-force closure semantics.

Volume integration uses `dx`. These forms use periodic boundaries; they also
apply to a stationary ideal conducting wall with vanishing tangential electric
field, compatible zero-tangential tests for `A`, and no flow, diffusive or
thermal flux through the wall. Nonzero boundary data require their surface
terms and work ledger.

With outward unit normal `N`, the boundary residual terms are

```
Rrho_boundary = (rho v - D grad rho).N r
Rv_boundary   = [rho v (v.N)/2 + p N - tau N + F.N v/2].w
RTs_boundary  = [3/2 n Ts v - n Ks grad Ts].N qs.
```

Induction has no integration-by-parts boundary term. A no-slip ideal wall
sets `v=0`, fixes `N cross A`, and fixes scalar Dirichlet values, so their
tests vanish at the boundary. The current is `J = J_wall + J0`, where `J0`
and its tests `K` belong to H0(curl). The prescribed lift `J_wall` has the
wall trace degrees of freedom and zero interior coefficients. Thus
`(N cross B).K` vanishes and the Ampere
residual has no surface term, even for nonzero tangential magnetic field.
Manufactured nonzero Dirichlet data supply the matching current trace.
These prescribed ports need their work ledger before a closed-system
energy statement applies.

For a stationary solve, the mixed Coulomb gauge adds
`integral grad(lambda).W dx` to induction and
`integral A.grad(q_lambda) dx = 0`. Multiplier `lambda` and test
`q_lambda` belong to `V_s` with zero boundary trace. A compatible continuous
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
eta(Te) = eta0 (Te/T0)^a_eta
nu(Ti)  = nu0 (Ti/T0)^a_nu
chi_parallel(Ts) = chi_parallel0_s (Ts/T0)^a_chi.
```

`I` is the Cartesian identity. Each perpendicular coefficient is constant.
The exponents default to the Spitzer and Braginskii scalings
`a_eta = a_nu = -3/2` and `a_chi = 5/2`. Benchmarks that prescribe a
temperature-independent viscosity or parallel conductivity, as the JOREK
cylinder runs do, set `a_nu = 0` or `a_chi = 0`; a run states its exponents
with its closure.
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
d_perpendicular, nu_parallel, equilibration_rate, epsilon_b`, the exponents
`eta_exponent, nu_exponent, chi_parallel_exponent` and `ohmic_fraction`
(see the energy exchange below). The caller must
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
primal and JVP leaves and retains the state-dependent mass coupling.

## Mixed current and energy exchange

Weak Ampère recovery makes `J_h` a discrete curl-curl of `A_h`. Its curl
is an additional recovered discrete derivative, for which this formulation gives
no convergence guarantee. Here `h` denotes mesh spacing. The diagnostic
interior `curl J_h` L2 error converges at approximately order
0.4–0.5 despite the prescribed current trace; no term of the formulation
contains `curl J`. Induction, Lorentz force and Ohmic
heating use only L2 pairings of `J`. The solved-current gate therefore uses
L2, with the same design order as `curl A` in L2 (one for the lowest-order
spaces). The curl error remains a diagnostic, reported over the whole mesh
and after excluding the two cell layers next to the boundary.

Stationary stability concerns the primal fields
`p=(rho,v,A,Te,Ti,lambda)`, with auxiliary current eliminated by the exact
algebraic recovery `mu0 M_J J = C^T B` on the free current coefficients after
subtracting the fixed lift's mass load.
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
by a second pointwise derivative of `A`. When `J` is an admissible induction
test, use `W = J` in induction and `K = A_dot` in Ampère. For zero voltage
and boundary power the discrete magnetic-energy derivative is then

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

`ohmic_fraction = f`, default one, multiplies the Ohmic term of `He`. With
`f < 1` the power `(1 - f) integral eta J.(J - Js) dx` leaves the system: it is
an explicit sink of the energy ledger, not a change of the induction or
Lorentz terms. Models without Ohmic heating of the temperature correspond to
`f = 0`.

The viscous force contribution is `-tau:grad w`, so with `w = v` it cancels
the ion heating `tau:grad v`. Both use the same stress, and its dissipation
is `2 nu S:S + 3 nu_parallel (P:S)^2`, which is nonnegative for nonnegative
coefficients. Pressure work cancels the two thermal compression terms.
For constant thermal tests the volume advection and conduction terms vanish,
while electron-ion equipartition cancels between species.

For `w=v`, the two advection contractions cancel algebraically at each
quadrature point, and the momentum mass is exactly
`d/dt integral rho |v|^2/2 dx`. No continuity test with `|v|^2/2`, integration
by parts or exact quadrature is required for this identity. It holds with
the rotation and screw metric weights as well. Constant thermal tests remove
thermal advection and conduction from the global internal-energy balance.
The price of this energy-compatible split is loss of discrete global momentum
conservation: a constant velocity test no longer gives a conservative flux
balance when the discrete continuity residual does not vanish pointwise.

The former conservative momentum form left the energy defect
`-integral |v|^2/2 [rho_dot + div(rho v)] dx`. Degree-one scalar tests do not
contain the piecewise quadratic `|v|^2/2`; a converged mass solve could not
remove this defect. The closed hydrodynamic subcase of `vrmhd_section_energy`
solves the production mass equation for random positive scalar fields and
no-slip random velocity, then independently differentiates physical energy.
It requires relative closure at `1e-12` before the magnetic ledger is tested.
Scalar tests remain unrestricted, retaining constants and natural no-flux
boundaries. The independent stationary manufactured source is derived from
the conservative strong balance minus half the continuity defect times `v`;
its prescribed volume loads are separate from the `Srho` port.

### Conducting wall current and its power

For a stationary conducting wall, let the subscript `t` denote tangential
components. With `v.N=0` and `B.N=0`, both vectors are tangent to the wall,
so `(v cross B)_t=0`. Ohm's law and `N cross (E - V h)=0` therefore give

```
eta (J - Js)_t = (V h)_t,
J_t = Js_t + (V h)_t/eta,    eta > 0.
```

The trace is prescribed boundary data. For ideal studies (`eta=0`), it is
instead supplied by the given equilibrium; division by eta is not used.
The wall resistivity, including its temperature dependence, is part of those
prescribed data. Perturbations at fixed data have zero current trace.
The same total `J = J_wall + J0` enters induction, Lorentz force and heating.
Ampere is tested only with `K` in H0(curl), just as induction uses `W` in
H0(curl). This removes the former unconstrained-current wall mismatch.

Write `F = Atilde_dot - v cross B + eta (J - Js) - V h` for the induction
balance in V/m, and use `(a,b)=integral a.b dx`. Testing Ampere with
`K=Atilde_dot` gives `dE_mag/dt=(J,Atilde_dot)`. Testing induction with `J0`
gives `(J0,F)=0`; consequently

```
dE_mag/dt = (J,v cross B) - (J,eta (J - Js)) + (J,V h) + P_wall,
P_wall = (J_wall,Atilde_dot) - (J_wall,v cross B)
         + (J_wall,eta (J - Js)) - (J_wall,V h).
```

`E_mag` is magnetic energy in joules and `P_wall` is the explicit power of
the prescribed affine current constraint, in watts (per unit orbit coordinate
in the section). This is a discrete lift pairing, not an added Ampere surface
residual or an identification with continuum Poynting flux. It vanishes for
homogeneous current data. It is evaluated from the four displayed integrals,
independently of the energy defect; contraction with the constrained induction
rows checks the same identity. No power is inferred by subtracting the other
ledger entries. The global ledger includes `P_wall` in external power, while
all exchanges use the total current and retain their local cancellation.

The random-state ledger keeps compact magnetic support, no-slip velocity and
natural scalar no-flux tests. It prescribes a nonzero current lift: the
current-source case supplies its trace through `Js`, while the voltage case
represents `V h/eta` with constant wall electron temperature in the current
space. Finite-element trace data are coefficient data; the analytic voltage
trace supplies horizontal edge line integrals and orbit nodal values. Interior
current coefficients alone solve Ampere. The test independently differentiates
physical energy and evaluates every power at the production quadrature points.
It requires nonzero wall power in all six symmetry/drive combinations.

The literature supports the homogeneous wall structure. In
[Gawlik and Gay-Balmaz's 2021 author manuscript](https://math.hawaii.edu/~egawlik/pdf/GaGB2021.pdf),
p. 1, equation (1.5) imposes zero velocity, zero normal magnetic field and zero
tangential curl of B for positive resistivity. Section 3, p. 7 defines the
auxiliary Nedelec space inside H0(curl); equations (3.13), (3.18) on p. 9 and
(3.27), (3.32) on p. 10 recover current using that space, without a surface
load. The ideal boundary statement below (1.6) does not impose the resistive
current condition. In the
[Hu, Ma and Xu preprint of the 2017 paper](https://arxiv.org/pdf/1410.1095),
p. 2, equations (1.2), (1.5)-(1.7) give `j=E+u cross B` in their normalization,
zero velocity, zero normal B and zero tangential E. Section 2.2, p. 6 places E
in H0(curl); for smooth traces Ohm's law then implies `j_t=0`. This does not
assert that their discrete pointwise product `u_h cross B_h` itself lies in
H(curl). The nonzero driven lift and its power above are the present extension,
not a formula verified in either paper. The related incompressible 2022 paper's
wall treatment has not been checked here.

## Leaf interface and verification scope

The generator is [gen_vrmhd_leaves.f90](../gen/gen_vrmhd_leaves.f90); committed
[Fortran leaves](../src/generated/milca_vrmhd_leaves.f90) need no symbolic
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
for mass and both external heat sources. `momentum_mass_source` and
`momentum_density_diffusion` return the signed residual couplings for
particles born at rest. Temperature leaves use `temperature`
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
