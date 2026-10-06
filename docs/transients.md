# Time-dependent extension

Status: planned secondary branch using the same operators as the primary
direct steady solver. Tiny kinetic and HM student dynamics enter the early
[laboratory gallery](student-laboratory.md); global kinetic/turbulent trajectory
certification remains later work. This extends the [shared evolution operator](../DESIGN.md#canonical-equation-interface).

## Equation and constraints

Use the same orbit, collision, moment, and source implementations for
`dot U=G(U,t)` and `G(U)=0`. In a Galerkin representation the transient driver
generally solves `W_h dot U_h=G_h(U_h,t)` with the relevant constraints.
Stationary potential elimination is a consumer-specific reduction of the
same field equations; it does not require identical solver unknowns.

Full transient Maxwell uses

$$
\partial_t B=-\nabla\times E,\qquad
\epsilon_0\partial_t E=\mu_0^{-1}\nabla\times B-j.
$$

Gauss's laws hold initially and must propagate consistently. Charge
continuity, source charge balance, boundary fluxes, and compatible discrete
divergence/curl are required. Prescribed time-varying external fields must
satisfy their Maxwell relations, including induced electric fields.

Finite element exterior calculus and the
[GEMPIC construction](https://arxiv.org/abs/1609.03053) provide references for
compatible field geometry. GEMPIC is a particle method; borrowing its
geometric identities does not make this deterministic harmonic solver a
GEMPIC implementation. Select the smallest compatible field representation
when self-consistent Maxwell becomes the next capability.

Each species model uses its same physical action in steady and transient
consumers; FK/GK/DK actions themselves differ. Reduced species require their
consistent dynamic polarization/magnetization and field constraints. Removing
electron gyromotion does not automatically remove fast Maxwell waves or
parallel streaming. State the selected field model before claiming a temporal
scale saving. Fixed reservoirs can exchange momentum/energy with the evolved
state and do not receive a transient distribution block.

Compute timescales from each frozen device case rather than reusing one
tokamak parameter set for TCV, SPARC, LHD, or reactor stellarators. Report the
certified horizon and requested widths for the selected full or reduced
model; long-time evolution remains a separate gate from steady existence.
[Species models](species-models.md) and [device envelope](devices.md) define
the corresponding physics scope.

## Physical timescales

The following are estimates, not promised run lengths or certified physical
inputs. Take `n_e=10^20 m^-3`, `B=2–6 T`, singly charged deuterium, and use
nonrelativistic bulk frequencies for illustration:

$$
\omega_{pe}=\sqrt{n_e e^2/(\epsilon_0m_e)},\qquad
\Omega_s=|q_s|B/m_s.
$$

| Process | Illustrative scale | Definition or qualification |
|---|---|---|
| Electron plasma oscillation | `1.11e-11 s` | Period `2pi/omega_pe`; inverse angular frequency is `1.77e-12 s` |
| Electron gyromotion | `5.95e-12–1.79e-11 s` | Period `2pi/Omega_e`; relativistic orbit period has an additional `gamma` factor |
| Deuterium gyromotion | `2.19e-8–6.56e-8 s` | Period `2pi/Omega_D` |
| Thermal electron transit over 10 m | `about 1e-7 s` | `length/sqrt(2 T_e/m_e)` for `T_e=10–30 keV` |
| Thermal deuterium transit over 10 m | `about 1e-5 s` | Same estimate at `T_i=10–30 keV` |
| Collision relaxation | Case-dependent | Depends on density, charge, energy, Coulomb logarithm, and the quantity relaxing |
| Heating and transport evolution | Milliseconds to seconds as a target | Depends on discharge, sources, boundaries, and instability dynamics |

These estimates follow from the displayed frequencies and thermal transit
definition. They distinguish periods from inverse angular frequencies.
One second contains roughly `10^11` electron plasma periods at this density:
about eleven decades of scale separation. There is no universal collision
time or fixed ordering valid for every case.

A direct stationary solve avoids traversing those periods. A transient
calculation must bound the components present in its initial data and forcing.
An implicit integrator can be stable across fast periods while inaccurately
damping them. The temporal certificate, rather than stability alone, decides
whether a requested interval is adequately resolved.

## Temporal representations

The first transient is a short homogeneous relaxation or prescribed-field
orbit problem with a known trajectory. It uses one temporal method and shares
the already validated physics. Later alternatives must solve a demonstrated
accuracy or cost problem:

- Fully resolved short transients for waves, plasma oscillations, and rapid
  instability growth.
- Implicit collision treatment and exponential magnetic rotations for stiff
  kinetic relaxation.
- Multirate or temporal spectral elements when slow observables can be
  enclosed without storing each fast oscillation separately.
- Space–time residual solves, `partial_t U-G(U,t)=0`, on bounded slabs when
  high-order temporal reconstruction proves worthwhile.
- Periodic-orbit problems with their own periodic boundary and phase
  constraints, rather than assumed attraction to a steady state.

For constant magnetic field and fixed momentum magnitude, gyromotion is a
known harmonic rotation. In a magnetic-axis basis the angular coefficients
acquire phase factors; in the chosen cylindrical basis apply equivalent
rotation matrices within each degree. For varying fields, electric
acceleration, and moving orbits, a frozen rotation is only a substep.
Commutators, splitting error, and field variation need bounds. Exact treatment
of this substep does not remove the need to resolve coupled plasma waves.

Coarse slow evolution is allowed only if a proved trajectory or observable
bound controls unresolved fast motion. Gyrokinetic or quasineutral reductions
are separately named models with separate model-error obligations if added;
they are not hidden effects of a large timestep.

## Structure and stability

For a closed compatible system, the intended structure is

$$
\dot U=J(U)DH(U)+K(U)DS(U),\qquad
J^*=-J,\quad K^*=K\succeq0.
$$

Energy conservation and entropy production also require the degeneracy
conditions `K DH=0` and `J DS=0`, appropriate domains, and no unaccounted
boundary/source work. Antisymmetry and positive semidefiniteness alone do not
prove both laws. Sources, radiation, walls, and imposed fields contribute
explicit balance terms. The
[metriplectic Landau discretization](https://arxiv.org/abs/1707.01801)
provides a reference for discrete invariants and entropy, with its own scope.

Preserving energy does not imply that distances between nonlinear solutions
are non-expansive. Entropy production does not by itself supply a quantitative
decay rate for perturbations. Plasma instabilities can amplify errors.
Long-time certification therefore needs a case-specific stability theorem,
which may use a coercive relative energy, hypocoercivity, a logarithmic norm,
or observable-adjoint estimates where their assumptions can be verified.

## Trajectory certificate

For a continuous reconstruction `ubar(t)` on `[0,T]`, bound the defect
`d(t)=dot ubar-G(ubar,t)` and initial uncertainty. A valid stability estimate
can take the form

$$
\|U^*(t)-\bar U(t)\|_X\le
e^{\int_0^t a(\tau)d\tau}r_0+
\int_0^t e^{\int_s^t a(\tau)d\tau}\eta(s)\,ds,
$$

where `eta` rigorously bounds the defect and `a` is a justified local growth
bound on the entire enclosed tube. Prove existence and continuation inside
that tube; a formal inequality alone does not establish either. Include all
spatial, momentum, angular, field, temporal, and roundoff complements. Bounds
at timestep endpoints alone do not enclose the trajectory between them.

Weak stability estimates can grow too quickly for useful long-time bounds.
Refinement, shorter slabs, validated stability improvements, or a narrower
observable theorem may help. Physical instability can make a tight enclosure
impossible at the desired cost. Report the certified horizon and widths;
do not promise second-long full-state enclosures for arbitrary reactor cases.

Temporal slab certificates must propagate the previous endpoint uncertainty.
Checkpoint restart and time-series compression preserve the reconstruction
or include a rigorous compression defect. Store selected observables,
checkpoints, temporal coefficients when useful, and certificate artifacts.
Avoid saving every 6D state at every step.

Chris&AI
