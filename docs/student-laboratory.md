# Steady kinetic laboratory and early turbulence gallery

Status: the first prescribed-field steady/free-streaming slice is implemented
([run and evidence](first-example.md)); other examples remain planned. Direct
steady kinetic science determines the architecture. Student dynamics and
toy turbulence exercise the same operators from the first visible release.
The interface has two small example groups, rather than separate solvers.

## Milestone 0 and shared equations

Build a deterministic CPU 1X1V laboratory, periodic in the slab coordinate
`x`, with a declared finite velocity interval and any whole-line tail scope.
The stationary residual is primary:

`R_f = v d_x f + (q/m) E d_v f - C[f] - S`,
`R_phi = -d_x^2 phi - rho/epsilon0`, `E=-d_x phi`.

For multiple populations `rho=sum_s q_s integral f_s dv + rho_background`.
Prescribed-field and self-consistent Poisson configurations are distinct.
Periodic Poisson needs exact net-charge compatibility and a mean-potential
gauge; do not silently remove incompatible charge. Use explicit SI-to-lab
normalization, species/velocity measures and source/reservoir balances.

The transient driver uses the same residual through `W_h dot U_h+R_h=0`
with algebraic field constraints, or an equivalent derived split. It is not
the method for finding every steady state. Keep a clean spatial/velocity
reconstruction, direct residual/JVP/VJP and a separate continuous certifier.
FortFEM supplies the smallest supported spatial representation; FortSym emits
all owned numerical leaves. CPU build/run/proof/plot needs no Python, GPU,
cluster, authenticated data or universal registry.

## First runnable steady increment

Use prescribed `E=0`, Krook `C=-nu(f-F)` with known `nu>0`, and
`f_exact(x,v)=F(v)(1+a cos(kx))`, `|a|<1`, `k=2 pi n/L`.
Define `S=(v d_x+nu)f_exact-nu F` symbolically. The residual is
`(v d_x+nu)f-nu F-S` on periodic x and bounded v. Solve this nontrivial
source-balanced state directly; the known formula is an independent oracle,
not the candidate solver. Its source/target make an open-system model.

Streaming is skew in the stated periodic pairing and `nu I` is coercive.
The Fourier resolvent `1/(nu+i k_n v)` establishes existence and the
continuous L2 residual bound `||f-fbar|| <= ||R_f(fbar)||/nu`.
An L-infinity bound at each v follows from the damped-characteristic inverse:
`sup_x |f-fbar| <= sup_x |R_f(fbar)|/nu`. Prove the applicable reconstruction
regularity/boundary assumptions and enclose residuals over whole cells before
plotting this band. Finite-domain velocity data needs no fictitious periodic
velocity wrap. Whole-line moments of a Gaussian require actual tail bounds.

The first increment has a numerical candidate, independent rigorous evaluator,
Lean-checked theorem instantiation, one bounded observable, two resolutions,
an invalid-certificate rejection, Fortran plots and a small timing/RSS record.
It is the first configuration of the shared engine. Do not postpone it behind
Landau, gyrofluid, GPU or a general symbolic compiler rewrite. The retained
bounded 0X1P FP proof case is an auxiliary coercive validator, not a separate
required first executable.

## Teaching sequence

Retain this presentation order as named examples; development follows the
steady spine and activates the smallest complete case rather than waiting
for the entire gallery. Each entry includes its equation, prediction,
numerical result, units, structure/error ledger and actual certificate status.

| Example | Prediction and purpose | Scope |
|---|---|---|
| 01 Free streaming | `f=f0(x-vt,v)`; phase-space shear and mass | Periodic x; exact flow and reconstruction error |
| 02 Constant electric field | `f=f0(x-vt+a t^2/2,v-a t)` | Prescribed uniform field; velocity support/inflow or tails |
| 03 Harmonic rotation | Exact rotation of `(x,v/omega)`; phase-area/invariants | Separate nonperiodic harmonic domain and time interval |
| 04 Poisson Fourier mode | `phi_k=rho_k/(epsilon0 k^2)`, `k!=0` | Compatibility, gauge, signs, field/energy enclosure |
| 05 Plasma oscillation | Linear frequency and initial-value response | Named linear limit, fit window and numerical/remainder errors |
| 06 Landau damping | Plasma-dispersion roots plus phase mixing | Initial-value/transient contributions and finite-amplitude scope |
| 07 Two-stream instability | Early linear growth, later trapping | Growth oracle applies before nonlinear saturation |
| 08 BGK equilibrium | Known trapped/passing distributions and consistent Poisson | Direct nonlinear residual, phase-space boundaries and constraints |
| 09 Krook relaxation/steady balance | Exact exponential and source-balanced state | Fixed reservoir; conservation/relative-entropy assumptions |
| 10 OU/FP | Gaussian mean/variance and stationary Gaussian | Velocity flux/tails; relative entropy and heat-bath exchange |
| 11 Lorentz scattering | Exact harmonic decay rates | 0X3V lab with normalized angular/radial measure |
| 12 Landau relaxation/steady balance | Equilibrium/invariants and independently benchmarked forced states | No general exact nonlinear relaxation is assumed |
| 13 HM/drift turbulence | Linear wave, bracket identities, energy/enstrophy | First early turbulence configuration |
| 14 Gyrofluid turbulence | Declared linear/FLR/free-energy targets | Actual closure and TIFF/Scott model overlap |

Uniform acceleration on periodic x is a prescribed electric drive; it cannot
come from a globally periodic scalar potential. A quadratic harmonic potential
likewise belongs to a separate domain with correct boundary/support treatment.
Never wrap velocity or harmonic x merely to reuse a periodic implementation.

For Landau damping, define thermal speed, Debye length and the dispersion
function convention before rigorous complex-root evaluation. A fitted root
is not the entire finite-time or nonlinear solution. Stationary `f=F(H)`
solves collisionless transport only when H is an exact invariant; charge,
Poisson, boundary and collision/source conditions remain. Use an actual
consistent BGK construction rather than claiming any F(H) is a plasma state.

Differentiate a smooth parameter such as amplitude, temperature or beam speed
from the first generated leaf, with an independent analytic/directional check.
Periodic k is an integer-mode choice, not a continuous parameter at fixed L;
varying L requires a chart/domain derivative. Later Poisson/BGK examples
demonstrate implicit state/adjoint equations with their branch assumptions.

## Turbulence door in the first visible release

Introduce the two-dimensional bracket
`[a,b]=d_x a d_y b-d_y a d_x b` from a shared FortSym definition.
Start with periodic HM `q=laplacian phi-phi` and
`d_t q+[phi,q]-kappa d_y phi=0`. With this convention its linear dispersion
is `omega=kappa k_y/(1+k_perp^2)`. The background-gradient term is essential:
the homogeneous equation with `kappa=0` has no such drift-wave frequency.
State the normalized sound-radius/time units and any zonal modification.

For the undriven ideal periodic model, declare
`E=1/2 integral (|grad phi|^2+phi^2)` and `Z=1/2 integral q^2`.
First verify a single linear Fourier wave and exact bracket identities; then
activate multimode nonlinear transfer and monitor the full discrete ledger.
Use alias-controlled products, the correct pairing and mean/gauge treatment.
Antisymmetry, cyclic integration and invariants are independent tests; none
alone proves a finite Poisson bracket's Jacobi identity. Full-time errors and
roundoff are measured separately from semi-discrete conservation.

The early gallery then grows through HW resistive drift waves, interchange/
blob dynamics, minimal ITG gyrofluid and electromagnetic drift-Alfvén models.
Derive complete equations, equilibrium-gradient/curvature drives, closure,
parallel/zonal response, FLR and boundary assumptions before implementation.
HW's large-response limit needs the matching zonal convention; an unqualified
`alpha->infinity` claim cannot select a different HM model automatically.
Each starts with a matching analytic dispersion/growth limit and an exact
energy/free-energy source, sink and boundary ledger. A small ITG model needs
its own validated closure; isothermal T3FF is not its thermal ITG oracle.

All gallery cases use common Laplace, bracket, parallel and closure actions,
with differing state/reduction/bracket definitions. Simple Fourier slab
actions can stay small while the main programme adds magnetic geometry.
Toy turbulence must not turn the primary API into `advance_one_timestep()`.
The first release includes the steady example and HM linear/nonlinear
demonstration, then expands the gallery incrementally; it does not need all
five turbulence models implemented before a useful release.

## Output and performance gates

Plot f(x,v), characteristic paths, fields, energy/entropy transfers and named
observable intervals. Analytic flows and linear waves get rigorous reference
evaluation and explicit solution-error bounds. A nonlinear movie without a
continuous trajectory theorem is labeled candidate; its conservation/error
ledger is useful evidence but not a fabricated interval certificate. Figures
distinguish candidate, reference and certified quantities. Where pointwise
bounds are unavailable, plot certified cell averages or bounded observables.

Fortran/fortplot produces snapshots and animations from local output. Keep
sources, recipes and numeric/certificate artifacts; generated graphics remain
untracked. Every integrated increment preserves the working certified default
and [demonstration contract](demonstrations.md). Reading the source should
let a student trace the selected equation through generated algebra, residual,
solver, certifier and plot without a framework hiding the mathematics.

A nominal `256x256` nodal phase grid has 65,536 state entries before species,
field, basis and work-vector costs; it is a scale target, not an instant-runtime
promise or a high-order DOF count. Measure laptop CPU runtime and attributed
memory, then a representative batch. Keep strict reference/optimized CPU
actions; add CUDA only when a real measured case warrants it. GPU availability
does not gate the first release. The gallery is also a derivative/performance
laboratory, with the [best-known frontier](performance.md) recorded by shape
and achieved accuracy.

## Progression beyond the slab

The primary next stages are steady 0X3V collision balance, 1X3V slab/cylinder
transport and classical coefficients, then NEO-2/SFINCS-level DK neoclassical
transport on real tokamak/stellarator inputs. Add gyromoments/GK, nonlinear
multispecies Landau/BB, global FK/Maxwell and deterministic neutral/atomic/wall
states under their own scoped certificates. Student transients and turbulence
remain a secondary early branch at every stage, retaining the same operators,
structural contracts, limits and physical observables.

Chris&AI
