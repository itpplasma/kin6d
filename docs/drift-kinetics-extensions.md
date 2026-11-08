# Drift-kinetic extensions: adaptive velocity space, finite orbits, nonlinear resonances

Status: plan (2026-11-08). It extends [drift-kinetics.md](drift-kinetics.md): its local tokamak
ansatz is the current thin-orbit model, and its K2 section (radial drift acting on the unknown) is
the global problem used here. The work proceeds in the steps below; a step starts only when the
gate of the step before it has passed, and a failed gate is recorded as a negative result.

## Equation

All modes solve the steady guiding-centre drift-kinetic equation in the perturbed field
B = B0 + delta B, with B0 axisymmetric and delta B the non-axisymmetric perturbation (toroidal mode
number n). psi is the poloidal flux per radian labelling the unperturbed surfaces, theta and phi
are the Boozer angles, and b = B/|B|. The unknown is f(psi, theta, phi; H, mu, sigma), with total
energy H = m v^2/2 + Z e Phi(psi), magnetic moment mu and sigma = sign(v_par). With f = F_M + f1
and the local Maxwellian F_M(psi, H),

    (v_par b + v_d) . grad f1 |_(H, mu)  -  C_L[f1]  =  -(v_d . grad psi) d_psi F_M |_H
    d_psi F_M |_H = F_M [ n'/n + Z e Phi'/T + (x^2 - 3/2) T'/T ]

v_par = sigma sqrt(2 (H - Z e Phi - mu |B|)/m) and the drift v_d (E x B, grad B, curvature) are
taken in the full field; v_d0 and v_d1 are its parts from B0 and delta B. C_L is the linearised
Fokker-Planck operator of all species, x = v/v_th the normalised speed and xi = v_par/v the pitch.
The equation is linear in f1 and contains delta B to all orders. The torque density is
<int d^3v p_phi_dot f1>, with <.> the flux-surface average and p_phi = m v_par B_phi/|B| + Z e psi the
canonical toroidal momentum (B_phi the covariant toroidal component, psi = A_phi); for thin orbits it reduces to Z e Gamma_psi with the radial particle
flux Gamma_psi. Expanding about F_M(psi, H) assumes the orbit width is small against the gradient
lengths.

With `orbit_width` the background is the canonical Maxwellian F_c(psi*, H), F_M with psi
replaced by psi* = p_phi/(Z e), which is constant along unperturbed orbits (NEO-RT and POTATO use the same background). It solves the collisionless
unperturbed equation exactly, so the equation becomes

    (v_par b + v_d) . grad f1 |_(H, mu)  -  C_L[f1]  =  C[F_c]  -  p_phi_dot d_p_phi F_c

with p_phi_dot = -d_phi H from delta B alone. The orbit-width ordering then no longer enters the
background, and the drive pairs with the torque <int d^3v p_phi_dot f1>. Expanded to first order in
the poloidal gyroradius over the gradient length it reproduces the drive of the thin-orbit
model, which is a test of step 3. Two caveats: the profiles n(psi*), T(psi*), Phi(psi*) are not the
measured profiles n(psi), T(psi), which follow from moments of F_c, so device inputs need this map;
and where F_c is not close to a local Maxwellian (pedestal), C is evaluated on F_c itself rather
than linearised about F_M.

The modes are subsets of the terms of the first equation:

| Mode | Terms kept | Structure |
|---|---|---|
| thin-orbit quasilinear (current) | streaming along b0, tangential drift (magnetic precession omega_B, E x B rotation Omega_tE), delta B to first order: the n = 0 response f10 and one harmonic e^{i n phi} | local per surface, unknowns (theta, xi, x) |
| `orbit_width` | + v_d0 . grad psi d_psi f1 | global in psi, harmonics separate |
| `nonlinear_orbits` | + delta B in the operator (v_par delta b, v_d1, delta|B| in v_par), including v_d1 . grad psi d_psi f1 | global in psi, harmonics k n coupled to (k +- 1) n, truncated at abs(k) <= K |
| both | all terms | full equation |

Resonances satisfy n (Omega_tE + omega_B) + l omega_b = 0 with omega_b the bounce or transit
frequency; l = 0 for trapped particles is the superbanana resonance. In `nonlinear_orbits`
resonant particles are trapped in islands of the resonant angle, which attenuates the
quasilinear torque. For an isolated resonance this attenuation is a universal function Theta(D)
of the dimensionless resonant diffusivity D (collisional diffusion across the island against its
nonlinear frequency), with Theta -> 1 for large D (Albert, PhD thesis, TU Graz 2020, Sec. 4.9, as
implemented in NEO-RT). Here Theta(D) is a limit to verify against; the solver computes the
attenuated distribution directly. Where islands of neighbouring resonances overlap (Chirikov
criterion), the island barriers vanish and transport returns towards the quasilinear value.

The unknowns grow by the number of radial points for `orbit_width` and by 2K + 1 for
`nonlinear_orbits`. The global systems are solved by GMRES preconditioned with the per-surface
factorisations, which works while the psi and harmonic couplings are weak. The cost risk is the
superbanana regime at very low collisionality, where layers at the island separatrix must be resolved.

Flux coordinates describe the field and label the outputs in every mode. With `orbit_width` the
phase-space coordinates are the invariants: psi* replaces psi, orbits are lines of constant
(psi*, H, mu), and the radial drift disappears from the operator, while collisions, which change
v_par at fixed position, and delta B carry the radial coupling. Moments map back to psi through
the orbit-weighted moment map. The limits are the magnetic axis, where orbit classes change, and
orbits that cross the last closed surface, which leave this solver and are handed to the edge
treatment.

## Two representations of one term list

For axisymmetric fields with perturbations, the full particle motion is integrable with the
actions H, p_phi and the gyro-action J_g = m mu/(Z e) (adiabatic), and the angles advance at the
canonical frequencies Omega(J) of real orbits (as in NEO-RT and POTATO). In these action-angle
variables the unperturbed operator is diagonal in the angle harmonics, i m.Omega(J); finite orbit
width is exact, finite Larmor radius enters through the gyro-action (Bessel factors of k_perp rho),
a perturbation H_m(J) e^{i(m.theta - omega t)} couples harmonics only through resonances
m.Omega = omega, and superbanana and resonance islands are pendulums in (J, resonant angle). The
same operator covers NTV (omega = 0), RF heating and Alfven-mode transport (omega != 0). Collisions
become the non-diagonal part.

This code implements the coordinate form only (flux or canonical coordinates, all
collisionalities, collisional boundary layers, stellarators); NEO-RT and POTATO are the
action-angle references at low collisionality, and a second implementation here is deferred. Separatrix layers, where Omega has its logarithmic
singularity and collisions dominate, stay in the coordinate form. Stellarators have no global
actions; at low collisionality the bounce-averaged second invariant J_par (as in KNOSOS) takes
their place. 
## Velocity-space mesh

The hp mesh (finite elements adapting both element size h and polynomial degree p) places element
edges on level sets of invariants of the motion, so that collisional layers lie on edges and
convergence stays spectral at low collisionality. The invariant, and with it the mesh, depends on
the mode:

| Mode | Layers | Fitted mesh |
|---|---|---|
| axisymmetric, one energy (prototype) | trapped-passing separatrix; lambda = (1 - xi^2) B_ref/|B| is constant along orbits | (theta, xi) elements bounded by lambda level curves, collapsing at the B maximum where the separatrix closes |
| thin-orbit quasilinear, all energies | separatrix and the resonance curves lambda_res(x) | one mesh per speed node, with projection between node meshes for energy scattering |
| `orbit_width` | orbits on level sets of p_phi at fixed (H, mu); the separatrix is a surface in (psi, theta, xi) | curved elements fitted to p_phi; the topology changes near the axis (potato and stagnation orbits) |
| `nonlinear_orbits`, isolated resonance | island separatrix in (p_phi, resonant angle) | local mesh per resonance, on top of the mesh of the mode without it |
| `nonlinear_orbits`, overlapping resonances | no invariant | unfitted, refined by dual-weighted residuals; convergence algebraic |

The direct sparse factorisation of the prototype is limited to two or three dimensions; global
modes use the iterative solver above. If hp fails its gate (step 1), the later steps use the
Legendre-Fourier solver of the current model with radial elements.

## Steps and gates

Every step lands as its own pull requests with tests in the repository's test suite on synthetic
data. Each gate lists the verification (manufactured solutions, conservation, analytic limits,
resolution ladders), the benchmarks against other codes, and the cost record (wall time,
core-hours, peak memory against the reference at matched accuracy). nu* is the normalised
collisionality. Reference codes: NEO-2 (field-line-following drift-kinetic solver), NEO-RT
(resonant transport in action-angle variables, thin orbits), POTATO (collisionless resonant
transport with finite orbit width), PERFECT (global Eulerian drift kinetics for pedestals,
Landreman et al. 2014), rabe (bootstrap coefficients at low collisionality), GORILLA and SIMPLE
(guiding-centre Monte Carlo). Device cases are ASDEX Upgrade discharge 30835 and the ITER
scenario TC24; their data stay outside the repository.

| Step | Content | Gate |
|---|---|---|
| 0 | thin-orbit quasilinear model with the full operator, device profiles | passed: torque profiles within 5% of NEO-2 with 10 Laguerre energy basis functions on both devices. Open: energy convergence of the AUG ion full operator at resonances |
| 1 | hp prototype: axisymmetric, monoenergetic, (theta, xi) mesh fitted to lambda | spectral convergence down to nu* ~ 1e-6; banana, plateau and Pfirsch-Schlueter limits of D11, D31, D33; axisymmetric bootstrap coefficient against the Sauter fit and NEO-2; faster or leaner than the Legendre solver at the same accuracy |
| 2 | hp for the thin-orbit quasilinear model with all energies (after a passed step 1) | agrees with the spectral full-operator kernel; 1/nu, sqrt(nu)-nu, superbanana-plateau and rapid-rotation limits; converges the AUG ion resonances where uniform speed elements did not |
| 3 | `orbit_width` on the canonical Maxwellian: radial coordinate psi* = p_phi/(Z e) at fixed (H, mu), so that unperturbed orbits stay at one radial node and only collisions and delta B couple radii; orbit classes (trapped, passing, potato) as in POTATO; the Eulerian form in psi (as PERFECT) as an independent check; outer boundary at the last closed surface absorbing outgoing orbits; GMRES with surface factors | orbit width -> 0 reproduces step 0 with the expected scaling; first-order expansion of F_c reproduces the thin-orbit drive; moments map n(psi*) to n(psi); manufactured solution; particle and p_phi conservation; GMRES iterations bounded under psi refinement; POTATO in its thin-orbit limit and at finite orbit width; PERFECT pedestal case |
| 4 | validity indicators per surface, species and energy: D; island width over resonance spacing; d(resonant frequency)/d p_phi, which sets the island width; orbit width over gradient length | D against NEO-RT and canonical frequencies against POTATO on shared cases |
| 5 | `nonlinear_orbits` for one isolated resonance, thin orbit | delta B -> 0 reproduces step 0 with torque proportional to delta B^2; attenuation against NEO-RT Theta(D) over D; Theta -> 1 in the quasilinear limit |
| 6 | `nonlinear_orbits` with `orbit_width` for isolated resonances | each switch off reproduces steps 3 and 5; resonance shifts against POTATO canonical frequencies |
| 7 | coupled harmonics, overlapping resonances | isolated limits reproduce steps 5 and 6; Monte Carlo with finite orbits and the perturbation in the orbits; cost below Monte Carlo at the same accuracy. Overlap is expected in the pedestal (sheared E_r, large delta B), on the canonical Maxwellian background |
| 8 | device curves | AUG 30835 and ITER TC24 with each switch alone and combined, against NEO-2, NEO-RT, POTATO and Monte Carlo where they apply, with the indicators of step 4 |

## NEO-2 feature parity

The steps above extend the physics beyond NEO-2. Replacing NEO-2 also needs its remaining
features, each checked against NEO-2 on the same inputs before it counts as covered:

| NEO-2 feature | Status |
|---|---|
| quasilinear NTV torque, all species, full linearised operator, both quasilinear sources (NEO-2-QL) | covered (step 0) |
| E_r from the radial force balance of one species with prescribed toroidal rotation (`isw_calc_er = 1`) | missing; input-side computation |
| non-axisymmetric heat flux (`isw_qflux_NA`) | missing; needs the K_2j coefficients next to K_1j |
| axisymmetric transport matrix with the full operator: bootstrap current, parallel conductivity, Ware pinch | the n = 0 solve exists inside the NTV kernel; not yet an output |
| relativistic collision operator (Braams-Karney) and the generalised Spitzer function for current drive | missing |
| misaligned-potential drive: radial E x B drift of the potential harmonic not aligned with the perturbed surfaces, Phi^MA_m = Phi_m + E_0r Delta_m (Delta_m the surface displacement), acting on d f0/dr, with the local parallel current over the surface (NEO-2 pull request 240, draft) | missing; one more right-hand side of the kernel, which already resolves theta. Check: a corrugation without misaligned potential gives zero parallel current at any collisionality with a momentum-conserving operator; collisionless limit against KiLCA |
| quasineutral Phi_1 on the surface from the species responses (planned NEO-2 stack after 240, SFINCS benchmark) | missing; a bordered unknown per angle in the multispecies solve |
| stellarator transport coefficients (NEO-2-PAR), including ripple trapping in local wells of the full |B| | K1 in [drift-kinetics.md](drift-kinetics.md); bootstrap against rabe (near-omnigenous fields, Lorentz, monoenergetic) and down to the Shaing-Callen limit belongs here |

## Complete scope

The goal is one code for steady neoclassical drift kinetics on closed flux surfaces, tokamaks and
stellarators, with every effect below as a switch that can be combined with the others and checked
against its limit. It supersedes the out-of-scope list of [drift-kinetics.md](drift-kinetics.md)
for sonic rotation, magnetic islands and kinetic-MHD feedback. Each effect counts as available only
after its check has passed.

| Area | Effect | Stage | Check |
|---|---|---|---|
| geometry | axisymmetric tokamak, up-down asymmetric shaping, finite beta from EFIT/VMEC/GVEC via Boozer files | 0 | NEO-2, NEO |
| geometry | non-axisymmetric perturbation of a tokamak, several n, plasma response from MARS, GPEC or MEPHIT as input | 0, 6 | NEO-2-QL |
| geometry | 3D stellarator field, ripple wells of the full |B|, tokamak TF ripple | K1 | DKES, MONKES, SFINCS, NEO-2-PAR; ripple-plateau and ripple-trapped nu limits |
| geometry | magnetic islands with helical flux coordinates: bootstrap and polarisation currents in the island | after 7 | island-width asymptotics (flattening, small-island threshold) |
| drifts | precession, E x B with E_r, magnetic shear term, tangential magnetic drift in stellarators | 0, K1 | SFINCS and KNOSOS tangential-drift cases |
| drifts | finite orbit width on the canonical Maxwellian, potato orbits near the axis | 3 | POTATO, PERFECT |
| drifts | perturbation in the orbits: superbanana regime, resonance islands, overlap | 5-7 | NEO-RT Theta(D), Monte Carlo |
| rotation | slow rotation (Mach number M << 1, current ordering) | 0 | NEO-2 |
| rotation | sonic rotation M ~ 1: rotating-frame Maxwellian, centrifugal and Coriolis drifts, centrifugal trapping, poloidal density asymmetry with Phi_1 | new | Hinton-Wong and Catto-Bernstein-Tessarotto theory; NEO (Belli-Candy) sonic option; heavy-impurity asymmetry |
| rotation | strong electric field, v_E ~ v_th: guiding-centre equations from the Lagrangian in the frame drifting with v_E (Littlejohn 1983; Cary and Brizard, Rev. Mod. Phys. 81, 693 (2009); Brizard 1995 for rotating tokamaks), with the centrifugal, Coriolis and polarisation-type drifts and the modified magnetic moment; poloidal E x B comparable to the poloidal projection of parallel streaming in the pedestal, which changes trapping | new, with sonic rotation | Kagan-Catto pedestal banana-regime theory; NEO sonic option; weak-field limit reproduces the current drifts |
| potential | prescribed E_r; E_r from force balance; ambipolar E_r with ion and electron roots | 0, K1, K2 | NEO-2 (force balance), SFINCS, DKES-based root finding |
| potential | Phi_1 on the surface from quasineutrality; misaligned-potential drive and local parallel current; delta Phi_n of the perturbation | parity | NEO-2 pull request 240, KiLCA, SFINCS Phi_1, EUTERPE/KNOSOS |
| potential | inductive E_par: Ware pinch, parallel conductivity | parity | Sauter fits, NEO-2 |
| perturbation | |B|-only drive; full delta H with field-line bending kappa.xi and b.delta B | 0, (b) of drift-kinetics.md | NEO-2-QL, MARS-K delta W_k |
| perturbation | rotating perturbation omega != 0 (rotating modes, locking) | (b) | analytic resonance shift, MARS-K |
| perturbation | kinetic delta W_k and torque fed back to the MHD response | new | MARS-K, MISK |
| collisions | Lorentz, momentum correction, full linearised Fokker-Planck, all species and impurities, pair Coulomb logarithms | 0 | NEO-2, NEO, SFINCS |
| collisions | relativistic Braams-Karney operator | parity | NEO-2-PAR |
| collisions | classical (gyro-) transport of impurities next to neoclassical | new | Braginskii, Buller-Helander stellarator impurity flux |
| background | Maxwellian; canonical Maxwellian; slowing-down fast ions and alphas | 0, 3, new | ASCOT, SIMPLE Monte Carlo |
| ordering | second order in the gyroradius over the scale length: neoclassical momentum flux and intrinsic rotation in up-down symmetric tokamaks, which vanish at first order | new | Parra-Catto higher-order theory; global finite-orbit-width limit of step 3 |
| gyrokinetic | finite Larmor radius in the steady linear response: gyroaverages J0, J1 (k_perp rho) of delta Phi, delta A_par, delta B_par, polarisation density in quasineutrality, FLR terms of the collision operator (which also give classical transport); radially nonlocal kernel where k_perp rho ~ 1. Fields closed by Poisson and Ampere in the kinetic layer strips of [response-unified.md](response-unified.md); no turbulence | R-K strips, new | KiLCA and KIM of KAMEL in the cylinder limit; drift-kinetic limit k_perp rho -> 0 |
| sources | neutral friction (charge exchange, ionisation), beam particle, momentum and heat sources, orbit loss through the last closed surface as the outer boundary of the global solve | K2, new | edge E_r and rotation analytics; Monte Carlo orbit loss |
| regimes | banana, plateau, Pfirsch-Schlueter; 1/nu, sqrt(nu)-nu, nu, superbanana plateau, superbanana, bounce and transit resonances, ripple plateau; nu* down to 1e-6 and up to the Pfirsch-Schlueter limit | 0-7 | analytic limits per regime |
| outputs | full transport matrix, particle, heat and momentum fluxes per species, non-ambipolar fluxes, NTV torque per species and n, offset rotation, bootstrap current, conductivity, poloidal rotation coefficient, parallel flows and local parallel current, impurity flux, current-drive efficiency (generalised Spitzer function by the adjoint) | 0, parity | NEO-2, NEO, SFINCS, Sauter fits |
| numerics | derivatives of every output with respect to profiles, geometry and perturbation (tangent and adjoint); certified bands (#274); hp adaptivity | 0, 1-2 | finite differences, interval bounds |

E_r-shear orbit squeezing needs no term of its own: it follows from Phi(psi) in the total energy
H once orbit width is retained.

Outside this scope: computing turbulence itself (it enters as a spectrum, below), time dependence beyond the
steady state, open field lines and the scrape-off layer, and stochastic magnetic fields.

## Design: one operator, composed terms

Every effect above is a term of the one equation, not a mode of the code. A run is described by
the list of terms it contains (streaming, drifts, each drive, each collision part, each closure),
and the operator is assembled from that list by one kernel. Consequences:

- no Boolean flag per effect inside numerical routines; a term that is absent is not assembled;
- one solver path for all combinations; a combination needs no new code, only its check;
- each term owns its coefficients, its source, its output moments and its limit test (removing
  it reproduces the run without it), so tests compose with the terms;
- geometry, collisions, velocity mesh and solver are separate interfaces; a new effect touches
  one of them, never copies a solver;
- presets are term lists, never flags; the `neo2-ql` preset reproduces the validated thin-orbit
  quasilinear model, and a term requested outside
  its validity is refused with the indicator that failed.

Sonic rotation, Phi_1 and the misaligned potential share one building block: a background whose
density and potential vary on the surface (centrifugal force, perturbation, or both), closed by the
same quasineutrality condition.

## Complete model and reduced presets

One input runs the most complete model the case admits: all terms whose validity indicators pass,
with the full operator. Named presets are fixed term lists that reproduce the models of established
codes, so a reduced run is the complete model with terms removed, on the same geometry, profiles,
mesh and solver:

| Preset | Terms (removed from the complete model) |
|---|---|
| `complete` | all admitted terms |
| `neo2-ql` | thin orbit, slow rotation, Phi_1 = 0, delta B to first order with the |B|-only drive, local Maxwellian |
| `neo2-par` | 3D |B|, thin orbit, no perturbation, Phi_1 = 0 |
| `neo-rt` | thin orbit, bounce-averaged resonant transport, model collision diffusion, nonlinear attenuation |
| `potato` | finite orbit width, canonical Maxwellian, collisionless resonant limit |
| `sfincs` | 3D, radially local, full operator, optional Phi_1 |
| `dkes`, `monkes` | 3D, monoenergetic, Lorentz operator (DKES with momentum correction) |
| `neo` | axisymmetric, full operator, sonic rotation with poloidal asymmetry |
| `perfect` | axisymmetric, global finite orbit width, full operator |
| `kilca` | cylinder, linear response with FLR, self-consistent fields in the resonant layer |

The standard product of a device case is the complete result together with the reduced presets
that were used before, so that the model error of each established approach is shown directly as
the difference to the complete model. A preset reproduces its code only after the code benchmark
of its stage has passed; until then the difference is labelled as numerical plus model error.

## Implementation order

Each item lands with its checks before the next starts; derivatives and cost records accompany
every item.

1. Close step 0: energy convergence of the AUG ion full operator at resonances.
2. Step 1: hp prototype and its go/no-go.
3. Rename to MILCA.
4. Composed-term architecture: the current kernel expressed as a term list with the `neo2-ql`
   preset, results unchanged within solver tolerance. All later terms are added to this.
5. Local tokamak parity on the existing kernel: axisymmetric transport matrix (bootstrap,
   conductivity, Ware pinch), heat flux K_2j, local parallel current, poloidal rotation coefficient,
   E_r from force balance.
6. Full perturbation drive in the local kernel: field-line bending and b.delta B, delta Phi_n, the
   misaligned potential with its zero-current test, rotating perturbations.
7. Step 2 (hp with all energies), if step 1 passed.
8. Validity indicators (step 4), needed to decide where the complete model must run.
9. Finite orbit width on the canonical Maxwellian (step 3).
10. Perturbation in the orbits, isolated resonances, thin orbit and with finite orbit width
    (steps 5 and 6).
11. AUG 30835 and TC24 curves of these effects separately and combined (step 8, first part).
12. Surface-varying background: Phi_1 from quasineutrality, sonic rotation, strong-E guiding-centre
    equations.
13. Stellarators (K1): 3D |B| with ripple wells, tangential drift, ambipolar roots, rabe.
14. Edge, scrape-off layer and divertor: the continuum model hierarchy below, with divertor heat
    loads as the first product.

Deferred until a concrete case needs them, and then only after the items above: a second
action-angle implementation (NEO-RT and POTATO remain the action-angle references), overlapping
resonances in the global solve (Monte Carlo covers them), second order in the gyroradius, edge
sources and orbit loss, fast-ion backgrounds (ASCOT, SIMPLE), the relativistic operator and
current drive (NEO-2-PAR), FLR and self-consistent fields (KAMEL, MEPHIT and
[response-unified.md](response-unified.md)), kinetic delta W feedback (MARS-K), magnetic islands,
the turbulence-spectrum operator, and presets other than
`neo2-ql`, `neo-rt`, `potato` and `sfincs`. hp beyond step 1 continues only if step 1 shows a
clear gain.

## Turbulence, time dependence and the edge

The turbulence and time-dependence parts are deferred; the edge is item 14 of the order above.

Turbulence enters as a given saturated spectrum (from a gyrokinetic code or a quasilinear model
such as TGLF or QuaLiKiz, with amplitudes, frequencies and cross-phases), not as a computed
turbulent state. Its term is a phase-space resolved quasilinear operator: diffusion and
convection in psi and energy with the resonance of each mode with parallel streaming, drifts and
rotation, gyroaveraged with J0(k_perp rho) and orbit-averaged. Compared with a scalar
diffusivity, this gives the energy, charge and mass dependence of turbulent impurity and fast-ion
transport (pinches, screening, FLR reduction for heavy ions) from the same spectrum. Beyond the
quasilinear limit, resonance broadening (decorrelation by the turbulence itself, the analogue of
Theta(D) for the superbanana) is the first correction; for trace species, test particles in
synthetic or imported turbulent fields (GORILLA) are exact for the given fields and serve as the
reference where the Kubo number is not small.

Time dependence: a perturbation that is periodic in time is a sum of frequencies, each a steady
solve of the same linear operator. Transients (mode slowing, ELMs, ramps) use implicit time
stepping of the same operator, (1/dt + L) f = rhs, with the same factorisations.

Edge, scrape-off layer and wall belong to the same framework (the long-range goal of the
[roadmap](../ROADMAP.md)): flux coordinates fail at the X point and on open field lines, and the
physics there needs full-f, a nonlinear collision operator, sheath and wall boundaries, neutrals,
atomic reactions and turbulence, but none of this needs a second code. The rules that keep it one
framework:

- An atlas of blocks, each with its own coordinate map and basis: Boozer or canonical coordinates
  with Fourier and Legendre bases in the core; locally field-aligned mapped blocks for the
  near-separatrix closed region, the scrape-off layer and the private flux region, matched at the
  separatrix with the X point as a block corner (as COGENT does); unstructured blocks only next to
  the wall. Terms are written against a metric and coordinate interface, so each term exists once
  and every block evaluates it; blocks exchange conservative fluxes at their faces.
- Full-f is the base operator, and the delta-f operators used now are its linearisation about the
  Maxwellian or canonical Maxwellian, obtained by differentiation (JVP), not written twice. The
  same holds for the nonlinear Landau operator and the linearised Fokker-Planck operator.
- Sheath and wall are boundary terms on block faces, with wall data (recycling, sputtering yields)
  from external tables or codes; neutrals are species with their own fluid or kinetic equation;
  atomic and impurity charge-state processes are reaction terms between species.
- Turbulence is either the spectrum term above or a resolved gyrokinetic region where an
  adjoint indicator asks for it.
- Model fidelity is chosen locally (delta-f or full-f, linearised or nonlinear collisions,
  fluid or kinetic neutrals) by the indicators, as the roadmap's adaptive model hierarchy states.
- Tests carry model-class labels and run by what a change affects, never as the full suite.

All edge models are continuum models in this framework, in a hierarchy whose levels are term
lists on the same block atlas and check each other:

| Level | Model | Benchmarks |
|---|---|---|
| E0 | field-line fluid scrape-off layer with sheath boundary | two-point model, Stangeby sheath analytics |
| E1 | 2D axisymmetric multi-fluid Braginskii with drifts and fluid neutrals | SOLPS-ITER, UEDGE |
| E2 | 3D fluid on field-aligned or FCI blocks: RMP edges, stellarator island divertors | EMC3-EIRENE, W7-X |
| E3 | deterministic kinetic neutrals (discrete-velocity Boltzmann/BGK), hybrid with fluid neutrals by Knudsen number | EIRENE |
| E4 | full-f drift kinetics in the edge blocks where fluid closures fail: nonlocal parallel heat flux, kinetic sheath | COGENT, analytic nonlocal heat-flux limits |
| E5 | turbulence as spectrum term or resolved region | GENE-X, experiment |

Divertor heat loads (convective, conductive, radiative and neutral), target profiles, impurity
sources and detachment are the products, compared with infrared heat-load measurements on AUG
and W7-X. Reactions are bookkept in species and reaction tables with a conservation ledger:
atomic and molecular rates (ADAS, AMJUEL/HYDHEL), sputtering yields and fusion reaction cross
sections. GORILLA Monte Carlo serves as one independent benchmark, not as the edge model.

## Working rule

Exploration uses focused tests and resolution ladders. Pre-registered acceptance rules, error
budgets (numerical, model, input) and certified bands (#274) apply when a result is claimed,
published or used for a decision, not to every development step. A step that shows no clear
advantage in accuracy, range or cost is reported as a negative result before the next step starts.
