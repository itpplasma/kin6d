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
<int d^3v p_phi_dot f1>, with <.> the flux-surface average and p_phi = m R v_phi - Z e psi the
canonical toroidal momentum; for thin orbits it reduces to Z e Gamma_psi with the radial particle
flux Gamma_psi. Expanding about F_M(psi, H) assumes the orbit width is small against the gradient
lengths; where it is not (pedestal), F_M is replaced by a function of the invariants
(p_phi, H, mu) or by a full-f formulation.

The modes are subsets of its terms:

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
| 1 | hp prototype: axisymmetric, monoenergetic, (theta, xi) mesh fitted to lambda | spectral convergence down to nu* ~ 1e-6; banana, plateau and Pfirsch-Schlueter limits of D11, D31, D33; bootstrap coefficient against rabe; faster or leaner than the Legendre solver at the same accuracy |
| 2 | hp for the thin-orbit quasilinear model with all energies (after a passed step 1) | agrees with the spectral full-operator kernel; 1/nu, sqrt(nu)-nu, superbanana-plateau and rapid-rotation limits; converges the AUG ion resonances where uniform speed elements did not |
| 3 | `orbit_width`: radial elements in psi, GMRES with surface factors, radial buffer boundaries, axis region excluded | orbit width -> 0 reproduces step 0 with the expected scaling; manufactured solution; particle and p_phi conservation; GMRES iterations bounded under psi refinement; POTATO in its thin-orbit limit and at finite orbit width; PERFECT pedestal case |
| 4 | validity indicators per surface, species and energy: D; island width over resonance spacing; d(resonant frequency)/d p_phi, which sets the island width; orbit width over gradient length | D against NEO-RT and canonical frequencies against POTATO on shared cases |
| 5 | `nonlinear_orbits` for one isolated resonance, thin orbit | delta B -> 0 reproduces step 0 with torque proportional to delta B^2; attenuation against NEO-RT Theta(D) over D; Theta -> 1 in the quasilinear limit |
| 6 | `nonlinear_orbits` with `orbit_width` for isolated resonances | each switch off reproduces steps 3 and 5; resonance shifts against POTATO canonical frequencies |
| 7 | coupled harmonics, overlapping resonances | isolated limits reproduce steps 5 and 6; Monte Carlo with finite orbits and the perturbation in the orbits; cost below Monte Carlo at the same accuracy. Overlap is expected in the pedestal (sheared E_r, large delta B), where the expansion about F_M is replaced as stated above |
| 8 | device curves | AUG 30835 and ITER TC24 with each switch alone and combined, against NEO-2, NEO-RT, POTATO and Monte Carlo where they apply, with the indicators of step 4 |

## Working rule

Exploration uses focused tests and resolution ladders. Pre-registered acceptance rules, error
budgets (numerical, model, input) and certified bands (#274) apply when a result is claimed,
published or used for a decision, not to every development step. A step that shows no clear
advantage in accuracy, range or cost is reported as a negative result before the next step starts.
