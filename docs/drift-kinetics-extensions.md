# Drift-kinetic extensions: adaptive velocity space, finite orbits, nonlinear resonances

Status: plan (2026-11-08). It extends [drift-kinetics.md](drift-kinetics.md), whose K1/K2 section
defines the local and global steady problems. Each stage below is a switch on the same steady
Eulerian operator, so the effect of every addition can be measured on its own or combined.

## Stages and switches

| Stage | Physics added | Switch | Contains as limit |
|---|---|---|---|
| hp velocity mesh | orbit-class-fitted spectral elements in pitch, adaptive by dual-weighted residuals | discretisation | spectral Legendre solver where it converges |
| K1 | local 3D transport, tangential drift, ambipolar E_r | geometry | axisymmetric local solver |
| K2 | radial drift on the unknown: finite orbit width in the unperturbed field, global in psi | `orbit_width` | K1 as orbit width -> 0 |
| perturbed characteristics | the perturbation in the guiding-centre motion: resonance islands, superbanana regime, nonlinear in delta B | `nonlinear_orbits` | quasilinear torque as delta B -> 0 |

`orbit_width` and `nonlinear_orbits` are independent: either alone, both, or neither (the
current quasilinear, thin-orbit model). The solution stays Eulerian and steady; where orbits are
complicated the operator coefficients come from integration along the actual characteristics
(as NEO-2 along field lines and NEO-RT in action-angle variables), not from time stepping.

## Verification and benchmarks

| Mode | Analytic limit | Code benchmark |
|---|---|---|
| hp mesh | 1/nu, nu-sqrt(nu), superbanana plateau, rapid rotation, Pfirsch-Schlueter, banana; manufactured solutions | spectral solver, DKES/MONKES/SFINCS, bootstrap against rabe down to nu* < 1e-6 |
| K2 | thin-orbit limit recovering K1 | POTATO thin-orbit limit and finite orbit width; PERFECT pedestal cases |
| perturbed characteristics | quasilinear limit (attenuation 1); isolated-resonance attenuation Theta(D) (Albert, dissertation 2020, Sec. 4.9) | NEO-RT nonlinear attenuation (thin orbit); GORILLA/SIMPLE Monte Carlo |
| combined | each switch off reproduces the stage below | Monte Carlo with finite orbits and the perturbation in the orbits |

Device curves: AUG #30835 and ITER TC24 with each switch separately and combined, against
NEO-2, NEO-RT, POTATO and Monte Carlo where they apply. Every result carries a validity
indicator per surface: superbanana onset (resonant diffusivity D of NEO-RT), orbit width over
gradient length, distance to rational surfaces, Mach number.

## Working rule

Exploration uses focused tests and resolution ladders. Pre-registered acceptance rules, error
budgets (numerical, model, input) and certified bands (#274) apply when a result is claimed,
published or used for a decision, not to every development step. Each benchmark records wall
time, core-hours and peak memory next to the accuracy, against the reference at matched accuracy.
A stage that shows no clear advantage in accuracy, range or cost is reported as a negative
result before the next stage starts.
