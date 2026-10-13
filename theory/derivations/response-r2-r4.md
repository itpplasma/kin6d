# T8: linearised rotating MHD response per n, (R2)-(R4)

- Scope: the per-n momentum equation (R2), Ohm's law (R3) and pressure
  closure (R4) of [response-unified.md](../../docs/response-unified.md):
  linearisation about a rigidly rotating axisymmetric state, Lagrangian and
  Eulerian forms, the ideal limit and its gauges, the static ideal closure of
  `kin6d_mhd_ideal_n`, the resistive and drift-MHD layer equations with their
  known limits, and the pressure-tensor interface to the DK solver.
- Status keys follow [theory-map.md](../../docs/theory-map.md). Every
  identity marked "checked" is evaluated by `test/check_response_r2_r4.py`
  on concrete fields (residuals 1e-16 to 1e-15); [K] marks input stated from
  background knowledge and not derived here.

## Assumptions

1. Per-n perturbation X(R, Z) exp(i n phi - i omega t), complex omega,
   cylindrical (R, phi, Z) right-handed; e_Z x e_R = e_phi.
2. Equilibrium B0 = sigma grad psi x grad phi + F(psi) grad phi, sigma = +-1
   (the code's `s` in `mhd_ideal_n.f90:26-31` is sigma = -s).
3. Flow V0 = Omega(psi) R^2 grad phi (rigid toroidal rotation per surface,
   no poloidal flow), steady, div(rho0 V0) = 0.
4. Generalised Ohm's law without electron inertia,
   E + V x B = eta J + s_2f (J x B - grad p_e)/(n e), V the ion (centre of
   mass) velocity, s_2f in {0, 1}; quasineutral, n0 constant across a layer.
5. Momentum rho (d_t + V.grad) V = J x B - div P + div(rho nu grad V);
   two-fluid runs need in addition the ion gyroviscous stress [K] (see R6).
6. Single-fluid ideal equilibrium: E0 + V0 x B0 = 0; two-fluid equilibrium:
   E0 + V0 x B0 = s_2f grad p_i0/(n e), E0 = -Phi0'(psi) grad psi.
7. Subsonic rotation is NOT needed for R2-R4 below; it is needed only where
   p0 = p0(psi) is used (R5) and is then O(M^2).
8. Isotropic Maxwellian background for R7 (p_par0 = p_perp0 = p0).

## R1. Rotating equilibrium

V0 x B0 = sigma Omega grad psi, so assumption 6 gives Phi0' = sigma Omega
(single fluid) and, with s_2f = 1, the three surface frequencies

    sigma Phi0' = Omega_E,  Omega_i = Omega_E + sigma p_i0'/(n e),  Omega_e = Omega_E - sigma p_e0'/(n e).   (T8.1)

Checked: E0 + V0 x B0 = 0 for sigma = +-1.

## R2. Lagrangian and Eulerian velocity

V0 is the vector field Omega(psi) d/dphi. For any per-n vector field xi the
Lie bracket gives, with physical components,

    v1   = d_t xi + V0.grad xi - xi.grad V0 = -i w xi - Omega' (xi.grad psi) R^2 grad phi,   (T8.2)
    Dv   = d_t xi + V0.grad xi             = (-i w + Omega e_Z x) xi,                       (T8.3)

w = omega - n Omega (Doppler frequency). v1 is the Eulerian velocity
perturbation, Dv = v1 + xi.grad V0 the Lagrangian one; Q_L = Q1 + xi.grad Q0
for every scalar. Checked (both).

## R3. Momentum (R2)

Frieman-Rotenberg form [Frieman and Rotenberg 1960], with V0.grad = Omega d_phi
acting on physical components as Omega (i n + e_Z x):

    rho0 (-i w + Omega e_Z x)^2 xi = J1 x B0 + J0 x B1 - div P1 + div(rho0 nu grad v1)
                                     + div(rho0 xi) a_c + rho0 xi.grad a_c,  a_c = -Omega^2 R e_R.   (T8.4)

The inertia operator is checked (rho0 (d_t + V0.grad)^2 xi); the force-side
rotation terms are the Frieman-Rotenberg operator for this flow (stated from
the paper, not re-derived). Expanded, the left side is -rho w^2 xi
(Doppler-shifted inertia), -2 i w rho Omega e_Z x xi (Coriolis) and
-rho Omega^2 (xi_R e_R + xi_phi e_phi); the e_phi part cancels against
rho xi.grad a_c. Consequences for the design form:

- For the static coil drive w = -n Omega the Coriolis and centrifugal terms
  are O(rho Omega^2), the same order as rho w^2 for n = 1-3. Keeping rho w^2
  without them is not a consistent ordering.
- Viscosity acts on v1 of (T8.2): div(rho nu grad v1) equals
  -i w div(rho nu grad xi) only for uniform w and Omega' = 0.

## R4. Ohm's law (R3), co-moving form and gauges

For arbitrary A1, Phi1, xi, with E1 = i omega A1 - grad Phi1, B1 = curl A1
and v1 from (T8.2), the exact identity (checked)

    E1 + v1 x B0 + V0 x B1 = i w G - Omega' G_phi grad psi - grad phi_c,
    G = A1 - xi x B0,  G_phi = R^2 grad phi . G,  phi_c = Phi1 - V0 . A1                      (T8.5)

follows from L_V0 A1 <-> B1 x V0 + grad(V0 . A1) and L_V0 B0 = 0. The
right side is gauge invariant (A1 -> A1 + grad chi, Phi1 -> Phi1 + i omega
chi). The linearised Ohm's law of assumption 4 is therefore

    i w G - Omega' G_phi grad psi - grad phi_c
      = eta J1 + s_2f [J1 x B0 + J0 x B1 - grad p_e1 - (n1/n0) grad p_i0]/(n0 e).             (T8.6)

Ideal limit (eta = s_2f = 0, w != 0 or not): every solution has the form
G = grad chi, phi_c = i w chi; in particular

- gauge A1 = xi x B0 (the DK gauge, drift-kinetics.md Eq. 2):
  Phi1 = -sigma Omega xi.grad psi. Checked. Hence the Lagrangian potential
  delta Phi_L = Phi1 + xi.grad Phi0 = 0 exactly for ideal single-fluid
  rotation: the DK drive term Z e delta Phi_L vanishes in this gauge and is
  driven only by G != 0 (resistive, two-fluid) or by Phi0' != sigma Omega.
- code gauge R A_phi = 0 (`mhd_ampere_n.f90:12-13`, A -> A - grad(a/(i n))):
  A1 = xi x B0 - (i sigma/n) grad(xi.grad psi), Phi1 = sigma (w/n) xi.grad psi,
  V0 . A1 = 0. Checked.

In the code gauge the unknowns are A_pol (2), xi (3) and phi_c (1) against
(T8.4) and (T8.6), three components each; phi_c is a genuine unknown, so
(R3) without it is overdetermined. phi_c carries the electrostatic response;
with w = 0 in a non-rotating plasma (T8.6) gives eta J1 = -grad phi_c (no
shielding), as expected.

Parallel projection. b0 . (V0 x B1) = -sigma Omega_i B1.grad psi/B and
b0 . (J0 x B1) = -p0' B1.grad psi/B; with (T8.1),

    b0 . E1 = eta J1_par + sigma Omega_e B1.grad psi/B - s_2f b0.grad p_e1/(n e)              (T8.7)

(sigma Omega_e = Phi0' - s_2f p_e0'/(n e)). The flux is frozen into the
electron fluid; the J0 x B1 term of (T8.6) is what turns Omega_i into
Omega_e. Status: hand derivation; enters the checked layer equations below.

## R5. Static ideal outer region and `kin6d_mhd_ideal_n`

Drop rho w^2, rotation terms (O(M^2), O(w^2 tau_A^2)), nu, eta, s_2f; fluid
P1 = p1 I. (T8.4) is grad p1 = J1 x B0 + J0 x B1 with J0 x B0 = grad p0. Then

1. b0 . (T8.4): B0.grad p1 = B0.(J0 x B1) = -p0' B1.grad psi, independent of
   the p1 closure (`mhd_ideal_n.f90:9-10`, RHS at 330, MDE `mde_solve` 268-296).
2. J1_perp = B0 x (grad p1 - J0 x B1)/B0^2 (header line 11, code 364).
3. div J1 = 0 (from curl of (R1)): B0.grad h = -div J1_perp, J1_par = h B0,
   and the code form of -div J1_perp (lines 12-16, code 373-381) is an exact
   identity for div B1 = div J0 = 0 (uses (a x b).curl c =
   b.(a.grad)c - a.(b.grad)c). Checked.
4. Ideal B1 = curl(xi x B0) gives B1.grad psi = B0.grad(xi.grad psi) (used at
   line 399). Checked. With 1, B0.grad(p1 + xi.grad p0) = 0, so off rational
   surfaces (n != 0) p1 = -xi.grad p0 and Gamma p0 div xi = 0: the Lagrangian
   pressure perturbation vanishes, Gamma does not enter the static closure,
   and xi need not be carried (as the design states).

Rotation enters the outer region only through the dropped inertia; the
frozen-in law (T8.5) is w-independent, so the same closure serves the
rotating lanes outside the strips. Status: derived; code identity checked.

## R6. Layer equations: resistive and drift-MHD

Slab near q = m/n: x = r - r_s, k_y = m/r_s, B0 = B_z e_z + B_y'(x) x e_y,
k_par = k_y B_y' x/B; perturbed flux psi (B1x = i k_y psi), E x B
displacement xi = xi_x of the E x B velocity, constant density, reduced
(incompressible, strong guide field) dynamics, p_e1 = -xi.grad p_e0 (no
parallel electron compression: SLAYER ds -> 0), mu0 = 1. Frequencies in the
frame of the E x B flow: gamma_E = -i(omega - omega_E), gamma_e = -i(omega -
omega_E - omega_*e), gamma_i = -i(omega - omega_E - omega_*i).

(T8.7) with the advection -i k_par Phi1 written through xi gives the Ohm
equation; the vorticity equation with ion polarisation in the ion
diamagnetic frame (gyroviscous cancellation) is [K]:

    gamma_e (psi - B_y' x xi) = eta psi'',
    rho gamma_E gamma_i xi'' - rho nu gamma_i xi'''' = -k_y^2 B_y' x psi''.                 (T8.8)

Single fluid: gamma_E = gamma_e = gamma_i = gamma = -i w.

Constant psi (psi = psi_s across the layer), x = delta X, xi = psi_s Y/(B_y'
delta): (T8.8) becomes parameter free (checked),

    RI:  Y'' - X^2 Y = -X,       delta^4 = rho eta gamma_E gamma_i/(k_y^2 B_y'^2 gamma_e),
    VR:  Y'''' + X^2 Y = X,      delta^6 = rho nu eta gamma_i/(k_y^2 B_y'^2 gamma_e),       (T8.9)

with Y -> 1/X (ideal outer solution) and

    Delta_layer = (gamma_e delta/eta) Dhat,  Dhat_RI = 2 pi Gamma(3/4)/Gamma(1/4) = 2.1236,
                                             Dhat_VR = 6^(2/3) pi Gamma(5/6)/Gamma(1/6) = 2.1036.   (T8.10)

Both constants are reproduced by a BVP solve to 2.4e-7 and 1.8e-8 (checked;
also the integral representation Y = (X/2) int_0^1 (1-t^2)^(-1/4)
exp(-X^2 t/2) dt). Limits (all checked):

- Delta_RI ~ gamma_e^(3/4) (gamma_E gamma_i)^(1/4) eta^(-3/4) and
  Delta_VR ~ gamma_e^(5/6) gamma_i^(1/6) (nu rho)^(1/6) eta^(-5/6): the RIi and
  VRi regimes of Cole and Fitzpatrick 2006, Eqs. (41) and (47) and Table I
  (constants 2.124, 2.104), with gamma_E = iQ, gamma_e = i(Q - Q_e),
  gamma_i = i(Q - Q_i) in their normalisation.
- Single fluid: Delta_RI ~ gamma^(5/4) S^(3/4). FKR tearing (Delta_layer =
  Delta'): gamma ~ S^(-3/5) Delta'^(4/5), width ~ S^(-2/5) [Furth et al. 1963].
- Constant-psi validity delta Delta_layer << 1 fails at gamma tau ~ S^(-1/3),
  where the width is S^(-1/3) independent of w: the non-constant-psi
  (resistive-Alfven) scale [Coppi et al. 1966]. VR width
  P_m^(1/6) S^(-1/3), also w independent.
- Penetration: Delta_layer -> 0 as gamma_e -> 0, i.e. shielding is lost when
  the electron fluid is at rest in the frame of the drive (omega = omega_E +
  omega_*e); for the static drive the relevant rotation is Omega_e of (T8.1).
- SLAYER (GPEC `slayer/delta.f`, `w_der`, PeOhmOnly): its layer ODE
  (C3 Y')' = G Y in the Fourier variable p reduces, for pr, c_beta, ds -> 0,
  to C3 = p^2/(p^2 + gamma_e), G = (gamma_E gamma_i/gamma_e) p^2, which is the
  Fourier transform (x -> i d/dp, d/dx -> i p) of (T8.8) with nu = 0.
  Checked symbolically against the source expression.

Dimensional widths. With k_par' = k_y B_y'/B = n s/(R r_s) and
eta_m = eta/mu0, delta_RI^4 = gamma eta_m/(k_par' v_A)^2 (single fluid), so
(delta/a)^4 = (gamma tau_A) S^(-1) (n s)^(-2) (R r_s/a^2)^2 for tau_A = a/v_A,
S = a^2/(eta_m tau_A). The design table's widths hold up to the geometric
factor (R r_s/a^2)^(1/2) ~ 2-3 for ITER, and the constant-psi RI row needs
w tau_A << (n s)^(2/3) S^(-1/3) (R r_s/a^2)^(-2/3) ~ 1e-3 at S ~ 1e8-1e9:
the quoted w tau_A ~ 1e-3-1e-2 is at or beyond the boundary, i.e. ITER layers
are non-constant-psi or two-fluid, not FKR constant-psi. Status: scalings
derived and checked; the geometric factor is a hand step.

## R7. Pressure closure (R4)

Exact for any divergence-free B (checked):

    delta B_L / B = (b.delta B + xi.grad B)/B = -div xi + b.grad xi_par - kappa.xi,          (T8.11)

kappa = b.grad b, xi_par = b.xi. Closures (Lagrangian, delta rho_L/rho =
-div xi):

- fluid: delta p_L = -Gamma p0 div xi, P1 = p1 I, p1 = delta p_L - xi.grad p0;
- CGL [Chew et al. 1956]: delta p_perp,L = p0(-div xi + delta B_L/B),
  delta p_par,L = p0(-3 div xi - 2 delta B_L/B), i.e.
  delta p_perp,L = -p0 (2 div xi - b.grad xi_par + kappa.xi),
  delta p_par,L = -p0 (div xi + 2 b.grad xi_par - 2 kappa.xi);
- kinetic: P1 = sum_a int d^3v m v v f_a1, gyrotropic to the DK order and
  with assumption 8
  P1 = delta p_perp (I - b b) + delta p_par b b (the (p_par0 - p_perp0)(b b1 + b1 b)
  term vanishes), delta p_perp = sum int (m v_perp^2/2) f1,
  delta p_par = sum int m v_par^2 f1 (Eulerian f1 at fixed x).

Energy form. For gyrotropic P1 (checked, with xi* in place of xi),

    grad xi* : P1 = -delta p_perp delta B_L*/B + delta p_par (b.grad xi_par* - kappa.xi*),

so, with delta H = mu delta B_L + m v_par^2 kappa.xi of drift-kinetics.md (2)
(delta Phi_L = 0 by R4) and P1 -> 0 at the plasma boundary,

    delta W_k = (1/2) int xi*.div P1 d^3x
              = (1/2) int d^3x [ sum_a int d^3v delta H* f_a1 - delta p_par b.grad xi_par* ].   (T8.12)

The second term is the parallel-displacement part; whether it cancels
against the parallel kinetic force balance (Kruskal-Oberman-type) for a
DK solution is not shown here (open). With the stated torque form of
drift-kinetics.md, T_phi,n = -n Im(first term) per unit volume [K form,
sign not verified].

Quantities the DK solver must return, per n, per flux surface psi_j, on the
surface theta grid, summed over species (per species for diagnostics):

1. delta p_par,n(theta) and delta p_perp,n(theta), complex, Eulerian, of the
   part of f1 that (R2) does not already carry: the full f1 when the kinetic
   species replace the fluid pressure, only the non-adiabatic part when the
   fluid p1 of (T8.4) is kept for that species (MARS-K split [Liu et al.
   2008], to be confirmed against its source; no double counting either way).
2. The surface-integrated bilinear w_k,n(psi_j) = (1/2) <sum int d^3v
   delta H* f1> (Re: energy, Im: torque) and, separately, the
   b.grad xi_par* term of (T8.12), so that the cancellation is measured.
3. A flag naming the representation (Eulerian or Lagrangian, conversion
   delta p_E = delta p_L - xi.grad p0 for both components).

Inputs from the response lane are the L7 export: xi (incl. xi_par),
delta B_L, kappa.xi, w(psi) and delta Phi_L (zero in ideal single-fluid
regions). Perturbative (MARS-K-type): one DK pass on the ideal xi gives
delta W_k; self-consistent: -div P1 from item 1 enters (T8.4) in the
fixed point. Status: interface derived; DK moments not yet computed.

## Verification

`python3 test/check_response_r2_r4.py` (SymPy, NumPy, SciPy, mpmath; about
a minute). Gates: identities below 1e-12 relative (observed 5e-17 to 1e-15)
for sigma = +-1; layer reductions exactly parameter free; Dhat to 1e-6;
FKR and constant-psi exponents exact rationals; SLAYER limit exact.

## Disagreements with response-unified.md (corrected there)

1. (R3) as written, i w (A - xi x B0) = ..., is the co-moving law only in a
   gauge with phi_c = 0 and without the flow-shear term; in the code gauge
   R A_phi = 0 it is overdetermined. Correct form (T8.6).
2. The ideal limit in the code gauge is A = xi x B0 - (i sigma/n)
   grad(xi.grad psi), gauge-equivalent to xi x B0 (same B1), with the scalar
   potential of R4.
3. The two-fluid term needs the full linearisation (T8.6) incl. J0 x B1;
   without it the frozen-in frame is Omega_i, not Omega_e. The s_2f row also
   needs the ion gyroviscous stress in (R2) [K].
4. (R2): Coriolis and centrifugal terms are the same order as rho w^2 for the
   static drive; viscosity acts on v1.
5. Width table: S^(-1/3) is the w-independent constant-psi boundary, and the
   ITER w tau_A range is at that boundary (constant-psi FKR not applicable);
   geometric factor (R r_s/a^2)^(1/2) omitted in the table.

## Open questions

- needs review: gyroviscous drift-MHD vorticity of (T8.8) [K]; derive from
  Braginskii gyroviscosity in toroidal form before R-B2's s_2f gate.
- needs review: the force-side Frieman-Rotenberg terms of (T8.4) are quoted,
  not derived; only the inertia operator is checked.
- needs review: xi_par term of (T8.12) and the MARS-K fluid/kinetic split.
- needs review: sign of T_phi,n versus Im delta W_k (follows the DK1 sign fix).
- needs review: ion-sound / semi-collisional terms (SLAYER ds, c_beta) and
  parallel electron compression are not derived.

## References (checked against Crossref 2026-10-13; Cole-Fitzpatrick also in Zotero)

- Frieman, Rotenberg, Rev. Mod. Phys. 32, 898 (1960), https://doi.org/10.1103/RevModPhys.32.898
- Chew, Goldberger, Low, Proc. R. Soc. A 236, 112 (1956), https://doi.org/10.1098/rspa.1956.0116
- Kruskal, Oberman, Phys. Fluids 1, 275 (1958), https://doi.org/10.1063/1.1705885
- Furth, Killeen, Rosenbluth, Phys. Fluids 6, 459 (1963), https://doi.org/10.1063/1.1706761
- Coppi, Greene, Johnson, Nucl. Fusion 6, 101 (1966), https://doi.org/10.1088/0029-5515/6/2/003
- Cole, Fitzpatrick, Phys. Plasmas 13, 032503 (2006), https://doi.org/10.1063/1.2178167
  (Eqs. 17, 41, 47 and Table I read from the Zotero PDF)
- Liu, Chu, Chapman et al., Phys. Plasmas 15, 112503 (2008), https://doi.org/10.1063/1.3008045
  (title checked; the fluid/kinetic split is [K])

Chris&AI
