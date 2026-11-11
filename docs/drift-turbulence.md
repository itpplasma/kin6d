# Drift-wave turbulence: spectra and transport coefficients

Status: design specification with the fluid part of stage L, stage HW and
the output contract implemented; the nonlinear benchmark follows.
Each further stage becomes active with the pull request that brings its tests.

## Consumer

Models of the plasma response to a resonant magnetic perturbation resolve a
thin layer around a rational surface. Inside that layer the perturbed electron
motion along the field competes with cross-field transport that turbulence
supplies. Present layer models take this transport as a constant diffusivity
or viscosity chosen by hand. The consumer of this branch is such a layer model:
it needs, at the plasma parameters of one rational surface,

- the fluctuation spectrum \(S(k_x,k_y,\omega)\) of the electrostatic
  potential (later also of the parallel vector potential), with the cross
  spectra that fix the fluxes;
- correlation times and lengths, since a layer narrower than a correlation
  length, or a response faster than a correlation time, does not see a
  diffusion;
- particle, heat and momentum fluxes with the diffusivities derived from
  them, each with a statistical standard error;
- the same quantities from a cheap quasilinear rule that is calibrated against
  the nonlinear runs where both apply.

Spectra are delivered in the frame without mean \(E\times B\) flow. The
consumer applies the Doppler shift \(\omega+k_yV_E\).

This is the "prescribed electromagnetic turbulence cross-spectrum" that
[TR0](../ROADMAP.md#tr0--macro-transport-core-and-cheap-integrated-modelling-cross-cutting-starts-before-r6)
asks for before a live provider, and it uses the closure contract of
[transport-multiscale.md](transport-multiscale.md).

## Stages

| Stage | Model | Geometry | Role |
|---|---|---|---|
| L | local linear roots of the fluid models below | one wave vector | linear oracle, quasilinear rule |
| HW | Hasegawa--Wakatani, original and modified | doubly periodic \((x,y)\) | verification of bracket, ledger, statistics and output; calibration case |
| DA | electromagnetic drift-Alfvén, isothermal electrons, gyrofluid ions | sheared slab, \(k_\parallel=k_yx/L_s\) | first model with electron inertia, \(A_\parallel\) and magnetic shear |
| K | kinetic local dispersion (gyrokinetic ions, drift-kinetic electrons) | one wave vector | quasilinear provider where no collisional fluid closure holds |

The order follows [S0](../ROADMAP.md#s0--stabilityresponse-consumer-cross-cutting-starts-with-r1r3):
the linear spectrum precedes the nonlinear run.

## Regime statement

The Hasegawa--Wakatani model describes resistive drift waves with a single
constant parallel wave number. Its control parameter is the adiabaticity
\(\alpha\) divided by the drive \(\kappa\); the model produces strong
turbulence for \(\alpha/\kappa\) of order one. At the top of an H-mode
pedestal the electrons are nearly collisionless, \(\alpha/\kappa\gg1\), the
electron skin depth is comparable to \(\rho_s\) and \(\nu_{*e}<1\). There the
parallel electron response is inductive and inertial, not resistive, and
trapped electrons are not detrapped by collisions. Stage HW is therefore
verification infrastructure and a calibration case, not a predictor of
pedestal-top transport. Stage DA contains the electromagnetic and inertial
response but no trapped particles, curvature, temperature-gradient dynamics or
microtearing. Any output of this branch names its model and carries the claim
label `candidate`.

## Stage HW: equations and ledger

Normalisation of Numata, Ball and Dewar, Phys. Plasmas 14, 102312 (2007):
lengths in \(\rho_s=c_s/\Omega_i\), time in \(1/\Omega_i\),
\(\phi=e\varphi/T_e\), \(n=n_1/n_0\), vorticity \(\zeta=\nabla^2\phi\),
\(\{a,b\}=\partial_xa\,\partial_yb-\partial_ya\,\partial_xb\),

\[
\partial_t\zeta+\{\phi,\zeta\}=\alpha(\phi-n)-D(-\nabla^2)^p\zeta,\qquad
\partial_tn+\{\phi,n\}=\alpha(\phi-n)-\kappa\,\partial_y\phi-D(-\nabla^2)^pn .
\]

Here \(\kappa=-\mathrm d\ln n_0/\mathrm dx=\rho_s/L_n\), so that a positive
flux \(\Gamma_n\) runs down the gradient. Numata et al. use \(p=2\). The
modified model applies the coupling \(\alpha(\phi-n)\) to the non-zonal parts
only (\(k_y\neq0\)). With \(n=\kappa n'\), \(\phi=\kappa\phi'\),
\(t=t'/\kappa\) the equations depend on \(\alpha/\kappa\) and
\(D/\kappa\) alone; the normalisation of Camargo, Biskamp and Scott,
Phys. Plasmas 2, 48 (1995), is the case \(\kappa=1\).

Energy and potential enstrophy per unit area, with \(\langle\cdot\rangle\) the
box average,

\[
E=\tfrac12\langle n^2+|\nabla\phi|^2\rangle,\qquad
W=\tfrac12\langle(n-\zeta)^2\rangle,
\]
\[
\dot E=\Gamma_n-D_\alpha-D_E,\qquad \dot W=\Gamma_n-D_W,\qquad
\Gamma_n=-\kappa\langle n\,\partial_y\phi\rangle,\qquad
D_\alpha=\alpha\langle(n-\phi)^2\rangle .
\]

\(D_E\) and \(D_W\) are the hyperviscous sinks. The discrete system keeps both
balances: the bracket is evaluated in a form for which
\(\langle f\{g,f\}\rangle=0\) holds after dealiasing, and every reported run
carries the residual of both balances.

Linear modes \(\exp[i(k_xx+k_yy-\omega t)]\), \(k_y\ne0\), obey

\[
\mu^2+b\mu+ib\omega_*=0,\qquad \mu=-i\omega+Dk^{2p},\qquad
b=\alpha\frac{1+k^2}{k^2},\qquad \omega_*=\frac{k_y\kappa}{1+k^2},
\]

which is Eq. (8) of Numata et al. for \(p=2\). For \(D=0\) and
\(\alpha\to\infty\), \(\omega\to\omega_*+i\omega_*^2k^2/[\alpha(1+k^2)]\).

## Stage L: fluid roots

`milca_drift_dispersion` returns, for one wave vector,

- `hw_linear_roots`: both Hasegawa--Wakatani roots in closed form, written so
  that the growth rate keeps its relative accuracy for \(\alpha\to\infty\),
  with the eigenvector \(n/\phi=1+k^2\mu/\alpha\); `hw_linear_roots_jvp` is
  the tangent in \(\alpha,\kappa,D,k_x,k_y\);
- `drift_alfven_roots`: the three roots of the homogeneous drift-Alfvén
  dispersion relation of stage DA for cold ions without the sound wave,

  \[
  \omega^2\left[\hat\beta(\omega-k_y)+\hat\mu k^2\omega+iCk^2\right]
  =k_\parallel^2\left[(1+k^2)\omega-k_y\right],
  \]

  in the units of that stage (lengths \(\rho_s\), time \(L_n/c_s\),
  parallel length \(L_\parallel\),
  \(\hat\beta=(\beta_e/2)(L_\parallel/L_n)^2\),
  \(\hat\mu=(m_e/m_i)(L_\parallel/L_n)^2\),
  \(C=0.51\,(\nu_eL_n/c_s)\hat\mu\)), with the residual of each call;
  `drift_alfven_root_tangent` applies the implicit rule at a simple root.

For \(\hat\beta,\hat\mu\to0\) the cubic reduces to the Hasegawa--Wakatani
relation with \(\alpha=k_\parallel^2/C\), \(\kappa=1\); for
\(k_\parallel\to\infty\) one root tends to the drift wave
\(k_y/(1+k^2)\); for \(k_y=0\), \(C=0\) the roots are \(0\) and the kinetic
or inertial Alfvén wave
\(\pm k_\parallel[(1+k^2)/(\hat\beta+\hat\mu k^2)]^{1/2}\). The formulas are
the leaves of `gen/gen_drift_fluid.f90`. `drift_dispersion_fluid_roots` checks
the roots against binary128 eigenvalues of the linearised equations, Eq. (8)
of Numata et al., these limits and central differences of the tangents; its
oracle with the density-gradient drive moved to the vorticity equation is
rejected.

## Periodic spectral grid

`milca_drift_spectral` holds a real field on a doubly periodic box as its
Fourier coefficients. Per direction the modes \(|m|\le(n-1)/3\) are retained,
which removes the aliases of quadratic products. `drift_bracket` evaluates

\[
\{\phi,f\}=\tfrac13\left[\phi_xf_y-\phi_yf_x
+\partial_x(\phi f_y-\phi_yf)+\partial_y(\phi_xf-\phi f_x)\right],
\]

the average of the advective and the two conservative forms, with the
pointwise products from the leaf `bracket_split_products`. Summation by parts
of the collocated derivative then gives
\(\langle\{\phi,f\}\rangle=\langle f\{\phi,f\}\rangle=\langle\phi\{\phi,f\}\rangle=0\)
term by term. The bracket is bilinear; under the box average its adjoint in
\(f\) is \(g\mapsto-\{\phi,g\}\) and its adjoint in \(\phi\) is
\(g\mapsto\{f,g\}\). `drift_spectral_bracket` checks a manufactured solution,
these identities to rounding and spectral convergence; a collocated advective
bracket on untruncated fields violates the quadratic invariant and is
rejected.

## Stage HW: implementation

`milca_hasegawa_wakatani` advances the coefficients of \(\zeta\) and \(n\).
One step is the Strang composition of the exact propagator
\(\exp[(M-\nu)\,\Delta t/2]\) of each mode, the bracket terms over
\(\Delta t\) with the third-order strong-stability-preserving Runge--Kutta
scheme, and the propagator again. Coupling, drive and hyperviscosity therefore
carry no time-step error for a single mode, and the scheme is second order for
the coupled system. After every step the mirror symmetry
\(\hat f(-k)=\hat f(k)^*\) of both fields is restored exactly: rounding leaves
a part of the opposite symmetry, which the linear instability would otherwise
amplify to order one after about \(37/\gamma_{\max}\). The propagator, the matrix \(M\) and the ledger terms of
one mode are leaves of `gen/gen_drift_fluid.f90`. `hw_ledger` returns \(E\),
\(W\), \(\Gamma_n\), \(D_\alpha\), \(D_E\), \(D_W\) and the zonal share of the
kinetic energy. `hw_rhs`, `hw_rhs_jvp` and `hw_rhs_vjp` give the right-hand
side with its tangent and adjoint under the box average.

`hasegawa_wakatani_oracles` checks: seeded eigenmodes of both roots against
the dispersion relation in the hydrodynamic, intermediate and adiabatic regime
for both models (error below \(10^{-11}\) after 40 steps); conservation of
\(E\) and \(W\) for \(\alpha=\kappa=D=0\) with third-order convergence of the
drift; closure of both balances in a driven dissipative nonlinear state with
second-order convergence; the zonal modes of the modified model; the tangent
against central differences and the adjoint identity. A run whose propagator
has the opposite sign of \(\kappa\) is rejected.

## Output contract

`milca_turbulence_statistics` turns equally spaced samples of the coefficients
\(\hat\phi(k_x,k_y)\), \(\hat n(k_x,k_y)\) of a stationary homogeneous run
into `turbulence_spectrum_t`; `milca_turbulence_output` writes it into the
result file of [interface-reproducibility.md](interface-reproducibility.md).
Every provider of this branch writes the same variables.

| Variable | Dimensions | Meaning |
|---|---|---|
| `kx`, `ky`, `omega` | axes | symmetric about zero; time dependence \(e^{-i\omega t}\) |
| `S_phi_phi` | (kx, ky, omega) | power spectral density; its sum times \(\Delta k_x\Delta k_y\Delta\omega\) is the variance of \(\phi\) |
| `S_n_phi_re`, `S_n_phi_im` | (kx, ky, omega) | cross spectral density \(\langle\hat n\hat\phi^*\rangle\); \(\Gamma_n=-\kappa\sum k_y\,\mathrm{Im}\,S_{n\phi}\,\Delta k_x\Delta k_y\Delta\omega\) |
| `lag_t`, `lag_x`, `lag_y`, `C_phi_t`, `C_phi_x`, `C_phi_y` | lags | cuts of the correlation function of \(\phi\) through the origin |
| `tau_c_euler` | scalar | first \(1/e\) crossing of the single-point autocorrelation |
| `tau_c_mode` | scalar | first \(1/e\) crossing of \(\sum_k\lvert\langle\hat\phi_k(t+s)\hat\phi_k^*(t)\rangle\rvert\): decay of the mode amplitudes without their phase rotation |
| `ell_x`, `ell_y` | scalars | first \(1/e\) crossings of the correlation cuts |
| `kubo` | scalar | `tau_c_euler` times the rms radial \(E\times B\) velocity over `ell_x` |
| `Gamma_n`, `D_n`, `D_alpha`, `D_E`, `D_W`, `E`, `W`, `zonal_fraction` | scalars | time averages of the ledger; \(D_n=\Gamma_n/\kappa^2\) |
| `unit_length`, `unit_time`, `unit_potential` | scalars | size of the normalised units in cm, s and statvolt (Gaussian cgs) |

Each scalar `X` has a companion `X_error`, one statistical standard error;
a negative value means "not available". Arrays are in the normalised units of
the model. Spectra are Welch estimates from Hann-windowed, half-overlapping
segments, normalised by the window power. Correlation functions are
transforms of the spectrum divided by the autocorrelation of the window, and
their lag range is a quarter of a segment. Errors of correlation scales come
from batches of contiguous segments, errors of time averages from the
integrated autocorrelation time of the series. `tau_c_mode` sums absolute
values of estimates and is therefore biased upward when few segments are
available. Global attributes name the model, the provider (`nonlinear` or a
quasilinear rule), the frame, the parallel structure and the derivative status.

`turbulence_statistics_oracles` recovers a synthetic Ornstein--Uhlenbeck mode:
Parseval identity to rounding, the Lorentzian spectrum, both correlation
times, both correlation lengths and the Kubo number within their standard
errors and within 5 %, and the standard error of a correlated series within
10 %. `hasegawa_wakatani_result_contract` runs the application
`hasegawa_wakatani` (namelist `&hasegawa_wakatani`) on a small saturated case
and reads the file back: Parseval identity, the flux recovered from the cross
spectrum, the mean balance and the unit and status attributes. A spectrum
normalised without the window power and a cross spectrum of the opposite sign
are rejected.

## Admission record

1. Use case: the consumer above; teaching use as the R3 gallery entry.
2. Model and ledger: the equations and balances above for stage HW; stages DA
   and K add theirs with their derivation before code.
3. Structural oracles: bracket identities, the two balances, the exact linear
   propagator of each Fourier mode.
4. Literature oracles: the dispersion relation above; saturated particle flux
   against adiabaticity from published two-dimensional runs. The published
   values differ between codes by tens of percent, so the acceptance band is
   stated per source.
5. Derivatives: the right-hand side and the linear roots carry JVP/VJP. A time
   average over a chaotic state has no usable derivative; it reports a
   statistical standard error and the status "derivative unavailable".
6. Certification: `candidate`. Linear roots carry their residual. Time
   averages carry a statistical interval from the autocorrelation time, which
   is never a deterministic enclosure.
7. Performance: one saturated \(256^2\) run of stage HW to \(10^3\) time
   units within minutes on one node; measured with the first benchmark.

## What the stages cannot represent

- HW: any parallel structure beyond one constant \(k_\parallel\); magnetic
  shear, Alfvén dynamics, magnetic flutter; ion temperature and finite Larmor
  radius; temperature gradients; curvature; trapped particles.
- DA: trapped-electron modes, toroidal drift resonances, curvature drive,
  electron-temperature-gradient turbulence, microtearing, shaping.
- L and K: saturation amplitudes (one declared constant per rule), zonal
  flows, radially non-local effects.
