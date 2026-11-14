# Edge thermal density limit: radiative fold, not an invented Greenwald theory

**Question:** For prescribed current density, edge heat transport, impurity
radiation and heating, at what **local edge electron density** does the hot
steady electron-temperature branch disappear?

This is an inexpensive 0D edge-shell energy model. It is **not** a derivation
of a universal Greenwald coefficient, nor a disruption predictor. The
empirical Greenwald line-averaged density is reported as a **separate
reference quantity**, not substituted for the model's local-edge threshold.

## Governing equation and units

At fixed edge density `n` and local current density `j_parallel`, write
the electron energy-balance residual in W/m³:

```math
C_e \dot T_e = F(T_e,n)
= P_{\rm ext} + \eta_{\parallel}(T_e)j_{\parallel}^2
- \frac{3}{2} n e\,\frac{\chi_e}{\delta^2}(T_e-T_{\rm sink})
- n^2\left[f_Z L_Z(T_e)+L_{\rm D,ff}(T_e)\right].
```

Here `T_e` and `T_sink` are in **eV**, `e` converts eV to joules,
`chi_e` is a constant prescribed electron heat diffusivity (m²/s),
`delta` is the shell's transport length (m), `P_ext` is deposition in
the shell (W/m³), `f_Z=n_Z/n_e` is a constant trace-impurity fraction.
The coefficient `C_e>0` is an effective electron heat capacity; its exact
value does not affect the steady fold or the sign of linear thermal growth.

Use the stated simple Spitzer-like closure

```math
\eta_\parallel = 5.2\times 10^{-5}
\frac{Z_{\rm eff}\ln\Lambda}{T_e^{3/2}}\quad[\Omega\,\mathrm{m}],
```

and the deuterium free-free cooling rate from published literature,

```math
L_{\rm D,ff} = 5.35\times 10^{-37}
\sqrt{T_e[\mathrm{keV}]}\quad[\mathrm{W\,m^3}].
```

The currently implemented smooth, tunable **surrogate** for impurity
radiation (not an atomic-physics calculation) is

```math
L_Z(T) = L_{\rm peak}\,
\frac{(a+b)(T/T_p)^a}{b+a(T/T_p)^{a+b}}.
```

This has peak exactly `L_peak` at `T_p` and low/high-temperature
asymptotic slopes `a` and `-b`. These values **must be fitted to
declared atomic cooling-rate measurements/tables** for predictive use.
Defaults are a transparent synthetic carbon-like illustration only; they
are not a universal physical fit, and uncertainties in emissivities,
charge-state balance and radial transport dominate the numerical error.

## How the solver actually finds a density threshold

Let `Q(T)=P_ext+eta(T) j^2`, `D(T)=(3/2)e chi_e (T-T_sink)/delta²`,
and `L(T)=f_Z L_Z(T)+L_D,ff(T)`. The steady equation is

```math
Q(T)-nD(T)-n^2L(T)=0.
```

Solving this quadratic without subtractive cancellation gives

```math
n_{\rm eq}(T) = \frac{2Q(T)}{D(T)+\sqrt{D(T)^2+4L(T)Q(T)}}.
```

A fixed-density electron-temperature perturbation is locally stable if
`F_T(T,n_eq)<0`. The solver tracks the **hot stable branch toward lower
temperature** and detects its first fold `F_T=0` by log-temperature
bracketing and bisection. The threshold is `n_eq(T_fold)`, in m⁻³.

No hot stable branch, invalid input, and no detected fold on the requested
temperature interval have distinct statuses; **no fold does not mean an
infinite density limit**. A fold is a local thermal instability of the
selected energy-balance model, not by itself a terminal tokamak disruption.
Current redistribution, magnetic islands, particle transport, atomic-rate
non-equilibrium, neutrals, H-mode physics, profile evolution and MHD modes
are outside this toy model.

For comparison only,

```math
n_{\rm GW} = 10^{20}\frac{I_p[\mathrm{MA}]}{\pi a^2[\mathrm m^2]}
\quad[\mathrm{m^{-3}}]
```

is the **empirical line-averaged** Greenwald reference for a circular plasma.
Never compare `n_edge/n_GW` as though both were the same observable:
a separately specified density profile is needed to convert edge density
to line-averaged density. The `j_parallel` used for the Ohmic term is a
**local** current density, not automatically `I_p/(pi*a²)`.

## Reproducible illustrative input (NOT a device prediction)

```text
I_p=1.0 MA, a=1 m; j_parallel=I_p/(pi a²) as an illustrative uniform-current assumption
P_ext=3e5 W/m³, chi_e=1 m²/s, delta=0.1 m, T_sink=2 eV
Z_eff=1, lnLambda=15, f_Z=0.005
synthetic L_peak=1e-31 W m³, T_peak=15 eV, a=2, b=1
scan 5...800 eV
```

Standalone gfortran test returns a thermal fold near
`T_fold=115 eV`, `n_edge=4.87e19 m^-3`. The independent empirical
reference at `1 MA, a=1 m` is `n_GW=3.18e19 m^-3` **line-averaged**.
The numerical proximity has no validation value because the radiative
coefficients and transport inputs are illustrative.

The solver has constant storage, no I/O and `bind(C)` for a WASM/JS host
(the 14-input ABI is documented in the source). Compile with the usual
MILCA dependencies once their WASM ports are available; no dependency
exclusion is implied. Actual browser compilation/performance remains
**unverified**.

## Scientific verification and next fidelity steps

- Test equilibrium density by substitution into the independent balance.
- Check `F_T=0` and loss of hot-branch stability across the fold.
- Vary power, impurity fraction and transport length: the predicted
  threshold **need not** follow the empirical Greenwald scaling.
- Fit `L_Z(T)` to an explicitly identified ADAS/experiment-derived dataset;
  record licensing, atomic regime and provenance. Then compare edge
  observables against independently published machine data.
- Add a spatial 1D edge heat equation and radiation-zone / thermo-resistive
  island coupling as independently verified higher-fidelity rungs.
  Reuse MILCA's radiation/transport and MHD/stability operators;
  do not add a separate density-limit universe.

## Public sources

- [NRL/Spitzer resistivity](https://www-internal.psfc.mit.edu/research/MFEFormulary/files/MFEFormulary.pdf),
  and F. F. Chen, *Introduction to Plasma Physics and Controlled Fusion*.
- Q. Teng, D. P. Brennan, L. Delgado-Aparicio, D. A. Gates et al.,
  [*A predictive model for the tokamak density limit*](https://collaborate.princeton.edu/en/publications/a-predictive-model-for-the-tokamak-density-limit/),
  Nuclear Fusion **56** (2016), radiative thermo-resistive mechanism.
- A. D. Maris et al., [*Correlation of the L-mode density limit with edge collisionality*](https://doi.org/10.1088/1741-4326/ad90f0),
  Nuclear Fusion **65** (2025): Greenwald does not fully predict all regimes.
- [Tokamak density limit review](https://link.springer.com/article/10.1007/s41614-024-00178-z)
  on heating, radiation and thermal collapse.

Implementation is new and based on public physical equations.
No non-public code of others was used or in context.
