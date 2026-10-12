# Electromagnetic wave at a plasma edge

Teaching case for the plasma-oscillation chapter: a plane electromagnetic wave
at normal incidence on a cold electron plasma with a sharp edge. Below the
electron plasma frequency the wave is reflected and the field decays inside
over the skin depth; above it the wave is transmitted with the cold-plasma
wavelength. It is the first kin6d case that couples particles to the full
(one-dimensional) Maxwell equations rather than to an electrostatic field.

Status: **candidate**. The dynamics is an untrusted PIC realization; the
reported reflection/transmission coefficients, skin depth and wavelength are
lock-in estimates compared with the linear cold-plasma oracle. No enclosure is
claimed. Files: `src/em_slab.f90`, generator `gen/gen_em_slab.f90` with the
emitted leaves `src/generated/kin6d_em_slab_leaves.f90`, tests
`test/test_em_slab.f90`, application `app/em_slab_demo.f90`.

## Model

Relativistic Vlasov--Maxwell for electrons in one spatial dimension \(x\),
with momentum components \(p_x, p_y\) and fields \(E_x, E_y, B_z\); ions are a
fixed neutralizing background of density \(n_0\) for \(x>x_e\) (zero in the
vacuum \(x<x_e\)). In SI,

\[
\partial_t B_z = -\partial_x E_y,\qquad
\partial_t E_y = -c^2\partial_x B_z - J_y/\varepsilon_0,\qquad
\partial_t E_x = -J_x/\varepsilon_0,\qquad
\partial_x E_x = \rho/\varepsilon_0,
\]

\[
\frac{d\mathbf p}{dt} = -e\,(\mathbf E + \mathbf v\times\mathbf B),\qquad
\mathbf v = \frac{\mathbf p}{\gamma m_e},\qquad \gamma=\sqrt{1+p^2/(m_ec)^2}.
\]

The ignorable transverse coordinate \(y\) is still integrated for each sheet,
\(\dot y = v_y\), because the animation shows the transverse electron motion.
The \(v\times B\) force and the longitudinal field \(E_x\) are retained; they
enter only at second order in the wave amplitude.

Normalization: \(t\,\omega_{pe}\), \(x\,\omega_{pe}/c\),
\(\hat E = eE/(m_ec\,\omega_{pe})\), \(\hat B = ecB/(m_ec\,\omega_{pe})\),
\(u = p/(m_ec)\), \(\hat J = J/(en_0c)\), \(\hat\rho=\rho/(en_0)\). Then

\[
\partial_t \hat B = -\partial_x \hat E_y,\quad
\partial_t \hat E_y = -\partial_x \hat B - \hat J_y,\quad
\partial_t \hat E_x = -\hat J_x,\quad
\partial_x \hat E_x = \hat\rho,\quad
\dot{\mathbf u} = -(\hat{\mathbf E} + \mathbf v\times \hat{\mathbf B}),
\]

field energy density \((\hat E_x^2+\hat E_y^2+\hat B^2)/2\) and Poynting flux
\(\hat E_y\hat B\), both in units of \(n_0m_ec^2\) and \(n_0m_ec^3\). The
dimensionless amplitude is \(a_0 = eE_0/(m_e\omega c)=\hat E_0\,\omega_{pe}/\omega\);
the default \(\hat E_0 = 5\times10^{-3}\) gives \(a_0=10^{-2}\) at
\(\omega=\omega_{pe}/2\) and \(2.5\times10^{-3}\) at \(\omega=2\omega_{pe}\).

## Linear cold-plasma theory (oracle)

Linearize about the cold uniform state with \(e^{i(kx-\omega t)}\). The
transverse electron momentum obeys \(-i\omega m_e v_y = -eE_y\), so
\(J_y = -en_0v_y = i\varepsilon_0\omega_{pe}^2E_y/\omega\) and the
susceptibility is \(\chi=-\omega_{pe}^2/\omega^2\). Eliminating \(B_z\) gives

\[
\omega^2 = \omega_{pe}^2 + c^2k^2,\qquad
N^2 = \frac{c^2k^2}{\omega^2} = 1-\frac{\omega_{pe}^2}{\omega^2}.
\]

At normal incidence \(E_y\) and \(B_z\) are transverse and continuous at the
edge (no surface current or charge at first order). With incident amplitude 1,
\(E_y=e^{ik_0x}+Re^{-ik_0x}\) and \(cB_z = e^{ik_0x}-Re^{-ik_0x}\) for
\(x<0\), \(E_y = Te^{ikx}\), \(cB_z=NTe^{ikx}\) for \(x>0\):

\[
1+R=T,\qquad 1-R=NT\quad\Longrightarrow\quad
R=\frac{1-N}{1+N},\qquad T=\frac{2}{1+N}.
\]

* \(\omega>\omega_{pe}\): \(N\) is real, \(R^2+NT^2=1\) (Poynting balance),
  \(k=\sqrt{\omega^2-\omega_{pe}^2}/c\), \(\lambda=2\pi/k\).
* \(\omega<\omega_{pe}\): \(N=i\kappa\), \(\kappa=\sqrt{\omega_{pe}^2/\omega^2-1}\),
  \(R=(1-i\kappa)/(1+i\kappa)\), \(|R|=1\), \(\arg R=-2\arctan\kappa\),
  \(|T|=2\omega/\omega_{pe}\). The field inside is evanescent,
  \(E_y\propto e^{-x/\delta}\) with the skin depth
  \(\delta = c/\sqrt{\omega_{pe}^2-\omega^2}\).

The steady vacuum envelope is \(|E_y|=|1+Re^{-2ik_0x}|\). `gen_em_slab`
proves with the FortSym zero test, before emitting: both matching conditions
and the Poynting balance for real \(N\); \(|R|=1\) and the two matching
conditions for \(N=i\kappa\); the \(\omega\) forms of \(R\) and
\(|T|^2=4\omega^2/\omega_{pe}^2\); the dispersion relation; and that the Boris
rotation preserves \(|u|\) for every rotation parameter. The oracle leaves
(`em_slab_fresnel`, `em_slab_evanescent`, the envelopes) are never called by
the solver.

## Discretization

* **Grid.** Yee staggering on \([x_L,x_R]\): \(E_y\) at nodes \(i\Delta x\);
  \(B_z\), \(E_x\), \(J_x\) at half nodes; \(E\) at integer and \(B\), \(J\),
  \(u\) at half-integer time levels; \(c\Delta t=0.95\,\Delta x\), which
  also satisfies the Yee--plasma bound \((c\Delta t/\Delta x)^2+(\omega_{pe}\Delta t/2)^2<1\).
* **Boundaries.** First-order Mur outgoing condition for \(E_y\) at both ends
  (exact for normally incident vacuum waves up to grid dispersion). The plasma
  reaches two cells before \(x_R\); `em_slab_run` stops with an error if a
  signal travelling at \(c\) from the source could reach \(x_R\), so the plasma
  is effectively semi-infinite and the right Mur condition is never exercised
  by a physical wave.
* **Source.** A soft sheet current \(\hat K(t)=-2\hat E_0\,s(t)\sin\omega t\)
  at \(x_s\) radiates \(\hat E_y=-\hat K/2\) to both sides, with
  \(s=\sin^2(\pi t/2t_r)\) for \(t<t_r\) and 1 afterwards. The left-going half
  leaves through the left Mur boundary; reflected waves pass the transparent
  source.
* **Particles.** Cold quiet start: \(p\) equally spaced sheets per cell, weight
  \(\Delta x/p\), so the CIC-deposited density is \(n_0\) inside and a one-cell
  linear ramp centred at \(x_e\). The ion background is the deposited initial
  electron charge, so the initial state is exactly neutral and field-free.
  Gather: CIC from nodes (\(E_y\)) and half nodes (\(E_x\), \(B_z^n\)).
  Push: relativistic Boris with \(B^n=\tfrac12(B^{n-1/2}+B^{n+1/2})\) obtained
  by two half Faraday steps. Deposit: \(J_y\) with CIC at the midpoint position;
  \(J_x\) by the one-dimensional charge-conserving rule
  \(J_{i+1/2}-J_{i-1/2} = -q w (S_i^{n+1}-S_i^n)/\Delta t\), so the discrete
  continuity equation and hence the discrete Gauss law hold to rounding.
* **Leaves.** Boris map, velocity, CIC weights, Faraday/Ampere updates, Mur,
  ramp/sheet source, energy densities and the vacuum characteristics
  \((\hat E_y\pm\hat B)/2\) are FortSym leaves with generated JVP/VJP.
* **Threads.** The push/deposit loop is OpenMP-parallel over static particle
  blocks with private current columns, summed in thread order: fixed thread
  count is bitwise reproducible; a different thread count changes rounding
  only.

## Measurements

Over a window of whole periods after the flat top has reached the probes,
the lock-in amplitude \(a=(2/T_w)\sum f(t_n)e^{i\omega t_n}\Delta t\) is taken
for \(E_y\) on every node and for the vacuum characteristics
\(F_\pm=(\hat E_y\pm\hat B)/2\) at the probe \(x_p\) (between source and edge,
\(B\) averaged to the node and time level). Then
\(R=(a_-/a_+)\,e^{2ik_0(x_p-x_e)}\); \(|T|\) is the mean \(|a(x)|/|a_+|\) on
\([1,10]\,c/\omega_{pe}\) inside; \(k\) is the least-squares slope of the
unwrapped phase on the same interval; \(\delta\) is the least-squares decay
length of \(\log|a(x)|\) on \([0.25,3]\,c/\omega_{pe}\).

Energy ledger (per area, right of the probe): \(\int F_+^2dt\) (incident) \(=\)
\(\int F_-^2dt\) (reflected) + \(\int \hat E_y\hat B\,dt\) at the last interior
node (outflow) + change of field energy + change of electron kinetic energy
\(\sum w(\gamma-1)\), with \(\gamma\) at integer levels taken as the mean of
the adjacent half levels. The reported residual is normalized by the incident
energy.

## Verification (CTest `em_slab_*`, defaults \(\Delta x=0.04\), \(p=8\))

| Test | Quantity | Measured | Linear theory | Tolerance |
|---|---|---|---|---|
| `em_slab_oracles` | Boris rotation, kick, JVP vs central differences, VJP pairing; oracle leaves vs complex matching; Gauss law after a driven run (\(a_0\approx0.4\)) | Gauss residual \(3.6\times10^{-14}\) | -- | rounding |
| `em_slab_vacuum` | phase speed between probes 4 apart (\(\omega=2\)) | \(1.000018\,c\) | \(c\) | \(2\times10^{-4}\) |
| | \(|a_-/a_+|\) without plasma | \(7.6\times10^{-4}\) | 0 | \(2\times10^{-3}\) |
| `em_slab_reflection` (\(\omega=0.5\)) | \(|R|\) | 0.999591 | 1 | \(5\times10^{-3}\) |
| | \(\arg R\) | \(-2.094844\) | \(-2.094395\) | \(5\times10^{-3}\) |
| | \(\delta\) | 1.154797 | 1.154701 | \(2\times10^{-3}\) rel. |
| | energy residual | \(-3.8\times10^{-5}\) | 0 | \(2\times10^{-3}\) |
| `em_slab_transmission` (\(\omega=2\)) | \(R\) | 0.072127 | 0.071797 | \(2\times10^{-3}\) |
| | \(T\) | 1.072370 | 1.071797 | \(2\times10^{-3}\) rel. |
| | \(\lambda\) | 3.627596 | 3.627599 | \(10^{-4}\) rel. |
| | energy residual | \(-8.3\times10^{-4}\) | 0 | \(2\times10^{-3}\) |
| `em_slab_refinement` | \(|R|-R\), \(T/T_{th}-1\) at \(\Delta x=0.08\to0.04\) | \(1.27\times10^{-3}\to3.3\times10^{-4}\), \(3.3\times10^{-3}\to5.3\times10^{-4}\) | -- | must decrease |

The discrete Gauss residual is \(3.4\times10^{-13}\) at the end of both
production runs. The residual vacuum \(|a_-/a_+|\) is the
\(O((k_0\Delta x)^2)\) staggering error of the characteristic split, not a
boundary reflection. The refinement indicator shows that the remaining
Fresnel deviations are dominated by discretization (about second order);
the finite amplitude enters at \(O(a_0^2)\).

At \(\omega=\omega_{pe}/2\) the second-harmonic \(v\times B\) force oscillates
at \(2\omega=\omega_{pe}\) and resonantly drives longitudinal Langmuir motion
near the surface. Its displacement grows secularly at \(O(a_0^2)\)
(\(\lesssim4\times10^{-2}\,c/\omega_{pe}\) for the edge sheet by
\(t\omega_{pe}=125\), \(\lesssim10^{-3}\) one skin depth inside) and leaves the
first-order quantities above unchanged within their tolerances; it is real
cold-plasma physics of this frequency choice, not a numerical instability.

## Differentiability and certification

The Boris map and all field-update, source, energy and oracle leaves carry
generated JVP/VJP; `em_slab_oracles` checks the Boris JVP against central
differences and the VJP pairing. Derivatives of whole trajectories or of the
lock-in observables are **unavailable**: CIC cell crossings make the deposit
piecewise smooth, and no discrete adjoint is provided. Smooth parameter
derivatives of the oracle (for example \(dR/d\omega\)) are available from the
generated leaves.

No certificate is claimed. The tolerances above are verification bounds for
this candidate realization; the independent indicator is the refinement
test. A rigorous version would need a continuous model-error bound for the
PIC sampling and the edge ramp, which is not planned for this teaching slice.

## Performance

Default demo (two runs, about 27000 sheets, 5600 steps in total): 11.6 s on 1
thread, 6.4 s on 2, 4.9 s on 4 and 2.8 s on 8 (shared workstation, Release).
The CTests run with `OMP_NUM_THREADS=4` in about 7 s together.

## Reproduction and lecture data

```sh
cmake --preset cpu && cmake --build --preset cpu
ctest --preset cpu -R em_slab
em_slab_demo [--config FILE] [--export DIR]   # namelist &em_slab
```

The application prints the comparison table, writes `output/em_slab_demo.nc`
and with `--export DIR` the lecture-scene data
`DIR/em-wave-plasma-reflection/` through the animation export: \(E_y/E_0\) on
the display grid \([-15,10]\,c/\omega_{pe}\) per frame for both runs, the
transverse and longitudinal displacements of sample electron sheets in
\([0,4]\,c/\omega_{pe}\), the linear-theory and measured steady envelopes, and
the measured and theoretical coefficients as declared parameters. Regenerate
the leaves with a FortSym build at the pinned revision
(`-DKIN6D_FORTSYM_BUILD=...`, target `gen_em_slab`, run from the source
directory).
