# Saturated Hasegawa-Wakatani benchmark

The original Hasegawa-Wakatani model ([docs/drift-turbulence.md](../../docs/drift-turbulence.md),
stage HW) on the doubly periodic box of side \(2\pi/0.15\) with
\(\kappa = 1\). Lengths are in \(\rho_s\), time in \(L_n/c_s\), and the particle
flux \(\Gamma_n = -\langle n\,\partial_y\phi\rangle\) in \(n_0 c_s (\rho_s/L_n)^2\).
Errors of MILCA values are standard errors of the time mean from the
autocorrelation time of the flux series; they are statistical intervals, not
enclosures.

## Inputs

- `cases.tsv`: one case per line (adiabaticity, drive, dissipation, order
  \(p\) of the damping \(d\,k^{2p}\), modified or original model, grid, time
  step, end of run, start of the average, sampling, seed, threads).
- `scan.sbatch`: one case per Slurm array task, with `HW_ROOT` holding
  `cases.tsv` and the executable `hasegawa_wakatani`.
- `reference.json`: published values with their location and read status.

## Results

Runs of `hasegawa_wakatani` built from the revision that added the
application (its solver and statistics are unchanged since), GNU Fortran
12.2, AMD EPYC 7452, \(d = 5\times10^{-8}\), \(p = 3\), \(\Delta t = 0.01\)
(0.005 at \(\alpha = 0.1\)), average over \(t \ge 300\) (250, 400 and 900 at
\(\alpha = 0.1\), 2 and 5):

| \(\alpha\) | MILCA \(256^2\) | hw2d \(512^2\) (std. dev. in time) |
|---:|---:|---:|
| 0.1 | 1.499 ± 0.035 | 1.51 (0.21) |
| 0.5 | 0.881 ± 0.011 | 0.89 (0.09) |
| 1.0 | 0.584 ± 0.006 (three seeds) | 0.60 (0.05); 0.60 ± 0.01 over 25 runs |
| 2.0 | 0.364 ± 0.006 | 0.38 (0.03) |
| 5.0 | 0.115 ± 0.002 | 0.13 (0.02) |

Every MILCA value lies within one hw2d time standard deviation of the hw2d
value. The two quantities differ: the hw2d bar measures the fluctuation of the
flux, the MILCA bar the uncertainty of its mean.

Resolution and dissipation at \(\alpha = 1\):

| case | \(\Gamma_n\) |
|---|---:|
| \(256^2\), \(\Delta t = 0.005\) | 0.590 ± 0.006 |
| \(512^2\), \(\Delta t = 0.005\) | 0.616 ± 0.008 |
| \(256^2\), \(d = 10^{-6}\) | 0.659 ± 0.005 |
| \(128^2\), \(d = 10^{-5}\) (the CTest case) | 0.668 ± 0.002 (twelve seeds) |
| \(128^2\), \(d = 5\times10^{-8}\) | 0.140 ± 0.0002 |
| \(256^2\), \(d = 5\times10^{-9}\) | 0.303 ± 0.005 |
| \(256^2\), \(p = 2\), \(d = 10^{-4}\) (Numata et al.) | 0.700 ± 0.008 |
| \(512^2\), \(p = 2\), \(d = 10^{-4}\) | 0.710 ± 0.011 |

The time step is converged within the statistical error. Between \(256^2\)
and \(512^2\) the flux changes by 0.032, four times the combined standard
error. A dissipation too weak for the grid (\(128^2\) with
\(d = 5\times10^{-8}\), \(256^2\) with \(d = 5\times10^{-9}\)) lowers the
flux by a factor of two to four; the cause was not diagnosed. The hw2d
README quotes \(d = 5\times10^{-8}\) in its text and \(5\times10^{-9}\) in
its table. Among the resolved cases the dissipation choice moves the flux by
15 to 20 %, comparable to the spread between published codes.

The modified model with the settings of Numata et al. (\(256^2\), \(p = 2\),
\(d = 10^{-4}\)) saturates at \(\Gamma_n = 0.025 \pm 0.001\) with 90 % of the
kinetic energy in the zonal flow, against 4 % in the original model.

## Acceptance tests

`hasegawa_wakatani_benchmark_flux` (label `slow`, about 11 minutes on one
core) runs the \(128^2\) case with \(d = 10^{-5}\) to \(t = 1000\) and
requires the flux within 0.04 of 0.668, a standard error below 5 % of the
flux, the mean energy balance within 1 % of the flux, the flux balanced by
the resistive sink within 10 % and the hyperviscous sink below 10 % of the
flux, as Numata et al. state for the saturated state, and a zonal share of
the kinetic energy below 0.2.

The reference is the mean of this case over the seeds 1 to 12 (0.659 to
0.685, standard deviation 0.008, standard error 0.002). The test runs seed 1;
another seed, thread count or compiler changes the chaotic trajectory, and
the band of five seed standard deviations covers that. It lies below a 10 %
error of the flux (0.067). The flux at this resolution and damping is 14 %
above the \(256^2\), \(d = 5\times10^{-8}\) value of the table above; the
test checks this case, not the converged flux.

Two variants have to be rejected. `hasegawa_wakatani_benchmark_zonal_mutant`
runs the modified model (flux 0.021). `hasegawa_wakatani_benchmark_coupling_mutant`
runs \(\alpha = 1.2\), a 20 % error in the parallel coupling, which gives
0.588 to 0.599 over three seeds: only the flux band rejects it, and the
earlier band (0.584 ± 0.10) accepted it.

## Cost

\(256^2\) to \(t = 1000\) (\(10^5\) steps) takes 2000 to 2400 s on four
threads and about 2100 s extrapolated on twelve; \(512^2\) to \(t = 600\)
takes 9600 s on twelve threads. Peak memory is 0.13 GiB at \(256^2\) and
0.17 GiB at \(512^2\).

## Sources

- R. Greif, HW2D, J. Open Source Softw. 8, 5959 (2023): values from the
  repository README at commit `623841f1`; the refereed text holds no numbers.
- R. Numata, R. Ball, R. L. Dewar, Phys. Plasmas 14, 102312 (2007): model,
  dispersion relation and the statements on the saturated state; it
  tabulates no flux.
- S. J. Camargo, D. Biskamp, B. D. Scott, Phys. Plasmas 2, 48 (1995):
  \(\Gamma_n = 0.73\) at \(\alpha = 1\) is quoted from the hw2d README and
  was not checked against the paper.
