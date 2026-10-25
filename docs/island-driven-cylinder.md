# Driven cylindrical island reference

This rung is separate from the reusable diagnostics and vacuum comparison.
Its NN2 gates remain unpromoted until the mesh, Fourier and current-regularisation
ladders establish the stated reference limit. Pressure is zero, so this case
cannot verify pressure flattening.

## Reference and conventions

[V] The reference assumes a circular, large-aspect-ratio cylinder, zero beta,
no flow, one (m,n)=(2,1) harmonic, a thin fully reconnected layer, constant-psi
matching and zero helical sheet current. Units are a=1 m, B_z=1 T and mu0=1.
For q=q0+c r^2, B_theta=r/(R0 q) and J=2q0/q^2, the outer equation is

    f'' + f'/r - m^2 f/r^2 - J' f/[r(1/q-n/m)] = 0.

This is Fitzpatrick's cylindrical tearing equation, Eq. (3.60) of the
[author's text](https://farside.ph.utexas.edu/teaching/plasma1/Fusionhtml/node55.html).
The matching convention is [f']_s=Delta' f_s+D f_b. Zero sheet current gives
f_s/f_b=-D/Delta', as in the static matching in
[Eq. (7.1)](https://farside.ph.utexas.edu/teaching/plasma1/Fusionhtml/node102.html).
The underlying paper is Fitzpatrick (1993),
[DOI 10.1088/0029-5515/33/7/I08](https://doi.org/10.1088/0029-5515/33/7/I08).
Its local Zotero body was checked for the distinct outer and inner problems.

[V] `island_outer_reference` integrates the axis-regular solution with FortNum
and transfers across r_s using equal finite Frobenius derivative parts. The
pole coefficient is A=2q0/r_s and
f=f_s(1+A x log|x|)+L x+O(x^2 log|x|). The cut is halved eight times, separately
from the ODE tolerance, and a zero-current-gradient reduction recovers r_s^m.
This is a numerical reference without a rigorous Frobenius remainder bound.

The optional radial samples are the signed edge-normalized solution f(r)/f_b.
The resonant scalar f_s/f_b is reported separately and is not a radial profile.
Each sampled field value is compared with the outer solution at the same radius,
whose boundary value is +1. Samples strictly inside the unresolved Frobenius cut,
apart from r_s itself, return a failure status rather than an interpolation.
Positive radii below the axis-regular seed cutoff rmin=1e-4 are also rejected;
the exact axis value is zero. Neither unresolved region is extrapolated.
The existing outer-reference tests check the edge and resonant normalization,
the independent zero-current-gradient profile f(r)/f_b=r^2, and the signed
Fourier projection of the imposed potential's analytic curl, including both
polarities. The latter is an extraction oracle, not a claim about the native
finite-aspect-ratio vacuum solution. These checks do not close the original
resonant response or model-matching gates below.

The local helical flux is

    H = (B_theta,s q'_s/q_s) x^2/2 - |f_s| cos(2 theta-phi-alpha).

[V] The O-to-X flux difference is 2|f_s|. Equating the separatrix flux to its
value at the X-point gives the full radial width

    W = 4 sqrt(|f_s| / |B_theta,s q'_s/q_s|).

This follows the constant-psi structure of
[Fitzpatrick, Eq. (8.55)](https://farside.ph.utexas.edu/teaching/plasma1/Fusionhtml/node113.html),
with the phase shifted to place the O-point at zero. The imposed toroidal
boundary datum is delta(R A_phi)=eps r^2 cos(2 theta-phi), so f_b=eps/R0.
Writing s_s=r_s q'_s/q_s gives

    W = 4 sqrt(2 |f_s/f_b|/s_s) sqrt(|eps|),

in the declared units. Positive response puts the O-point at helical phase
zero, while negative response puts it at pi. The Fourier convention is
B=B_0+2 Re(B_n exp(i n phi)); conversion to real blocks uses 2 Re(A_n) and
-2 Im(A_n). The boundary gauge is Lambda=i eps (R-R0-iZ)^2/2.

## Stable profile and gates

[V] The new stable case selects q0=0.5, c=2.5, hence r_s=sqrt(0.6) and
s_s=1.5. The outer cut ladder gives Delta'=-4.245974/m, D=4.361632/m and
f_s/f_b=1.027240. The prefactor is about 4.681 m in the declared eps units.
The statement of stability is restricted to the tested 2/1 outer problem.
Other rational surfaces exist, so toroidal sidebands and Fourier refinement
must still be checked.

The stable NN2 case uses R0/a=20, g_perp/g_parallel=1e-8,
g_parallel shape=(q0/q)^2 and fixed current -2 pi/[R0 q(a)]. The amplitude
ladder is 5e-4, 1.25e-4 and 3.125e-5. It retains the O/X classification and
phase gates, an exponent tolerance of 0.5 +/- 0.05 and the necessary prefactor
gate: tracing, duration and scan indicators plus W/R0. The O(a/R0) term is an
asymptotic scale, not a proven bound. No mesh allowance is invented from one
run. Comparing nr=8/16 and N=1/2 is required to quantify mesh/mode errors
before a final prefactor claim. A failed prefactor is retained while all three
widths are collected. No failed test is marked WILL_FAIL.

## Original sign reversal

[V] The original q0=1.2,c=1.8 case has Delta'=+4.524741/m and response
-0.4940385. Its original O/X expectation remains unchanged in its separate
registered test. The response study samples the actual radial-field harmonic
with the stated Fourier and gauge conventions before classifying roots.
On P3 nr=4 at g_perp=1e-4, the resonant flux response is +0.44972, close to
the vacuum r_s^2=0.44444. The opposite O/X signs therefore reflect opposite
resonant fluxes, rather than an absolute-value or residue-sign correction.

[V] NN2's current in `src/mhd/nonnested_coupled.f90` includes perpendicular
potential current -g_perp grad_perp Phi and does not impose the reference's
zero-sheet-current matching condition. The outer reference has neither term.
[H] A broad regularised layer explains the vacuum-like response: the leading
balance gives a radial smoothing scale proportional to
(g_perp R0^2/s_s^2)^(1/4), about 0.5 m for the original parameters. This
estimate is not a measured layer width or a completed causal isolation.
The response sweep keeps the original profile and varies only g_perp;
its converged mesh evidence is required before assigning the sign reversal
to that mechanism. Tearing instability and unresolved mesh/modes remain
additional possibilities. The stable profile is a new case, not a reversal
of the failed original expectation.
