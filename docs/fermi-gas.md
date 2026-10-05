# Finite-temperature ideal Fermi gas

This QF0 slice solves the chemical potential of a homogeneous, nonrelativistic
three-dimensional ideal electron gas with spin degeneracy two. Temperature and
density are supplied inputs. The production checker receives neither an exact
chemical potential nor an exact EOS. Its certified objects are momentum
integrals, the implicit density root, thermodynamic outputs and their equation
derivatives. Correlations, ion pressure, relativity, finite-wavelength Lindhard
response and physical-data uncertainty are outside this model certificate.

## Normalization and outputs

For density `n`, let `pF = hbar*(3*pi**2*n)**(1/3)`,
`EF = pF**2/(2*m)`, `theta = kB*T/EF`, `r = p/pF` and
`mu = chemical_potential/EF`. All supplied constants are caller data; the module
does not silently choose a species or a constants release. The occupation is

```text
f(r;theta,mu) = 1/(1+exp((r*r-mu)/theta)), theta > 0,
I2 = integral_0^infinity r*r*f dr,
I4 = integral_0^infinity r**4*f dr,
3*I2 = 1.
```

The integral convention agrees with the normalized Fermi–Dirac integral after
the substitution `u=r*r/theta`; its gamma normalization is stated in
[NIST DLMF 25.12.14](https://dlmf.nist.gov/25.12.E14).

The six normalized outputs, in derivative-array order, are
`[mu, pressure, energy, susceptibility, screening_ratio, bulk_modulus]`:

| Field | Definition |
| --- | --- |
| pressure | `P/(n*EF) = 2*I4` |
| energy | `U/(n*EF) = 3*I4`, with `U` the energy density |
| susceptibility | `chi = (EF/n)*(dn/dchemical_potential)_T = 3*I2_mu` |
| screening_ratio | `ks/kTF0 = sqrt(2*chi/3)` |
| bulk_modulus | `B/(n*EF) = 1/chi` |

The isothermal compressibility is `kappa_T = chi/(n*EF)`. It is different from
the bulk modulus. Static screening uses
`ks**2 = charge_magnitude**2*n*chi/(epsilon0*EF)`, for an infinitesimal
long-wavelength test charge. This does not certify a finite charge cloud or a
strongly coupled screening closure.

## Candidate and strict certificate

`initialize_fermi_workspace(work,max_cells,ok)` preallocates the velocity
edges, Gauss nodes and weights. At least 64 cells are required. The candidate
uses upstream FortNum Gauss–Legendre quadrature and safeguarded Newton, starting
with 128 nodes and refining when the independent checker exposes inadequate
candidate accuracy. These are untrusted floating-point approximations. No
Gauss-quadrature estimate is used as a certified error.

The retained FortSym generator proves the differential identities of the
Fermi occupation and emits their primal/JVP/VJP and interval expressions. For
positive exponent `z`, it evaluates `exp(-z)/(1+exp(-z))`; for negative `z`,
it evaluates `1/(1+exp(z))`. Both charts and their derivatives agree at zero.
The interval occupation uses monotonic endpoint evaluation and intersection
with `[0,1]`. Extreme derived overflow or unsupported arithmetic fails closed.

For fixed positive temperature, density is continuous and strictly increasing
in `mu`, since `I2_mu = integral r*r*f*(1-f)/theta dr > 0`. It tends to zero
and infinity at the two chemical-potential limits, giving existence and
uniqueness. A candidate bracket is
`[-theta*(4+2*log(1+theta)), 4+4*theta]`. The lower endpoint has density below
one by the Maxwell–Boltzmann upper bound and
`theta**(3/2)/(1+theta)**2 <= 1`. At the upper endpoint, `f(r)>1/2` on
`0<=r<=2`, so the normalized density exceeds four. These analytic facts help
the candidate; strict sign tests are still required for the final bracket.

The checker integrates independently on the exact stored binary64 edges from
zero to `R=sqrt(4+36*theta)`. Generic integration belongs to FortNum's
`composite_midpoint_batch`: on each cell it outwardly encloses the midpoint
value and adds the exact midpoint remainder bound `width**3*M2/24`.
Here the callback supplies bounds on the second momentum derivative throughout
the entire cell. FortSym obtains these by twice applying
`D_r = partial_r + f_r*partial_f`, with
`f_r = -2*r*f*(1-f)/theta`, to each of the eight moment integrands. Interval
dependency can widen these bounds but cannot invalidate them. There is no
sampling-based assertion of a whole-cell curvature bound.

For `r>=R>0`, `f <= exp((mu-r*r)/theta)`. Writing
`a=exp((mu-R*R)/theta)`, integration by parts and the Gaussian Mills bound give

```text
T2 = a*(theta*R/2 + theta**2/(4*R)),
T4 = a*(theta*R**3/2 + 3*theta**2*R/4 + 3*theta**3/(8*R)),
T6 = a*(theta*R**5/2 + 5*theta**2*R**3/4
        + 15*theta**3*R/8 + 15*theta**4/(16*R)).
```

These upper bounds enclose the missing second, fourth and sixth occupation
moments. The underlying recurrence is
`Gm = theta*R**(m-1)*exp(-R*R/theta)/2 + theta*(m-1)*G(m-2)/2`,
starting from `G0 <= theta*exp(-R*R/theta)/(2*R)`.
Since `|1-2*f|<=1` and `0<=f*(1-f)<=f`, the absolute tail bounds for
`[I2,I4,I2_mu,I2_theta,I4_mu,I4_theta,I2_mumu,I2_mutheta]` are

```text
[T2, T4, T2/theta, (T4+abs(mu)*T2)/theta**2,
 T4/theta, (T6+abs(mu)*T4)/theta**2,
 T2/theta**2, (T4+abs(mu)*T2)/theta**3+T2/theta**2].
```

Positive moments receive `[0,T]`; signed derivatives receive `[-T,T]`.
All displayed tail arithmetic is emitted and evaluated outwardly. The proof
does not depend on the candidate bracket making the tail small.
The owning quadrature implementation is pinned at FortNum
`015bb8b1490d99d55493b3936e0dd6adf77bc293`; the retained generator uses FortSym
`a6142712af998b727c245fb6cfbb3745b296ec24`, including its verified logarithm hook.

The final chemical-potential halfwidth is `0.005*theta`. Outward density
integrals must prove a strictly negative residual at its lower endpoint and
a strictly positive residual at its upper endpoint. The full EOS is then
evaluated over that chemical-potential interval. Each of its six output
intervals must lie within `0.01*max(1,abs(candidate))` of the candidate. This
is a stated output accuracy budget, not an arbitrary padding rule. The
exported intervals also include their floating-point candidate centers by a
hull, which preserves both enclosure and this acceptance budget. The
occupation amplification bound is `1/(4*theta)`; the nominal chemical-potential
budget therefore contributes at most 0.00125 to an occupation radius.
Increasing the quadrature budget cannot make an unresolved root or an
inaccurate candidate into a certificate. Exhausted budgets return `ok=false`.

## Derivative and SI input contracts

At fixed density, the implicit function theorem gives
`mu_theta = -I2_theta/I2_mu`. The other five derivatives follow from the same
EOS DAG and the integrated moment derivatives. A strictly positive interval
lower bound on `I2_mu` is required. `fermi_gas_jvp` and `fermi_gas_vjp` expose
these equation derivatives with Euclidean pairing on the six outputs. They
do not differentiate the adaptive quadrature schedule or the Newton iteration
map. Their intervals are rigorous but have no separate 1% width acceptance
gate; consumers must inspect widths after their observable amplification.
No temperature-curve interpolation or second-temperature-derivative certificate
is supplied. The domain is positive finite temperature; zero temperature is
a distinct limiting model, and adaptation is a scheduling event.

`solve_fermi_gas(work,theta,state,ok,certify,theta_interval)` supports an optional
positive finite interval containing the candidate temperature. Its strict
sign and moment tests hold for every temperature in that interval. The default
is the exact supplied binary64 point. `solve_fermi_gas_si` uses this facility:
its emitted normalization computes a temperature interval from the exact
supplied binary64 `n,T,m,hbar,kB,charge_magnitude,epsilon0`, then certifies the
normalized EOS over that interval before applying the emitted SI conversion.
Thus normalization rounding is propagated. Uncertain physical constants or
uncertain measured `n,T` are not included by pretending their central values
are exact physical data. An assembled SI-input tangent API is not yet exposed;
the normalization and SI-map leaves have generated JVP/VJP for future composition.

Certificates seal their prepared input, output and derivative packets. Altered
public centers or bounds, fabricated flags and candidate-only packets cannot
be used by the occupation or derivative queries. `certify=.false.` returns a
numerical candidate for matched solve-only measurements. Queries requiring
an enclosure reject it until a complete certified solve succeeds.

## Independent evidence and limits

`test_fermi_gas` has an independent binary128 composite-Simpson integrator and
bisection root oracle, with a second resolution check. It does not call the
production moment leaves for expected values. Cold `theta=.05`, crossover
`theta=1` and classical `theta=10` cases check EOS, occupation and implicit
derivative enclosure; an independent dual contraction checks the VJP. Exact
half occupation and chart derivatives, large-exponent stability, SI momentum,
thermodynamic identities and corrupted-packet cases supplement these oracles.
This is independent numerical evidence, not a formal proof of the oracle.

The leading zero-temperature values are
`[mu,P/(nEF),U/(nEF),chi,ks/kTF0,B/(nEF)] = [1,.4,.6,1.5,1,2/3]`.
The Maxwell–Boltzmann limit has `P/(nEF)=theta`, `U/(nEF)=1.5*theta`,
`chi=1/theta` and `ks/kTF0=sqrt(2/(3*theta))`. These are limiting references;
no uniform asymptotic remainder is asserted for a finite temperature. The
actual finite-temperature integral/root certificate remains authoritative.

`benchmark_fermi_gas [timing.csv]` retains nine raw warm timing repetitions
at the same supplied states and selected candidate/strict grids. Setup and
validation are outside timed sections. Its baseline performs exactly the
same floating-point candidate solve as the complete path. It also writes one
illustrative SI electron-pressure row for `rho=1e6 kg/m**3`, `Ye=.5`,
`T=1e6 K`, using explicitly supplied constants. `pF/(m*c)≈.0801` supports the
nonrelativistic use; this is not a hydrostatic star or a relativistic core.
The example excludes ion/Coulomb pressure and constants uncertainty. Run it
in a scratch output directory; generated CSV/media are not committed.

An isolated GNU Fortran `-O2 -fno-fast-math -ffp-contract=off` measurement on an
AMD Ryzen 9 5950X retained these nine-trial warm CPU-time ranges. Both paths use
128 candidate nodes; setup/refinement is excluded from this table.

| theta | strict cells | complete median (range), ms | candidate, microseconds |
| --- | --- | --- | --- |
| .05 | 256 | 10.312 (9.893–10.803) | 18.15–25.9 |
| 1 | 64 | 2.478 (2.476–2.484) | 17.5–17.6 |
| 10 | 64 | 2.577 (2.478–2.588) | 19.25–20.1 |

The complete-to-candidate factors are roughly 385–587, 141 and 129. The strict
certificate therefore has a large relative overhead despite its millisecond
absolute cost. Batch sharing and a repaired adaptation schedule are measured
here; no claim of a low-overhead certified workhorse or an optimized Fermi
quadrature method is made. Raw repetition CSVs, compiler flags and input
selection belong to the execution receipt. A future optimized quadrature must
preserve the independent integral/tail/root argument rather than replace it
with unchecked high-order quadrature.

The strict arithmetic assumes the supported FortNum outward interval runtime,
finite endpoints, and builds without reassociation or contraction. Every
callback must bound its derivative throughout its supplied cell. GNU Fortran
hosted callbacks need trampoline storage. The controller build uses the
compiler's heap-trampoline option where supported; the isolated measurements
above used its default stack-trampoline ABI. This is a portability distinction,
not additional enclosure slack.

## Relativistic extension

Campaign case Q8. `solve_fermi_gas(..., relativity=x, relativity_interval=X)`
uses the relativistic kinetic energy with `x = pF/(m*c)`. It keeps the
normalization above (`r = p/pF`, `theta = kB*T/EF`, `mu` in units of
`EF = pF**2/(2m)`). In these units

```text
eps(r) = E/EF = 2r**2/(s+1),   s = sqrt(1 + x**2 r**2),
P/(n EF) = integral r**3 eps'(r) f dr,
U/(n EF) = 3 integral r**2 eps f dr,
3 integral r**2 f dr = 1.
```

The six outputs keep their meaning. The relations `chi = 3*I2_mu` and
`B = 1/chi` follow from Gibbs–Duhem at fixed T. `U = 1.5P` holds only
nonrelativistically. In general `1.5P <= U < 3P`, because the ratio of
kinetic energy to `p*v` is `s/(s+1)`.

If `x` is absent or zero, the existing nonrelativistic leaves run bit for bit.
`solve_fermi_gas_si(..., speed=c)` derives `x` and its interval from the SI
inputs.

FortSym generates the relativistic leaves in `gen/gen_relativistic_fermi.f90`
(`kin6d_relativistic_fermi_leaves`). They carry `s` as a chart variable with
`ds/dr = x**2 r/s`. The chart identities are proved on `s**2 = 1 + x**2 r**2`,
both the occupation derivatives and `d eps/dr = 2r/s`.

There are 11 moments, compared with 8 in the nonrelativistic case: density,
pressure, energy, their `mu`/`theta` derivatives, `I2_mumu` and
`I2_mutheta`. The whole-cell second r-derivatives are emitted as strict
interval leaves.

**Tail bound.** `eps'' = 2/s**3 > 0`, so `eps` is convex. For `r >= R` this
gives `f <= a*exp(-kappa(r-R))` with `a = exp((mu - eps(R))/theta)` and
`kappa = eps'(R)/theta`, and

```text
integral_R^oo r**m f <= a * sum_j m!/(m-j)! R**(m-j) / kappa**(j+1).
```

Each moment integrand is then bounded with `eps <= r**2`, `eps' <= 2r`,
`f(1-f) <= f` and `|1-2f| <= 1`.

**Radius.** `R` solves `eps(R) = 4 eps(1) + 36 theta` exactly through
`r(eps) = sqrt(eps(1 + x**2 eps/4))`.

**Candidate bracket.** The upper end is `eps(2) + 4 theta`. The lower end is
extended until the candidate density residual changes sign.

**Graded mesh.** Degenerate relativistic states use a three-segment mesh with
the `1/4 : 1/2 : 1/4` cell split. It covers the window
`1 +- 20 theta sqrt(1+x**2)` around the cold Fermi surface. Candidate
Gauss–Legendre nodes use the same split. The mesh only places cells: the
midpoint/curvature certificate holds for any edges.

**Results.** At T=1e7 K, the relativistic path certifies every point from
`x = 0.01` to `x = 10`, with `theta` down to 3.4e-5, in at most 2048 cells and
0.1 s per point. Root band, accuracy budget, implicit theta-derivatives and
sealed packets are unchanged.

**Independent oracle.** `test_relativistic_fermi_gas` uses binary128
composite Simpson. It integrates in rapidity `r = sinh(t)/x` when the state is
nondegenerate, which removes the `r = +-i/x` branch-point limit. It checks:

- containment of all six outputs and their theta-derivatives for
  `(theta, x) = (.05, 1), (1, 1), (10, 3), (.002, 5)`;
- `x = 1e-6` against the nonrelativistic certificate, to 1e-9;
- the cold ultrarelativistic limit `P = (hbar c/4)(3 pi**2)**(1/3) n**(4/3)`,
  which is `1/(2x)` normalized, at `x = 50`, to 1e-3;
- the exact cold Chandrasekhar curve
  `[x(2x**2-3)s + 3 asinh x]/(4x**5)`, to 1e-4 at `x = 50` and 1e-6 at
  `x = 0.1`;
- the classical limit `P = n*kB*T` at `theta = 300`, to 1e-3;
- a centered finite difference of `dP/dtheta`;
- rejection of an altered relativity parameter;
- an SI state at `rho = 1e9 kg/m**3`.

**Demo.** `relativistic_fermi_demo` writes
`output/relativistic-fermi-{pressure,ratios,occupation}.png` and the CSV. The
nonrelativistic uniform-mesh certificate fails closed below `theta` of about
1e-3, at `x >= 2.4` on that curve. Those points are omitted from its curve.
Applying the graded mesh to the nonrelativistic path would change its
committed numbers; that decision is left open.

**Scope.** These are electrons only. Ion pressure, Coulomb/exchange-correlation
corrections, pair creation and hydrostatic structure are outside the model.
