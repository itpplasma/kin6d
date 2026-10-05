# Finite-N ideal fermion ensembles (Q3)

`kin6d_fermion_ensemble` is a microscopic statistical reference for
noninteracting spin-1/2 fermions at fixed particle number. It shows where the
same-spin exchange hole, Fermi pressure and degenerate occupations come from
without any Coulomb interaction. It does not claim that ideal electrons form a
self-consistent plasma. The relativistic homogeneous EOS (Q8) extends
[fermi-gas.md](fermi-gas.md#relativistic-extension).

## Model

Each spin species is an independent canonical ensemble with N particles. The
two species share one spectrum.

- Periodic box with side L: orbitals `k = 2*pi*n/L`, `n` in `Z**3`, energy
  `|n|**2` in units of `eu = (hbar**2/2m)(2*pi/L)**2`. All orbitals belong to
  the model. The ones above the listed cutoff enter every enclosure through a
  rigorous remainder.
- Isotropic harmonic trap: levels `k*hbar*omega` above the zero point, with
  degeneracy `(k+1)(k+2)/2`.
- Closed spectra (`closed=.true.`) are finite orbital models with no remainder.
  The brute-force oracles use them.

`b = eu/(kB*T)` is the inverse temperature. The box is normalized to the
thermodynamic-limit Fermi energy of the same total density,
`EF/eu = (3N/(8*pi))**(2/3)`, with `N` the number of electrons summed over
both spins, and `Theta = kB*T/EF`.

## Canonical sums and cancellation handling

`Z_N` is the elementary symmetric polynomial of the orbital weights
`y = exp(-b*(eps-c))`. The code builds it with the subtraction-free recursion
`Z_k <- Z_k + y*Z_(k-1)`, one orbital at a time, and carries the first and
second b-derivatives along (jets). `E = -dlnZ/db` and `Var(E) = d2lnZ/db2`
follow without any numerical differencing. The shift `c` is an exact binary64
estimate of `F/N`, taken from the grand-canonical Chernoff bound. It changes
only the scaling: every reported quantity is independent of it.

- Occupations: `n_a = y_a * Z_(N-1)[shell a minus one orbital] / Z_N`, from
  prefix/suffix shell products.
- Same-shell pair occupations: two orbitals removed.
- Cross-shell pair occupations: the exact identity
  `<n_a n_b> = (y_a*n_b - y_b*n_a)/(y_a - y_b)`.

The textbook alternating recursion
`Z_N = (1/N) sum_k (-1)**(k+1) Z_1(k b) Z_(N-k)` is provided only as a
binary128 cross-check (`power_sum_recursion`), together with its cancellation
factor `max sum|terms|/|k Z_k|`. Measured factors:

| Theta | Cancellation factor | Consequence |
| --- | --- | --- |
| 1 to 4 | 1.1 to 2.6 | Agrees with the production sums to about 1e-13 |
| 0.25, N=33 | 3e27 | Loses most of binary128 |
| 0.05 | 1e34 to 5e36 | lnZ wrong by hundreds |

Production code therefore never uses the alternating recursion.

## Certification

`solve_canonical(..., certify)` repeats the same recursion in outward FortNum
intervals. It adds a rigorous remainder for the omitted orbitals.

The bounds on the omitted orbitals are `tau_m >= sum |eps-c|**m y`:

- Box: unit cubes around the lattice points, followed by the Gaussian moment
  recurrence with a Mills start. This needs `b*a**2 >= 2` with
  `a = sqrt(cutoff+1) - sqrt(3)`, and the condition is checked in intervals.
- Trap: a geometric ratio bound.

Listed shells with candidate occupation below `1e-18` are moved into the same
remainder.

With these bounds, the omitted sums satisfy

- `e_j <= tau0**j/j!`
- `|e_j'| <= tau1*tau0**(j-1)/(j-1)!`
- `|e_j''| <= tau2*tau0**(j-1)/(j-1)! + tau1**2*tau0**(j-2)/(j-2)!`

Together these enclose `lnZ`, `E`, `Var(E)`, every shell occupation, and the
total occupation of the omitted orbitals.

`box_eos` maps the two species to normalized outputs. Each output carries an
interval and is accepted with a 1% budget on `max(1,|value|)`:

```text
e = E/(N EF)
p = P/(n EF) = 2e/3                                   exact: levels scale as L**-2
de/dTheta = b**2 Var(E)/N = C_V/(N kB)
dp/dTheta
B/(n EF) = n dP/dn|_T/(n EF) = 5p/3 - (2/3) Theta dp/dTheta
```

The status of each output is as follows.

| Output | Status |
| --- | --- |
| lnZ, E, Var(E), occupations, EOS | `certified-model`: an enclosure of the named finite-N ideal model, not of a physical plasma |
| Pair correlation and density matrix (below) | Floating-point candidates |
| Pair-correlation truncation term | A bound in exact arithmetic on the shells omitted at `tolerance` |

Measured EOS enclosure widths for N=33 per spin:

| Theta | e, p | Derivatives |
| --- | --- | --- |
| 0.05 | 3e-9 | 4e-6 |
| 1 | 2e-11 | 7e-9 |
| 10 | 1.8e-8 | 5e-6 |

The widest bounds at the extremes come from interval dependency in the
signed jets.

## Observables

`same_spin_pair` returns the sphere-averaged same-spin pair correlation:

```text
g(r) = [N(N-1) - sum_{i/=j} <n_i n_j> sin(qr)/(qr)] / N**2
```

Here q is the orbital momentum difference. Pauli exclusion gives `g(0)=0` at
every temperature; this is exactly the sum rule `sum_{i/=j}<n_i n_j> = N(N-1)`,
and the test checks it. In a closed-shell ground state `g -> 1` at large r. At
T>0 a canonical offset of `1-1/N` remains, which is a real finite-N effect.

The other observables are:

- `one_body_density`: `rho(r)/n`.
- `opposite_spin_pair`: equals 1 to rounding, because the two species are
  independent and plane-wave densities are uniform.
- `solve_ground`: the T -> 0 limit. An open shell is filled uniformly over its
  determinants.
- `solve_boltzmann`: the classical Maxwell–Boltzmann preparation in the same
  box, `Z = Z1**N/N!`. It is recorded as a different preparation, and its
  occupations can exceed one.
- `grand_canonical_density`: the thermodynamic-limit reference hole
  `1 - rho**2`, computed from the certified homogeneous Fermi gas.

## Independent oracles (`test_fermion_ensemble`)

- **Brute-force enumeration** in binary128 of all `C(19,N)` box determinants
  (N=2..4, two temperatures) and all `C(10,3)` trap determinants. The checks
  cover lnZ, E, Var(E), occupations, pair g(r) to 1e-12, and containment.
- **Ground-state asymptotics**: `lnZ + b*E0 = 96*exp(-b)` for 19 particles at
  b=20.
- **Remainder check**: the open trap enclosure contains a closed trap with a
  much larger cutoff.
- **Cold thermodynamic-limit hole** `1 - 9[(sin x - x cos x)/x**3]**2` with
  `x = kF*r`. The maximum deviation falls from 3.0e-3 at N=33 to 2.1e-3 at
  N=251 and 8.6e-5 at N=1021, all per spin. This is convergence, not equality.
- **Grand-canonical Fermi gas** at Theta=1. The energy difference
  `e_N - e_TL` is -4.6e-3, -1.2e-3 and -6.1e-4 for N = 33, 123 and 251 per
  spin, roughly `1/N`. Canonical occupations of N=123 match Fermi–Dirac to
  8e-4.
- **Classical limits**: e=3Theta/2 within 1% at Theta=10. The Boltzmann box
  gives e=3Theta/2 and p=Theta to 1e-9. Fermion pressure exceeds the
  classical pressure.
- **Finite differences** of independent solves for:
  - `de/dTheta`;
  - `P=2E/(3V)` from `dlnZ/db`;
  - the bulk modulus, by rescaling the box at fixed `kB*T`.

## Evidence

`fermion_ensemble_demo` writes the following to `output/`:

- `fermion-exchange-hole.{csv,png}`: N=123 at Theta = 0, 0.25, 1, 2, compared
  with the grand-canonical, T=0 analytic, Boltzmann and classical g=1 curves.
- `fermion-occupations.{csv,png}`: canonical occupations compared with
  Fermi–Dirac and Maxwell–Boltzmann.
- `fermion-eos-scan.{csv,png}` and `fermion-eos-limit.csv`: the pressure and
  heat-capacity scan for Theta from 0.05 to 10.
- `fermion-size-scan.csv`: finite-N differences and the recursion
  cancellation.

The full demo takes about 25 s.

## Limits and open items

- **Exponent range**: the binary64 range with one constant shift is about
  `exp(0.11 N/Theta)` cold and `exp(N/e)` hot. Cold, large-N states return
  `candidate-range-failed`; per-degree scaling would remove this limit.
- **Cost**: N=123 at Theta=10 takes 14 s, because the interval recursion is
  `O(orbitals * N)`.
- **Pair-correlation cost** is `O(occupied orbitals**2)`, so the demo stops at
  Theta=2.
- **Not modeled**: interaction, spin polarization beyond fixed `N_up`,
  `N_down`, and trap pair correlations.
