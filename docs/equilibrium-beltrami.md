# Taylor closure on the non-nested complex

[V] `kin6d_nonnested_beltrami` solves a single fixed-boundary region with constant
mu and `curl B - mu B = f`. Pressure is constant. This is one regional closure
of the equilibrium problem, with no interface motion or traction balance.
The magnetic state uses the NN1 cylindrical FEEC complex times Fourier phi,
`B = curl A + h`, so divergence vanishes by the discrete complex identity.
The initialized field must have constant pressure. Nonconstant transported
pressure and nonzero p-prime are rejected. A simply connected meridional
section is required. Extra holes and their
additional periods are rejected.

## Equations and data

[V] `beltrami_prepare(field, flux, data, nphi, system, status [,map,nfp])` fixes
the toroidally averaged section flux, in magnetic-field times area units. For
a zero-normal wall this is the same flux through every section.
The harmonic is `h = F/R e_phi`, with `F = (flux - flux(curl A_boundary)) / integral(dR dZ/R)`.
The boundary circulation is integrated in the same FEEC trace convention.
Zero-boundary variations add zero toroidal flux. Under a 3D map these are reference quantities, and B is
pushed forward by the contravariant Piola map. The flux is unchanged.

[V] `data` supplies the tangential potential trace in the existing real Fourier
V1 layout: mean, cos(phi), sin(phi), cos(2phi), sin(2phi), etc. Its surface curl
prescribes normal B, exactly as in NN1. Interior entries are ignored. At n=0,
the scalar `a=R A_phi` fixes normal B and the poloidal trace fixes a gauge
and boundary circulation. At n>0, `a=0` is the
exact gauge and the poloidal tangential trace supplies the normal field.
Compatible traces can be obtained with `nn_mapped_data`. The harmonic
coefficient compensates the prescribed mean boundary circulation, retaining
the requested total flux.

[V] With C the discrete curl, M2 the magnetic mass, and P the V1-V2 pairing, the
physical equation is

    (C^T M2 C - mu P C) A = - C^T M2 h + mu P h + load(f).

[V] The n=0 gradient kernel is removed by a Coulomb multiplier with coupling
`M1 grad_0`, using zero-boundary scalar functions. This gives a symmetric
indefinite saddle system after lifting the boundary trace. No penalty or
spectral shift is added. FortFEM owns P and the exact block-diagonal Fourier
masses. Their element integrands come from retained FortSym products, and
assembly uses OpenMP. For a mapped domain, the magnetic Hodge tensor changes
while P is invariant under paired covariant/contravariant Piola maps.

[V] The optional `nfp` declares the field periodicity, defaulting to one (the
full torus). Only multiples of nfp are free. Incompatible map harmonics,
traces and source harmonics are rejected, with a rounding tolerance for the
interpolated map. Inertia covers this selected periodic space, and the report
records nfp. The SPEC comparison selects nfp=2 explicitly.

[V] `beltrami_solve` returns a candidate, equation residual, divergence, inertia,
pivot diagnostics, magnetic energy, helicity and current. An optional full-V1
weak source supports manufactured solutions. The optional `dmu` output solves
the differentiated equation with the same factor, including the harmonic load.
Geometry and trace derivatives and continuous error enclosures are unavailable.
For scalar mu, `dmu*d_mu` is the JVP and `dot(A_bar,dmu)` is the VJP.
The reported bound status says `unavailable`. A discrete residual is not a
field-error certificate.

[V] An opt-in matrix-free path (`beltrami_prepare(..., mr=...)`) solves the
same saddle system by block-preconditioned MINRES and the helicity constraint
by Newton-MR; it has no inertia diagnostics. See
[nonnested-newton-mr.md](nonnested-newton-mr.md#taylorbeltrami-closure-minres-on-the-saddle-system).

## Resonances and constraints

[V] The solve explicitly requests native symmetric-indefinite LDLT, with AMD when
available. It explicitly normalizes assembly symmetry within a relative 1e-12
budget, as NN1 does. It rejects uncertain inertia, zero pivots, or a minimum
pivot at most 1000 times FortSparse's rounding threshold. A rejected solve does
not publish a new B. These are numerical singularity diagnostics, not a proved
spectral gap or condition-number enclosure. No positive-definiteness or
energy-decrease test is used for the Taylor operator. If the equation
residual fails, FortNum shift-invert Lanczos checks the same unshifted factor
and reports a near-zero eigenvalue or an inconclusive diagnostic.

[V] For a closed wall, the helicity convention is the full-torus relative helicity
with reference h: `H = 2pi (A^T P C A + 2 A^T P h)`. The harmonic reference may
itself be mapped. This quantity is gauge invariant for zero tangential A.
For nonzero boundary data the reported quadratic is a fixed-gauge diagnostic,
and the helicity/current constraint interface rejects such open-wall data.

[V] `beltrami_constraint(system,'helicity',H,[mu_lo,mu_hi],tol,field,report)` uses
bracketed bisection. The bracket must lie on one side of zero and its endpoint
inertias must agree. Each iterate checks the same branch. The current option
uses `mu = mu0 I_tor / flux` for the homogeneous Beltrami equation, without an
iterative root solve. Both options propagate a resonant linear-solve status. Helicity convergence
means `abs(H-H_target) <= tol*max(1,abs(H_target))` in the declared units.
The interface does not construct homogeneous zero-flux eigenstates with a
helicity-selected amplitude at a singular mu; those need a bordered
eigenproblem and are rejected by this linear-solve path.

[V] The SPEC analogues were checked against local revision
b23a2408fc3f1dfbc14896e16d73a57f7d1d9e88:

- [V] [src/inputlist.f90:65](https://github.com/PrincetonUniversity/SPEC/blob/b23a2408/src/inputlist.f90#L65)
  and [src/mp00ac.f90:616](https://github.com/PrincetonUniversity/SPEC/blob/b23a2408/src/mp00ac.f90#L616):
  Lconstraint=0 fixes mu and flux in a plasma volume.
- [V] [src/mp00ac.f90:668](https://github.com/PrincetonUniversity/SPEC/blob/b23a2408/src/mp00ac.f90#L668):
  Lconstraint=2 uses `lABintegral-helicity` as its scalar constraint. Its
  normalization and gauge are SPEC's, and its input number is not copied into
  this relative-helicity convention without conversion.
- [V] [src/preset.f90:541](https://github.com/PrincetonUniversity/SPEC/blob/b23a2408/src/preset.f90#L541)
  sets the innermost-volume mu from `Ivolume/(tflux*phiedge)` for Lconstraint=3.
  [src/inputlist.f90:123](https://github.com/PrincetonUniversity/SPEC/blob/b23a2408/src/inputlist.f90#L123)
  defines `Ivolume=mu0 I`. Outer volumes additionally constrain surface currents.
  This implementation supplies only the single-volume current analogue.

## Verification cases

[V] `test_nonnested_beltrami` has P2/P3 manufactured-field and analytic-field
ladders, a nonzero-mode NN1 vacuum reduction, a mu-tangent check, helicity and
current recovery, and an inertia-bisection test that approaches a resonance
and requires rejection without changing B.

[V] The order ladders use affine square-section meshes and physical sample
points throughout the section, together with volume L2 errors. The torus
comparisons use curved P2 meshes. No optimal P3 order is claimed for the
existing graded polar P3 geometry, whose V2 best-projection error limits
convergence in the investigated curved fixtures.

[V] The manufactured field is the curl of an independently defined smooth
potential with modes 0 and 1. Its source is calculated analytically as
`curl B_m - mu B_m`. The analytic Lundquist field is solved in the declared
periodic-cylinder limit with axial period `2pi R0`, restricted to the exact square section
`-1 < x,y < 1` with its analytic nonzero normal trace:

    B_x = J1(mu r) sin(theta), B_y = -J1(mu r) cos(theta),
    B_axial = J0(mu r).

[V] The square-section flux is independently computed from a Stokes boundary
integral. A separate circular conducting-wall case, with
`flux=2pi J1(mu)/mu`, checks helicity/current recovery.

[V] The (x, axial, y) component ordering follows the cylindrical complex. This is
an exact cylinder solution, with zero cylinder model error. Interpreting it as
a finite-aspect-ratio torus would introduce an O(a/R0) model error, which is
not included in the numerical error. No such toroidal-exact claim is made.

[V] `benchmarks/beltrami/axis.sp` and `helical.sp` are synthetic public inputs for
SPEC, with R0=6, minor radius 1, mu=1.1 and toroidal flux 1. The helical case
has nfp=2 and a rotating elliptic section of amplitude 0.1. The optional native
SPEC comparison and reproduction commands are described in that directory.

[V] The native affine-section gate gives observed volume L2 orders 1.879/2.923
for the manufactured field at P2/P3, and 1.834/2.832 for Lundquist. These
are observed orders across two mesh doublings. The mean/nonzero-mode vacuum
coefficient differences are 4.53e-13 (P2) and 5.69e-12 (P3). The circular
constraint test recovers mu=1.1, and its mu derivative agrees with an independent
central difference to 2.59e-11. The resonance probe reports a nearby saddle
operator eigenvalue 4.16e-7 with eigensolver residual 2.88e-14 and rejects the
candidate. These diagnostics carry no continuous spectral-gap certificate.

[V] Against SPEC's finest retained solution, the eight-ring P2 relative vector
errors are 1.240e-3 (axisymmetric) and 1.216e-3 (helical), and the magnitude
errors are 7.263e-4 and 7.113e-4. The independently sampled element-curl
residuals are 9.861e-2 and 9.824e-2, while the algebraic equation residuals
are below 4.7e-14. The element residual excludes facet jumps. SPEC's finest
native maximum curvilinear-component residuals are below 2.8e-8 and 4.9e-8.
The SPEC M6-to-M8 vector differences are 2.40e-9 and 5.50e-9. The separate
mesh and mode differences are retained by the benchmark reducer, rather than
being added into a certified error bound.
