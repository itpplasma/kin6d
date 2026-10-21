# S0 physical-mass spectrum: plasma-volume slice

`kin6d_s0_volume` assembles the existing fixed-pressure (`gamma=0`)
[displacement Hessian](ideal-xi.md) and physical kinetic inertia in a supplied
Galerkin space. `kin6d_s0_spectrum` consumes
[FortNum](https://github.com/lazy-fortran/fortnum) eigen iterations and
[FortSparse](https://github.com/lazy-fortran/fortsparse) native factors.
Reports say **candidate / bound unavailable**. For exp(-i omega t),

\[
K_{ij}=\int L_4(\nabla\xi_j,\nabla\xi_i^*)\,dV,\qquad
M_{ij}=\int\rho\,\xi_j\cdot\xi_i^*\,dV,\qquad Kx=\omega^2Mx.
\]

K is the Hessian, twice the complex modal delta-W quadratic form.
Inertia uses the same physical-volume and Fourier convention; their common
factor cancels. Negative omega-squared corresponds to gamma=sqrt(-omega-squared).
Density is physical input, rather than a Euclidean coefficient norm.

`s0_assemble_volume` takes connectivity `dof(nb,ne)`, physical weights
`weight(nq,ne)`, fields `b(3,nq,ne)`, `p(nq,ne),rho(nq,ne)`, basis values
`value(3,nb,nq,ne)` and physical gradients `gradient(3,3,nb,nq,ne)`.
All components use the same physical orthonormal frame. Cylindrical callers
can use `xi_gradient_n`, including its basis-rotation terms. DOF zero means a
homogeneous essential constraint; other indices are 1..ndof. Weights include
the physical Jacobian and integration period. Density, weights and permeability
must be positive and finite. Vacuum requires a separate operator.

One Hermitian element triangle is evaluated and reflected; quadratic diagonals
are explicitly real, as required by the native LDLH contract. Table arguments
are contiguous dummies, so any repacking occurs at entry rather than at each
quadrature contraction. The OpenMP
element loop writes disjoint triplets and allocates no memory. Storage/work
are O(ne nb^2) / O(ne nq nb^2), plus caller-owned basis tables. No dense global
matrix is formed. Returned full CSC matrices have the same pattern, including
stored zeros. The caller owns boundary conditions, gauge elimination and an
independent admissible basis. Mass has generated FP64/FP32 primal/JVP/VJP leaves;
stiffness uses KS3's leaves. The sparse/eigen ABI currently uses FP64 complex.

`s0_prepare(pencil,K,M,sigma,report)` requires full Hermitian matrices
with the same pattern and a positive-definite M. It retains factors of M and
K-sigma M. Sigma has units of omega-squared: pass omega_target**2 to target a
complex frequency. Real shifts use LDL^H and expose below/at/above counts with
FortSparse's reliability status. Sylvester's law relates negative pivots to
eigenvalues below sigma when M is positive definite. Floating-point counts
are not interval certificates. Singular shifts cannot be inverted. Complex
shifts use general supernodal LU and leave inertia unavailable (counts -1).

`s0_solve` accepts a caller-supplied start, tolerance and
`eigen_options_t` (EIGEN_NEAREST for gap targeting). Real/complex shifts use
Lanczos/Krylov-Schur respectively. Returned eigenvectors use the supplied DOF
basis and M normalization. The caller can export them with its geometry,
basis, inputs and revision provenance. Residuals are norms of the original
equation Kx-omega-squared Mx. Stopping uses independent Hermitian column-sum
bounds: ||r|| <= tol (||K||_1+|omega-squared| ||M||_1) ||x||.

`gen/gen_s0.f90` retains exact inertia/force and rigid-cylinder identities;
native FortSym zero tests reject UNKNOWN and probes. The wave regression uses
the exact P1 generalized sine spectrum, continuous shear-Alfven limit and
second-order eigenvalue convergence. Its discrete tolerance follows the
returned residual and M >= rho h/3 I, plus rounding. Continuum checks cover
the local omega_A^2=k_parallel^2 v_A^2 multiplication operator and a torsional
cylindrical annulus with Bz=1+r, p+Bz^2/(2mu)=constant. Its cross term is a
boundary term that vanishes under essential constraints. Cell spectra are
continuum approximations. The rigid long-wavelength uniform-current cylinder
oracle gamma^2=2q(1-q)/rho (unit radius, mu=1) does not supply missing boundary
terms to the volume assembler.

## Solved GS state to a fixed-wall m=1 kink

`kin6d_s0_gs_space` supplies a production radial-P1 displacement adapter for a
circular meridional disk revolved about the symmetry axis. `s0_gs_tables`
samples **the core state's field evaluator**, also used by D7, once in a batch;
it never substitutes analytic fields for the GS solution. It explicitly
reorders the core's `(R,Z,phi)` field to S0's `(R,phi,Z)` frame.
`s0_gs_assemble` passes these tables to the existing volume assembler.
Unsupported charts, a disk reaching R=0, n=0 and nonpositive density are
rejected. Status 1 means invalid input; 2 means a failed core sample.

The displacement is

\[
\xi=e^{i\theta+in\phi}
 [u(r)e_r+i v(r)e_\theta+i w(r)e_\phi],\qquad
R=R_0+r\cos\theta,\quad Z=r\sin\theta.
\tag{S2.1-1}
\]

Here `e_theta` lies in the meridional plane. Three continuous radial P1
amplitudes are independent away from the axis. At r=0 one combined basis
enforces **u(0)=v(0), w(0)=0**; its value is a regular transverse translation.
All three amplitudes vanish at the wall r=a. This is a conforming trial space
for a fixed conducting wall, with no plasma/vacuum interface, vacuum inertia
or missing external-kink surface term. The axis identification removes the
duplicate coefficient. There is no vector-potential gauge variable: physical
parallel displacements retain their kinetic mass and approach the gamma=0
neutral subspace. The mass matrix is positive definite; the zero-shift inertia
test includes every radial branch in this admitted m=1 space.

`gen/gen_s0_gs_space.f90` emits the physical values, R/Z derivatives and polar
volume Jacobian, with JVP/VJP leaves, into
`src/generated/kin6d_s0_gs_space_leaves.f90`. The adapter reuses
`xi_gradient_n` for the toroidal derivative and frame rotation, and FortFEM's
generated line basis for radial P1 values and derivatives. Element loops
use OpenMP and fixed-size local scratch. Four Gauss points per radial cell and
a periodic angular trapezoidal rule supply the quadrature. Retained table
storage is approximately `1200 * 4 * Ntheta * Nr` bytes; at `(32,512)` this is
75 MiB, plus setup and sparse-factor storage. The retained workstation target
is under 60 s per case and 384 MiB peak RSS; the Release focused run measured
249 MiB peak RSS and at most 2.16 s per new physics case with eight threads.
No optimized-throughput claim is made. The local geometry/basis derivatives
are verified; a composed GS-to-mode
parameter sensitivity is not yet exposed or claimed.

### Matched cylinder and independent oracle

The GS implementation solves a torus, not a straight periodic cylinder:
`initialize_gs` requires R>0 and its operator is Delta-star. The cylindrical
chart maps `(R,Z,phi)` into Cartesian space; it does not change that equation.
The synthetic family therefore solves the **forward GS equation from a zero
seed**, with boundary psi=0 on r=a=2 and constants

\[
\mu_0=\rho=1,\quad p'=-3/R_0,\quad FF'=R_0,\quad
F_{\rm edge}=R_0\sqrt{a^2+1/4}.
\tag{S2.1-2}
\]

The large-aspect-ratio cylinder limit is

\[
B_r=0,\quad B_\theta=r,\quad B_z=\sqrt{1/4+r^2},\quad
p=\tfrac32(a^2-r^2),\quad k=n/R_0=-0.4.
\tag{S2.1-3}
\]

Its radial force balance is exactly
`p_r + Btheta*(Btheta_r+Btheta/r) + Bz*Bz_r = 0`.
The perturbation pressure law is the existing **gamma_adiabatic=0** law.
With `tau_A=a*sqrt(mu0*rho)/Btheta(a)=1`, the independent reference gives
`omega^2=-0.1949752661` and `gamma*tau_A=0.4415600368`. This is a fixed-wall,
pressure-carrying ideal m=1 kink, not the uniform-current vacuum external kink.

`test/support/s0_cylinder_oracle.f90` independently assembles the Eulerian
force form, using Q=curl(xi cross B), the analytic current and pressure
gradient. It imports none of S0's Hessian, mass, geometry or basis leaves.
Its source is [Bernstein et al. (1958)](https://doi.org/10.1098/rspa.1958.0023),
Eqs. (2.24), (2.30), (3.8) and (3.13), with gravity absent and gamma=0:

\[
2\delta W=\int [|Q|^2/\mu_0-\xi^*\cdot(J\times Q)
                 +(\nabla\cdot\xi^*)(\xi\cdot\nabla p)]\,dV.
\tag{S2.1-4}
\]

It uses continuous P2 u and a discontinuous P2 amplitude along the unit
perpendicular direction `(Bz*e_theta-Btheta*e_z)/|B|`, eliminating the parallel
neutral motion by `w=-Btheta*v/Bz`, with the **physical** inertia of all three
components. This unit direction keeps the mass well scaled as Bz tends to
zero. A dense LAPACK generalized eigensolve is confined to this test
oracle. In the cylindrical energy form only u has radial derivatives;
tangential wall traces have no independent meaning in its closed energy
space. The production all-zero-wall H1 spaces approach this same normal-wall
problem through a tangential boundary layer. This explains the observed
first-order eigenvalue and half-order profile convergence, rather than an
assumed P1 second-order eigenvalue rate.

An analytic limit checks the oracle itself: setting Btheta=0 and constant Bz
gives `omega^2=k^2*Bz^2/(mu0*rho)=0.04` on divergence-free transverse
displacements. This follows directly from (S2.1-4); the P2 pair represents
`v=(r*u)'` exactly. A second, **unstable** analytic limit sets Btheta=r,
Bz -> 0 and |k| -> infinity with |k Bz| -> 0. For Bz=0, m=1 and gamma=0,
integration of (S2.1-4) by parts with u(a)=0 gives

\[
K=\int_0^a r[(krw-ru'+u)^2+w^2-3u^2]dr,\qquad
M=\int_0^a r(u^2+w^2)dr.
\tag{S2.1-5}
\]

Thus K+3M is nonnegative. Taking smooth u supported away from the axis and
wall, and `w=(r*u'-u)/(k*r)`, proves the lowest omega-squared tends to -3,
or **gamma -> sqrt(3)**. The same oracle checks k=-40,-80,-160 with Bz=1e-8:
omega-squared is -2.9908574, -2.9977096, -2.9994298, approaching this limit at
second order; the final relative growth error is below 2e-4. Bz is kept
positive to retain the regular-axis perpendicular frame, and |k Bz|<=1.6e-6
is a separate limit parameter. There is no claimed analytic growth value for
the finite-k reference (S2.1-3).

### Candidate acceptance gates

All cases use registered CTest tests, labels `stability;s0;mhd;fast`, and
timeouts below or equal to 60 s. On the Release CPU run with eight OpenMP
threads (plus a one-thread radial check agreeing at the printed precision):

| Gate | Measured result |
| --- | --- |
| Independent P2 oracle, 32/64/128 cells | omega-squared changes 3.3e-8, 2.0e-9; one negative eigenvalue |
| GS P3 forward residual | below 4e-13 on the retained ladder |
| Radial P1, 128/256/512 cells, R0=160 | omega-squared -0.187248486, -0.190908754, -0.192704057; observed order 1.028 |
| Physical profile distance, optimal complex phase | 0.12940, 0.09127, 0.06445; about half order |
| R0/a=10/20/40/80, Nr=512, n=-8/-16/-32/-64 | omega-squared -0.161902150, -0.185270761, -0.191208408, -0.192704057; orders 1.977, 1.989 |
| Original residual `norm(Kx-omega2*Mx)/(norm1(K)*norm2(x))` | below 1e-11 (measured below 1e-15) |
| Zero-shift LDLH inertia | one negative, zero zero, all remaining positive on every level |
| Finest relative gamma error against the cylinder | 0.005842 |

The measured Richardson gamma estimates, each padded by 25%, are below 0.006
for radial approximation and 0.0018 for aspect ratio. Separate GS refinement
and doubled angular quadrature gates each allow 2e-6; the independent radial
oracle allows 1e-7. Their sum, 0.0078041, is the final empirical acceptance
budget, with margin below one percent. These are refinement estimates and
acceptance thresholds, **not rigorous error bounds**. The state, trial-space,
aspect-ratio and eigen-residual contributions are kept separate.

This demonstrates the cylinder limit of a single-m admissible Ritz space.
It does not resolve toroidal coupling to other poloidal harmonics or count
instabilities in other m sectors. No device stability conclusion, free-boundary
threshold, continuous spectral certificate or eigenvalue sensitivity follows.
Reports retain **candidate / bound unavailable**.

Outstanding scope: general core-state displacement spaces and gauge
complements, vacuum/surface assembly, coupled-poloidal toroidal spectra,
toroidal GLISS/DCON comparison and continuous spectral certification.
FortSparse's complex CSC matvec allocates thread scratch above its parallel
threshold; a reusable-workspace action is needed before claiming an
allocation-free spectral path. OpenMP directives alone imply no scaling result.

Complete the global KS3 lift, residual-dependent chart terms and D6
boundary terms before admitting free-boundary thresholds or general toroidal
stability. AE1
can use the targeting/eigenvector contract when those physical operators and
their approximation-error gates are available.

## AE1 local toroidal continuum entry point

[V] `s0_toroidal_continuum_gap` tests the cold, large-aspect-ratio circular
field-line model

\[
-D_\parallel^2\Phi=\frac{\omega^2}{v_A^2}
 (1+2\epsilon\cos\theta)\Phi,\qquad
D_\parallel=R^{-1}(\partial_\zeta+q^{-1}\partial_\theta).
\]

Retaining adjacent harmonics exp[i(m theta - n zeta)] gives
K = rho v_A^2 diag(k_m^2,k_{m+1}^2),
M = rho [[1,epsilon],[epsilon,1]], with
k_m=(m/q-n)/R. At q=(m+1/2)/n their cylindrical frequencies cross at
omega0=v_A/(2qR). Toroidicity opens the local gap
omega_minus^2=omega0^2/(1+epsilon),
omega_plus^2=omega0^2/(1-epsilon).
The characteristic determinant is
(omega0^2-lambda)^2-epsilon^2 lambda^2.
`gen_s0` derives the blocks from the energy and inertia Hessians, proves the
angular Laurent constant, line bending, both edge roots and characteristic
polynomial exactly, and emits FP64/FP32 blocks with JVP/VJP.

[V] The test independently integrates the Fourier functional on a flat-q
annulus (q=3/2, m=1,2, n=1, R=10, r=0.5..1.5, rho=1.3,
v_A^2=2.1). P0 radial cells average epsilon=r/R exactly. Shifts inside the
gap target each edge separately; inertia must count precisely one lower
branch per cell. Discrete eigenvalue gates follow the original residual and
M >= rho h (1-r_outer/R) I. Refinement checks first-order convergence of the
sampled continuum support edges to the inner-surface analytic edges; it is
not convergence of an isolated eigenfunction. Frequency errors must be below
0.3%. An independent finite-difference check supplements JVP/VJP duality.

[V] This is a reduced local continuum model, not the complete toroidal
ideal-MHD operator. It neglects pressure/acoustic coupling, Shafranov shift,
other harmonics and second-order aspect-ratio terms; their model error is
unavailable. It supplies no radial envelope or discrete global TAE mode.
The leading-order physical gap width is O(epsilon omega0); higher-order
accuracy of these rational edge expressions is claimed only for the retained
model. Shared-energy GL2 toroidal qualification remains separate.

[V] `s0_external_kink_rigid_growth` projects the existing free-boundary
Schur energy onto the two degenerate real phases of the rigid helical
trial. It integrates rho*dV on the actual curved plasma mesh, rather than
using the free-boundary solver's coefficient metric, then solves that
physical-mass pencil with S0. The independent FortSym-checked cylinder
oracle is lambda=[q^2-1+f(ka)(q-1)^2]/rho,
f=K1(ka)/(K1(ka)+ka K0(ka)), ka=0.1. The existing seven exact kink
identities retain its equilibrium, vacuum/surface split and ka->0
Shafranov limit. Three toroidal aspect ratios approach the cylinder;
reported orders describe that limit at fixed mesh and ka, not spatial
convergence. This is a prescribed rigid-subspace growth estimate, not the
full external-kink eigenmode or an internal-kink result.
