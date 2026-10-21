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

One Hermitian element triangle is evaluated and reflected. The OpenMP
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

Limits: no core-state/FE displacement adapter, general gauge complement,
vacuum/surface assembly, internal/external kink growth-rate convergence,
toroidal GLISS/DCON comparison or continuous spectral certificate is supplied.
FortSparse's complex CSC matvec allocates thread scratch above its parallel
threshold; a reusable-workspace action is needed before claiming an
allocation-free spectral path. OpenMP directives alone imply no scaling result.

Complete the global KS3 lift, residual-dependent chart terms and D6
boundary terms before admitting kink thresholds or toroidal stability. AE1
can use the targeting/eigenvector contract when those physical operators and
their approximation-error gates are available.

## AE1 local toroidal continuum entry point

`s0_toroidal_continuum_gap` tests the cold, large-aspect-ratio circular
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

The test independently integrates the Fourier functional on a flat-q
annulus (q=3/2, m=1,2, n=1, R=10, r=0.5..1.5, rho=1.3,
v_A^2=2.1). P0 radial cells average epsilon=r/R exactly. Shifts inside the
gap target each edge separately; inertia must count precisely one lower
branch per cell. Discrete eigenvalue gates follow the original residual and
M >= rho h (1-r_outer/R) I. Refinement checks first-order convergence of the
sampled continuum support edges to the inner-surface analytic edges; it is
not convergence of an isolated eigenfunction. Frequency errors must be below
0.3%. An independent finite-difference check supplements JVP/VJP duality.

This is a reduced local continuum model, not the complete toroidal
ideal-MHD operator. It neglects pressure/acoustic coupling, Shafranov shift,
other harmonics and second-order aspect-ratio terms; their model error is
unavailable. It supplies no radial envelope or discrete global TAE mode.
The leading-order physical gap width is O(epsilon omega0); higher-order
accuracy of these rational edge expressions is claimed only for the retained
model. Shared-energy GL2 toroidal qualification remains separate.

`s0_external_kink_rigid_growth` projects the existing free-boundary
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
