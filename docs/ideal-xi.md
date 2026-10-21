# Displacement energy: first D5 increment

[V] `kin6d_ideal_xi` evaluates the plasma-volume second variation of the
map lane's existing fixed-pressure (`gamma=0`) frozen-in energy. It accepts
physical displacement gradients, without constructing a global matrix.
It does not yet supply a global stability eigenproblem or a continuous bound.

[V] Both map and xi generators use `kin6d_symbolic_ideal:frozen_energy`:

\[
w(F)=\frac{|F B|^2}{2\mu_0\det F}-p\det F,\qquad \det F>0.
\]

Here the reference coordinates are physical Cartesian coordinates, so the
reference map is the identity, B is the reference magnetic field, and p is
frozen pressure. This is the map lane's pressure constraint; an adiabatic
compressional energy requires a separate, derived model increment.

[V] For A = grad xi and C = grad eta, the local Hessian is

\[
\begin{split}
L_4(A,C)={}&\frac{(AB)\cdot(CB)}{\mu_0}
-\frac{(B\cdot AB)\operatorname{tr}C
       +(B\cdot CB)\operatorname{tr}A}{\mu_0}\\
&+\left(\frac{B^2}{2\mu_0}-p\right)
       \operatorname{tr}A\operatorname{tr}C
+\left(\frac{B^2}{2\mu_0}+p\right)\operatorname{tr}(AC).
\end{split}
\]

[V] `gen/gen_ideal_xi.f90` differentiates the retained energy twice, sets
F=I, and proves all 81 coefficients against this independently written xi
formula using FortSym's exact native zero test. A further exact test checks
the contraction. UNKNOWN and probe evidence fail. Production leaves are
emitted from the differentiated energy, together with JVP/VJP in B, p, mu0
and both gradients, in FP64 and FP32. FP64 is the default consumer kind;
the two schedules share one implementation include.

[V] `xi_gradient_n(n,R,xi,dr,dz,A,status)` constructs the physical gradient
for amplitudes proportional to exp(+i n phi), with local orthonormal
components ordered **(R, phi, Z)** and gradient columns **(R, phi, Z)**.
The toroidal column is

\[
A_{:,\phi}=R^{-1}(in\xi_R-\xi_\phi,
                         in\xi_\phi+\xi_R,
                         in\xi_Z)^T.
\]

[V] The two basis-rotation terms are essential. The public constructor
rejects R<=0 and nonfinite inputs or jets; regularity at the cylindrical
axis is not implemented here.
These are components in a local orthonormal frame of a physical tensor;
the magnetic field order differs from the core evaluator's (R,Z,phi).
A future core sampler must reorder it explicitly.

[V] `xi_bilinear(B,p,mu0,A,C)` extends L4 sesquilinearly, conjugating C.
`xi_volume_energy` reduces caller-supplied physical volume weights times
L4(A,A)/2 with OpenMP, fixed-size stack scratch, and no allocation in its
quadrature loop. For one complex Fourier amplitude its output is

\[
\delta W_{n,\mathrm{complex}}=\tfrac12\int L_4(A,A^*)\,dV.
\]

[V] Weights include the toroidal integration period. For the physical
perturbation Re[xi_n exp(i n phi)] with nonzero n, the real-mode energy is
one half of that complex-amplitude result. For n=0 with real coefficients,
the result is already the physical energy. Callers must use the same
convention for inertia. The volume consumer rejects negative weights,
nonpositive mu0 and nonfinite inputs. Reports always say `candidate`, bound `unavailable`;
no equilibrium residual, quadrature or discretization bound is implied.

[V] The batch work is O(Nq), with fixed-size scratch per thread and no
global Hessian storage. This is the performance/memory contract for the
first correctness increment; no optimized throughput is claimed.

[V] The native test compares all coefficients with the existing emitted map
Hessian, checks an independent tensor contraction, Hermitian pairing,
JVP/VJP duality and central differences, both kinds, Cartesian finite
differences of the cylindrical frame, and the n=1 translation null jet.
The screw-pinch rigid helical translation gives the raw map plasma Hessian
q^2 (mu0=1, q=k Bz), matching the plasma piece used by the existing external
kink test. An axial shear-Alfven jet checks k^2 Bz^2/mu0 and the factor 1/2.
These analytic limits take priority over future code comparisons.

[H] Next increment: evaluate these jets on the shared axisymmetric core
state and actual per-n FE trial/test basis, transform reference derivatives
into physical derivatives, and compare the assembled action with the map
Hessian pulled back to the same displacement basis. Retain the Fourier
pairing and residual-dependent chart acceleration terms explicitly. Then
handle the gauge complement and the n=0 prescribed-q Schur reduction;
include D6 vacuum/surface terms before making an external-kink threshold
claim. The existing strip-current response operator is not yet identified
with this displacement operator. Rigorous spectral and state bounds remain
unavailable until their stability and discretization premises are checked.
