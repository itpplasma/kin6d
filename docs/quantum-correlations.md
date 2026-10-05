# Finite-basis fermionic correlations

This helper supplies Wick exchange and beyond-Wick cumulants for the small
many-body examples. It does not evolve reduced density matrices or implement
BBGKY. The finite Fock parent retains the many-body state and its correlations.

For orthonormal spin orbitals, the shared convention is

\[
\gamma_{ij}=\langle a_j^\dagger a_i\rangle,\qquad
\Gamma_{ij,kl}=\langle a_k^\dagger a_l^\dagger a_j a_i\rangle.
\]

Gamma is antisymmetric in i,j and in k,l, with
\(\Gamma_{ij,kl}^*=\Gamma_{kl,ij}\). A normalized fixed-N state satisfies

\[
\operatorname{Tr}\gamma=N,\qquad
\sum_a\Gamma_{ia,ka}=(N-1)\gamma_{ik},\qquad
\operatorname{Tr}\Gamma=N(N-1),\qquad 0\leq\gamma\leq I.
\]

Here Tr Gamma means the full ordered-index trace
\(\sum_{ij}\Gamma_{ij,ij}\). For stored ordered pairs i<j, its value is twice
the pair-matrix trace. The anti-commutation relations establish these identities
for the exact fixed-N parent. Floating candidates and interval reconstructions
retain separate numerical-error contracts.

Wick exchange and the two-body cumulant are

\[
W_{ij,kl}=\gamma_{ik}\gamma_{jl}-\gamma_{il}\gamma_{jk},\qquad
K_{ij,kl}=\Gamma_{ij,kl}-W_{ij,kl}.
\]

The exchange term belongs to Wick even when K=0. Same-spin exchange holes are
therefore distinct from beyond-Slater correlation. These algebraic helpers
retain supplied Gamma; they never replace it with W. Spatial spin-density and
pair-density reconstructions require the declared orthonormal basis and their
own generated contraction/enclosure schedule. A finite-basis calculation alone
does not certify a continuum basis-truncation error.

A Slater determinant has idempotent gamma and K=0. A fractional-occupation Wick
state is generally quasifree/grandcanonical. If Tr gamma=N, then

\[
\sum_aW_{ia,ka}=N\gamma_{ik}-(\gamma^2)_{ik},\qquad
\operatorname{Tr}W=N^2-\operatorname{Tr}(\gamma^2).
\]

It meets the fixed-N partial contraction only when gamma is idempotent. An
ensemble with independent fractional occupations can be a valid grandcanonical
state without having a sharp particle number. The helper does not confuse that
ensemble with a fixed-N correlated state.

`kin6d_quantum_correlations` re-exports generated `quantum_wick_pair` and
`quantum_cumulant_pair`, each with `_jvp`, `_vjp` and `_interval` variants.
Coordinates are paired real/imaginary values. Supply x=gamma(i,k), y=gamma(j,l),
z=gamma(i,l), w=gamma(j,k); the cumulant additionally takes gr+i*gi=Gamma(i,j,k,l).
These scalar leaves also work directly on ordered-pair storage without expanding
Gamma into an O(M^4) tensor. Their independent-coordinate derivatives use the
real Euclidean pairing, extended to the real Frobenius pairing when accumulated
over arrays. Hermitian/antisymmetric tangent constraints are a caller choice.
The polynomial maps are smooth everywhere. Generated strict leaves enclose
exact algebraic values for the supplied interval inputs, inheriting FortNum's
outward arithmetic contract. They do not certify an underlying state, time
propagator, model closure error or sensitivity trajectory. Leaf callers must
validate input domains and reject nonfinite/empty output enclosures. Scalar
leaves and `quantum_pauli` allocate no workspace.

`quantum_pauli(gamma_real,gamma_imag,status,ok)` checks the exact supplied binary64
matrix. Non-Hermitian data or a real diagonal outside [0,1] violate necessary
conditions. For Hermitian inputs, it uses outward Gershgorin row bounds with
\(|z|\leq|\Re z|+|\Im z|\). If every row lies in [0,1], status is certified;
otherwise it is undetermined. Exact diagonal rows permit endpoint occupations
0 and 1. NaN, bad shapes, or nonfinite arithmetic fail with ok=false. This is a
sufficient finite-matrix Pauli test, not an eigenvalue solver or full
N-representability certificate. It does not establish positivity of an uncertain
physical RDM merely because its floating candidate passes. Fermionic operator
positivity is also distinct from positivity of a Wigner quasidistribution.

Independent tests build CAR operator matrices from occupation-bit parity in the
N=2, M=4 sector. Real and complex Slater states agree with Wick. The coherent
(|12>+|34>)/sqrt(2) state has gamma=.5I and a nonzero pair-coherence cumulant;
its incoherent mixture shares gamma and loses that coherence. Tests verify the
fixed-N trace/contraction and the failure of fractional Wick fixed-N identities.
Separate native complex-product derivatives, centered perturbations and
real-pairing identities check primal/JVP/VJP behavior. Binary128 exact dyadic
corner oracles check zero-crossing interval containment. Pauli tests include
valid diagonal and complex matrices, a rotated Slater with undetermined
Gershgorin status, necessary violations, malformed inputs and overflow.

Provenance: generator `gen/gen_quantum_correlations.f90`, FortSym
`a6142712af998b727c245fb6cfbb3745b296ec24`, FortNum
`1b1349948a6250ed5a38655d9b8042c894370d42`. Private checked GNU Fortran tests pass;
the controller owns registered CMake/CTest integration and release evidence.
The earlier unintegrated first-level RDM/BBGKY packet remains deferred outside
Git after the user's narrower milestone instruction.
