# Fermionic finite-basis microscopic reference

The reference stores the complete fixed-N fermionic state in a supplied set of
M spin orbitals. Neither N=2 nor a two-level sector is built into the engine.
It supplies many-body Hamiltonian actions and exact CAR scheduling, together
with a genuine periodic 3D jellium model. It applies no correlation closure.
The current coefficient specialization is a real Hermitian Hamiltonian; states
and density matrices use separate real and imaginary arrays.

## Basis and Hamiltonian convention

Orbital 1 is the least significant occupancy bit. Determinants are

    |S> = a_i1^dagger ... a_iN^dagger |0>, i1 < ... < iN.

Masks are sorted numerically. Creation and annihilation phases are exactly
(-1) to the number of occupied lower orbitals. Operator words are supplied in
written left-to-right order and applied right-to-left. Repeated creation of an
occupied orbital and annihilation of an empty orbital give zero.

The Hamiltonian is

    H = sum_ij h_ij a_i^dagger a_j
        + (1/4) sum_ijkl v_ij,kl a_i^dagger a_j^dagger a_l a_k.

Here h is real symmetric, and v is antisymmetric separately in ij and kl, with
v_ij,kl=v_kl,ij. The input stores only the ordered-pair matrix i<j,k<l in
lexicographic pair order. Its Hamiltonian coefficient is one: the four
antisymmetric permutations exactly cancel the displayed 1/4. This is a
representation of the general interaction tensor, not an interaction cutoff.

The integer-mask reference backend accepts 1<=M<=62 and 0<=N<=M. This is a
backend limit, not a physical limit on future orbital/RDM implementations.
Combinatorial dimension and explicitly supplied state/payload byte budgets
are checked before enumeration and allocation. Integer count intermediates
are bounded before the next multiplication. The basis payload is 8 D bytes
plus the two integer pair-index arrays. Hamiltonian payload is 24(M^2+P^2)
bytes for binary64 candidates and two-endpoint interval packets, P=M(M-1)/2.
These payload counts exclude compiler object metadata and caller-owned arrays.
A dense many-body output is optional and caller-owned; its memory is additional.

## Periodic 3D jellium

A spin orbital is V^(-1/2) exp(i k dot x) times an up/down spin state in a cubic
periodic cell, V=L^3 and k=2 pi n/L. The constructor accepts arbitrary distinct
integer 3-vectors n and spin labels +1/-1. Duplicate momentum with different
spin is permitted. L, hbar and mass must be positive; Coulomb coupling g>=0.
Integer mode components are bounded by 10^6 so every squared transfer norm is
exactly representable as binary64 and integer products cannot overflow.

The diagonal kinetic energy is hbar^2 |k|^2/(2m). The unscreened interaction
Fourier coefficient is 4 pi g/(V |q|^2) for q!=0; the q=0 mode is removed by the
uniform neutralizing background. Momentum is conserved exactly, and the
antisymmetrized pair matrix is direct minus exchange, with their respective
spin Kronecker deltas. No Yukawa force, Debye length, softening, or 1D particle
line is inserted. The Hamiltonian energy reference excludes an additional
additive background constant. The plane-wave Hamiltonian and neutral-background
zero-mode convention are described in the
[Georgetown jellium lecture](https://site.physics.georgetown.edu/~jkf/grad_quant2/html/lecture37.html).

The model is the Galerkin projection onto the supplied finite orbital list.
Strict interval leaves enclose its exact physical coefficients, interpreting
binary64 L,hbar,m,g as exact real parameters and pi as the mathematical
constant. Candidate coefficient arithmetic and these physical coefficient
packets are distinct. Removing modes changes the model; cutoff comparisons
are experiments unless a separate remainder theorem is supplied.

Finite-basis matrix elements/actions are certified-discrete when their strict
queries succeed. Continuum momentum-cutoff, thermodynamic-limit, omitted
radiation/relativity, physical-data and model-reduction bounds are unavailable.
Time evolution and its finite-basis error control belong to the separate
quantum-evolution consumer; this module does not infer trajectory bands from
energy or norm conservation.

## Actions, blocks and correlations

`fock_matvec` applies H to a complex vector without assembling the many-body
matrix. A shared CAR transition iterator scans occupied annihilation indices.
No per-action allocation, global kinetic Jacobian or derivative tape is used.
The reference stores the pair coefficient matrix; it is not a scalable
production electronic-structure solver or a full four-index tensor backend.

A selected determinant block can be assembled without a full D-by-D matrix.
Its `invariant` flag is true only when every outgoing term either lands within
the selected set or has an exact-zero coefficient packet. This is a sufficient
structural zero-leakage proof; cancellation-dependent invariant sectors may
remain unproved. Distinct valid indices are required. A returned submatrix
with `invariant=false` is a projection and cannot claim closed parent dynamics.
The total-momentum-zero, one-up/one-down jellium sector therefore has only
7/19/27 states for the corresponding symmetric spatial-mode lists, while the
general parent basis remains available.

The direct `assemble_jellium_pair_sector` constructor accepts an arbitrary unique
3D integer mode list of length K and caller-preallocated candidate/strict K-by-K
matrices. Its physical basis is |n up,-n down>, so no canonical determinant sign
is applied. It composes the same generated kinetic and Coulomb leaves: diagonal
entries are two kinetic energies, and off-diagonal entries are the nonzero
momentum-transfer Coulomb coefficient. Opposite spins remove exchange; q=0
remains absent. This is exactly the total-momentum-zero, one-up/one-down sector
of the parent Hamiltonian on those up modes and their opposite down modes.
Full-parent block tests on 7/19/27 spatial-mode lists verify the schedule.

This constructor allocates no full determinant basis or orbital-pair tensor and
has no 62-bit reference-backend restriction. Its explicit output payload budget
is 48 K^2 bytes for two binary64 matrices and two two-endpoint interval matrices.
The optional `certify=false` benchmark path computes identical candidate matrices
and returns entire interval packets with an explicit unavailable status, avoiding
strict arithmetic. Strict construction is the default.
The byte check uses overflow-safe division before numerical work; unique modes,
finite positive L/hbar/mass, nonnegative coupling and bounded integer momenta
are required. Successful strict packets certify the same finite-basis physical
matrix, not a continuum cutoff. The neutral-background energy convention does
not imply a positive-semidefinite Hamiltonian. Composed parameter derivatives
for this constructor are unavailable; the common generated physics leaves
retain primal/JVP/VJP formulas.

A canonical determinant may differ by a sign from a physically ordered
|k up,-k down> pair. Relative-coordinate Fourier amplitudes must include that
orbital-order sign. `fock_pair_action` instead computes A_second A_first psi
for arbitrary complex weighted annihilators and returns amplitudes in the
same-M, fixed-(N-2) residual basis. It iterates occupied annihilators and skips
only exact zero candidate amplitudes; its interval version skips only exact
zero interval packets. This permits direct pair-density queries without
constructing a full two-body density matrix. Its scalar numerical
triple-product accumulation has generated primal/JVP/VJP/interval leaves.

For arbitrary complex psi, `fock_rdms` returns

    gamma_ij = <a_j^dagger a_i>,
    Gamma_ij,kl = <a_k^dagger a_l^dagger a_j a_i>.

Gamma is again stored in ordered-pair form. Full ordered-tensor trace is twice
the pair-matrix trace. No normalization is silently imposed: Tr gamma=N||psi||^2
and Tr Gamma=N(N-1)||psi||^2. For a normalized fixed-N state,
sum_a Gamma_ia,ka=(N-1)gamma_ik. A Slater determinant obeys the Wick identity;
a general correlated state retains its complete two-body correlations.
`fock_operator_expectation` computes Tr(rho O) for a complete density matrix,
and `fock_density_rhs` computes -i[H,rho]/hbar without discarding any entries.
These candidate actions do not assert positivity or N-representability of
arbitrary supplied density data. Density normalization belongs to the caller.

## Observable and derivative status

`fock_rdms_interval` encloses strict arithmetic at the exact binary64 candidate
psi plus an optional, caller-proved complex Hilbert-norm state radius epsilon.
Each CAR word has operator norm at most one, so

    |<psi,O psi>-<c,O c>| <= (2||c||+epsilon) epsilon.

The candidate norm and enlargement are evaluated outward. This supplies each
real and imaginary RDM component separately; it is not a sum-free bound for
an arbitrary subsequent contraction. Passing epsilon=0 gives candidate-only
strict contraction intervals. A downstream density/pair observable adds its
own operator-norm amplification and state error once. For A_second A_first,
its squared-norm observable has operator norm at most
||first_weights||_2^2 ||second_weights||_2^2.

The strict pair action also accepts interval coefficient and state coordinates.
It encloses that coordinate box, and does not itself claim a Hilbert error
radius or normalization theorem. Other complete-density actions and their
sensitivities are candidate/unavailable until a requested checker is supplied.

Scalar continuous physics and accumulation leaves have retained FortSym
primal/JVP/VJP definitions. Composed matrix-free JVP/VJP actions are implemented
for symmetric real h/v and complex state coordinates. The pairing is real
Euclidean for real/imaginary state components and full Frobenius for symmetric
coefficient matrices; the VJP projects onto the symmetric matrix space.
Discrete mode/spin/basis selection is fixed under differentiation. The generic
pair-query composed derivative interface is unavailable; its smooth numerical
leaf derivatives are retained. No continuous physical trajectory adjoint or
certified sensitivity is inferred from these discrete actions.

## Independent evidence and performance target

The registered test checks exhaustive CAR identities, an independent ordered-list
insertion/deletion oracle for the full four-index Hamiltonian including 1/4,
noninteracting vacuum/one/three/full sectors, complex Slater Wick identities,
correlated-state RDMs and fixed-N contraction, strict real128 containment,
directional JVP finite differences and VJP pairing, full-density commutator
trace, genuine 3D jellium momentum/spin selection, budget failures, generic
weighted pair annihilation/Pauli exclusion, and selected-sector physical
coefficients with mixed canonical signs and zero-leakage status.

The immediate workload target is the budgeted 7/19/27-state invariant-sector
teaching calculation. Larger full-basis actions are a verification reference;
the target is measured total solve-plus-requested-observable-bound cost, not a
speed ranking. Every consumer reports its actual state/observable error status
and distinguishes finite-basis numerical control from unavailable cutoff error.
