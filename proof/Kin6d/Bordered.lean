import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Tactic.Module
import Mathlib.Algebra.Order.Star.Real

set_option autoImplicit false

/-!
# Bordered saddle-point matrices (FB-D1: Lemmas B1 and B2)

Discrete linear-algebra lemmas behind the gauged free-boundary Newton system
(`theory/derivations/free-boundary-gauge-vertical.md`, §4). Let `H : n × n`,
`M : m × n`, `K = [H Mᵀ; M 0]`, and let the columns of `Z : n × k` form a basis
of `ker M` (`IsKernelBasis`).

* `lemma_B1`: `K` is nonsingular iff `rank M = card m` and `Zᵀ H Z` is
  nonsingular (no symmetry needed).
* `lemma_B2`: for symmetric `H` and `rank M = card m = r`,
  `In(K) = In(Zᵀ H Z) + (r, r, 0)`.

Inertia is defined variationally (`posIndex`: largest dimension of a subspace on
which the quadratic form is positive definite; `negIndex A = posIndex (-A)`;
`zeroIndex A = card - rank A`). For symmetric matrices these equal the numbers of
positive, negative and zero eigenvalues (`posIndex_eq_card_pos`,
`negIndex_eq_card_neg`, `zeroIndex_eq_card_zero`), so `lemma_B2` is the
eigenvalue-count statement of the derivation. The proof of `lemma_B2` uses an
explicit congruence `T = [Z G; -GᵀHZ, σI - ½GᵀHG]` with `MG = I`, for which
`Tᵀ K T = diag(Zᵀ H Z, 2σ I)`, and the uniqueness half of Sylvester's law
(`posIndex_eq_of_split`).
-/

namespace Kin6d.Bordered

open Matrix

variable {n m k : Type*} [Fintype n] [Fintype m] [Fintype k]
  [DecidableEq n] [DecidableEq m] [DecidableEq k]

/-- The bordered (KKT) matrix `[H Mᵀ; M 0]`. -/
def bordered (H : Matrix n n ℝ) (M : Matrix m n ℝ) : Matrix (n ⊕ m) (n ⊕ m) ℝ :=
  fromBlocks H Mᵀ M 0

/-- The columns of `Z` form a basis of `ker M`. -/
structure IsKernelBasis (M : Matrix m n ℝ) (Z : Matrix n k ℝ) : Prop where
  mul_eq_zero : M * Z = 0
  injective : Function.Injective Z.mulVec
  spans : ∀ x, M *ᵥ x = 0 → ∃ z, Z *ᵥ z = x

omit [DecidableEq n] [DecidableEq m] [DecidableEq k] in
theorem mulVec_injective_iff {p q : Type*} [Fintype q] (A : Matrix p q ℝ) :
    Function.Injective A.mulVec ↔ ∀ v, A *ᵥ v = 0 → v = 0 := by
  rw [← Matrix.coe_mulVecLin]
  exact injective_iff_map_eq_zero _

omit [DecidableEq n] [DecidableEq m] [DecidableEq k] in
/-- Full row rank is injectivity of `Mᵀ`. -/
theorem rank_eq_card_iff_injective (M : Matrix m n ℝ) :
    M.rank = Fintype.card m ↔ Function.Injective Mᵀ.mulVec := by
  have h := LinearMap.finrank_range_add_finrank_ker Mᵀ.mulVecLin
  rw [Module.finrank_fintype_fun_eq_card] at h
  have hr : M.rank = Module.finrank ℝ (LinearMap.range Mᵀ.mulVecLin) := by
    rw [← rank_transpose M]; rfl
  rw [hr, ← Matrix.coe_mulVecLin, ← LinearMap.ker_eq_bot, ← Submodule.finrank_eq_zero]
  omega

omit [DecidableEq n] [DecidableEq m] [DecidableEq k] in
/-- `range Mᵀ = (ker M)^⊥ = ker Zᵀ` when `Z` is a kernel basis. -/
theorem exists_transpose_of_orth {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    (hZ : IsKernelBasis M Z) (w : n → ℝ) (hw : Zᵀ *ᵥ w = 0) : ∃ y, Mᵀ *ᵥ y = w := by
  classical
  -- dimension count: finrank ker M = card k
  have hker : LinearMap.ker M.mulVecLin = LinearMap.range Z.mulVecLin := by
    ext x
    simp only [LinearMap.mem_ker, LinearMap.mem_range, Matrix.mulVecLin_apply]
    constructor
    · exact hZ.spans x
    · rintro ⟨z, rfl⟩
      rw [Matrix.mulVec_mulVec, hZ.mul_eq_zero, Matrix.zero_mulVec]
  have hZinj : LinearMap.ker Z.mulVecLin = ⊥ := by
    rw [LinearMap.ker_eq_bot]; exact hZ.injective
  have hk : Module.finrank ℝ (LinearMap.range Z.mulVecLin) = Fintype.card k := by
    rw [LinearMap.finrank_range_of_inj hZ.injective, Module.finrank_fintype_fun_eq_card]
  have h1 := LinearMap.finrank_range_add_finrank_ker M.mulVecLin
  have h2 := LinearMap.finrank_range_add_finrank_ker Zᵀ.mulVecLin
  rw [Module.finrank_fintype_fun_eq_card, hker, hk] at h1
  rw [Module.finrank_fintype_fun_eq_card] at h2
  have hrM : Module.finrank ℝ (LinearMap.range Mᵀ.mulVecLin) =
      Module.finrank ℝ (LinearMap.range M.mulVecLin) := rank_transpose M
  have hrZ : Module.finrank ℝ (LinearMap.range Zᵀ.mulVecLin) =
      Module.finrank ℝ (LinearMap.range Z.mulVecLin) := rank_transpose Z
  have hle : LinearMap.range Mᵀ.mulVecLin ≤ LinearMap.ker Zᵀ.mulVecLin := by
    rintro _ ⟨y, rfl⟩
    simp only [LinearMap.mem_ker, Matrix.mulVecLin_apply, Matrix.mulVec_mulVec]
    rw [← Matrix.transpose_mul, hZ.mul_eq_zero, Matrix.transpose_zero, Matrix.zero_mulVec]
  have heq : LinearMap.range Mᵀ.mulVecLin = LinearMap.ker Zᵀ.mulVecLin :=
    Submodule.eq_of_le_of_finrank_eq hle (by omega)
  have : w ∈ LinearMap.range Mᵀ.mulVecLin := by
    rw [heq]; exact hw
  obtain ⟨y, hy⟩ := this
  exact ⟨y, hy⟩

/-- **Lemma B1.** `K = [H Mᵀ; M 0]` is nonsingular iff `rank M = card m` and
`Zᵀ H Z` is nonsingular, for any basis `Z` of `ker M`. `H` need not be symmetric. -/
theorem lemma_B1 (H : Matrix n n ℝ) (M : Matrix m n ℝ) {Z : Matrix n k ℝ}
    (hZ : IsKernelBasis M Z) :
    IsUnit (bordered H M).det ↔ M.rank = Fintype.card m ∧ IsUnit (Zᵀ * H * Z).det := by
  rw [← Matrix.isUnit_iff_isUnit_det, ← Matrix.isUnit_iff_isUnit_det,
    ← Matrix.mulVec_injective_iff_isUnit, ← Matrix.mulVec_injective_iff_isUnit,
    rank_eq_card_iff_injective]
  have hMZ : ∀ z, M *ᵥ (Z *ᵥ z) = 0 := fun z => by
    rw [Matrix.mulVec_mulVec, hZ.mul_eq_zero, Matrix.zero_mulVec]
  have hZM : ∀ y, Zᵀ *ᵥ (Mᵀ *ᵥ y) = 0 := fun y => by
    rw [Matrix.mulVec_mulVec, ← Matrix.transpose_mul, hZ.mul_eq_zero,
      Matrix.transpose_zero, Matrix.zero_mulVec]
  -- `K (x, y) = (H x + Mᵀ y, M x)`
  have hK : ∀ v : n ⊕ m → ℝ, bordered H M *ᵥ v =
      Sum.elim (H *ᵥ (v ∘ Sum.inl) + Mᵀ *ᵥ (v ∘ Sum.inr)) (M *ᵥ (v ∘ Sum.inl)) := by
    intro v
    rw [bordered, fromBlocks_mulVec, Matrix.zero_mulVec, add_zero]
  simp only [mulVec_injective_iff]
  constructor
  · intro hKinj
    constructor
    · intro y hy
      have := hKinj (Sum.elim 0 y) (by
        rw [hK]; ext (i | i) <;> simp [hy])
      simpa using congrArg (fun v => v ∘ Sum.inr) this
    · intro z hz
      have hz' : Zᵀ *ᵥ (H *ᵥ (Z *ᵥ z)) = 0 := by
        simpa [Matrix.mulVec_mulVec, Matrix.mul_assoc] using hz
      obtain ⟨y, hy⟩ := exists_transpose_of_orth hZ _ (by
        rw [Matrix.mulVec_neg, hz', neg_zero] : Zᵀ *ᵥ (-(H *ᵥ (Z *ᵥ z))) = 0)
      have h0 := hKinj (Sum.elim (Z *ᵥ z) y) (by
        rw [hK]; ext (i | i)
        · simp [hy]
        · simp [hMZ])
      have hx : Z *ᵥ z = 0 := by simpa using congrArg (fun v => v ∘ Sum.inl) h0
      exact (mulVec_injective_iff Z).1 hZ.injective z hx
  · rintro ⟨hMinj, hSinj⟩ v hv
    rw [hK] at hv
    have hv1 : H *ᵥ (v ∘ Sum.inl) + Mᵀ *ᵥ (v ∘ Sum.inr) = 0 := by
      simpa using congrArg (fun w => w ∘ Sum.inl) hv
    have hv2 : M *ᵥ (v ∘ Sum.inl) = 0 := by
      simpa using congrArg (fun w => w ∘ Sum.inr) hv
    obtain ⟨z, hz⟩ := hZ.spans _ hv2
    have hz0 : z = 0 := by
      apply hSinj z
      have := congrArg (fun w => Zᵀ *ᵥ w) hv1
      simp only [Matrix.mulVec_add, hZM, add_zero, Matrix.mulVec_zero] at this
      rw [← hz] at this
      simpa [Matrix.mulVec_mulVec, Matrix.mul_assoc] using this
    have hx : v ∘ Sum.inl = 0 := by rw [← hz, hz0, Matrix.mulVec_zero]
    rw [hx, Matrix.mulVec_zero, zero_add] at hv1
    have hy : v ∘ Sum.inr = 0 := hMinj _ hv1
    ext (i | i)
    · exact congrFun hx i
    · exact congrFun hy i

/-! ## Inertia (Lemma B2)

`posIndex A` is the largest dimension of a subspace on which the quadratic form
`x ↦ xᵀ A x` is positive definite, `negIndex A = posIndex (-A)` and
`zeroIndex A = card - rank A`. For symmetric `A` these are the numbers of positive,
negative and zero eigenvalues (`posIndex_eq_card_pos`, `negIndex_eq_card_neg`,
`zeroIndex_eq_card_zero`). -/

section Inertia

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- `xᵀ A x > 0` for every nonzero `x ∈ V`. -/
def IsPosOn (A : Matrix ι ι ℝ) (V : Submodule ℝ (ι → ℝ)) : Prop :=
  ∀ x ∈ V, x ≠ 0 → 0 < x ⬝ᵥ (A *ᵥ x)

/-- `xᵀ A x ≤ 0` for every `x ∈ V`. -/
def IsNonposOn (A : Matrix ι ι ℝ) (V : Submodule ℝ (ι → ℝ)) : Prop :=
  ∀ x ∈ V, x ⬝ᵥ (A *ᵥ x) ≤ 0

/-- Positive index of inertia (variational definition). -/
noncomputable def posIndex (A : Matrix ι ι ℝ) : ℕ :=
  sSup {d | ∃ V : Submodule ℝ (ι → ℝ), IsPosOn A V ∧ Module.finrank ℝ V = d}

/-- Negative index of inertia. -/
noncomputable def negIndex (A : Matrix ι ι ℝ) : ℕ := posIndex (-A)

/-- Zero index of inertia (nullity). -/
noncomputable def zeroIndex (A : Matrix ι ι ℝ) : ℕ := Fintype.card ι - A.rank

omit [DecidableEq ι] in
theorem finrank_le_posIndex {A : Matrix ι ι ℝ} {V : Submodule ℝ (ι → ℝ)}
    (h : IsPosOn A V) : Module.finrank ℝ V ≤ posIndex A :=
  le_csSup ⟨Fintype.card ι, by
    rintro _ ⟨W, -, rfl⟩
    exact (Submodule.finrank_le W).trans_eq (Module.finrank_fintype_fun_eq_card ℝ)⟩
    ⟨V, h, rfl⟩

omit [DecidableEq ι] in
theorem finrank_add_finrank_le {A : Matrix ι ι ℝ} {P N : Submodule ℝ (ι → ℝ)}
    (hP : IsPosOn A P) (hN : IsNonposOn A N) :
    Module.finrank ℝ P + Module.finrank ℝ N ≤ Fintype.card ι := by
  have hinf : P ⊓ N = ⊥ := by
    rw [eq_bot_iff]
    intro x hx
    rw [Submodule.mem_bot]
    by_contra h0
    exact absurd (hN x hx.2) (not_le.2 (hP x hx.1 h0))
  have h1 := Submodule.finrank_sup_add_finrank_inf_eq P N
  rw [hinf, finrank_bot, add_zero] at h1
  have h2 := Submodule.finrank_le (P ⊔ N)
  rw [Module.finrank_fintype_fun_eq_card] at h2
  omega

omit [DecidableEq ι] in
/-- A complementary pair (positive definite, negative semidefinite) determines the
positive index; this is the uniqueness half of Sylvester's law of inertia. -/
theorem posIndex_eq_of_split {A : Matrix ι ι ℝ} {P N : Submodule ℝ (ι → ℝ)}
    (hP : IsPosOn A P) (hN : IsNonposOn A N)
    (hd : Module.finrank ℝ P + Module.finrank ℝ N = Fintype.card ι) :
    posIndex A = Module.finrank ℝ P := by
  apply le_antisymm
  · unfold posIndex
    refine csSup_le ⟨0, ⊥, fun x hx h => absurd ((Submodule.mem_bot ℝ).1 hx) h,
      finrank_bot ℝ _⟩ ?_
    rintro _ ⟨V, hV, rfl⟩
    have := finrank_add_finrank_le hV hN
    omega
  · exact finrank_le_posIndex hP

omit [DecidableEq ι] in
theorem quad_mulVec {κ : Type*} [Fintype κ] (A : Matrix ι ι ℝ) (T : Matrix ι κ ℝ)
    (u : κ → ℝ) : (T *ᵥ u) ⬝ᵥ (A *ᵥ (T *ᵥ u)) = u ⬝ᵥ ((Tᵀ * A * T) *ᵥ u) := by
  rw [Matrix.mulVec_mulVec, dotProduct_mulVec, dotProduct_mulVec, ← vecMul_transpose,
    vecMul_vecMul, Matrix.mul_assoc]

omit [DecidableEq ι] in
theorem injective_of_transpose_mul {κ : Type*} [Fintype κ] [DecidableEq κ]
    {B : Matrix ι κ ℝ} (h : Bᵀ * B = 1) : Function.Injective B.mulVec := by
  rw [mulVec_injective_iff]
  intro u hu
  have := congrArg (fun v => Bᵀ *ᵥ v) hu
  simpa [Matrix.mulVec_mulVec, h] using this

/-- The image under `U` of the coordinate subspace supported on `{i | s i}`. -/
noncomputable def coordSub (U : Matrix ι ι ℝ) (s : ι → Prop) [DecidablePred s] :
    Submodule ℝ (ι → ℝ) :=
  LinearMap.range (U * (1 : Matrix ι ι ℝ).submatrix id (Subtype.val : {i // s i} → ι)).mulVecLin

theorem finrank_coordSub {U : Matrix ι ι ℝ} (hU : Uᵀ * U = 1) (s : ι → Prop)
    [DecidablePred s] : Module.finrank ℝ (coordSub U s) = Fintype.card {i // s i} := by
  set E := (1 : Matrix ι ι ℝ).submatrix id (Subtype.val : {i // s i} → ι)
  have hE : Eᵀ * E = 1 := by
    rw [Matrix.transpose_submatrix, Matrix.transpose_one, ← Matrix.submatrix_mul _ _ _ _ _
      Function.bijective_id, Matrix.one_mul, Matrix.submatrix_one _ Subtype.val_injective]
  have hUE : (U * E)ᵀ * (U * E) = 1 := by
    rw [Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc Uᵀ, hU, Matrix.one_mul, hE]
  rw [coordSub, LinearMap.finrank_range_of_inj (injective_of_transpose_mul hUE),
    Module.finrank_fintype_fun_eq_card]

theorem mem_coordSub {U : Matrix ι ι ℝ} {s : ι → Prop} [DecidablePred s] {x : ι → ℝ}
    (hx : x ∈ coordSub U s) : ∃ e, x = U *ᵥ e ∧ ∀ i, ¬ s i → e i = 0 := by
  obtain ⟨w, rfl⟩ := hx
  refine ⟨(1 : Matrix ι ι ℝ).submatrix id (Subtype.val : {i // s i} → ι) *ᵥ w, ?_, ?_⟩
  · simp [Matrix.mulVec_mulVec]
  · intro i hi
    simp only [mulVec, dotProduct, submatrix_apply, id_eq, one_apply]
    refine Finset.sum_eq_zero fun j _ => ?_
    have : i ≠ j.val := fun h => hi (h ▸ j.property)
    simp [this]

/-- Quadratic form of a diagonalised matrix. -/
theorem quad_diag {A U : Matrix ι ι ℝ} {d : ι → ℝ} (hD : Uᵀ * A * U = diagonal d)
    (e : ι → ℝ) : (U *ᵥ e) ⬝ᵥ (A *ᵥ (U *ᵥ e)) = ∑ i, d i * e i ^ 2 := by
  rw [quad_mulVec, hD]
  simp only [dotProduct, mulVec_diagonal]
  exact Finset.sum_congr rfl fun i _ => by ring

theorem isPosOn_coordSub {A U : Matrix ι ι ℝ} {d : ι → ℝ} (hD : Uᵀ * A * U = diagonal d)
    (s : ι → Prop) [DecidablePred s] (hs : ∀ i, s i → 0 < d i) :
    IsPosOn A (coordSub U s) := by
  intro x hx hx0
  obtain ⟨e, rfl, he⟩ := mem_coordSub hx
  rw [quad_diag hD]
  have hne : e ≠ 0 := by rintro rfl; exact hx0 (Matrix.mulVec_zero U)
  obtain ⟨i, hi⟩ := Function.ne_iff.1 hne
  have hsi : s i := by by_contra h; exact hi (he i h)
  refine Finset.sum_pos' (fun j _ => ?_) ⟨i, Finset.mem_univ _, ?_⟩
  · by_cases hj : s j
    · exact mul_nonneg (hs j hj).le (sq_nonneg _)
    · simp [he j hj]
  · exact mul_pos (hs i hsi) (lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hi)))

theorem isNonposOn_coordSub {A U : Matrix ι ι ℝ} {d : ι → ℝ} (hD : Uᵀ * A * U = diagonal d)
    (s : ι → Prop) [DecidablePred s] (hs : ∀ i, s i → d i ≤ 0) :
    IsNonposOn A (coordSub U s) := by
  intro x hx
  obtain ⟨e, rfl, he⟩ := mem_coordSub hx
  rw [quad_diag hD]
  refine Finset.sum_nonpos fun j _ => ?_
  by_cases hj : s j
  · exact mul_nonpos_of_nonpos_of_nonneg (hs j hj) (sq_nonneg _)
  · simp [he j hj]

/-- Sylvester's law for a diagonalised matrix: the positive index is the number of
positive diagonal entries. -/
theorem posIndex_eq_of_diag {A U : Matrix ι ι ℝ} {d : ι → ℝ} (hU : Uᵀ * U = 1)
    (hD : Uᵀ * A * U = diagonal d) : posIndex A = Fintype.card {i // 0 < d i} := by
  rw [← finrank_coordSub hU (fun i => 0 < d i)]
  apply posIndex_eq_of_split (isPosOn_coordSub hD _ fun _ h => h)
    (isNonposOn_coordSub hD (fun i => ¬ 0 < d i) fun _ h => not_lt.1 h)
  rw [finrank_coordSub hU, finrank_coordSub hU, Fintype.card_subtype_compl]
  have := Fintype.card_subtype_le (fun i => 0 < d i)
  omega

theorem exists_split_of_diag {A U : Matrix ι ι ℝ} {d : ι → ℝ} (hU : Uᵀ * U = 1)
    (hD : Uᵀ * A * U = diagonal d) : ∃ P N : Submodule ℝ (ι → ℝ), IsPosOn A P ∧
      IsNonposOn A N ∧ Module.finrank ℝ P + Module.finrank ℝ N = Fintype.card ι := by
  refine ⟨_, _, isPosOn_coordSub hD (fun i => 0 < d i) fun _ h => h,
    isNonposOn_coordSub hD (fun i => ¬ 0 < d i) fun _ h => not_lt.1 h, ?_⟩
  rw [finrank_coordSub hU, finrank_coordSub hU, Fintype.card_subtype_compl]
  have := Fintype.card_subtype_le (fun i => 0 < d i)
  omega

theorem exists_diag {A : Matrix ι ι ℝ} (hA : A.IsHermitian) :
    ∃ U : Matrix ι ι ℝ, Uᵀ * U = 1 ∧ Uᵀ * A * U = diagonal hA.eigenvalues := by
  refine ⟨(hA.eigenvectorUnitary : Matrix ι ι ℝ), ?_, ?_⟩
  · have := Unitary.coe_star_mul_self hA.eigenvectorUnitary
    simpa [Matrix.star_eq_conjTranspose] using this
  · have := hA.conjStarAlgAut_star_eigenvectorUnitary
    rw [Unitary.conjStarAlgAut_star_apply] at this
    simpa [Matrix.star_eq_conjTranspose, RCLike.ofReal_real_eq_id] using this

/-- For symmetric `A`, `posIndex A` is the number of positive eigenvalues. -/
theorem posIndex_eq_card_pos {A : Matrix ι ι ℝ} (hA : A.IsHermitian) :
    posIndex A = Fintype.card {i // 0 < hA.eigenvalues i} := by
  obtain ⟨U, hU, hD⟩ := exists_diag hA
  exact posIndex_eq_of_diag hU hD

/-- For symmetric `A`, `negIndex A` is the number of negative eigenvalues. -/
theorem negIndex_eq_card_neg {A : Matrix ι ι ℝ} (hA : A.IsHermitian) :
    negIndex A = Fintype.card {i // hA.eigenvalues i < 0} := by
  obtain ⟨U, hU, hD⟩ := exists_diag hA
  have hD' : Uᵀ * (-A) * U = diagonal (-hA.eigenvalues) := by
    rw [Matrix.mul_neg, Matrix.neg_mul, hD, Matrix.diagonal_neg]; rfl
  rw [negIndex, posIndex_eq_of_diag hU hD']
  simp [Pi.neg_apply]

/-- For symmetric `A`, `zeroIndex A` is the number of zero eigenvalues. -/
theorem zeroIndex_eq_card_zero {A : Matrix ι ι ℝ} (hA : A.IsHermitian) :
    zeroIndex A = Fintype.card {i // hA.eigenvalues i = 0} := by
  rw [zeroIndex, hA.rank_eq_card_non_zero_eigs]
  have := Fintype.card_subtype_compl (fun i => hA.eigenvalues i = 0)
  have h2 := Fintype.card_subtype_le (fun i => hA.eigenvalues i = 0)
  simp only [ne_eq]
  omega

theorem posIndex_add_negIndex {A : Matrix ι ι ℝ} (hA : A.IsHermitian) :
    posIndex A + negIndex A = A.rank := by
  rw [posIndex_eq_card_pos hA, negIndex_eq_card_neg hA, hA.rank_eq_card_non_zero_eigs,
    ← Fintype.card_subtype_or_disjoint]
  · exact Fintype.card_congr (Equiv.subtypeEquivRight fun i => by
      constructor
      · rintro (h | h) <;> [exact h.ne'; exact h.ne]
      · intro h; exact lt_or_gt_of_ne h |>.symm)
  · intro p hp hq i hi
    exact absurd ((hp i hi).trans (hq i hi)) (lt_irrefl _)

end Inertia

/-! ### The bordered inertia count -/

section B2

omit [Fintype m] [DecidableEq n] [DecidableEq m] [DecidableEq k] in
theorem ker_eq_range {M : Matrix m n ℝ} {Z : Matrix n k ℝ} (hZ : IsKernelBasis M Z) :
    LinearMap.ker M.mulVecLin = LinearMap.range Z.mulVecLin := by
  ext x
  simp only [LinearMap.mem_ker, LinearMap.mem_range, Matrix.mulVecLin_apply]
  constructor
  · exact hZ.spans x
  · rintro ⟨z, rfl⟩
    rw [Matrix.mulVec_mulVec, hZ.mul_eq_zero, Matrix.zero_mulVec]

omit [Fintype m] [DecidableEq n] [DecidableEq m] [DecidableEq k] in
/-- Rank–nullity with a kernel basis: `card n = card k + rank M`. -/
theorem card_eq_of_kernelBasis {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    (hZ : IsKernelBasis M Z) : Fintype.card n = Fintype.card k + M.rank := by
  have h1 := LinearMap.finrank_range_add_finrank_ker M.mulVecLin
  rw [Module.finrank_fintype_fun_eq_card, ker_eq_range hZ,
    LinearMap.finrank_range_of_inj (f := Z.mulVecLin) hZ.injective,
    Module.finrank_fintype_fun_eq_card] at h1
  have : M.rank = Module.finrank ℝ (LinearMap.range M.mulVecLin) := rfl
  omega

omit [DecidableEq n] [DecidableEq k] in
/-- A full-row-rank matrix has a right inverse. -/
theorem exists_mul_eq_one {M : Matrix m n ℝ} (hM : M.rank = Fintype.card m) :
    ∃ G : Matrix n m ℝ, M * G = 1 := by
  have htop : LinearMap.range M.mulVecLin = ⊤ := by
    apply Submodule.eq_top_of_finrank_eq
    rw [Module.finrank_fintype_fun_eq_card]
    exact hM
  have hs : ∀ j : m, ∃ g : n → ℝ, M *ᵥ g = Pi.single j 1 := fun j => by
    have : Pi.single j (1 : ℝ) ∈ LinearMap.range M.mulVecLin := htop ▸ Submodule.mem_top
    obtain ⟨g, hg⟩ := this
    exact ⟨g, hg⟩
  choose g hg using hs
  refine ⟨Matrix.of fun i j => g j i, ?_⟩
  ext i j
  have := congrFun (hg j) i
  simp only [mulVec, dotProduct] at this
  simp [Matrix.mul_apply, this, Matrix.one_apply, Pi.single_apply]

/-- The congruence `[Z G; -GᵀHZ, σI - ½GᵀHG]` that block-diagonalises `K`. -/
noncomputable def liftMat (H : Matrix n n ℝ) (Z : Matrix n k ℝ) (G : Matrix n m ℝ) (σ : ℝ) :
    Matrix (n ⊕ m) (k ⊕ m) ℝ :=
  fromBlocks Z G (-(Gᵀ * H * Z)) (σ • 1 - (1 / 2 : ℝ) • (Gᵀ * H * G))

omit [Fintype k] [DecidableEq n] [DecidableEq k] in
theorem liftMat_congr {H : Matrix n n ℝ} {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    {G : Matrix n m ℝ} (σ : ℝ) (hH : Hᵀ = H) (hMZ : M * Z = 0) (hMG : M * G = 1) :
    (liftMat H Z G σ)ᵀ * bordered H M * liftMat H Z G σ =
      fromBlocks (Zᵀ * H * Z) 0 0 ((2 * σ) • (1 : Matrix m m ℝ)) := by
  have hZt0 : Zᵀ * Mᵀ = 0 := by rw [← transpose_mul, hMZ, transpose_zero]
  have hGt0 : Gᵀ * Mᵀ = 1 := by rw [← transpose_mul, hMG, transpose_one]
  simp only [liftMat, bordered, fromBlocks_transpose, fromBlocks_multiply]
  rw [fromBlocks_inj]
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
    simp only [Matrix.add_mul, Matrix.mul_assoc, hZt0, hGt0, hMZ, hMG, smul_zero,
      Matrix.mul_zero, Matrix.zero_mul, add_zero, transpose_neg, transpose_mul,
      transpose_transpose, hH, transpose_sub, transpose_smul, transpose_one, Matrix.neg_mul,
      Matrix.mul_neg, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul, Matrix.mul_smul,
      Matrix.one_mul, Matrix.mul_one]
  all_goals module

omit [DecidableEq n] [DecidableEq k] in
theorem quad_fromBlocks (S : Matrix k k ℝ) (c : ℝ) (z : k → ℝ) (w : m → ℝ) :
    Sum.elim z w ⬝ᵥ (fromBlocks S 0 0 (c • (1 : Matrix m m ℝ)) *ᵥ Sum.elim z w) =
      z ⬝ᵥ (S *ᵥ z) + c * (w ⬝ᵥ w) := by
  rw [fromBlocks_mulVec]
  simp [sumElim_dotProduct_sumElim, Matrix.smul_mulVec, dotProduct_smul]

/-- `(z, w) ↦ T (z, w)` for `z ∈ V`. -/
noncomputable def liftLin (V : Submodule ℝ (k → ℝ)) (T : Matrix (n ⊕ m) (k ⊕ m) ℝ) :
    V × (m → ℝ) →ₗ[ℝ] (n ⊕ m → ℝ) :=
  T.mulVecLin ∘ₗ (LinearEquiv.sumArrowLequivProdArrow k m ℝ ℝ).symm.toLinearMap ∘ₗ
    V.subtype.prodMap LinearMap.id

omit [Fintype n] [DecidableEq n] [DecidableEq m] [DecidableEq k] in
theorem liftLin_apply (V : Submodule ℝ (k → ℝ)) (T : Matrix (n ⊕ m) (k ⊕ m) ℝ)
    (p : V × (m → ℝ)) : liftLin V T p = T *ᵥ Sum.elim (p.1 : k → ℝ) p.2 := by
  simp only [liftLin, LinearMap.coe_comp, Function.comp_apply, Matrix.mulVecLin_apply]
  congr 1

omit [DecidableEq n] [DecidableEq k] in
theorem liftMat_injective (H : Matrix n n ℝ) {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    {G : Matrix n m ℝ} (hZ : IsKernelBasis M Z) (hMG : M * G = 1) (σ : ℝ) :
    Function.Injective (liftMat H Z G σ).mulVec := by
  rw [mulVec_injective_iff]
  intro u hu
  rw [liftMat, fromBlocks_mulVec] at hu
  have h1 : Z *ᵥ (u ∘ Sum.inl) + G *ᵥ (u ∘ Sum.inr) = 0 := by
    simpa using congrArg (fun v => v ∘ Sum.inl) hu
  have h2 : u ∘ Sum.inr = 0 := by
    have := congrArg (fun v => M *ᵥ v) h1
    simpa [Matrix.mulVec_add, Matrix.mulVec_mulVec, hZ.mul_eq_zero, hMG] using this
  rw [h2, Matrix.mulVec_zero, add_zero] at h1
  have h3 : u ∘ Sum.inl = 0 := (mulVec_injective_iff Z).1 hZ.injective _ h1
  ext (i | i)
  · exact congrFun h3 i
  · exact congrFun h2 i

omit [DecidableEq n] [DecidableEq k] in
theorem finrank_range_liftLin (H : Matrix n n ℝ) {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    {G : Matrix n m ℝ} (hZ : IsKernelBasis M Z) (hMG : M * G = 1) (σ : ℝ)
    (V : Submodule ℝ (k → ℝ)) :
    Module.finrank ℝ (LinearMap.range (liftLin V (liftMat H Z G σ))) =
      Module.finrank ℝ V + Fintype.card m := by
  have hinj : Function.Injective (liftLin V (liftMat H Z G σ)) := by
    intro p q hpq
    rw [liftLin_apply, liftLin_apply] at hpq
    have := liftMat_injective H hZ hMG σ hpq
    have h1 := congrArg (fun v => v ∘ Sum.inl) this
    have h2 := congrArg (fun v => v ∘ Sum.inr) this
    simp only [Sum.elim_comp_inl, Sum.elim_comp_inr] at h1 h2
    exact Prod.ext (Subtype.ext h1) h2
  rw [LinearMap.finrank_range_of_inj hinj, Module.finrank_prod,
    Module.finrank_fintype_fun_eq_card]

theorem dotProduct_self_nonneg' {ι : Type*} [Fintype ι] (w : ι → ℝ) : 0 ≤ w ⬝ᵥ w := by
  simpa using dotProduct_star_self_nonneg w

omit [DecidableEq n] [DecidableEq k] in
theorem isPosOn_lift {H : Matrix n n ℝ} {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    {G : Matrix n m ℝ} (hH : Hᵀ = H) (hZ : IsKernelBasis M Z) (hMG : M * G = 1)
    {P : Submodule ℝ (k → ℝ)} (hP : IsPosOn (Zᵀ * H * Z) P) :
    IsPosOn (bordered H M) (LinearMap.range (liftLin P (liftMat H Z G 1))) := by
  rintro _ ⟨p, rfl⟩ hx
  rw [liftLin_apply] at hx ⊢
  rw [quad_mulVec, liftMat_congr 1 hH hZ.mul_eq_zero hMG, quad_fromBlocks]
  have hw := dotProduct_self_nonneg' p.2
  by_cases h1 : (p.1 : k → ℝ) = 0
  · have h2 : p.2 ≠ 0 := by
      rintro h
      apply hx
      rw [h1, h]
      convert Matrix.mulVec_zero (liftMat H Z G 1)
      ext (i | i) <;> rfl
    have h3 : 0 < p.2 ⬝ᵥ p.2 :=
      lt_of_le_of_ne hw (Ne.symm (fun h => h2 (dotProduct_self_eq_zero.1 h)))
    rw [h1]
    simp only [Matrix.mulVec_zero, dotProduct_zero, zero_add]
    linarith
  · have := hP _ p.1.2 h1
    nlinarith

omit [DecidableEq n] [DecidableEq k] in
theorem isNonposOn_lift {H : Matrix n n ℝ} {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    {G : Matrix n m ℝ} (hH : Hᵀ = H) (hZ : IsKernelBasis M Z) (hMG : M * G = 1)
    {N : Submodule ℝ (k → ℝ)} (hN : IsNonposOn (Zᵀ * H * Z) N) :
    IsNonposOn (bordered H M) (LinearMap.range (liftLin N (liftMat H Z G (-1)))) := by
  rintro _ ⟨p, rfl⟩
  rw [liftLin_apply, quad_mulVec, liftMat_congr (-1) hH hZ.mul_eq_zero hMG, quad_fromBlocks]
  have := hN _ p.1.2
  have hw := dotProduct_self_nonneg' p.2
  nlinarith

omit [Fintype k] [DecidableEq n] [DecidableEq k] in
theorem isHermitian_kernelForm {H : Matrix n n ℝ} (hH : Hᵀ = H) (Z : Matrix n k ℝ) :
    (Zᵀ * H * Z).IsHermitian := by
  simp [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial, transpose_mul, hH,
    Matrix.mul_assoc]

omit [DecidableEq n] in
theorem posIndex_bordered {H : Matrix n n ℝ} {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    (hH : Hᵀ = H) (hM : M.rank = Fintype.card m) (hZ : IsKernelBasis M Z) :
    posIndex (bordered H M) = posIndex (Zᵀ * H * Z) + Fintype.card m := by
  obtain ⟨U, hU, hD⟩ := exists_diag (isHermitian_kernelForm hH Z)
  obtain ⟨P, N, hP, hN, hd⟩ := exists_split_of_diag hU hD
  obtain ⟨G, hMG⟩ := exists_mul_eq_one hM
  have hcard := card_eq_of_kernelBasis hZ
  rw [posIndex_eq_of_split hP hN hd,
    posIndex_eq_of_split (isPosOn_lift hH hZ hMG hP) (isNonposOn_lift hH hZ hMG hN),
    finrank_range_liftLin H hZ hMG]
  rw [finrank_range_liftLin H hZ hMG, finrank_range_liftLin H hZ hMG, Fintype.card_sum]
  omega

omit [DecidableEq n] in
theorem negIndex_bordered {H : Matrix n n ℝ} {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    (hH : Hᵀ = H) (hM : M.rank = Fintype.card m) (hZ : IsKernelBasis M Z) :
    negIndex (bordered H M) = negIndex (Zᵀ * H * Z) + Fintype.card m := by
  have hneg : -bordered H M = bordered (-H) (-M) := by
    simp [bordered, fromBlocks_neg]
  have hH' : (-H)ᵀ = -H := by rw [transpose_neg, hH]
  have hM' : (-M).rank = Fintype.card m := by
    rw [rank_eq_card_iff_injective] at hM ⊢
    intro x y hxy
    apply hM
    simpa [transpose_neg, Matrix.neg_mulVec] using hxy
  have hZ' : IsKernelBasis (-M) Z :=
    ⟨by rw [Matrix.neg_mul, hZ.mul_eq_zero, neg_zero], hZ.injective,
      fun x hx => hZ.spans x (by simpa [Matrix.neg_mulVec] using hx)⟩
  rw [negIndex, negIndex, hneg, posIndex_bordered hH' hM' hZ', Matrix.mul_neg, Matrix.neg_mul]

/-- **Lemma B2 (inertia).** For symmetric `H`, `M : m × n` of full row rank
`r = card m` and a kernel basis `Z`, `In(K) = In(Zᵀ H Z) + (r, r, 0)`. -/
theorem lemma_B2 {H : Matrix n n ℝ} {M : Matrix m n ℝ} {Z : Matrix n k ℝ}
    (hH : Hᵀ = H) (hM : M.rank = Fintype.card m) (hZ : IsKernelBasis M Z) :
    posIndex (bordered H M) = posIndex (Zᵀ * H * Z) + Fintype.card m ∧
      negIndex (bordered H M) = negIndex (Zᵀ * H * Z) + Fintype.card m ∧
      zeroIndex (bordered H M) = zeroIndex (Zᵀ * H * Z) := by
  have hp := posIndex_bordered hH hM hZ
  have hn := negIndex_bordered hH hM hZ
  refine ⟨hp, hn, ?_⟩
  have hHh : H.IsHermitian := by
    simp [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial, hH]
  have hK : (bordered H M).IsHermitian :=
    IsHermitian.fromBlocks hHh (by simp [Matrix.conjTranspose_eq_transpose_of_trivial])
      isHermitian_zero
  have hrK := posIndex_add_negIndex hK
  have hrS := posIndex_add_negIndex (isHermitian_kernelForm hH Z)
  have hcard := card_eq_of_kernelBasis hZ
  rw [zeroIndex, zeroIndex, Fintype.card_sum]
  omega

end B2

end Kin6d.Bordered
