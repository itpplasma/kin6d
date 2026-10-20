import Kin6d

set_option autoImplicit false

namespace Kin6d.Test.Bordered

open Matrix Kin6d.Bordered

/-- One constraint `x₀ = 0` on `ℝ²`; `Z = e₁` spans its kernel. -/
def M : Matrix (Fin 1) (Fin 2) ℝ := !![1, 0]
def Z : Matrix (Fin 2) (Fin 1) ℝ := !![0; 1]
/-- Indefinite Hessian; the constrained block `Zᵀ H₁ Z = [-1]` is nonsingular. -/
def H₁ : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, -1]
/-- Hessian whose constrained block `Zᵀ H₀ Z = [0]` is singular. -/
def H₀ : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, 0]

theorem hZ : IsKernelBasis M Z where
  mul_eq_zero := by
    ext i j; fin_cases i; fin_cases j; simp [M, Z, Matrix.mul_apply, Fin.sum_univ_two]
  injective := by
    intro a b h
    have := congrFun h 1
    funext i; fin_cases i
    simpa [Z, mulVec, dotProduct] using this
  spans := by
    intro x hx
    have h0 := congrFun hx 0
    simp [M, mulVec, dotProduct, Fin.sum_univ_two] at h0
    refine ⟨![x 1], ?_⟩
    funext i; fin_cases i <;> simp [Z, mulVec, dotProduct, h0]

theorem hM : M.rank = Fintype.card (Fin 1) := by
  rw [rank_eq_card_iff_injective]
  intro a b h
  have := congrFun h 0
  funext i; fin_cases i
  simpa [M, mulVec, dotProduct, Fin.sum_univ_two] using this

theorem S₁ : Zᵀ * H₁ * Z = diagonal ![-1] := by
  ext i j; fin_cases i; fin_cases j; simp [Z, H₁, Matrix.mul_apply, Fin.sum_univ_two]

theorem S₀ : Zᵀ * H₀ * Z = 0 := by
  ext i j; fin_cases i; fin_cases j; simp [Z, H₀, Matrix.mul_apply, Fin.sum_univ_two]

/-- Lemma B1 applied: `K₁` is nonsingular. -/
theorem K₁_nonsingular : IsUnit (bordered H₁ M).det :=
  (lemma_B1 H₁ M hZ).2 ⟨hM, by simp [S₁]⟩

/-- Independent oracle: an explicit inverse of `K₁`. -/
theorem K₁_inverse :
    bordered H₁ M * fromBlocks !![0, 0; 0, -1] !![1; 0] !![1, 0] !![-1] = 1 := by
  ext (i | i) (j | j) <;> fin_cases i <;> fin_cases j <;>
    simp [bordered, H₁, M, Matrix.mul_apply, Fintype.sum_sum_type, Fin.sum_univ_two]

/-- Lemma B1 applied: `K₀` is singular. -/
theorem K₀_singular : ¬ IsUnit (bordered H₀ M).det := fun h =>
  by simpa [S₀] using ((lemma_B1 H₀ M hZ).1 h).2

/-- Independent oracle: `(Z e, 0)` is a nonzero kernel vector of `K₀`. -/
theorem K₀_kernel : bordered H₀ M *ᵥ Sum.elim ![0, 1] ![0] = 0 := by
  ext (i | i) <;> fin_cases i <;>
    simp [bordered, H₀, M, mulVec, dotProduct, Fintype.sum_sum_type, Fin.sum_univ_two]

/-- Lemma B2 applied: `In(K₁) = In([-1]) + (1, 1, 0) = (1, 2, ·)`. -/
theorem K₁_inertia : posIndex (bordered H₁ M) = 1 ∧ negIndex (bordered H₁ M) = 2 := by
  have hU : (1 : Matrix (Fin 1) (Fin 1) ℝ)ᵀ * 1 = 1 := by simp
  have hp : posIndex (Zᵀ * H₁ * Z) = 0 := by
    rw [posIndex_eq_of_diag hU (by simp [S₁] : (1 : Matrix (Fin 1) (Fin 1) ℝ)ᵀ *
      (Zᵀ * H₁ * Z) * 1 = diagonal ![-1])]
    exact Fintype.card_eq_zero_iff.2 ⟨fun ⟨i, hi⟩ => by fin_cases i; norm_num at hi⟩
  have hn : negIndex (Zᵀ * H₁ * Z) = 1 := by
    rw [negIndex, posIndex_eq_of_diag hU (by
      ext i j; fin_cases i; fin_cases j; simp [S₁] : (1 : Matrix (Fin 1) (Fin 1) ℝ)ᵀ *
      (-(Zᵀ * H₁ * Z)) * 1 = diagonal ![1])]
    simp
  obtain ⟨h1, h2, -⟩ := lemma_B2 (H := H₁) (by ext i j; fin_cases i <;> fin_cases j <;> simp [H₁]) hM hZ
  exact ⟨by rw [h1, hp]; rfl, by rw [h2, hn]; rfl⟩

#print axioms Kin6d.Bordered.lemma_B1
#print axioms Kin6d.Bordered.lemma_B2
#print axioms Kin6d.Bordered.posIndex_eq_card_pos
#print axioms Kin6d.Bordered.negIndex_eq_card_neg
#print axioms Kin6d.Bordered.zeroIndex_eq_card_zero

end Kin6d.Test.Bordered
