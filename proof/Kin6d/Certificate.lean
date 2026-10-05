import Kin6d.ManufacturedKrook

/-!
# Exact-rational data and a checked continuous theorem application

A certificate describes precisely the polynomial/cosine profiles from
`ManufacturedKrook`, with real parameters obtained by the canonical ℚ → ℝ
embedding. It does not assert that an external array reconstructs this profile.
-/

set_option autoImplicit false

namespace Kin6d

structure Certificate where
  exactAmplitude : ℚ
  candidateAmplitude : ℚ
  residualBound : ℚ
  radius : ℚ
  deriving Repr

def Certificate.Valid (c : Certificate) : Prop :=
  (3 / 2) * |c.candidateAmplitude - c.exactAmplitude| ≤ c.residualBound ∧
  c.residualBound ≤ c.radius ∧
  |c.exactAmplitude| ≤ 1 ∧ |c.candidateAmplitude| ≤ 1

instance (c : Certificate) : Decidable c.Valid := inferInstanceAs (Decidable (_ ∧ _ ∧ _ ∧ _))

/-- The numerical predicate only accepts data that instantiate the proved
continuous enclosure. All arithmetic uses exact rationals. -/
theorem Certificate.sound (c : Certificate) (hc : c.Valid) :
    ∀ v ∈ Set.Icc (-1 : ℝ) 1, ∀ x : ℝ,
      |profile (c.candidateAmplitude : ℝ) v x - profile (c.exactAmplitude : ℝ) v x| ≤
        (c.radius : ℝ) := by
  apply manufactured_krook_enclosure
    (c.exactAmplitude : ℝ) (c.candidateAmplitude : ℝ) (c.residualBound : ℝ)
  · have hq : ((3 / 2 : ℚ) : ℝ) * |(c.candidateAmplitude : ℝ) - (c.exactAmplitude : ℝ)| ≤ (c.residualBound : ℝ) := by
      exact_mod_cast hc.1
    norm_num at hq
    exact hq
  · exact_mod_cast hc.2.1

/-- The exact profile is a continuous solution of the stated source-balanced
problem; existence is discharged, rather than left as an input premise. -/
theorem Certificate.existence (c : Certificate) :
    ∀ v : ℝ, ∃ f : ℝ → ℝ,
      Differentiable ℝ f ∧ Function.Periodic f (2 * Real.pi) ∧
      ∀ x, v * deriv f x + f x = forcing (c.exactAmplitude : ℝ) v x := by
  intro v
  exact ⟨profile (c.exactAmplitude : ℝ) v,
    profile_differentiable _ _, profile_periodic _ _, profile_exact _ _⟩

/-- Uniqueness among differentiable periodic solutions for every fixed v.
The proof uses the continuous maximum principle with a zero residual. -/
theorem Certificate.uniqueness (c : Certificate) (v : ℝ) (g : ℝ → ℝ)
    (hg : Differentiable ℝ g) (hp : Function.Periodic g (2 * Real.pi))
    (heq : ∀ x, v * deriv g x + g x = forcing (c.exactAmplitude : ℝ) v x) :
    g = profile (c.exactAmplitude : ℝ) v := by
  have he := periodic_transport_error
    (f := profile (c.exactAmplitude : ℝ) v) (g := g)
    (s := forcing (c.exactAmplitude : ℝ) v)
    (L := 2 * Real.pi) (v := v) (nu := 1) (R := 0)
    (by positivity) (by norm_num) (profile_differentiable _ _) hg
    (profile_periodic _ _) hp
    (fun x => by simpa using profile_exact (c.exactAmplitude : ℝ) v x)
    (fun x => by simp only [one_mul, heq x, sub_self, abs_zero, le_refl])
  funext x
  have hz : |g x - profile (c.exactAmplitude : ℝ) v x| = 0 :=
    le_antisymm (by simpa using he x) (abs_nonneg _)
  exact sub_eq_zero.mp (abs_eq_zero.mp hz)

/-- A physical nonnegativity fact, independent of residual certification. -/
theorem Certificate.exact_nonnegative (c : Certificate) (hc : c.Valid)
    {v : ℝ} (hv : v ∈ Set.Icc (-1) 1) (x : ℝ) :
    0 ≤ profile (c.exactAmplitude : ℝ) v x := by
  have ha : |(c.exactAmplitude : ℝ)| ≤ 1 := by exact_mod_cast hc.2.2.1
  have hprod : |(c.exactAmplitude : ℝ) * Real.cos x| ≤ 1 := by
    rw [abs_mul]
    calc
      |(c.exactAmplitude : ℝ)| * |Real.cos x| ≤ 1 * 1 := by gcongr; exact Real.abs_cos_le_one x
      _ = 1 := by norm_num
  unfold profile
  apply mul_nonneg (weight_bounds hv).1
  linarith [neg_le_of_abs_le hprod]

end Kin6d
