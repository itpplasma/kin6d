import Kin6d

set_option autoImplicit false

namespace Kin6d.Test

/-- A coarse explicitly reconstructed candidate. -/
def coarse : Certificate := ⟨1 / 2, 3 / 8, 3 / 16, 3 / 16⟩

/-- A finer explicitly reconstructed candidate with a narrower band. -/
def fine : Certificate := ⟨1 / 2, 7 / 16, 3 / 32, 3 / 32⟩

theorem coarse_valid : coarse.Valid := by norm_num [coarse, Certificate.Valid]
theorem fine_valid : fine.Valid := by norm_num [fine, Certificate.Valid]

/-- These are instantiated continuous theorem applications checked by Lean's
kernel, including all rational numerical premises. -/
theorem coarse_continuous : ∀ v ∈ Set.Icc (-1 : ℝ) 1, ∀ x : ℝ,
    |profile (3 / 8) v x - profile (1 / 2) v x| ≤ 3 / 16 := by
  simpa [coarse] using coarse.sound coarse_valid

theorem fine_continuous : ∀ v ∈ Set.Icc (-1 : ℝ) 1, ∀ x : ℝ,
    |profile (7 / 16) v x - profile (1 / 2) v x| ≤ 3 / 32 := by
  simpa [fine] using fine.sound fine_valid

/-- Independent physical counterexample for a corrupted radius: the actual
continuous solution difference already exceeds it at (x,v)=(0,0). -/
theorem undersized_radius_false :
    ¬ |profile (3 / 8) 0 0 - profile (1 / 2) 0 0| ≤ (1 / 32 : ℝ) := by
  norm_num [profile, weight]

/-- The continuous differential residual also exceeds the corrupted residual
bound at (x,v)=(0,0); no floating-point computation supplies this oracle. -/
theorem undersized_residual_false :
    ¬ |(0 : ℝ) * deriv (profile (3 / 8) 0) 0 + profile (3 / 8) 0 0 -
        forcing (1 / 2) 0 0| ≤ 1 / 32 := by
  norm_num [profile, weight, forcing]

#print axioms Kin6d.periodic_transport_abs_bound
#print axioms Kin6d.Certificate.sound
#print axioms coarse_continuous
#print axioms fine_continuous
#print axioms Kin6d.poisson_equation
#print axioms Kin6d.potential_gauge
#print axioms Kin6d.electric_energy

end Kin6d.Test
