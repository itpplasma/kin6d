import Kin6d.ManufacturedKrook
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-!
# Unit Fourier-mode Poisson component

This component theorem uses normalized epsilon0=1 and k=1. It proves the
continuous Poisson equation, electric-field sign, exact charge compatibility
and mean-potential gauge. It supplies no nonlinear coupled solution-error bound.
-/

set_option autoImplicit false

namespace Kin6d

open Real intervalIntegral

noncomputable def potential (a x : ℝ) : ℝ := a * cos x
noncomputable def electric (a x : ℝ) : ℝ := a * sin x
noncomputable def charge (a x : ℝ) : ℝ := a * cos x

theorem potential_hasDerivAt (a x : ℝ) :
    HasDerivAt (potential a) (-a * sin x) x := by
  have hd := (Real.hasDerivAt_cos x).const_mul a
  have he : -a * sin x = a * -sin x := by ring
  rw [he]
  exact hd

theorem potential_deriv (a : ℝ) : deriv (potential a) = fun x => -a * sin x := by
  funext x
  exact (potential_hasDerivAt a x).deriv

/-- The field is the negative potential gradient. -/
theorem electric_sign (a x : ℝ) : electric a x = -deriv (potential a) x := by
  rw [potential_deriv]
  unfold electric
  ring

/-- The continuous normalized periodic Poisson equation. -/
theorem poisson_equation (a x : ℝ) :
    -deriv (deriv (potential a)) x = charge a x := by
  rw [potential_deriv]
  have hd : HasDerivAt (fun x => -a * sin x) (-a * cos x) x :=
    (Real.hasDerivAt_sin x).const_mul (-a)
  rw [hd.deriv]
  unfold charge
  ring

theorem potential_periodic (a : ℝ) :
    Function.Periodic (potential a) (2 * Real.pi) := by
  intro x
  simp only [potential, cos_add_two_pi]

/-- Zero net charge is proved exactly; no incompatible zero mode is removed. -/
theorem charge_neutral (a : ℝ) : ∫ x in (0 : ℝ)..(2 * Real.pi), charge a x = 0 := by
  simp only [charge, integral_const_mul, integral_cos, sin_two_pi, sin_zero,
    sub_self, mul_zero]

/-- The potential obeys the prescribed zero-mean gauge. -/
theorem potential_gauge (a : ℝ) : ∫ x in (0 : ℝ)..(2 * Real.pi), potential a x = 0 := by
  exact charge_neutral a

/-- The electric energy is nonnegative and its exact integral is Real.pi*a². -/
theorem electric_energy (a : ℝ) :
    ∫ x in (0 : ℝ)..(2 * Real.pi), (electric a x) ^ 2 = Real.pi * a ^ 2 := by
  have hid : (fun x => (electric a x) ^ 2) = fun x => a ^ 2 * (sin x) ^ 2 := by
    funext x
    unfold electric
    ring
  rw [hid, integral_const_mul, integral_sin_sq]
  simp only [sin_zero, cos_zero, sin_two_pi, cos_two_pi, zero_mul,
    sub_zero]
  ring

end Kin6d
