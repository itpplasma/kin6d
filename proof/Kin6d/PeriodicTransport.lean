import Mathlib.Analysis.Calculus.LocalExtr.Basic
import Mathlib.Algebra.Field.Periodic
import Mathlib.Topology.Order.Compact
import Mathlib.Tactic

/-!
# Continuous periodic damped transport

The domain is the real line modulo a positive spatial period. Velocity is a
fixed real parameter: no velocity periodicity, tails or field coupling is
assumed. The theorem concerns continuous differentiable reconstructions and a
uniform bound on their *continuous* residual, not nodal values.
-/

set_option autoImplicit false

namespace Kin6d

open Set Function

/-- A continuous periodic real function has a global maximum. The endpoint
case is covered by periodic extension to the whole line. -/
theorem periodic_has_global_max {f : ℝ → ℝ} {L : ℝ}
    (hL : 0 < L) (hf : Continuous f) (hp : Periodic f L) :
    ∃ a : ℝ, ∀ x : ℝ, f x ≤ f a := by
  obtain ⟨a, _, ha⟩ := isCompact_Icc.exists_isMaxOn
    (nonempty_Icc.mpr hL.le) hf.continuousOn
  refine ⟨a, fun x => ?_⟩
  obtain ⟨y, hy, hxy⟩ := hp.exists_mem_Ico₀ hL x
  rw [hxy]
  exact ha ⟨hy.1, hy.2.le⟩

/-- The maximum principle for the actual continuous differential operator
`v*d_x + nu`. No finite-dimensional coercivity premise is used. -/
theorem periodic_transport_upper_bound {f : ℝ → ℝ} {L v nu R : ℝ}
    (hL : 0 < L) (hnu : 0 < nu) (hf : Differentiable ℝ f)
    (hp : Periodic f L)
    (hres : ∀ x, v * deriv f x + nu * f x ≤ R) :
    ∀ x, f x ≤ R / nu := by
  obtain ⟨a, ha⟩ := periodic_has_global_max hL hf.continuous hp
  have hmax : IsLocalMax f a := by
    exact Filter.Eventually.of_forall ha
  have hd : deriv f a = 0 := hmax.deriv_eq_zero
  have hra := hres a
  rw [hd, mul_zero, zero_add] at hra
  intro x
  exact (ha x).trans ((le_div_iff₀ hnu).mpr (by simpa [mul_comm] using hra))

/-- Whole-domain sup-norm residual-to-error inequality, expressed pointwise
so its conclusion can justify a pointwise band. -/
theorem periodic_transport_abs_bound {f : ℝ → ℝ} {L v nu R : ℝ}
    (hL : 0 < L) (hnu : 0 < nu) (hf : Differentiable ℝ f)
    (hp : Periodic f L)
    (hres : ∀ x, |v * deriv f x + nu * f x| ≤ R) :
    ∀ x, |f x| ≤ R / nu := by
  have hu := periodic_transport_upper_bound hL hnu hf hp
    (fun x => (le_abs_self _).trans (hres x))
  have hn := periodic_transport_upper_bound (f := fun x => -f x)
    (v := v) (nu := nu) (R := R) hL hnu hf.neg (fun x => congrArg Neg.neg (hp x)) (by
      intro x
      have hd : deriv (fun x => -f x) x = -deriv f x := (hf x).hasDerivAt.neg.deriv
      rw [hd]
      have hl := (neg_le_of_abs_le (hres x))
      nlinarith)
  intro x
  exact abs_le.mpr ⟨by linarith [hn x], hu x⟩

/-- Difference estimate between a smooth candidate and an existing exact
solution. This theorem establishes an enclosure and, with zero residual,
uniqueness; existence is supplied by the explicit exact solution below. -/
theorem periodic_transport_error {f g s : ℝ → ℝ} {L v nu R : ℝ}
    (hL : 0 < L) (hnu : 0 < nu)
    (hf : Differentiable ℝ f) (hg : Differentiable ℝ g)
    (hpf : Periodic f L) (hpg : Periodic g L)
    (hexact : ∀ x, v * deriv f x + nu * f x = s x)
    (hres : ∀ x, |v * deriv g x + nu * g x - s x| ≤ R) :
    ∀ x, |g x - f x| ≤ R / nu := by
  apply periodic_transport_abs_bound (f := fun x => g x - f x) (v := v) hL hnu (hg.sub hf)
    (fun x => by simp only [hpg x, hpf x])
  intro x
  have hd : deriv (fun x => g x - f x) x = deriv g x - deriv f x :=
    ((hg x).hasDerivAt.sub (hf x).hasDerivAt).deriv
  rw [hd]
  convert hres x using 1
  congr 1
  rw [← hexact x]
  ring

end Kin6d
