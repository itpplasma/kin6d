import Kin6d.PeriodicTransport
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv

/-!
# A fully instantiated continuous bounded-velocity Krook reference

The interval is v ∈ [-1,1], x has period 2π, and nu = 1. The reservoir weight
is the normalized polynomial 3/4*(1-v²). The equation is
v*d_x f + f = weight(v) + S, where S is manufactured independently from the
exact amplitude. Only spatial differentiation is present. This is E=0 Krook,
not the nonlinear self-consistent Vlasov–Poisson or an FP certificate.
-/

set_option autoImplicit false

namespace Kin6d

open Function Set Real

noncomputable def weight (v : ℝ) : ℝ := (3 / 4) * (1 - v ^ 2)

noncomputable def profile (a v x : ℝ) : ℝ := weight v * (1 + a * cos x)

/-- Total forcing, including the fixed reservoir term. -/
noncomputable def forcing (a v x : ℝ) : ℝ :=
  -v * weight v * a * sin x + profile a v x

theorem profile_hasDerivAt (a v x : ℝ) :
    HasDerivAt (profile a v) (-weight v * a * sin x) x := by
  have hd := ((hasDerivAt_const x (1 : ℝ)).add
    ((Real.hasDerivAt_cos x).const_mul a)).const_mul (weight v)
  have heq : -weight v * a * sin x = weight v * (0 + a * -sin x) := by ring
  rw [heq]
  exact hd

theorem profile_differentiable (a v : ℝ) : Differentiable ℝ (profile a v) :=
  fun x => (profile_hasDerivAt a v x).differentiableAt

theorem profile_periodic (a v : ℝ) : Periodic (profile a v) (2 * π) := by
  intro x
  simp only [profile, Real.cos_add_two_pi]

theorem profile_exact (a v x : ℝ) :
    v * deriv (profile a v) x + profile a v x = forcing a v x := by
  rw [(profile_hasDerivAt a v x).deriv]
  unfold forcing
  ring

theorem weight_bounds {v : ℝ} (hv : v ∈ Icc (-1) 1) :
    0 ≤ weight v ∧ weight v ≤ 3 / 4 := by
  have hv2 : v ^ 2 ≤ 1 := by nlinarith [hv.1, hv.2]
  unfold weight
  constructor <;> nlinarith [sq_nonneg v]

/-- A global continuous residual bound. Trigonometric inequalities apply to
all real x; this estimate never samples a grid. -/
theorem profile_residual_bound (a b x : ℝ) {v : ℝ}
    (hv : v ∈ Icc (-1) 1) :
    |v * deriv (profile b v) x + profile b v x - forcing a v x| ≤
      (3 / 2) * |b - a| := by
  have hid : v * deriv (profile b v) x + profile b v x - forcing a v x =
      weight v * (b - a) * (cos x - v * sin x) := by
    rw [(profile_hasDerivAt b v x).deriv]
    unfold forcing profile
    ring
  rw [hid, abs_mul, abs_mul, abs_of_nonneg (weight_bounds hv).1]
  have hvabs : |v| ≤ 1 := abs_le.mpr hv
  have htrig : |cos x - v * sin x| ≤ 2 := by
    calc
      |cos x - v * sin x| ≤ |cos x| + |v * sin x| := by simpa only [sub_eq_add_neg, abs_neg] using abs_add_le (cos x) (-(v * sin x))
      _ = |cos x| + |v| * |sin x| := by rw [abs_mul]
      _ ≤ 1 + 1 * 1 := by gcongr <;> first | exact Real.abs_cos_le_one x | exact hvabs | exact Real.abs_sin_le_one x
      _ = 2 := by norm_num
  calc
    weight v * |b - a| * |cos x - v * sin x| ≤ (3 / 4) * |b - a| * 2 := by
      gcongr
      exact (weight_bounds hv).2
    _ = (3 / 2) * |b - a| := by ring

/-- Every reconstructed smooth polynomial/trigonometric candidate in this
family has a continuous solution-error enclosure at every x and v. -/
theorem manufactured_krook_enclosure (a b R r : ℝ)
    (hR : (3 / 2) * |b - a| ≤ R) (hr : R ≤ r) :
    ∀ v ∈ Icc (-1 : ℝ) 1, ∀ x,
      |profile b v x - profile a v x| ≤ r := by
  intro v hv x
  have he := periodic_transport_error
    (f := profile a v) (g := profile b v) (s := forcing a v)
    (L := 2 * π) (v := v) (nu := 1) (R := R)
    (by positivity) (by norm_num)
    (profile_differentiable a v) (profile_differentiable b v)
    (profile_periodic a v) (profile_periodic b v)
    (fun y => by simpa using profile_exact a v y)
    (fun y => by simpa using (profile_residual_bound a b y hv).trans hR)
  simpa using (he x).trans (by simpa using hr)

end Kin6d
