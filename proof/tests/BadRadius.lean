import Kin6d

set_option autoImplicit false

namespace Kin6d.Test

def corrupt : Certificate := ⟨1 / 2, 3 / 8, 3 / 16, 1 / 32⟩

-- This purported certificate must be rejected by the Lean kernel.
theorem corrupt_valid : corrupt.Valid := by
  norm_num [corrupt, Certificate.Valid]

end Kin6d.Test
