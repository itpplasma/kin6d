"""Plot registered manufactured-field projection errors in normalized units.

Regenerate data with test_vrmhd_mms (no argument, mixed1, mixed2); retain the
eight natural-norm errors at each printed symmetry and refinement level.
The analytic fields and strong sources are in gen/gen_vrmhd_mms.f90.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--data", type=Path,
                    default=Path(__file__).parent / "data/vrmhd-mixed-section.csv")
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
data = np.genfromtxt(args.data, delimiter=",", names=True, dtype=None, encoding="utf-8")
fig, axes = plt.subplots(1, 3, figsize=(10, 3.5), sharey=True)
for kind, (ax, title) in enumerate(zip(axes, ["Rotation", "Translation", "Screw"]), 1):
    for family, label, marker in [("equal1", "Existing: velocity P₁", "o"),
                                  ("mixed1", "Mixed: velocity P₂", "s")]:
        sample = data[(data["symmetry"] == kind) & (data["family"] == family)]
        h = 1 / sample["cells"]
        error = np.sqrt(sum(sample[name]**2 for name in ["vx_h1", "vy_h1", "vz_h1"]))
        ax.loglog(h, error, marker=marker, label=label)
        ax.loglog(h, error[-1] * (h/h[-1])**(1 if family == "equal1" else 2),
                  color="0.55", linestyle=":", linewidth=1)
    ax.set_title(title)
    ax.set_xticks([1/16, 1/8, 1/4], ["1/16", "1/8", "1/4"])
    ax.minorticks_off()
    ax.set_xlabel("Cell width / section width [1]")
    ax.grid(alpha=0.25)
axes[0].set_ylabel("Velocity H¹ error [normalized units]")
axes[0].legend(fontsize=8)
fig.suptitle("Projection of exact manufactured velocity; thermodynamic degree 1")
fig.text(0.5, 0.01, "Dotted references: first and second order. Fields and coordinates are nondimensional.",
         ha="center", fontsize=9)
fig.tight_layout(rect=(0, 0.06, 1, 0.94))
fig.savefig(args.output, dpi=180)
