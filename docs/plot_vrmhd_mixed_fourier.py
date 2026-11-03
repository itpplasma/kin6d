"""Plot the retained equal/mixed Fourier comparison against analytic waves."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = Path(__file__).parent / "data/vrmhd-mixed-fourier.csv"
    with data.open() as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout="constrained")
    for ax, degree, branch, gate in zip(axes, [1, 2], ["slow", "fast"], [1.8, 3.8]):
        for family, color in [("equal", "#b25435"), ("mixed", "#236e9b")]:
            selected = sorted(
                (r for r in rows if r["family"] == family
                 and int(r["scalar_degree"]) == degree and r["branch"] == branch),
                key=lambda r: int(r["cells_per_period"]),
            )
            cells = np.array([int(r["cells_per_period"]) for r in selected])
            errors = np.array([float(r["relative_frequency_error"]) for r in selected])
            order = np.log2(errors[-2] / errors[-1])
            label = f"{family.capitalize()}: final order {order:.3f}"
            ax.loglog(cells, errors, "o-", color=color, label=label)
        ax.set(title=f"Distorted mesh, n=2, {branch}, scalar degree {degree}",
               xlabel="Cells per period [1]", ylabel="Relative frequency error [1]")
        ax.set_xticks(cells)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.text(.04, .05, f"Original order threshold: {gate}", transform=ax.transAxes)
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=.2, which="both")
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
