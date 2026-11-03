"""Plot nonlinear Fourier convergence for the nondimensional manufactured state."""
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
    with (Path(__file__).parent / "data/vrmhd-nonlinear-fourier.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), layout="constrained")
    for ax, symmetry in zip(axes, ["rotation", "translation", "helical"]):
        for family, color in [("equal", "#b25435"), ("mixed", "#236e9b")]:
            for degree, style in [(1, "o-"), (2, "s--")]:
                subset = sorted((r for r in rows if r["symmetry"] == symmetry
                                 and r["family"] == family
                                 and int(r["scalar_degree"]) == degree),
                                key=lambda r: int(r["cells"]))
                cells = np.array([int(r["cells"]) for r in subset])
                errors = np.array([float(r["vx"]) for r in subset])
                order = np.log2(errors[-2] / errors[-1])
                ax.loglog(cells, errors, style, color=color,
                          label=f"{family}, scalar degree {degree}: order {order:.2f}")
        ax.set(title=symmetry.capitalize(), xlabel="Cells per section side [1]",
               ylabel="Velocity x-component H1 error [1]")
        ax.set_xticks([2, 4, 8, 16])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(alpha=.2, which="both")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.23),
                  frameon=False, fontsize=8)
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
