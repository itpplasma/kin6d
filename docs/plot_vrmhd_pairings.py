#!/usr/bin/env python3
"""Plot the Fortran pairing study; matplotlib is an optional documentation tool."""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parents[1]
                        / "docs/data/vrmhd-pairings.csv")
    parser.add_argument("--output", type=Path, default=Path("vrmhd-pairings.png"))
    args = parser.parse_args()
    with args.data.open() as stream:
        rows = list(csv.DictReader(stream))
    expected = {(mesh, n, ks, ka) for mesh in ("structured", "perturbed")
                for n in (4, 8) for ks in range(1, 5) for ka in range(1, 5)
                if ks >= ka - 1}
    observed = {(r["mesh"], int(r["cells_per_period"]), int(r["ks"]), int(r["ka"]))
                for r in rows}
    if observed != expected or len(rows) != len(expected):
        raise ValueError("The table must contain all 52 distinct study cases")
    for r in rows:
        if int(r["unmatched"]) != sum(int(r[c]) for c in ("magnetic", "velocity", "mixed")):
            raise ValueError("Energy classifications do not sum to the raw audit")
    maximum = max(int(r["unmatched"]) for r in rows)
    errors = [max(float(r[c]) for c in ("alfven_error", "fast_error", "slow_error"))
              for r in rows if int(r["cells_per_period"]) == 8]
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), layout="constrained")
    cmap = plt.get_cmap("cividis").copy()
    cmap.set_bad("#eeeeee")
    for row_index, mesh in enumerate(("structured", "perturbed")):
        for column, n in enumerate((4, 8, 8)):
            matrix = np.full((4, 4), np.nan)
            text = {}
            for r in rows:
                if r["mesh"] != mesh or int(r["cells_per_period"]) != n:
                    continue
                ks, ka = int(r["ks"]), int(r["ka"])
                value = (int(r["unmatched"]) if column < 2 else
                         max(float(r[c]) for c in ("alfven_error", "fast_error", "slow_error")))
                matrix[ks-1, ka-1] = value
                text[ks-1, ka-1] = f"{value:d}" if column < 2 else f"{value:.1e}"
            norm = Normalize(0, maximum) if column < 2 else LogNorm(min(errors), max(errors))
            ax = axes[row_index, column]
            artist = ax.imshow(matrix, origin="lower", cmap=cmap, norm=norm)
            for (i, j), label in text.items():
                ax.text(j, i, label, ha="center", va="center", fontsize=10,
                        color="white" if norm(matrix[i, j]) < 0.5 else "black")
            ax.set(xticks=range(4), xticklabels=range(1, 5), yticks=range(4),
                   yticklabels=range(1, 5), xlabel=r"Magnetic edge order $k_A$",
                   ylabel=r"Scalar / velocity degree $k_s$",
                   title=f"{mesh.capitalize()}, {n} cells per period")
            if column == 2:
                fig.colorbar(artist, ax=ax, label="Maximum relative frequency error [1]")
            elif column == 1:
                fig.colorbar(artist, ax=axes[row_index, :2],
                             label="Unmatched eigenvalues per Bloch block [count]")
    fig.suptitle("Uniform ideal MHD pairing study at wave vector (1, 0, 0) m⁻¹", fontsize=15)
    fig.supxlabel("Structured block: one cell. Perturbed block: four cells. "
                  "Fixed null allowance: three per block.\n"
                  "Gray cells: excluded degree pairs. "
                  "Errors use the nearest positive frequency to each analytic branch; "
                  "small errors do not certify a clean spectrum.", fontsize=10)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    columns = ("static_magnetic", "static_velocity")
    maximum = max(int(r[c]) for r in rows for c in columns)
    fig, axes = plt.subplots(2, 2, figsize=(9, 8), layout="constrained")
    norm = Normalize(0, maximum)
    for row_index, mesh in enumerate(("structured", "perturbed")):
        for column, key in enumerate(columns):
            matrix = np.full((4, 4), np.nan)
            for r in rows:
                if r["mesh"] == mesh and int(r["cells_per_period"]) == 8:
                    matrix[int(r["ks"])-1, int(r["ka"])-1] = int(r[key])
            ax = axes[row_index, column]
            artist = ax.imshow(matrix, origin="lower", cmap=cmap, norm=norm)
            for i, j in zip(*np.where(np.isfinite(matrix))):
                ax.text(j, i, f"{matrix[i, j]:.0f}", ha="center", va="center",
                        color="white" if norm(matrix[i, j]) < 0.5 else "black")
            ax.set(xticks=range(4), xticklabels=range(1, 5), yticks=range(4),
                   yticklabels=range(1, 5), xlabel=r"Magnetic edge order $k_A$",
                   ylabel=r"Scalar / velocity degree $k_s$",
                   title=f"{mesh.capitalize()}: {key.removeprefix('static_')} nulls")
    fig.colorbar(artist, ax=axes, label="Energy-bearing static nulls per Bloch block [count]")
    fig.suptitle("Energy split of the static coupling nullspace, 8 cells per period")
    fig.supxlabel("SVD quotient by the three analytic nulls; energy fraction ≥ 90%.\n"
                  "Structured block: one cell. Perturbed block: four cells.\n"
                  "Gray cells: excluded degree pairs. "
                  "SVD counts are separate from roundoff-sensitive QZ audit counts.", fontsize=10)
    fig.savefig(args.output.with_stem(args.output.stem + "-energy"), dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
