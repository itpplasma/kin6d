"""Compare direct factors on the same fixed-step candidate island trajectory."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
with (root / "data/vrmhd-island-preconditioner.csv").open() as stream:
    trajectory = list(csv.DictReader(stream))
with (root / "data/vrmhd-island-preconditioner-time.csv").open() as stream:
    timing = list(csv.DictReader(stream))
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout="constrained")
methods = [("scalar", "Scalar direct", "#0072B2", "-"),
           ("supernodal", "Supernodal direct", "#D55E00", "o")]
for x, (method, label, color, style) in enumerate(methods):
    data = [row for row in trajectory if row["method"] == method]
    axes[0].plot([float(row["time"]) for row in data],
                 [float(row["width"]) for row in data], style,
                 color=color, label=label, markersize=4, markerfacecolor="none")
    seconds = np.array([float(row["seconds"]) for row in timing
                        if row["method"] == method])
    axes[1].bar(x, seconds.mean(), width=.55, color=color, alpha=.65)
    axes[1].plot(x + np.linspace(-.07, .07, len(seconds)), seconds,
                 "o", color="black", markersize=4)
axes[0].set(xlabel="Time [1]", ylabel="Full cut width [1]",
            title="Matched candidate trajectory")
axes[0].legend(frameon=False)
axes[1].set(xticks=[0, 1], xticklabels=[row[1] for row in methods],
            ylabel="Whole-run elapsed time [s]", title="Whole-run cost, four-thread setting")
for ax in axes:
    ax.grid(axis="y", alpha=.2)
fig.savefig("island-preconditioner.png", dpi=180)
