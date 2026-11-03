"""Render independent current-wall accounting and time-refinement checks."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
with (root / "data/vrmhd-island-ledgers.csv").open() as stream:
    accounting = list(csv.DictReader(stream))
with (root / "data/vrmhd-island-time.csv").open() as stream:
    time = list(csv.DictReader(stream))
fig, axes = plt.subplots(1, 2, figsize=(9, 3.7), layout="constrained")
x = np.arange(len(accounting))
for key, label, color, marker, shift in [
    ("energy_rate", "Wall power omitted", "#D55E00", "s", -.07),
    ("balance_defect", "Wall power included", "#0072B2", "o", .07),
]:
    axes[0].semilogy(x + shift, [abs(float(row[key])) for row in accounting],
                     marker, color=color, label=label)
axes[0].set(xticks=x, xticklabels=[row["family"] for row in accounting],
            ylabel="Energy-rate balance defect [1]",
            title="Independent instantaneous balance")
axes[0].legend(loc="upper center", bbox_to_anchor=(.5, -.16), frameon=False)
dt = np.array([float(row["dt"]) for row in time])
for key, label, color, marker in [
    ("energy_defect", "Total energy", "#0072B2", "o"),
    ("mechanical_defect", "Mechanical energy", "#D55E00", "s"),
]:
    axes[1].loglog(dt, [abs(float(row[key])) for row in time],
                   marker + "-", color=color, label=label)
reference = abs(float(time[-1]["energy_defect"])) * (dt / dt[-1]) ** 2
axes[1].loglog(dt, reference, ":", color=".4", label=r"$\Delta t^2$ reference")
axes[1].set(xlabel=r"Time step $\Delta t$ [1]",
            ylabel="Integrated energy defect [1]",
            title="Mixed family, final time 0.5 [1]")
axes[1].legend(loc="upper center", bbox_to_anchor=(.5, -.2), frameon=False)
for ax in axes:
    ax.grid(axis="y", alpha=.2)
fig.savefig("island-ledgers.png", dpi=180)
