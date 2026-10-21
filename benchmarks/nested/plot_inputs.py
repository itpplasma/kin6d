"""Plot the measured input-preparation ladder, independently of solver costs."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
with Path(__file__).with_name("boundary-inputs.tsv").open() as stream:
    rows = list(csv.DictReader(stream, delimiter="\t"))
for case in ("landreman-integer", "issan-integer"):
    data = [r for r in rows if r["case"] == case]
    modes = [int(r["modes"]) for r in data]
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5), constrained_layout=True)
    axes[0].semilogy(modes, [float(r["boundary_sampled_error_m"]) for r in data],
                     "o-", color="#0072B2", label="Sampled boundary error")
    axes[0].axhline(0.00125, ls="--", color="#D55E00", label="Requested 1% radius")
    axes[0].set_ylabel("Boundary position error (m)")
    axes[0].legend(fontsize=8)
    axes[1].plot(modes, [float(r["export_wall_s"]) for r in data],
                 "s-", color="#009E73")
    axes[1].set_ylim(bottom=0)
    axes[1].set_ylabel("Export wall time (s), eight threads")
    for axis in axes:
        axis.set_xlabel("Maximum poloidal and toroidal mode")
        axis.set_xticks(modes)
        axis.grid(alpha=0.25)
    fig.suptitle(case.replace("-", " ").title() + ": input preparation")
    fig.text(0.5, -0.03, "One measurement per point; sampled indicators; equilibrium solve excluded",
             ha="center", fontsize=8)
    fig.savefig(args.output / f"{case}-inputs.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
