"""Native nondimensional initial-induction errors; see vrmhd-island-benchmark.md."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path)
args = parser.parse_args()
with Path(__file__).with_name("data").joinpath("vrmhd-island-induction.csv").open() as f:
    rows = list(csv.DictReader(f))
fig, ax = plt.subplots(figsize=(6.2, 4.2))
for family, color, marker, style in [
    ("equal", "#D55E00", "s", "--"),
    ("mixed", "#0072B2", "o", "-"),
]:
    data = [row for row in rows if row["family"] == family]
    ax.loglog([int(row["nx"]) for row in data],
              [100 * float(row["relative_error"]) for row in data],
              color=color, marker=marker, linestyle=style, label=family)
ax.axhline(2, color="0.35", linestyle=":", label="2% error budget")
ax.set_xticks([8, 16], labels=["8", "16"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("Cells per direction [1]")
ax.set_ylabel("Initial flux derivative error [%]")
ax.set_title("Independent resistive induction check")
ax.grid(True, which="both", alpha=0.2)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.23), ncol=3,
          frameon=False)
fig.subplots_adjust(left=0.15, right=0.97, top=0.88, bottom=0.29)
fig.savefig(args.output, dpi=180)
