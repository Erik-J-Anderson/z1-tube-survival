import csv
from pathlib import Path

import matplotlib.pyplot as plt


input_path = Path("affine_smoke_test.csv")
output_path = Path("affine_smoke_test.png")

s_fraction = []
uncorrected = []
corrected = []

with input_path.open(newline="", encoding="utf-8") as input_file:
    reader = csv.DictReader(input_file)

    for row in reader:
        s_fraction.append(float(row["s_fraction"]))
        uncorrected.append(float(row["uncorrected"]))
        corrected.append(float(row["corrected"]))

fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)

ax.plot(
    s_fraction,
    uncorrected,
    linewidth=2.0,
    label="Uncorrected",
)

ax.plot(
    s_fraction,
    corrected,
    linewidth=2.0,
    linestyle="--",
    label="Affine corrected",
)

ax.set_xlabel(r"Normalized contour position, $s$")
ax.set_ylabel("Segment survival")
ax.set_title("Affine-correction smoke test")
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-0.05, 1.05)
ax.grid(alpha=0.25)
ax.legend(frameon=False)

fig.savefig(output_path, dpi=200)

print(f"Wrote {output_path}")
