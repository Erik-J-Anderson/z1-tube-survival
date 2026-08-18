#!/usr/bin/env python3

import csv
import math
import sys

import matplotlib.pyplot as plt


FRAME_TIME = float(sys.argv[1])
FOCUS_DIAMETER = float(sys.argv[2])


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def f(row, key):
    return float(row[key])


def select_diameter(rows, diameter):
    return [
        row for row in rows
        if math.isclose(
            float(row["tube_diameter"]),
            diameter,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
    ]


tube_raw = select_diameter(
    read_csv("history_smoke_tube_raw.csv"),
    FOCUS_DIAMETER,
)

tube_affine = select_diameter(
    read_csv("history_smoke_tube_affine.csv"),
    FOCUS_DIAMETER,
)

tube_history = select_diameter(
    read_csv("history_smoke_tube_history.csv"),
    FOCUS_DIAMETER,
)

tube_permanent = select_diameter(
    read_csv("history_smoke_tube_permanent_escape.csv"),
    FOCUS_DIAMETER,
)

segment_raw = select_diameter(
    read_csv("history_smoke_segment_raw.csv"),
    FOCUS_DIAMETER,
)

segment_affine = select_diameter(
    read_csv("history_smoke_segment_affine.csv"),
    FOCUS_DIAMETER,
)

segment_history = select_diameter(
    read_csv("history_smoke_segment_history.csv"),
    FOCUS_DIAMETER,
)

segment_permanent = select_diameter(
    read_csv("history_smoke_segment_permanent_escape.csv"),
    FOCUS_DIAMETER,
)

for rows in (
    tube_raw,
    tube_affine,
    tube_history,
    tube_permanent,
):
    rows.sort(key=lambda row: f(row, "lag_time"))


mid_index = 50

mid_raw = [
    row for row in segment_raw
    if int(row["segment_index"]) == mid_index
]

mid_affine = [
    row for row in segment_affine
    if int(row["segment_index"]) == mid_index
]

mid_history = [
    row for row in segment_history
    if int(row["segment_index"]) == mid_index
]

mid_permanent = [
    row for row in segment_permanent
    if int(row["segment_index"]) == mid_index
]

for rows in (
    mid_raw,
    mid_affine,
    mid_history,
    mid_permanent,
):
    rows.sort(key=lambda row: f(row, "lag_time"))


fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

ax = axes[0]
ax.plot(
    [f(r, "lag_time") for r in tube_raw],
    [f(r, "mu") for r in tube_raw],
    marker="o",
    label="Raw instantaneous",
)
ax.plot(
    [f(r, "lag_time") for r in tube_affine],
    [f(r, "mu") for r in tube_affine],
    marker="o",
    linestyle="--",
    label="Affine instantaneous",
)
ax.plot(
    [f(r, "lag_time") for r in tube_history],
    [f(r, "mu") for r in tube_history],
    marker="o",
    linestyle="-.",
    label="Affine + permanent kill",
)
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel(r"Tube survival $\mu(t)$")
ax.set_title(f"Three-way tube comparison, a={FOCUS_DIAMETER:g}")
ax.set_ylim(-0.02, 1.02)
ax.legend(fontsize=8)


ax = axes[1]
ax.plot(
    [f(r, "lag_time") for r in mid_raw],
    [f(r, "phi") for r in mid_raw],
    marker="o",
    label="Raw instantaneous",
)
ax.plot(
    [f(r, "lag_time") for r in mid_affine],
    [f(r, "phi") for r in mid_affine],
    marker="o",
    linestyle="--",
    label="Affine instantaneous",
)
ax.plot(
    [f(r, "lag_time") for r in mid_history],
    [f(r, "phi") for r in mid_history],
    marker="o",
    linestyle="-.",
    label="Affine + permanent kill",
)
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel(r"Mid-chain survival $\phi(s=0.5,t)$")
ax.set_title("Mid-chain comparison")
ax.set_ylim(-0.02, 1.02)
ax.legend(fontsize=8)


ax = axes[2]
ax.plot(
    [f(r, "lag_time") for r in tube_permanent],
    [f(r, "mu") for r in tube_permanent],
    marker="o",
)
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel("Not-yet-permanently-escaped fraction")
ax.set_title(r"Pure absorbing boundary: $d_{\max}<a$")
ax.set_ylim(-0.02, 1.02)

fig.tight_layout()
fig.savefig(
    "history_escape_three_way_comparison.png",
    dpi=200,
)
plt.close(fig)


# Contour comparison at the largest requested lag, where the
# history effect should be easiest to see.
available_times = sorted(
    {f(row, "lag_time") for row in segment_affine}
)

selected_time = available_times[-1]

fig, ax = plt.subplots(figsize=(8, 6))

for rows, label, linestyle in (
    (segment_raw, "Raw instantaneous", "-"),
    (segment_affine, "Affine instantaneous", "--"),
    (segment_history, "Affine + permanent kill", "-."),
):
    selected = [
        row for row in rows
        if math.isclose(
            f(row, "lag_time"),
            selected_time,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
    ]

    selected.sort(
        key=lambda row: f(row, "s_fraction")
    )

    ax.plot(
        [f(r, "s_fraction") for r in selected],
        [f(r, "phi") for r in selected],
        linestyle=linestyle,
        label=label,
    )

ax.set_xlabel("Normalized contour position, s")
ax.set_ylabel(r"Segment survival $\phi(s,t)$")
ax.set_title(
    rf"Three-way contour comparison, $\Delta t={selected_time:g}$"
)
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-0.02, 1.02)
ax.legend()
fig.tight_layout()
fig.savefig(
    "history_escape_contour_comparison.png",
    dpi=200,
)
plt.close(fig)

print("Wrote history_escape_three_way_comparison.png")
print("Wrote history_escape_contour_comparison.png")
