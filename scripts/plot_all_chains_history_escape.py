#!/usr/bin/env python3

import csv
import math
import sys

import matplotlib.pyplot as plt

FOCUS_DIAMETER = float(sys.argv[1])


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


def sorted_tube(path):
    rows = select_diameter(read_csv(path), FOCUS_DIAMETER)
    rows.sort(key=lambda row: f(row, "lag_time"))
    return rows


def sorted_mid(path):
    rows = select_diameter(read_csv(path), FOCUS_DIAMETER)
    rows = [row for row in rows if int(row["segment_index"]) == 50]
    rows.sort(key=lambda row: f(row, "lag_time"))
    return rows


tube_raw = sorted_tube("all_chains_history_tube_raw.csv")
tube_affine = sorted_tube("all_chains_history_tube_affine.csv")
tube_history = sorted_tube("all_chains_history_tube_history.csv")
tube_permanent = sorted_tube("all_chains_history_tube_permanent_escape.csv")

mid_raw = sorted_mid("all_chains_history_segment_raw.csv")
mid_affine = sorted_mid("all_chains_history_segment_affine.csv")
mid_history = sorted_mid("all_chains_history_segment_history.csv")
mid_permanent = sorted_mid("all_chains_history_segment_permanent_escape.csv")

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

ax = axes[0, 0]
for rows, label, linestyle in (
    (tube_raw, "Raw instantaneous", "-"),
    (tube_affine, "Affine instantaneous", "--"),
    (tube_history, "Affine + permanent kill", "-."),
):
    ax.plot(
        [f(r, "lag_time") for r in rows],
        [f(r, "mu") for r in rows],
        marker="o",
        linestyle=linestyle,
        label=label,
    )
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel(r"Tube survival $\mu(t)$")
ax.set_title(f"All chains: three-way tube comparison, a={FOCUS_DIAMETER:g}")
ax.set_ylim(-0.02, 1.02)
ax.legend(fontsize=8)

ax = axes[0, 1]
for rows, label, linestyle in (
    (mid_raw, "Raw instantaneous", "-"),
    (mid_affine, "Affine instantaneous", "--"),
    (mid_history, "Affine + permanent kill", "-."),
):
    ax.plot(
        [f(r, "lag_time") for r in rows],
        [f(r, "phi") for r in rows],
        marker="o",
        linestyle=linestyle,
        label=label,
    )
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel(r"Mid-chain survival $\phi(s=0.5,t)$")
ax.set_title("All chains: mid-chain comparison")
ax.set_ylim(-0.02, 1.02)
ax.legend(fontsize=8)

ax = axes[1, 0]
ax.plot(
    [f(r, "lag_time") for r in tube_permanent],
    [f(r, "mu") for r in tube_permanent],
    marker="o",
    label="Tube-averaged",
)
ax.plot(
    [f(r, "lag_time") for r in mid_permanent],
    [f(r, "phi") for r in mid_permanent],
    marker="o",
    linestyle="--",
    label="Mid-chain",
)
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel("Not-yet-permanently-escaped fraction")
ax.set_title(r"Outer absorbing boundary: $d_{max}<a$")
ax.set_ylim(-0.02, 1.02)
ax.legend(fontsize=8)

ax = axes[1, 1]
penalty_tube = [
    f(a, "mu") - f(h, "mu")
    for a, h in zip(tube_affine, tube_history)
]
penalty_mid = [
    f(a, "phi") - f(h, "phi")
    for a, h in zip(mid_affine, mid_history)
]
ax.plot(
    [f(r, "lag_time") for r in tube_affine],
    penalty_tube,
    marker="o",
    label="Tube-averaged",
)
ax.plot(
    [f(r, "lag_time") for r in mid_affine],
    penalty_mid,
    marker="o",
    linestyle="--",
    label="Mid-chain",
)
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel("Affine instantaneous - history survival")
ax.set_title("Recrossings removed by permanent escape")
ax.set_ylim(bottom=0.0)
ax.legend(fontsize=8)

fig.suptitle("All-chain affine + permanent-escape diagnostic")
fig.tight_layout()
fig.savefig("all_chains_history_escape_dashboard.png", dpi=200)
plt.close(fig)


# Contour-resolved three-way comparison at the largest requested lag.
segment_raw = select_diameter(
    read_csv("all_chains_history_segment_raw.csv"), FOCUS_DIAMETER)
segment_affine = select_diameter(
    read_csv("all_chains_history_segment_affine.csv"), FOCUS_DIAMETER)
segment_history = select_diameter(
    read_csv("all_chains_history_segment_history.csv"), FOCUS_DIAMETER)
segment_permanent = select_diameter(
    read_csv("all_chains_history_segment_permanent_escape.csv"), FOCUS_DIAMETER)

selected_time = max(f(row, "lag_time") for row in segment_affine)

fig, ax = plt.subplots(figsize=(9, 6))
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
    selected.sort(key=lambda row: f(row, "s_fraction"))
    ax.plot(
        [f(r, "s_fraction") for r in selected],
        [f(r, "phi") for r in selected],
        linestyle=linestyle,
        label=label,
    )
ax.set_xlabel("Normalized contour position, s")
ax.set_ylabel(r"Segment survival $\phi(s,t)$")
ax.set_title(
    rf"All-chain three-way contour comparison, $\Delta t={selected_time:g}$"
)
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-0.02, 1.02)
ax.legend()
fig.tight_layout()
fig.savefig("all_chains_history_contour_comparison.png", dpi=200)
plt.close(fig)


# Pure outer-boundary survival profile at the largest requested lag.
selected = [
    row for row in segment_permanent
    if math.isclose(
        f(row, "lag_time"),
        selected_time,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
]
selected.sort(key=lambda row: f(row, "s_fraction"))

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(
    [f(r, "s_fraction") for r in selected],
    [f(r, "phi") for r in selected],
)
ax.set_xlabel("Normalized contour position, s")
ax.set_ylabel("Not-yet-permanently-escaped fraction")
ax.set_title(
    rf"All-chain permanent-escape profile, $\Delta t={selected_time:g}$"
)
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-0.02, 1.02)
fig.tight_layout()
fig.savefig("all_chains_permanent_escape_profile.png", dpi=200)
plt.close(fig)

print("Wrote all_chains_history_escape_dashboard.png")
print("Wrote all_chains_history_contour_comparison.png")
print("Wrote all_chains_permanent_escape_profile.png")
