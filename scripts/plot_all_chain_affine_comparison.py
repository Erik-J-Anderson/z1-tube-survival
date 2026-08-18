#!/usr/bin/env python3
import csv
import math
import os
import sys

import matplotlib.pyplot as plt

FRAME_TIME = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
FOCUS_DIAMETER = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def f(row, key):
    return float(row[key])


def nearest_diameter(rows, requested):
    ds = sorted({float(r["tube_diameter"]) for r in rows})
    return min(ds, key=lambda x: abs(x - requested))


def filter_diameter(rows, diameter):
    return [r for r in rows if math.isclose(float(r["tube_diameter"]), diameter,
                                             rel_tol=0.0, abs_tol=1e-12)]


tube_raw = read_csv("all_chains_tube_uncorrected.csv")
tube_aff = read_csv("all_chains_tube_affine.csv")
seg_raw = read_csv("all_chains_segment_uncorrected.csv")
seg_aff = read_csv("all_chains_segment_affine.csv")

focus = nearest_diameter(tube_raw, FOCUS_DIAMETER)
print(f"Plotting tube diameter a = {focus:g}")

tr = filter_diameter(tube_raw, focus)
ta = filter_diameter(tube_aff, focus)
sr = [r for r in filter_diameter(seg_raw, focus) if int(r["segment_index"]) == 50]
sa = [r for r in filter_diameter(seg_aff, focus) if int(r["segment_index"]) == 50]

tr.sort(key=lambda r: f(r, "lag_time"))
ta.sort(key=lambda r: f(r, "lag_time"))
sr.sort(key=lambda r: f(r, "lag_time"))
sa.sort(key=lambda r: f(r, "lag_time"))

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

ax = axes[0, 0]
ax.plot([f(r, "lag_time") for r in tr], [f(r, "mu") for r in tr],
        marker="o", label="Uncorrected")
ax.plot([f(r, "lag_time") for r in ta], [f(r, "mu") for r in ta],
        marker="o", linestyle="--", label="Affine corrected (fixed center)")
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel(r"Tube survival $\mu(t)$")
ax.set_title(f"All chains: tube survival, a={focus:g}")
ax.set_ylim(-0.02, 1.02)
ax.legend()

ax = axes[0, 1]
ax.plot([f(r, "lag_time") for r in sr], [f(r, "phi") for r in sr],
        marker="o", label="Uncorrected")
ax.plot([f(r, "lag_time") for r in sa], [f(r, "phi") for r in sa],
        marker="o", linestyle="--", label="Affine corrected (fixed center)")
ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
ax.set_ylabel(r"Center survival $\phi(s=0.5,t)$")
ax.set_title("Mid-chain segment survival")
ax.set_ylim(-0.02, 1.02)
ax.legend()

end_path = "all_chains_end_retraction.csv"
if os.path.exists(end_path):
    end_rows = filter_diameter(read_csv(end_path), focus)

    # Summary quantities are repeated for every s.  Keep s_index=0 once per lag.
    summaries = [r for r in end_rows if int(r["s_index"]) == 0]
    summaries.sort(key=lambda r: int(r["lag_frame"]))

    ax = axes[1, 0]
    ax.plot([int(r["lag_frame"]) * FRAME_TIME for r in summaries],
            [f(r, "mean_max_depth_end0") for r in summaries],
            marker="o", label="End 0")
    ax.plot([int(r["lag_frame"]) * FRAME_TIME for r in summaries],
            [f(r, "mean_max_depth_end1") for r in summaries],
            marker="o", label="End 1")
    ax.plot([int(r["lag_frame"]) * FRAME_TIME for r in summaries],
            [f(r, "mean_max_depth_combined") for r in summaries],
            marker="o", linestyle="--", label="Combined")
    ax.set_xlabel(r"Lag time [$\tau_{LJ}$]")
    ax.set_ylabel("Mean maximum penetration depth")
    ax.set_title("End retraction (current uncorrected geometry)")
    ax.set_ylim(bottom=0.0)
    ax.legend()

    ax = axes[1, 1]
    available_lags = sorted({int(r["lag_frame"]) for r in end_rows})
    requested = [10, 20, 40, 80]
    plot_lags = [lag for lag in requested if lag in available_lags]
    if not plot_lags:
        plot_lags = available_lags[-min(4, len(available_lags)):]

    for lag in plot_lags:
        rows = [r for r in end_rows if int(r["lag_frame"]) == lag]
        rows.sort(key=lambda r: f(r, "s_fraction"))
        ax.plot([f(r, "s_fraction") for r in rows],
                [f(r, "combined_reached_probability") for r in rows],
                label=fr"$\Delta t={lag * FRAME_TIME:g}$")
    ax.set_xlabel("Normalized contour penetration, s")
    ax.set_ylabel(r"$P(\ell_{max} \geq s)$")
    ax.set_title("End-retraction first-passage profiles")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.02, 1.02)
    ax.legend()
else:
    axes[1, 0].text(0.5, 0.5, "End retraction skipped",
                    ha="center", va="center", transform=axes[1, 0].transAxes)
    axes[1, 1].text(0.5, 0.5, "End retraction skipped",
                    ha="center", va="center", transform=axes[1, 1].transAxes)

fig.suptitle("All-chain deforming-box diagnostic")
fig.tight_layout()
fig.savefig("all_chains_affine_end_retraction_dashboard.png", dpi=200)
plt.close(fig)

# Also make a contour-profile comparison at selected lag times.
fig, ax = plt.subplots(figsize=(8, 6))
available_times = sorted({float(r["lag_time"]) for r in filter_diameter(seg_raw, focus)})
# Prefer roughly quarter/half/full-cycle-relevant entries if present in this dataset.
preferred_frames = [10, 20, 40, 80]
preferred_times = [FRAME_TIME * lag for lag in preferred_frames]
plot_times = [t for t in preferred_times if any(math.isclose(t, a, abs_tol=1e-12)
                                                for a in available_times)]
if not plot_times:
    plot_times = available_times[-min(4, len(available_times)):]

for lag_time in plot_times:
    raw_rows = [r for r in filter_diameter(seg_raw, focus)
                if math.isclose(f(r, "lag_time"), lag_time, abs_tol=1e-12)]
    aff_rows = [r for r in filter_diameter(seg_aff, focus)
                if math.isclose(f(r, "lag_time"), lag_time, abs_tol=1e-12)]
    raw_rows.sort(key=lambda r: f(r, "s_fraction"))
    aff_rows.sort(key=lambda r: f(r, "s_fraction"))
    ax.plot([f(r, "s_fraction") for r in raw_rows],
            [f(r, "phi") for r in raw_rows],
            label=fr"raw $\Delta t={lag_time:g}$")
    ax.plot([f(r, "s_fraction") for r in aff_rows],
            [f(r, "phi") for r in aff_rows],
            linestyle="--",
            label=fr"affine, fixed center $\Delta t={lag_time:g}$")

ax.set_xlabel("Normalized contour position, s")
ax.set_ylabel(r"Segment survival $\phi(s,t)$")
ax.set_title(f"All-chain contour profiles, a={focus:g}")
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-0.02, 1.02)
ax.legend(ncol=2, fontsize=8)
fig.tight_layout()
fig.savefig("all_chains_segment_profile_comparison.png", dpi=200)
plt.close(fig)

print("Wrote all_chains_affine_end_retraction_dashboard.png")
print("Wrote all_chains_segment_profile_comparison.png")
