# -*- coding: utf-8 -*-
"""
RF vs FR history-dependent tube-survival comparison.

Reads the output from:
    scripts/run_directional_history_comparison.sh

Produces:
    1. Full tube survival: RF vs FR vs longitudinal
    2. Tube-level RF - FR differences
    3. Full segment-survival RF - FR heatmap
    4. Full segment survival at largest lag
    5. Full RF - FR difference at largest lag
    6. Transverse survival at largest lag
    7. Permanent-escape survival at largest lag
    8. Longitudinal survival at largest lag

@author: ejanders
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# User settings
# ============================================================

project_dir = (
    r"C:\Users\ejanders\source\repos"
    r"\Parallel_Segment_Survival_Affine_Correction"
    r"\time_origin_affine_patch"
    r"\Parallel_Segment_Survival_TimeOrigin_Affine_FULL"
    r"\full_project"
)

output_dir = os.path.join(
    project_dir,
    "directional_history_output"
)

os.chdir(output_dir)


# ============================================================
# Helper functions
# ============================================================

def load_segment_pair(rf_filename, fr_filename):
    """
    Load RF and FR segment-survival files and align them by
    tube diameter, lag time, contour position, and sample count.
    """

    rf = pd.read_csv(rf_filename)
    fr = pd.read_csv(fr_filename)

    merged = rf.merge(
        fr,
        on=[
            "tube_diameter",
            "lag_time",
            "segment_index",
            "s_fraction",
            "sample_count"
        ],
        suffixes=("_rf", "_fr")
    )

    merged["delta_phi"] = (
        merged["phi_rf"]
        -
        merged["phi_fr"]
    )

    return merged


def plot_pair_at_lag(
    data,
    lag,
    ylabel,
    title,
    rf_label="Reference → Future",
    fr_label="Future → Reference"
):
    """
    Plot RF and FR contour profiles at one lag.
    """

    profile = data[
        np.isclose(data["lag_time"], lag)
    ].copy()

    profile = profile.sort_values("s_fraction")

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        profile["s_fraction"],
        profile["phi_rf"],
        label=rf_label
    )

    ax.plot(
        profile["s_fraction"],
        profile["phi_fr"],
        label=fr_label
    )

    ax.set_xlabel(
        r"Normalized contour position $s$"
    )

    ax.set_ylabel(ylabel)

    ax.set_ylim(0.0, 1.02)

    ax.set_title(title)

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()

    return profile


def plot_delta_at_lag(
    data,
    lag,
    ylabel,
    title
):
    """
    Plot RF - FR contour-survival difference at one lag.
    """

    profile = data[
        np.isclose(data["lag_time"], lag)
    ].copy()

    profile = profile.sort_values("s_fraction")

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        profile["s_fraction"],
        profile["delta_phi"],
        marker="o",
        markersize=3
    )

    ax.axhline(
        0.0,
        linestyle="--",
        linewidth=1
    )

    limit = np.max(
        np.abs(profile["delta_phi"].to_numpy())
    )

    if limit > 0.0:
        ax.set_ylim(
            -1.10 * limit,
            1.10 * limit
        )

    ax.set_xlabel(
        r"Normalized contour position $s$"
    )

    ax.set_ylabel(ylabel)

    ax.set_title(title)

    ax.grid(alpha=0.3)

    fig.tight_layout()


# ============================================================
# Load tube-level summary
# ============================================================

tube = pd.read_csv(
    "directional_history_tube_comparison.csv"
)

print("\nTube comparison columns:")
print(tube.columns.tolist())

print("\nTube comparison:")
print(tube)


# ============================================================
# 1. Full tube survival
# ============================================================

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(
    tube["lag_time"],
    tube["rf_full"],
    marker="o",
    label="Reference → Future"
)

ax.plot(
    tube["lag_time"],
    tube["fr_full"],
    marker="s",
    label="Future → Reference"
)

ax.plot(
    tube["lag_time"],
    tube["longitudinal"],
    marker="^",
    linestyle="--",
    label="Longitudinal"
)

ax.set_xlabel(
    r"Lag time [$\tau_{\mathrm{LJ}}$]"
)

ax.set_ylabel(
    r"Tube survival $\mu(t)$"
)

ax.set_ylim(0.0, 1.02)

ax.set_title(
    "Full tube survival"
)

ax.legend()
ax.grid(alpha=0.3)

fig.tight_layout()


# ============================================================
# 2. Tube-level directional differences
# ============================================================

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(
    tube["lag_time"],
    tube["delta_full"],
    marker="o",
    label="Full"
)

ax.plot(
    tube["lag_time"],
    tube["delta_transverse"],
    marker="s",
    label="Transverse"
)

ax.plot(
    tube["lag_time"],
    tube["delta_permanent"],
    marker="^",
    label="Permanent escape"
)

ax.axhline(
    0.0,
    linestyle="--",
    linewidth=1
)

ax.set_xlabel(
    r"Lag time [$\tau_{\mathrm{LJ}}$]"
)

ax.set_ylabel(
    r"$\Delta\mu"
    r"=\mu_{\mathrm{RF}}-\mu_{\mathrm{FR}}$"
)

ax.set_title(
    "Directional difference in tube survival"
)

ax.legend()
ax.grid(alpha=0.3)

fig.tight_layout()


# ============================================================
# Load segment-level fields
# ============================================================

full = load_segment_pair(
    "directional_history_segment_rf_full.csv",
    "directional_history_segment_fr_full.csv"
)

transverse = load_segment_pair(
    "directional_history_segment_rf_transverse.csv",
    "directional_history_segment_fr_transverse.csv"
)

permanent = load_segment_pair(
    "directional_history_segment_rf_permanent.csv",
    "directional_history_segment_fr_permanent.csv"
)

longitudinal = pd.read_csv(
    "directional_history_segment_longitudinal.csv"
)


# ============================================================
# Determine largest lag
# ============================================================

largest_lag = full["lag_time"].max()

print(
    "\nLargest lag time:",
    largest_lag
)


# ============================================================
# 3. Full RF - FR heatmap over s and t
# ============================================================

heatmap = full.pivot(
    index="lag_time",
    columns="s_fraction",
    values="delta_phi"
)

heatmap = heatmap.sort_index()
heatmap = heatmap.sort_index(axis=1)

fig, ax = plt.subplots(figsize=(8, 5))

max_abs_delta = np.max(
    np.abs(heatmap.to_numpy())
)

image = ax.imshow(
    heatmap.to_numpy(),
    aspect="auto",
    origin="lower",
    extent=[
        heatmap.columns.min(),
        heatmap.columns.max(),
        heatmap.index.min(),
        heatmap.index.max()
    ],
    vmin=-max_abs_delta,
    vmax=max_abs_delta
)

cbar = fig.colorbar(
    image,
    ax=ax
)

cbar.set_label(
    r"$\Delta\phi_{\mathrm{full}}(s,t)$"
)

ax.set_xlabel(
    r"Normalized contour position $s$"
)

ax.set_ylabel(
    r"Lag time [$\tau_{\mathrm{LJ}}$]"
)

ax.set_title(
    r"$\phi_{\mathrm{RF}}(s,t)"
    r"-\phi_{\mathrm{FR}}(s,t)$"
)

fig.tight_layout()


# ============================================================
# 4. Full survival profile at largest lag
# ============================================================

full_profile = plot_pair_at_lag(
    data=full,
    lag=largest_lag,
    ylabel=r"Full segment survival $\phi(s,t)$",
    title=(
        rf"Full survival, "
        rf"$t={largest_lag:g}\,\tau_{{\mathrm{{LJ}}}}$"
    )
)


# ============================================================
# 5. Full directional difference at largest lag
# ============================================================

plot_delta_at_lag(
    data=full,
    lag=largest_lag,
    ylabel=(
        r"$\phi_{\mathrm{RF}}"
        r"-\phi_{\mathrm{FR}}$"
    ),
    title=(
        rf"Full-survival directional difference, "
        rf"$t={largest_lag:g}\,\tau_{{\mathrm{{LJ}}}}$"
    )
)


# ============================================================
# 6. Transverse survival at largest lag
# ============================================================

transverse_profile = plot_pair_at_lag(
    data=transverse,
    lag=largest_lag,
    ylabel=(
        r"Transverse survival "
        r"$\phi_{\perp}(s,t)$"
    ),
    title=(
        rf"Transverse survival, "
        rf"$t={largest_lag:g}\,\tau_{{\mathrm{{LJ}}}}$"
    )
)


# ============================================================
# 7. Permanent-escape survival at largest lag
# ============================================================

permanent_profile = plot_pair_at_lag(
    data=permanent,
    lag=largest_lag,
    ylabel=(
        r"Permanent-escape survival "
        r"$\phi_{\mathrm{perm}}(s,t)$"
    ),
    title=(
        rf"Permanent-escape survival, "
        rf"$t={largest_lag:g}\,\tau_{{\mathrm{{LJ}}}}$"
    )
)


# ============================================================
# 8. Shared longitudinal survival at largest lag
# ============================================================

long_profile = longitudinal[
    np.isclose(
        longitudinal["lag_time"],
        largest_lag
    )
].copy()

long_profile = long_profile.sort_values(
    "s_fraction"
)

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(
    long_profile["s_fraction"],
    long_profile["phi"]
)

ax.set_xlabel(
    r"Normalized contour position $s$"
)

ax.set_ylabel(
    r"Longitudinal survival "
    r"$\phi_{\parallel}(s,t)$"
)

ax.set_ylim(0.0, 1.02)

ax.set_title(
    rf"Longitudinal survival, "
    rf"$t={largest_lag:g}\,\tau_{{\mathrm{{LJ}}}}$"
)

ax.grid(alpha=0.3)

fig.tight_layout()


# ============================================================
# 9. Transverse directional difference at largest lag
# ============================================================

plot_delta_at_lag(
    data=transverse,
    lag=largest_lag,
    ylabel=(
        r"$\phi_{\perp,\mathrm{RF}}"
        r"-\phi_{\perp,\mathrm{FR}}$"
    ),
    title=(
        rf"Transverse directional difference, "
        rf"$t={largest_lag:g}\,\tau_{{\mathrm{{LJ}}}}$"
    )
)


# ============================================================
# 10. Permanent-escape directional difference
# ============================================================

plot_delta_at_lag(
    data=permanent,
    lag=largest_lag,
    ylabel=(
        r"$\phi_{\mathrm{perm,RF}}"
        r"-\phi_{\mathrm{perm,FR}}$"
    ),
    title=(
        rf"Permanent-escape directional difference, "
        rf"$t={largest_lag:g}\,\tau_{{\mathrm{{LJ}}}}$"
    )
)


# ============================================================
# Numerical summary
# ============================================================

print("\n============================================")
print("DIRECTIONAL DIFFERENCE SUMMARY")
print("============================================")

print(
    "Max |delta phi_full|       :",
    np.max(
        np.abs(full["delta_phi"])
    )
)

print(
    "Max |delta phi_transverse| :",
    np.max(
        np.abs(transverse["delta_phi"])
    )
)

print(
    "Max |delta phi_permanent|  :",
    np.max(
        np.abs(permanent["delta_phi"])
    )
)

print(
    "Max |delta mu_full|         :",
    np.max(
        np.abs(tube["delta_full"])
    )
)

print(
    "Max |delta mu_transverse|   :",
    np.max(
        np.abs(tube["delta_transverse"])
    )
)

print(
    "Max |delta mu_permanent|    :",
    np.max(
        np.abs(tube["delta_permanent"])
    )
)


# ============================================================
# Show all figures
# ============================================================

plt.show()