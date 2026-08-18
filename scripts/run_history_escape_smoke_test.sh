#!/usr/bin/env bash
set -euo pipefail

CXX=/usr/bin/g++

INPUT_FILE="${1:-Z1+SP.dat}"
FRAME_TIME="${2:-20}"
FOCUS_DIAMETER="${3:-5.0}"

BOX_CENTER_Z="${BOX_CENTER_Z:-}"
CHAIN_INDEX="${CHAIN_INDEX:-0}"

if [[ -z "$BOX_CENTER_Z" ]]; then
    echo "ERROR: BOX_CENTER_Z is required." >&2
    echo "Example:" >&2
    echo "  BOX_CENTER_Z=39.22658609485464 bash $0 Z1+SP.dat 20 5.0" >&2
    exit 2
fi

for f in \
    Geometry_Utils.cpp Geometry_Utils.hpp \
    Parse_Z1_File.cpp Parse_Z1_File.hpp \
    Parser_Utils.cpp Parser_Utils.hpp \
    Primitive_Path_Trajectory.cpp Primitive_Path_Trajectory.hpp \
    Trajectory_Time.cpp Trajectory_Time.hpp \
    Tube_Survival.cpp Tube_Survival.hpp \
    Survival_IO.cpp Survival_IO.hpp \
    tests/History_Dependent_Escape_Smoke_Test.cpp
do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: run this script from the patched project root. Missing: $f" >&2
        exit 1
    fi
done

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: input file not found: $INPUT_FILE" >&2
    exit 1
fi

UNIT_EXE="history_dependent_escape_unit_test"
REAL_CPP="Serial_History_Escape_Three_Way.cpp"
REAL_EXE="Serial_History_Escape_Three_Way"
PLOT_FILE="plot_history_escape_three_way.py"

echo "Using compiler:"
"$CXX" --version | head -n 1
echo
echo "Input file       : $INPUT_FILE"
echo "Frame time       : $FRAME_TIME tau_LJ"
echo "Focus diameter   : $FOCUS_DIAMETER"
echo "Chain index      : $CHAIN_INDEX"
echo "Fixed z center   : $BOX_CENTER_Z"
echo

echo "1) Compiling fast synthetic history test..."
"$CXX" \
    -std=c++20 \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -I. \
    tests/History_Dependent_Escape_Smoke_Test.cpp \
    Geometry_Utils.cpp \
    Primitive_Path_Trajectory.cpp \
    Tube_Survival.cpp \
    -o "$UNIT_EXE"

echo "2) Running fast synthetic history test..."
./"$UNIT_EXE"
echo

cat > "$REAL_CPP" <<'CPP'
#include "Geometry_Utils.hpp"
#include "Parse_Z1_File.hpp"
#include "Survival_IO.hpp"
#include "Trajectory_Time.hpp"
#include "Tube_Survival.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{

std::vector<std::size_t> KeepValidLags(
    const std::vector<std::size_t>& requested,
    std::size_t num_frames)
{
    std::vector<std::size_t> result;

    for (const std::size_t lag : requested) {
        if (lag < num_frames) {
            result.push_back(lag);
        }
    }

    if (result.empty()) {
        throw std::runtime_error(
            "No requested lag is valid for this trajectory.");
    }

    return result;
}


double GetBoxLengthZ(const Box& box)
{
    const double xz = box.matrix.value[0][2];
    const double yz = box.matrix.value[1][2];
    const double zz = box.matrix.value[2][2];

    return std::sqrt(
        xz * xz +
        yz * yz +
        zz * zz
    );
}


void SetFixedZBoxCenter(
    std::vector<Box>& frame_boxes,
    double z_center)
{
    if (frame_boxes.empty()) {
        throw std::invalid_argument(
            "SetFixedZBoxCenter: no frame boxes supplied.");
    }

    for (Box& box : frame_boxes)
    {
        const double Lz =
            GetBoxLengthZ(box);

        if (!(Lz > 0.0)) {
            throw std::runtime_error(
                "SetFixedZBoxCenter: encountered non-positive Lz.");
        }

        box.origin.z =
            z_center - 0.5 * Lz;
    }
}


void CheckFixedZCenterMapping(
    const std::vector<Box>& frame_boxes,
    double z_center)
{
    if (frame_boxes.empty()) {
        throw std::invalid_argument(
            "CheckFixedZCenterMapping: no frame boxes supplied.");
    }

    const Vec3 physical_center{
        0.0,
        0.0,
        z_center
    };

    const Box& destination_box =
        frame_boxes.front();

    constexpr double tolerance =
        1.0e-10;

    for (const Box& source_box : frame_boxes)
    {
        const Vec3 mapped_center =
            geometry::MapPositionBetweenBoxes(
                physical_center,
                source_box,
                destination_box
            );

        if (std::abs(
                mapped_center.z -
                z_center) > tolerance)
        {
            throw std::runtime_error(
                "Fixed-z-center affine regression check failed.");
        }
    }
}


void CheckLagZero(
    const TubeSurvivalFunction& result,
    const char* label)
{
    if (result.lag_times.empty() ||
        result.lag_times.front() != 0.0)
    {
        throw std::runtime_error(
            std::string(label) +
            ": lag-zero entry missing.");
    }

    const std::size_t num_lags =
        result.lag_times.size();

    for (std::size_t d = 0;
         d < result.tube_diameters.size();
         ++d)
    {
        const double mu0 =
            result.survival[
                d * num_lags
            ];

        if (std::abs(mu0 - 1.0) >
            1.0e-12)
        {
            throw std::runtime_error(
                std::string(label) +
                ": survival(0) != 1.");
        }
    }
}


void CheckHistorySubset(
    const SegmentSurvivalFunction& affine,
    const SegmentSurvivalFunction& history)
{
    if (affine.lag_times != history.lag_times ||
        affine.tube_diameters != history.tube_diameters ||
        affine.sample_counts != history.sample_counts ||
        affine.survival.size() != history.survival.size())
    {
        throw std::runtime_error(
            "Affine/history dimensions or origin cohorts do not match.");
    }

    constexpr double tolerance =
        1.0e-12;

    for (std::size_t i = 0;
         i < affine.survival.size();
         ++i)
    {
        if (history.survival[i] >
            affine.survival[i] +
                tolerance)
        {
            throw std::runtime_error(
                "History-dependent survival exceeded "
                "instantaneous affine survival.");
        }
    }
}


void CheckPermanentContainsHistory(
    const SegmentSurvivalFunction& history,
    const SegmentSurvivalFunction& permanent)
{
    if (history.lag_times != permanent.lag_times ||
        history.tube_diameters != permanent.tube_diameters ||
        history.sample_counts != permanent.sample_counts ||
        history.survival.size() != permanent.survival.size())
    {
        throw std::runtime_error(
            "History/permanent dimensions do not match.");
    }

    constexpr double tolerance =
        1.0e-12;

    for (std::size_t i = 0;
         i < history.survival.size();
         ++i)
    {
        if (history.survival[i] >
            permanent.survival[i] +
                tolerance)
        {
            throw std::runtime_error(
                "Inner-tube history survival exceeded "
                "pure permanent-escape survival.");
        }
    }
}

} // namespace


int main(int argc, char* argv[])
{
    try
    {
        if (argc != 6)
        {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <Z1+SP.dat> <frame_time> "
                << "<tube_diameter> <chain_index> <z_center>\n";

            return 2;
        }

        const std::string filename =
            argv[1];

        const long frame_time =
            std::stol(argv[2]);

        const double tube_diameter =
            std::stod(argv[3]);

        const std::size_t chain_index =
            static_cast<std::size_t>(
                std::stoull(argv[4])
            );

        const double z_center =
            std::stod(argv[5]);

        if (frame_time <= 0) {
            throw std::invalid_argument(
                "frame_time must be positive.");
        }

        if (tube_diameter <= 0.0) {
            throw std::invalid_argument(
                "tube_diameter must be positive.");
        }

        PrimitivePathTrajectory trajectory =
            parse_z1_file(filename);

        if (trajectory.chains.empty() ||
            trajectory.frame_boxes.empty())
        {
            throw std::runtime_error(
                "Parser returned no chains or no frame boxes.");
        }

        if (chain_index >=
            trajectory.chains.size())
        {
            throw std::out_of_range(
                "Requested CHAIN_INDEX exceeds parsed chain count.");
        }

        SetFixedZBoxCenter(
            trajectory.frame_boxes,
            z_center
        );

        CheckFixedZCenterMapping(
            trajectory.frame_boxes,
            z_center
        );

        const std::size_t num_frames =
            trajectory.frame_boxes.size();

        AssignUniformTimesteps(
            trajectory.chains,
            0,
            frame_time
        );

        const ChainTrajectory& chain =
            trajectory.chains[chain_index];

        if (chain.frame_offsets.size() !=
            num_frames + 1)
        {
            throw std::runtime_error(
                "Selected chain frame count does not match box count.");
        }

        // Keep the real-history smoke test intentionally modest.
        // lag 80 = 1600 tau_LJ when FRAME_TIME=20, i.e. close to
        // one full 1750-tau_LJ oscillation.
        const std::vector<std::size_t> lags =
            KeepValidLags(
                {0, 1, 2, 5, 10, 20, 40, 80},
                num_frames
            );

        const std::vector<double> diameters{
            tube_diameter
        };

        const SegmentSurvivalOptions affine_options{
            .apply_affine_correction = true
        };

        std::cout
            << "Parsed chains    : "
            << trajectory.chains.size()
            << '\n'
            << "Parsed frames    : "
            << num_frames
            << '\n'
            << "Selected chain   : "
            << chain_index
            << " (ID "
            << chain.chain_id
            << ")\n"
            << "Tube diameter    : "
            << tube_diameter
            << '\n'
            << "Fixed z center   : "
            << std::setprecision(15)
            << z_center
            << '\n'
            << "First-frame zlo  : "
            << trajectory.frame_boxes
                .front().origin.z
            << "\n\n";

        std::cout
            << "Computing raw instantaneous survival...\n";

        const SegmentSurvivalFunction raw =
            ComputeSegmentSurvivalFunction(
                chain,
                lags,
                diameters
            );

        std::cout
            << "Computing affine instantaneous survival...\n";

        const SegmentSurvivalFunction affine =
            ComputeSegmentSurvivalFunction(
                chain,
                trajectory.frame_boxes,
                lags,
                diameters,
                affine_options
            );

        std::cout
            << "Computing affine history-dependent survival...\n";

        const HistoryDependentSurvivalResult history_result =
            ComputeHistoryDependentSurvivalFunction(
                chain,
                trajectory.frame_boxes,
                lags,
                diameters,
                affine_options,
                false
            );

        const SegmentSurvivalFunction& history =
            history_result.stephanou_survival;

        const SegmentSurvivalFunction& permanent =
            history_result.permanent_escape_survival;

        CheckHistorySubset(
            affine,
            history
        );

        CheckPermanentContainsHistory(
            history,
            permanent
        );

        if (raw.sample_counts !=
            affine.sample_counts ||
            affine.sample_counts !=
            history.sample_counts)
        {
            throw std::runtime_error(
                "The three estimators used different lag-wise origin cohorts.");
        }

        const TubeSurvivalFunction raw_tube =
            ComputeTubeSurvivalFunction(raw);

        const TubeSurvivalFunction affine_tube =
            ComputeTubeSurvivalFunction(affine);

        const TubeSurvivalFunction history_tube =
            ComputeTubeSurvivalFunction(history);

        const TubeSurvivalFunction permanent_tube =
            ComputeTubeSurvivalFunction(permanent);

        CheckLagZero(
            raw_tube,
            "raw tube survival"
        );

        CheckLagZero(
            affine_tube,
            "affine instantaneous tube survival"
        );

        CheckLagZero(
            history_tube,
            "history-dependent tube survival"
        );

        CheckLagZero(
            permanent_tube,
            "permanent-escape survival"
        );

        WriteSegmentSurvivalFunction(
            raw,
            "history_smoke_segment_raw.csv"
        );

        WriteSegmentSurvivalFunction(
            affine,
            "history_smoke_segment_affine.csv"
        );

        WriteSegmentSurvivalFunction(
            history,
            "history_smoke_segment_history.csv"
        );

        WriteSegmentSurvivalFunction(
            permanent,
            "history_smoke_segment_permanent_escape.csv"
        );

        WriteTubeSurvivalFunction(
            raw_tube,
            "history_smoke_tube_raw.csv"
        );

        WriteTubeSurvivalFunction(
            affine_tube,
            "history_smoke_tube_affine.csv"
        );

        WriteTubeSurvivalFunction(
            history_tube,
            "history_smoke_tube_history.csv"
        );

        WriteTubeSurvivalFunction(
            permanent_tube,
            "history_smoke_tube_permanent_escape.csv"
        );

        std::cout
            << "\nlag  time      raw      affine   history  permanent\n";

        for (std::size_t lag_index = 0;
             lag_index < lags.size();
             ++lag_index)
        {
            std::cout
                << std::setw(3)
                << lags[lag_index]
                << "  "
                << std::setw(7)
                << raw_tube.lag_times[lag_index]
                << "  "
                << std::fixed
                << std::setprecision(6)
                << raw_tube.survival[lag_index]
                << "  "
                << affine_tube.survival[lag_index]
                << "  "
                << history_tube.survival[lag_index]
                << "  "
                << permanent_tube.survival[lag_index]
                << '\n';
        }

        std::cout
            << "\nALL REAL HISTORY-ESCAPE SMOKE CHECKS PASSED\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "\nFAILED: "
            << error.what()
            << '\n';

        return 1;
    }
}
CPP

cat > "$PLOT_FILE" <<'PY'
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
PY

echo "3) Compiling real-trajectory three-way smoke test..."
"$CXX" \
    -std=c++20 \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -I. \
    "$REAL_CPP" \
    Parse_Z1_File.cpp \
    Parser_Utils.cpp \
    Primitive_Path_Trajectory.cpp \
    Geometry_Utils.cpp \
    Trajectory_Time.cpp \
    Tube_Survival.cpp \
    Survival_IO.cpp \
    -o "$REAL_EXE"

echo "4) Running real-trajectory three-way smoke test..."
./"$REAL_EXE" \
    "$INPUT_FILE" \
    "$FRAME_TIME" \
    "$FOCUS_DIAMETER" \
    "$CHAIN_INDEX" \
    "$BOX_CENTER_Z"

echo
echo "5) Plotting..."
python3 "$PLOT_FILE" \
    "$FRAME_TIME" \
    "$FOCUS_DIAMETER"

echo
echo "Done. Main outputs:"
echo "  history_escape_three_way_comparison.png"
echo "  history_escape_contour_comparison.png"
echo "  history_smoke_tube_raw.csv"
echo "  history_smoke_tube_affine.csv"
echo "  history_smoke_tube_history.csv"
echo "  history_smoke_tube_permanent_escape.csv"
