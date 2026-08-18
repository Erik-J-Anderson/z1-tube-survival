#!/usr/bin/env bash
set -euo pipefail

CXX=/usr/bin/g++

INPUT_FILE="${1:-Z1+SP.dat}"
FRAME_TIME="${2:-20}"
FOCUS_DIAMETER="${3:-5.0}"
BOX_CENTER_Z="${BOX_CENTER_Z:-}"
MAX_CHAINS="${MAX_CHAINS:-0}"

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
    Survival_IO.cpp Survival_IO.hpp
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

CPP_FILE="Serial_All_Chains_History_Escape.cpp"
EXE="Serial_All_Chains_History_Escape"
PLOT_FILE="plot_all_chains_history_escape.py"

cat > "$CPP_FILE" <<'CPP'
#include "Geometry_Utils.hpp"
#include "Parse_Z1_File.hpp"
#include "Survival_IO.hpp"
#include "Trajectory_Time.hpp"
#include "Tube_Survival.hpp"

#include <algorithm>
#include <chrono>
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
        throw std::runtime_error("No requested lag is valid for this trajectory.");
    }
    return result;
}


double GetBoxLengthZ(const Box& box)
{
    const double xz = box.matrix.value[0][2];
    const double yz = box.matrix.value[1][2];
    const double zz = box.matrix.value[2][2];
    return std::sqrt(xz * xz + yz * yz + zz * zz);
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
        const double Lz = GetBoxLengthZ(box);
        if (!(Lz > 0.0)) {
            throw std::runtime_error(
                "SetFixedZBoxCenter: encountered non-positive Lz.");
        }
        box.origin.z = z_center - 0.5 * Lz;
    }
}


void CheckFixedZCenterMapping(
    const std::vector<Box>& frame_boxes,
    double z_center)
{
    const Vec3 physical_center{0.0, 0.0, z_center};
    const Box& destination_box = frame_boxes.front();
    constexpr double tolerance = 1.0e-10;

    for (const Box& source_box : frame_boxes)
    {
        const Vec3 mapped = geometry::MapPositionBetweenBoxes(
            physical_center,
            source_box,
            destination_box);

        if (std::abs(mapped.z - z_center) > tolerance) {
            throw std::runtime_error(
                "Fixed-z-center affine regression check failed.");
        }
    }
}


void InitializeSegmentAccumulator(
    SegmentSurvivalFunction& dst,
    const SegmentSurvivalFunction& src)
{
    dst = src;
    std::fill(dst.survival.begin(), dst.survival.end(), 0.0);
    std::fill(dst.sample_counts.begin(), dst.sample_counts.end(), 0);
}


void AccumulateSegmentWeighted(
    SegmentSurvivalFunction& dst,
    const SegmentSurvivalFunction& src)
{
    const std::size_t num_lags = src.lag_times.size();
    const std::size_t num_diameters = src.tube_diameters.size();

    if (dst.lag_times != src.lag_times ||
        dst.tube_diameters != src.tube_diameters ||
        dst.survival.size() != src.survival.size() ||
        dst.sample_counts.size() != src.sample_counts.size())
    {
        throw std::runtime_error(
            "Segment accumulator dimensions changed between chains.");
    }

    for (std::size_t lag_index = 0;
         lag_index < num_lags;
         ++lag_index)
    {
        const std::uint64_t count = src.sample_counts[lag_index];
        const double weight = static_cast<double>(count);
        dst.sample_counts[lag_index] += count;

        for (std::size_t d = 0; d < num_diameters; ++d)
        {
            for (std::size_t s = 0; s < NUM_SAMPLE_POINTS; ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS +
                    lag_index * NUM_SAMPLE_POINTS +
                    s;

                dst.survival[index] +=
                    src.survival[index] * weight;
            }
        }
    }
}


void FinalizeSegmentAccumulator(
    SegmentSurvivalFunction& result)
{
    const std::size_t num_lags = result.lag_times.size();
    const std::size_t num_diameters = result.tube_diameters.size();

    for (std::size_t lag_index = 0;
         lag_index < num_lags;
         ++lag_index)
    {
        const std::uint64_t count = result.sample_counts[lag_index];
        if (count == 0) {
            continue;
        }

        const double denom = static_cast<double>(count);

        for (std::size_t d = 0; d < num_diameters; ++d)
        {
            for (std::size_t s = 0; s < NUM_SAMPLE_POINTS; ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS +
                    lag_index * NUM_SAMPLE_POINTS +
                    s;

                result.survival[index] /= denom;
            }
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

    constexpr double tolerance = 1.0e-12;
    for (std::size_t i = 0; i < affine.survival.size(); ++i)
    {
        if (history.survival[i] > affine.survival[i] + tolerance) {
            throw std::runtime_error(
                "History-dependent survival exceeded instantaneous affine survival.");
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

    constexpr double tolerance = 1.0e-12;
    for (std::size_t i = 0; i < history.survival.size(); ++i)
    {
        if (history.survival[i] > permanent.survival[i] + tolerance) {
            throw std::runtime_error(
                "History occupancy exceeded pure permanent-escape survival.");
        }
    }
}


void CheckLagZero(
    const TubeSurvivalFunction& result,
    const char* label)
{
    if (result.lag_times.empty() || result.lag_times.front() != 0.0) {
        throw std::runtime_error(
            std::string(label) + ": lag-zero entry missing.");
    }

    const std::size_t num_lags = result.lag_times.size();
    for (std::size_t d = 0;
         d < result.tube_diameters.size();
         ++d)
    {
        const double value = result.survival[d * num_lags];
        if (std::abs(value - 1.0) > 1.0e-12) {
            throw std::runtime_error(
                std::string(label) + ": survival(0) != 1.");
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
                << "Usage: " << argv[0]
                << " <Z1+SP.dat> <frame_time> <tube_diameter> "
                << "<z_center> <max_chains=0_for_all>\n";
            return 2;
        }

        const std::string filename = argv[1];
        const long frame_time = std::stol(argv[2]);
        const double tube_diameter = std::stod(argv[3]);
        const double z_center = std::stod(argv[4]);
        const std::size_t max_chains =
            static_cast<std::size_t>(std::stoull(argv[5]));

        if (frame_time <= 0) {
            throw std::invalid_argument("frame_time must be positive.");
        }
        if (tube_diameter <= 0.0) {
            throw std::invalid_argument("tube_diameter must be positive.");
        }

        std::cout
            << "==============================================\n"
            << "ALL-CHAIN HISTORY-DEPENDENT TUBE COMPARISON\n"
            << "==============================================\n\n"
            << "Parsing " << filename << " ...\n";

        const auto parse_start = std::chrono::steady_clock::now();
        PrimitivePathTrajectory trajectory = parse_z1_file(filename);
        const auto parse_stop = std::chrono::steady_clock::now();

        if (trajectory.chains.empty() || trajectory.frame_boxes.empty()) {
            throw std::runtime_error(
                "Parser returned no chains or no frame boxes.");
        }

        SetFixedZBoxCenter(trajectory.frame_boxes, z_center);
        CheckFixedZCenterMapping(trajectory.frame_boxes, z_center);

        const std::size_t num_frames = trajectory.frame_boxes.size();
        AssignUniformTimesteps(trajectory.chains, 0, frame_time);

        for (const ChainTrajectory& chain : trajectory.chains)
        {
            if (chain.frame_offsets.size() != num_frames + 1) {
                throw std::runtime_error(
                    "Chain frame count does not match box count.");
            }
        }

        const std::size_t chains_to_run =
            (max_chains == 0)
                ? trajectory.chains.size()
                : std::min(max_chains, trajectory.chains.size());

        // Same output lags used by the one-chain smoke test.
        // With frame_time=20, lag 80 = 1600 tau_LJ, close to one
        // 1750-tau_LJ oscillation period.
        const std::vector<std::size_t> lags = KeepValidLags(
            {0, 1, 2, 5, 10, 20, 40, 80},
            num_frames);

        const std::vector<double> diameters{tube_diameter};

        const SegmentSurvivalOptions affine_options{
            .apply_affine_correction = true
        };

        std::cout
            << "Parsed chains      : " << trajectory.chains.size() << '\n'
            << "Chains to analyze  : " << chains_to_run << '\n'
            << "Frames/chain       : " << num_frames << '\n'
            << "Frame time         : " << frame_time << '\n'
            << "Tube diameter      : " << tube_diameter << '\n'
            << "Fixed z center     : " << std::setprecision(15)
            << z_center << '\n'
            << "First-frame zlo    : "
            << trajectory.frame_boxes.front().origin.z << '\n'
            << "Requested lags     : ";

        for (const auto lag : lags) {
            std::cout << lag << ' ';
        }
        std::cout << "\n\n";

        SegmentSurvivalFunction ensemble_raw;
        SegmentSurvivalFunction ensemble_affine;
        SegmentSurvivalFunction ensemble_history;
        SegmentSurvivalFunction ensemble_permanent;

        bool initialized = false;
        const auto work_start = std::chrono::steady_clock::now();

        for (std::size_t chain_index = 0;
             chain_index < chains_to_run;
             ++chain_index)
        {
            const ChainTrajectory& chain = trajectory.chains[chain_index];

            const SegmentSurvivalFunction raw =
                ComputeSegmentSurvivalFunction(
                    chain,
                    lags,
                    diameters);

            const SegmentSurvivalFunction affine =
                ComputeSegmentSurvivalFunction(
                    chain,
                    trajectory.frame_boxes,
                    lags,
                    diameters,
                    affine_options);

            const HistoryDependentSurvivalResult history_result =
                ComputeHistoryDependentSurvivalFunction(
                    chain,
                    trajectory.frame_boxes,
                    lags,
                    diameters,
                    affine_options,
                    false);

            const SegmentSurvivalFunction& history =
                history_result.stephanou_survival;
            const SegmentSurvivalFunction& permanent =
                history_result.permanent_escape_survival;

            CheckHistorySubset(affine, history);
            CheckPermanentContainsHistory(history, permanent);

            if (raw.sample_counts != affine.sample_counts ||
                affine.sample_counts != history.sample_counts ||
                history.sample_counts != permanent.sample_counts)
            {
                throw std::runtime_error(
                    "Estimators used different lag-wise origin cohorts.");
            }

            if (!initialized)
            {
                InitializeSegmentAccumulator(ensemble_raw, raw);
                InitializeSegmentAccumulator(ensemble_affine, affine);
                InitializeSegmentAccumulator(ensemble_history, history);
                InitializeSegmentAccumulator(ensemble_permanent, permanent);
                initialized = true;
            }

            AccumulateSegmentWeighted(ensemble_raw, raw);
            AccumulateSegmentWeighted(ensemble_affine, affine);
            AccumulateSegmentWeighted(ensemble_history, history);
            AccumulateSegmentWeighted(ensemble_permanent, permanent);

            const std::size_t done = chain_index + 1;
            if (done == 1 || done % 25 == 0 || done == chains_to_run)
            {
                const auto now = std::chrono::steady_clock::now();
                const std::chrono::duration<double> elapsed = now - work_start;

                std::cout
                    << "Finished chain " << done << " / " << chains_to_run
                    << "  (elapsed " << std::fixed << std::setprecision(1)
                    << elapsed.count() << " s)\n";
            }
        }

        if (!initialized) {
            throw std::runtime_error("No chains were analyzed.");
        }

        FinalizeSegmentAccumulator(ensemble_raw);
        FinalizeSegmentAccumulator(ensemble_affine);
        FinalizeSegmentAccumulator(ensemble_history);
        FinalizeSegmentAccumulator(ensemble_permanent);

        CheckHistorySubset(ensemble_affine, ensemble_history);
        CheckPermanentContainsHistory(ensemble_history, ensemble_permanent);

        const TubeSurvivalFunction tube_raw =
            ComputeTubeSurvivalFunction(ensemble_raw);
        const TubeSurvivalFunction tube_affine =
            ComputeTubeSurvivalFunction(ensemble_affine);
        const TubeSurvivalFunction tube_history =
            ComputeTubeSurvivalFunction(ensemble_history);
        const TubeSurvivalFunction tube_permanent =
            ComputeTubeSurvivalFunction(ensemble_permanent);

        CheckLagZero(tube_raw, "raw tube survival");
        CheckLagZero(tube_affine, "affine instantaneous tube survival");
        CheckLagZero(tube_history, "history-dependent tube survival");
        CheckLagZero(tube_permanent, "permanent-escape survival");

        WriteSegmentSurvivalFunction(
            ensemble_raw,
            "all_chains_history_segment_raw.csv");
        WriteSegmentSurvivalFunction(
            ensemble_affine,
            "all_chains_history_segment_affine.csv");
        WriteSegmentSurvivalFunction(
            ensemble_history,
            "all_chains_history_segment_history.csv");
        WriteSegmentSurvivalFunction(
            ensemble_permanent,
            "all_chains_history_segment_permanent_escape.csv");

        WriteTubeSurvivalFunction(
            tube_raw,
            "all_chains_history_tube_raw.csv");
        WriteTubeSurvivalFunction(
            tube_affine,
            "all_chains_history_tube_affine.csv");
        WriteTubeSurvivalFunction(
            tube_history,
            "all_chains_history_tube_history.csv");
        WriteTubeSurvivalFunction(
            tube_permanent,
            "all_chains_history_tube_permanent_escape.csv");

        std::cout
            << "\nlag   time        raw       affine    history   permanent\n";

        for (std::size_t i = 0; i < lags.size(); ++i)
        {
            std::cout
                << std::setw(3) << lags[i] << "  "
                << std::setw(8) << tube_raw.lag_times[i] << "  "
                << std::fixed << std::setprecision(6)
                << tube_raw.survival[i] << "  "
                << tube_affine.survival[i] << "  "
                << tube_history.survival[i] << "  "
                << tube_permanent.survival[i] << '\n';
        }

        const std::chrono::duration<double> parse_elapsed =
            parse_stop - parse_start;
        const std::chrono::duration<double> total_elapsed =
            std::chrono::steady_clock::now() - work_start;

        std::cout
            << "\nParsing time      : " << parse_elapsed.count() << " s\n"
            << "Calculation time  : " << total_elapsed.count() << " s\n"
            << "Total origins lag0: " << tube_raw.sample_counts.front() << "\n"
            << "\nALL-CHAIN HISTORY COMPARISON COMPLETED SUCCESSFULLY\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "\nFAILED: " << error.what() << '\n';
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
PY

echo "Using compiler:"
"$CXX" --version | head -n 1
echo
echo "Input file       : $INPUT_FILE"
echo "Frame time       : $FRAME_TIME tau_LJ"
echo "Focus diameter   : $FOCUS_DIAMETER"
echo "Fixed z center   : $BOX_CENTER_Z"
echo "Max chains       : $MAX_CHAINS (0 = all)"
echo

echo "Compiling all-chain history comparison..."
"$CXX" \
    -std=c++20 \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -I. \
    "$CPP_FILE" \
    Parse_Z1_File.cpp \
    Parser_Utils.cpp \
    Primitive_Path_Trajectory.cpp \
    Geometry_Utils.cpp \
    Trajectory_Time.cpp \
    Tube_Survival.cpp \
    Survival_IO.cpp \
    -o "$EXE"

echo
echo "Running all-chain history comparison..."
./"$EXE" \
    "$INPUT_FILE" \
    "$FRAME_TIME" \
    "$FOCUS_DIAMETER" \
    "$BOX_CENTER_Z" \
    "$MAX_CHAINS"

echo
echo "Plotting..."
python3 "$PLOT_FILE" "$FOCUS_DIAMETER"

echo
echo "Done. Main outputs:"
echo "  all_chains_history_escape_dashboard.png"
echo "  all_chains_history_contour_comparison.png"
echo "  all_chains_permanent_escape_profile.png"
echo "  all_chains_history_tube_raw.csv"
echo "  all_chains_history_tube_affine.csv"
echo "  all_chains_history_tube_history.csv"
echo "  all_chains_history_tube_permanent_escape.csv"
