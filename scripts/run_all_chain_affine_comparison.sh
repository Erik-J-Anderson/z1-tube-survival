#!/usr/bin/env bash
set -euo pipefail

CXX=/usr/bin/g++
INPUT_FILE="${1:-Z1+SP.dat}"
FRAME_TIME="${2:-20}"
FOCUS_DIAMETER="${3:-5.0}"
RUN_END_RETRACTION="${RUN_END_RETRACTION:-1}"
BOX_CENTER_Z="${BOX_CENTER_Z:-}"

if [[ -z "$BOX_CENTER_Z" ]]; then
    echo "ERROR: BOX_CENTER_Z is required for the midpoint-corrected affine test." >&2
    echo "Example:" >&2
    echo "  BOX_CENTER_Z=39.22658609485464 bash $0 Z1+SP.dat 20 5.0" >&2
    exit 2
fi

CPP_FILE="Serial_All_Chains_Affine_Comparison.cpp"
EXE="Serial_All_Chains_Affine_Comparison"
PLOT_FILE="plot_all_chain_affine_comparison.py"

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
        echo "ERROR: run this script from the project root. Missing: $f" >&2
        exit 1
    fi
 done

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: input file not found: $INPUT_FILE" >&2
    exit 1
fi

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
        dst.survival.size() != src.survival.size()) {
        throw std::runtime_error("Segment accumulator dimensions changed between chains.");
    }

    for (std::size_t lag_index = 0; lag_index < num_lags; ++lag_index)
    {
        const std::uint64_t count = src.sample_counts[lag_index];
        dst.sample_counts[lag_index] += count;
        const double weight = static_cast<double>(count);

        for (std::size_t d = 0; d < num_diameters; ++d)
        {
            for (std::size_t s = 0; s < NUM_SAMPLE_POINTS; ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS +
                    lag_index * NUM_SAMPLE_POINTS + s;
                dst.survival[index] += src.survival[index] * weight;
            }
        }
    }
}

void FinalizeSegmentAccumulator(SegmentSurvivalFunction& result)
{
    const std::size_t num_lags = result.lag_times.size();
    const std::size_t num_diameters = result.tube_diameters.size();

    for (std::size_t lag_index = 0; lag_index < num_lags; ++lag_index)
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
                    lag_index * NUM_SAMPLE_POINTS + s;
                result.survival[index] /= denom;
            }
        }
    }
}

void InitializeEndAccumulator(
    EndRetractionFunction& dst,
    const EndRetractionFunction& src)
{
    dst = src;
    std::fill(dst.end0_reached.begin(), dst.end0_reached.end(), 0.0);
    std::fill(dst.end1_reached.begin(), dst.end1_reached.end(), 0.0);
    std::fill(dst.combined_reached.begin(), dst.combined_reached.end(), 0.0);
    std::fill(dst.mean_max_depth_end0.begin(), dst.mean_max_depth_end0.end(), 0.0);
    std::fill(dst.mean_max_depth_end1.begin(), dst.mean_max_depth_end1.end(), 0.0);
    std::fill(dst.mean_max_depth_combined.begin(), dst.mean_max_depth_combined.end(), 0.0);
    std::fill(dst.sample_counts.begin(), dst.sample_counts.end(), 0);
}

void AccumulateEndWeighted(
    EndRetractionFunction& dst,
    const EndRetractionFunction& src)
{
    const std::size_t num_diameters = src.tube_diameters.size();
    const std::size_t num_lags = src.lag_frames.size();
    const std::size_t num_s = src.s_fraction.size();

    if (dst.tube_diameters != src.tube_diameters ||
        dst.lag_frames != src.lag_frames ||
        dst.s_fraction != src.s_fraction) {
        throw std::runtime_error("End-retraction accumulator dimensions changed between chains.");
    }

    for (std::size_t d = 0; d < num_diameters; ++d)
    {
        for (std::size_t lag_index = 0; lag_index < num_lags; ++lag_index)
        {
            const std::size_t summary = d * num_lags + lag_index;
            const std::uint64_t count = src.sample_counts[summary];
            const double weight = static_cast<double>(count);

            dst.sample_counts[summary] += count;
            dst.mean_max_depth_end0[summary] +=
                src.mean_max_depth_end0[summary] * weight;
            dst.mean_max_depth_end1[summary] +=
                src.mean_max_depth_end1[summary] * weight;
            dst.mean_max_depth_combined[summary] +=
                src.mean_max_depth_combined[summary] * weight;

            for (std::size_t s = 0; s < num_s; ++s)
            {
                const std::size_t field =
                    (d * num_lags + lag_index) * num_s + s;
                dst.end0_reached[field] += src.end0_reached[field] * weight;
                dst.end1_reached[field] += src.end1_reached[field] * weight;
                dst.combined_reached[field] += src.combined_reached[field] * weight;
            }
        }
    }
}

void FinalizeEndAccumulator(EndRetractionFunction& result)
{
    const std::size_t num_diameters = result.tube_diameters.size();
    const std::size_t num_lags = result.lag_frames.size();
    const std::size_t num_s = result.s_fraction.size();

    for (std::size_t d = 0; d < num_diameters; ++d)
    {
        for (std::size_t lag_index = 0; lag_index < num_lags; ++lag_index)
        {
            const std::size_t summary = d * num_lags + lag_index;
            const std::uint64_t count = result.sample_counts[summary];
            if (count == 0) {
                continue;
            }
            const double denom = static_cast<double>(count);

            result.mean_max_depth_end0[summary] /= denom;
            result.mean_max_depth_end1[summary] /= denom;
            result.mean_max_depth_combined[summary] /= denom;

            for (std::size_t s = 0; s < num_s; ++s)
            {
                const std::size_t field =
                    (d * num_lags + lag_index) * num_s + s;
                result.end0_reached[field] /= denom;
                result.end1_reached[field] /= denom;
                result.combined_reached[field] /= denom;
            }
        }
    }
}

void CheckLagZero(const TubeSurvivalFunction& result, const char* label)
{
    if (result.lag_times.empty() || result.lag_times.front() != 0.0) {
        throw std::runtime_error(std::string(label) + ": lag-zero entry missing.");
    }

    const std::size_t num_lags = result.lag_times.size();
    for (std::size_t d = 0; d < result.tube_diameters.size(); ++d)
    {
        const double mu0 = result.survival[d * num_lags];
        if (std::abs(mu0 - 1.0) > 1.0e-12) {
            throw std::runtime_error(std::string(label) + ": mu(0) != 1.");
        }
    }
}


double GetBoxLengthZ(const Box& box)
{
    // Z1+ is orthorhombic for this dataset.  Use the length of
    // the third box-vector column so this also remains valid if
    // the storage ever contains nonzero xz/yz components.
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

        // The Z1+ 3-value box record stores only Lx Ly Lz, so the
        // parser necessarily sets origin=(0,0,0).  For this LAMMPS
        // fix-deform z wiggle trajectory, the physical box center
        // is fixed and the lower z boundary moves as
        //
        //     zlo(t) = z_center - Lz(t)/2.
        //
        // Supply that missing origin information before the
        // source-box -> time-origin-box affine pullback.
        box.origin.z = z_center - 0.5 * Lz;
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

    const Box& destination_box = frame_boxes.front();

    constexpr double tolerance = 1.0e-10;

    for (const Box& source_box : frame_boxes)
    {
        const Vec3 mapped_center =
            geometry::MapPositionBetweenBoxes(
                physical_center,
                source_box,
                destination_box);

        if (std::abs(mapped_center.z - z_center) > tolerance) {
            throw std::runtime_error(
                "Fixed-z-center affine regression check failed.");
        }
    }
}

} // namespace

int main(int argc, char* argv[])
{
    try
    {
        if (argc != 5)
        {
            std::cerr
                << "Usage: " << argv[0]
                << " <Z1+SP.dat> <frame_time> <run_end_retraction> <z_center>\n";
            return 2;
        }

        const std::string filename = argv[1];
        const long frame_time = std::stol(argv[2]);
        const bool run_end_retraction =
            (std::stoi(argv[3]) != 0);
        const double z_center = std::stod(argv[4]);

        if (frame_time <= 0) {
            throw std::invalid_argument("frame_time must be positive.");
        }

        std::cout
            << "=============================================\n"
            << "ALL-CHAIN AFFINE / END-RETRACTION COMPARISON\n"
            << "=============================================\n\n"
            << "Parsing " << filename << " ...\n";

        const auto parse_start = std::chrono::steady_clock::now();
        PrimitivePathTrajectory trajectory = parse_z1_file(filename);
        const auto parse_stop = std::chrono::steady_clock::now();

        if (trajectory.chains.empty() || trajectory.frame_boxes.empty()) {
            throw std::runtime_error("Parser returned no chains or no frame boxes.");
        }

        SetFixedZBoxCenter(
            trajectory.frame_boxes,
            z_center);

        CheckFixedZCenterMapping(
            trajectory.frame_boxes,
            z_center);

        const std::size_t num_frames = trajectory.frame_boxes.size();
        for (const ChainTrajectory& chain : trajectory.chains) {
            if (chain.frame_offsets.size() != num_frames + 1) {
                throw std::runtime_error("Chain frame count does not match box count.");
            }
        }

        AssignUniformTimesteps(trajectory.chains, 0, frame_time);

        // Tube/segment comparison: extends through several parts of the cycle.
        const std::vector<std::size_t> survival_lags = KeepValidLags(
            {0, 1, 2, 5, 10, 20, 40, 80, 160, 320}, num_frames);

        // Several tube diameters for robustness.  Plotting defaults to a=5.
        const std::vector<double> survival_diameters{
            3.0, 5.0, 7.0, 9.0
        };

        // The history scan is much more expensive, so keep its maximum lag
        // smaller for this all-chain diagnostic.
        const std::vector<std::size_t> end_lags = KeepValidLags(
            {0, 5, 10, 20, 40, 80}, num_frames);
        const std::vector<double> end_diameters{5.0};

        std::cout
            << "Chains       : " << trajectory.chains.size() << '\n'
            << "Frames/chain : " << num_frames << '\n'
            << "Frame time   : " << frame_time << '\n'
            << "Fixed z center: " << std::setprecision(15) << z_center << '\n'
            << "First-frame zlo: " << trajectory.frame_boxes.front().origin.z << '\n'
            << "Survival lags: ";
        for (const auto lag : survival_lags) std::cout << lag << ' ';
        std::cout << "\nEnd lags     : ";
        for (const auto lag : end_lags) std::cout << lag << ' ';
        std::cout << "\n\n";

        SegmentSurvivalFunction ensemble_raw;
        SegmentSurvivalFunction ensemble_affine;
        EndRetractionFunction ensemble_end;
        bool raw_initialized = false;
        bool affine_initialized = false;
        bool end_initialized = false;

        const SegmentSurvivalOptions affine_options{
            .apply_affine_correction = true
        };

        const auto work_start = std::chrono::steady_clock::now();

        for (std::size_t chain_index = 0;
             chain_index < trajectory.chains.size();
             ++chain_index)
        {
            const ChainTrajectory& chain = trajectory.chains[chain_index];

            const SegmentSurvivalFunction raw =
                ComputeSegmentSurvivalFunction(
                    chain,
                    survival_lags,
                    survival_diameters);

            const SegmentSurvivalFunction affine =
                ComputeSegmentSurvivalFunction(
                    chain,
                    trajectory.frame_boxes,
                    survival_lags,
                    survival_diameters,
                    affine_options);

            if (!raw_initialized) {
                InitializeSegmentAccumulator(ensemble_raw, raw);
                raw_initialized = true;
            }
            if (!affine_initialized) {
                InitializeSegmentAccumulator(ensemble_affine, affine);
                affine_initialized = true;
            }

            AccumulateSegmentWeighted(ensemble_raw, raw);
            AccumulateSegmentWeighted(ensemble_affine, affine);

            if (run_end_retraction)
            {
                const EndRetractionFunction end =
                    ComputeEndRetractionFunction(
                        chain,
                        end_lags,
                        end_diameters);

                if (!end_initialized) {
                    InitializeEndAccumulator(ensemble_end, end);
                    end_initialized = true;
                }
                AccumulateEndWeighted(ensemble_end, end);
            }

            const std::size_t done = chain_index + 1;
            if (done == 1 || done % 25 == 0 || done == trajectory.chains.size())
            {
                const auto now = std::chrono::steady_clock::now();
                const std::chrono::duration<double> elapsed = now - work_start;
                std::cout
                    << "Finished chain " << done << " / "
                    << trajectory.chains.size()
                    << "  (elapsed " << std::fixed << std::setprecision(1)
                    << elapsed.count() << " s)\n";
            }
        }

        FinalizeSegmentAccumulator(ensemble_raw);
        FinalizeSegmentAccumulator(ensemble_affine);
        if (run_end_retraction) {
            FinalizeEndAccumulator(ensemble_end);
        }

        const TubeSurvivalFunction tube_raw =
            ComputeTubeSurvivalFunction(ensemble_raw);
        const TubeSurvivalFunction tube_affine =
            ComputeTubeSurvivalFunction(ensemble_affine);

        CheckLagZero(tube_raw, "uncorrected tube survival");
        CheckLagZero(tube_affine, "affine-corrected tube survival");

        if (tube_raw.sample_counts != tube_affine.sample_counts) {
            throw std::runtime_error(
                "Raw and affine calculations used different total origin counts.");
        }

        WriteSegmentSurvivalFunction(
            ensemble_raw, "all_chains_segment_uncorrected.csv");
        WriteSegmentSurvivalFunction(
            ensemble_affine, "all_chains_segment_affine.csv");
        WriteTubeSurvivalFunction(
            tube_raw, "all_chains_tube_uncorrected.csv");
        WriteTubeSurvivalFunction(
            tube_affine, "all_chains_tube_affine.csv");

        if (run_end_retraction) {
            WriteEndRetractionFunction(
                ensemble_end, "all_chains_end_retraction.csv");
        }

        double max_delta = 0.0;
        for (std::size_t i = 0; i < tube_raw.survival.size(); ++i) {
            max_delta = std::max(
                max_delta,
                std::abs(tube_affine.survival[i] - tube_raw.survival[i]));
        }

        const std::chrono::duration<double> parse_elapsed = parse_stop - parse_start;
        const std::chrono::duration<double> total_elapsed =
            std::chrono::steady_clock::now() - work_start;

        std::cout
            << "\nParsing time                         : "
            << parse_elapsed.count() << " s\n"
            << "Calculation time                     : "
            << total_elapsed.count() << " s\n"
            << "Max |affine - raw| tube survival     : "
            << max_delta << "\n"
            << "Total origins at lag 0               : "
            << tube_raw.sample_counts.front() << "\n\n"
            << "Wrote:\n"
            << "  all_chains_segment_uncorrected.csv\n"
            << "  all_chains_segment_affine.csv\n"
            << "  all_chains_tube_uncorrected.csv\n"
            << "  all_chains_tube_affine.csv\n";

        if (run_end_retraction) {
            std::cout << "  all_chains_end_retraction.csv\n";
        }

        std::cout
            << "\nNOTE: tube/segment affine correction uses the supplied fixed z center.\n"
            << "      The current production end-retraction routine is still\n"
            << "      uncorrected for affine box deformation.  This test plots\n"
            << "      it as a separate history-dependent diagnostic; it is not\n"
            << "      being labeled as affine-corrected.\n\n"
            << "ALL-CHAIN COMPARISON COMPLETED SUCCESSFULLY\n";

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
PY

echo "Using compiler:"
"$CXX" --version | head -n 1
echo
echo "Input file      : $INPUT_FILE"
echo "Frame time      : $FRAME_TIME tau_LJ"
echo "Focus diameter  : $FOCUS_DIAMETER"
echo "End retraction  : $RUN_END_RETRACTION"
echo "Fixed z center  : $BOX_CENTER_Z"
echo

echo "Compiling all-chain comparison test..."
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
echo "Running all-chain comparison..."
./"$EXE" "$INPUT_FILE" "$FRAME_TIME" "$RUN_END_RETRACTION" "$BOX_CENTER_Z"

echo
echo "Plotting..."
python3 "$PLOT_FILE" "$FRAME_TIME" "$FOCUS_DIAMETER"

echo
echo "Done. Main outputs:"
echo "  all_chains_tube_uncorrected.csv"
echo "  all_chains_tube_affine.csv"
echo "  all_chains_segment_uncorrected.csv"
echo "  all_chains_segment_affine.csv"
if [[ "$RUN_END_RETRACTION" != "0" ]]; then
    echo "  all_chains_end_retraction.csv"
fi
echo "  all_chains_affine_end_retraction_dashboard.png"
echo "  all_chains_segment_profile_comparison.png"
