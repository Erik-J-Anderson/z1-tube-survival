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
