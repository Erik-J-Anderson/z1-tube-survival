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
                history_result.reference_to_future.transverse_survival;
            const SegmentSurvivalFunction& permanent =
                history_result.reference_to_future.permanent_escape_survival;

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
