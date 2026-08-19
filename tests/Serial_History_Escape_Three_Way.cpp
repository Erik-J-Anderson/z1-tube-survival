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
            history_result.reference_to_future.transverse_survival;

        const SegmentSurvivalFunction& permanent =
            history_result.reference_to_future.permanent_escape_survival;

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
