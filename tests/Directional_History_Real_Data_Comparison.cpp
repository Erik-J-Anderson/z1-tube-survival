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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

std::vector<std::size_t> ParseLagList(const std::string& text)
{
    std::vector<std::size_t> lags;
    std::stringstream stream(text);
    std::string token;

    while (std::getline(stream, token, ','))
    {
        if (token.empty()) {
            continue;
        }

        const auto value =
            static_cast<std::size_t>(std::stoull(token));

        lags.push_back(value);
    }

    if (lags.empty()) {
        throw std::invalid_argument("No lag frames were supplied.");
    }

    std::sort(lags.begin(), lags.end());
    lags.erase(std::unique(lags.begin(), lags.end()), lags.end());

    return lags;
}


std::vector<std::size_t> KeepValidLags(
    const std::vector<std::size_t>& requested,
    std::size_t num_frames)
{
    std::vector<std::size_t> result;

    for (const std::size_t lag : requested)
    {
        if (lag < num_frames) {
            result.push_back(lag);
        }
    }

    if (result.empty()) {
        throw std::runtime_error(
            "No requested lag is valid for this trajectory.");
    }

    if (result.front() != 0) {
        result.insert(result.begin(), 0);
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
        const double Lz = GetBoxLengthZ(box);

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
        const Vec3 mapped =
            geometry::MapPositionBetweenBoxes(
                physical_center,
                source_box,
                destination_box
            );

        if (std::abs(
                mapped.z -
                z_center) > tolerance)
        {
            throw std::runtime_error(
                "Fixed-z-center affine regression check failed.");
        }
    }
}


void CheckCompatible(
    const SegmentSurvivalFunction& a,
    const SegmentSurvivalFunction& b,
    const char* label)
{
    if (
        a.lag_times != b.lag_times ||
        a.tube_diameters != b.tube_diameters ||
        a.sample_counts != b.sample_counts ||
        a.survival.size() != b.survival.size()
    )
    {
        throw std::runtime_error(
            std::string(label) +
            ": dimensions or lag-wise origin cohorts differ.");
    }
}


void CheckLessEqual(
    const SegmentSurvivalFunction& lhs,
    const SegmentSurvivalFunction& rhs,
    const char* label)
{
    CheckCompatible(lhs, rhs, label);

    constexpr double tolerance =
        1.0e-12;

    for (std::size_t i = 0;
         i < lhs.survival.size();
         ++i)
    {
        if (lhs.survival[i] >
            rhs.survival[i] + tolerance)
        {
            throw std::runtime_error(
                std::string(label) +
                ": expected lhs <= rhs.");
        }
    }
}


void CheckHistoryInvariants(
    const HistoryDependentSurvivalResult& result)
{
    const auto& rf =
        result.reference_to_future;

    const auto& fr =
        result.future_to_reference;

    const auto& longitudinal =
        result.longitudinal_survival;

    CheckCompatible(
        rf.transverse_survival,
        fr.transverse_survival,
        "RF/FR transverse");

    CheckCompatible(
        rf.permanent_escape_survival,
        fr.permanent_escape_survival,
        "RF/FR permanent");

    CheckCompatible(
        rf.full_survival,
        fr.full_survival,
        "RF/FR full");

    CheckCompatible(
        rf.full_survival,
        longitudinal,
        "full/longitudinal");

    CheckLessEqual(
        rf.transverse_survival,
        rf.permanent_escape_survival,
        "RF transverse <= RF permanent");

    CheckLessEqual(
        fr.transverse_survival,
        fr.permanent_escape_survival,
        "FR transverse <= FR permanent");

    CheckLessEqual(
        rf.full_survival,
        rf.transverse_survival,
        "RF full <= RF transverse");

    CheckLessEqual(
        fr.full_survival,
        fr.transverse_survival,
        "FR full <= FR transverse");

    CheckLessEqual(
        rf.full_survival,
        longitudinal,
        "RF full <= longitudinal");

    CheckLessEqual(
        fr.full_survival,
        longitudinal,
        "FR full <= longitudinal");
}


void InitializeAccumulator(
    SegmentSurvivalFunction& dst,
    const SegmentSurvivalFunction& src)
{
    dst = src;

    std::fill(
        dst.survival.begin(),
        dst.survival.end(),
        0.0
    );

    std::fill(
        dst.sample_counts.begin(),
        dst.sample_counts.end(),
        0
    );
}


void AccumulateWeighted(
    SegmentSurvivalFunction& dst,
    const SegmentSurvivalFunction& src)
{
    if (
        dst.lag_times != src.lag_times ||
        dst.tube_diameters != src.tube_diameters ||
        dst.survival.size() != src.survival.size() ||
        dst.sample_counts.size() != src.sample_counts.size()
    )
    {
        throw std::runtime_error(
            "Accumulator dimensions changed between chains.");
    }

    const std::size_t num_lags =
        src.lag_times.size();

    const std::size_t num_diameters =
        src.tube_diameters.size();

    for (std::size_t lag_index = 0;
         lag_index < num_lags;
         ++lag_index)
    {
        const std::uint64_t count =
            src.sample_counts[lag_index];

        const double weight =
            static_cast<double>(count);

        dst.sample_counts[lag_index] +=
            count;

        for (std::size_t d = 0;
             d < num_diameters;
             ++d)
        {
            for (std::size_t s = 0;
                 s < NUM_SAMPLE_POINTS;
                 ++s)
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


void FinalizeAccumulator(
    SegmentSurvivalFunction& result)
{
    const std::size_t num_lags =
        result.lag_times.size();

    const std::size_t num_diameters =
        result.tube_diameters.size();

    for (std::size_t lag_index = 0;
         lag_index < num_lags;
         ++lag_index)
    {
        const std::uint64_t count =
            result.sample_counts[lag_index];

        if (count == 0) {
            continue;
        }

        const double denominator =
            static_cast<double>(count);

        for (std::size_t d = 0;
             d < num_diameters;
             ++d)
        {
            for (std::size_t s = 0;
                 s < NUM_SAMPLE_POINTS;
                 ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS +
                    lag_index * NUM_SAMPLE_POINTS +
                    s;

                result.survival[index] /=
                    denominator;
            }
        }
    }
}


struct EnsembleHistoryFields
{
    SegmentSurvivalFunction rf_transverse;
    SegmentSurvivalFunction rf_permanent;
    SegmentSurvivalFunction rf_full;

    SegmentSurvivalFunction fr_transverse;
    SegmentSurvivalFunction fr_permanent;
    SegmentSurvivalFunction fr_full;

    SegmentSurvivalFunction longitudinal;

    bool initialized{false};
};


void InitializeEnsemble(
    EnsembleHistoryFields& ensemble,
    const HistoryDependentSurvivalResult& result)
{
    InitializeAccumulator(
        ensemble.rf_transverse,
        result.reference_to_future.transverse_survival);

    InitializeAccumulator(
        ensemble.rf_permanent,
        result.reference_to_future.permanent_escape_survival);

    InitializeAccumulator(
        ensemble.rf_full,
        result.reference_to_future.full_survival);

    InitializeAccumulator(
        ensemble.fr_transverse,
        result.future_to_reference.transverse_survival);

    InitializeAccumulator(
        ensemble.fr_permanent,
        result.future_to_reference.permanent_escape_survival);

    InitializeAccumulator(
        ensemble.fr_full,
        result.future_to_reference.full_survival);

    InitializeAccumulator(
        ensemble.longitudinal,
        result.longitudinal_survival);

    ensemble.initialized = true;
}


void AccumulateEnsemble(
    EnsembleHistoryFields& ensemble,
    const HistoryDependentSurvivalResult& result)
{
    AccumulateWeighted(
        ensemble.rf_transverse,
        result.reference_to_future.transverse_survival);

    AccumulateWeighted(
        ensemble.rf_permanent,
        result.reference_to_future.permanent_escape_survival);

    AccumulateWeighted(
        ensemble.rf_full,
        result.reference_to_future.full_survival);

    AccumulateWeighted(
        ensemble.fr_transverse,
        result.future_to_reference.transverse_survival);

    AccumulateWeighted(
        ensemble.fr_permanent,
        result.future_to_reference.permanent_escape_survival);

    AccumulateWeighted(
        ensemble.fr_full,
        result.future_to_reference.full_survival);

    AccumulateWeighted(
        ensemble.longitudinal,
        result.longitudinal_survival);
}


void FinalizeEnsemble(
    EnsembleHistoryFields& ensemble)
{
    FinalizeAccumulator(ensemble.rf_transverse);
    FinalizeAccumulator(ensemble.rf_permanent);
    FinalizeAccumulator(ensemble.rf_full);

    FinalizeAccumulator(ensemble.fr_transverse);
    FinalizeAccumulator(ensemble.fr_permanent);
    FinalizeAccumulator(ensemble.fr_full);

    FinalizeAccumulator(ensemble.longitudinal);
}


double MaxAbsDifference(
    const SegmentSurvivalFunction& a,
    const SegmentSurvivalFunction& b)
{
    CheckCompatible(
        a,
        b,
        "MaxAbsDifference");

    double maximum = 0.0;

    for (std::size_t i = 0;
         i < a.survival.size();
         ++i)
    {
        maximum =
            std::max(
                maximum,
                std::abs(
                    a.survival[i] -
                    b.survival[i]
                )
            );
    }

    return maximum;
}


double MaxAbsDifference(
    const TubeSurvivalFunction& a,
    const TubeSurvivalFunction& b)
{
    if (
        a.lag_times != b.lag_times ||
        a.tube_diameters != b.tube_diameters ||
        a.sample_counts != b.sample_counts ||
        a.survival.size() != b.survival.size()
    )
    {
        throw std::runtime_error(
            "Tube fields are incompatible.");
    }

    double maximum = 0.0;

    for (std::size_t i = 0;
         i < a.survival.size();
         ++i)
    {
        maximum =
            std::max(
                maximum,
                std::abs(
                    a.survival[i] -
                    b.survival[i]
                )
            );
    }

    return maximum;
}


void CheckLagZero(
    const TubeSurvivalFunction& result,
    const char* label)
{
    if (
        result.lag_times.empty() ||
        result.lag_times.front() != 0.0
    )
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
        const double value =
            result.survival[
                d * num_lags
            ];

        if (std::abs(value - 1.0) >
            1.0e-12)
        {
            throw std::runtime_error(
                std::string(label) +
                ": survival(0) != 1.");
        }
    }
}


std::string Out(
    const std::filesystem::path& directory,
    const std::string& filename)
{
    return (
        directory /
        filename
    ).string();
}


void WriteTubeComparisonSummary(
    const std::filesystem::path& output_directory,
    const TubeSurvivalFunction& rf_full,
    const TubeSurvivalFunction& fr_full,
    const TubeSurvivalFunction& rf_transverse,
    const TubeSurvivalFunction& fr_transverse,
    const TubeSurvivalFunction& rf_permanent,
    const TubeSurvivalFunction& fr_permanent,
    const TubeSurvivalFunction& longitudinal)
{
    std::ofstream file(
        output_directory /
        "directional_history_tube_comparison.csv"
    );

    if (!file) {
        throw std::runtime_error(
            "Could not open directional tube comparison CSV.");
    }

    file << std::setprecision(17);

    file
        << "tube_diameter,"
        << "lag_time,"
        << "rf_full,"
        << "fr_full,"
        << "delta_full,"
        << "rf_transverse,"
        << "fr_transverse,"
        << "delta_transverse,"
        << "rf_permanent,"
        << "fr_permanent,"
        << "delta_permanent,"
        << "longitudinal,"
        << "sample_count\n";

    const std::size_t num_lags =
        rf_full.lag_times.size();

    const std::size_t num_diameters =
        rf_full.tube_diameters.size();

    for (std::size_t d = 0;
         d < num_diameters;
         ++d)
    {
        for (std::size_t lag = 0;
             lag < num_lags;
             ++lag)
        {
            const std::size_t index =
                d * num_lags + lag;

            file
                << rf_full.tube_diameters[d]
                << ','
                << rf_full.lag_times[lag]
                << ','
                << rf_full.survival[index]
                << ','
                << fr_full.survival[index]
                << ','
                << (
                    rf_full.survival[index] -
                    fr_full.survival[index]
                )
                << ','
                << rf_transverse.survival[index]
                << ','
                << fr_transverse.survival[index]
                << ','
                << (
                    rf_transverse.survival[index] -
                    fr_transverse.survival[index]
                )
                << ','
                << rf_permanent.survival[index]
                << ','
                << fr_permanent.survival[index]
                << ','
                << (
                    rf_permanent.survival[index] -
                    fr_permanent.survival[index]
                )
                << ','
                << longitudinal.survival[index]
                << ','
                << rf_full.sample_counts[lag]
                << '\n';
        }
    }
}

} // namespace


int main(int argc, char* argv[])
{
    try
    {
        if (argc != 8)
        {
            std::cerr
                << "Usage: "
                << argv[0]
                << " <Z1+SP.dat>"
                << " <frame_time>"
                << " <tube_diameter>"
                << " <z_center>"
                << " <max_chains=0_for_all>"
                << " <lags_csv>"
                << " <output_directory>\n";

            return 2;
        }

        const std::string filename =
            argv[1];

        const long frame_time =
            std::stol(argv[2]);

        const double tube_diameter =
            std::stod(argv[3]);

        const double z_center =
            std::stod(argv[4]);

        const std::size_t max_chains =
            static_cast<std::size_t>(
                std::stoull(argv[5])
            );

        const std::vector<std::size_t>
            requested_lags =
                ParseLagList(argv[6]);

        const std::filesystem::path
            output_directory =
                argv[7];

        if (frame_time <= 0) {
            throw std::invalid_argument(
                "frame_time must be positive.");
        }

        if (tube_diameter <= 0.0) {
            throw std::invalid_argument(
                "tube_diameter must be positive.");
        }

        std::filesystem::create_directories(
            output_directory
        );

        std::cout
            << "=============================================\n"
            << "REAL-DATA DIRECTIONAL HISTORY COMPARISON\n"
            << "=============================================\n\n"
            << "Parsing "
            << filename
            << " ...\n";

        const auto parse_start =
            std::chrono::steady_clock::now();

        PrimitivePathTrajectory trajectory =
            parse_z1_file(filename);

        const auto parse_stop =
            std::chrono::steady_clock::now();

        if (
            trajectory.chains.empty() ||
            trajectory.frame_boxes.empty()
        )
        {
            throw std::runtime_error(
                "Parser returned no chains or no frame boxes.");
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

        for (const ChainTrajectory& chain :
             trajectory.chains)
        {
            if (
                chain.frame_offsets.size() !=
                num_frames + 1
            )
            {
                throw std::runtime_error(
                    "Chain frame count does not match box count.");
            }
        }

        const std::size_t chains_to_run =
            (max_chains == 0)
                ? trajectory.chains.size()
                : std::min(
                    max_chains,
                    trajectory.chains.size()
                );

        const std::vector<std::size_t> lags =
            KeepValidLags(
                requested_lags,
                num_frames
            );

        const std::vector<double> diameters{
            tube_diameter
        };

        const SegmentSurvivalOptions options{
            .apply_affine_correction = true
        };

        std::cout
            << "Parsed chains      : "
            << trajectory.chains.size()
            << '\n'
            << "Chains to analyze  : "
            << chains_to_run
            << '\n'
            << "Frames/chain       : "
            << num_frames
            << '\n'
            << "Frame time         : "
            << frame_time
            << '\n'
            << "Tube diameter      : "
            << tube_diameter
            << '\n'
            << "Fixed z center     : "
            << std::setprecision(15)
            << z_center
            << '\n'
            << "Output directory   : "
            << output_directory.string()
            << '\n'
            << "Requested lags     : ";

        for (const auto lag : lags) {
            std::cout << lag << ' ';
        }

        std::cout
            << "\n\n"
            << "Computing both directed definitions with\n"
            << "the SAME affine pullback, transverse history,\n"
            << "longitudinal first-passage history, and time origins.\n\n";

        EnsembleHistoryFields ensemble;

        const auto work_start =
            std::chrono::steady_clock::now();

        for (std::size_t chain_index = 0;
             chain_index < chains_to_run;
             ++chain_index)
        {
            const ChainTrajectory& chain =
                trajectory.chains[chain_index];

            const HistoryDependentSurvivalResult result =
                ComputeHistoryDependentSurvivalFunction(
                    chain,
                    trajectory.frame_boxes,
                    lags,
                    diameters,
                    options,
                    false
                );

            CheckHistoryInvariants(
                result
            );

            if (!ensemble.initialized) {
                InitializeEnsemble(
                    ensemble,
                    result
                );
            }

            AccumulateEnsemble(
                ensemble,
                result
            );

            const std::size_t done =
                chain_index + 1;

            if (
                done == 1 ||
                done % 25 == 0 ||
                done == chains_to_run
            )
            {
                const auto now =
                    std::chrono::steady_clock::now();

                const std::chrono::duration<double>
                    elapsed =
                        now - work_start;

                std::cout
                    << "Finished chain "
                    << done
                    << " / "
                    << chains_to_run
                    << "  (elapsed "
                    << std::fixed
                    << std::setprecision(1)
                    << elapsed.count()
                    << " s)\n";
            }
        }

        if (!ensemble.initialized) {
            throw std::runtime_error(
                "No chains were analyzed.");
        }

        FinalizeEnsemble(
            ensemble
        );

        HistoryDependentSurvivalResult ensemble_result;

        ensemble_result.reference_to_future.transverse_survival =
            ensemble.rf_transverse;

        ensemble_result.reference_to_future.permanent_escape_survival =
            ensemble.rf_permanent;

        ensemble_result.reference_to_future.full_survival =
            ensemble.rf_full;

        ensemble_result.future_to_reference.transverse_survival =
            ensemble.fr_transverse;

        ensemble_result.future_to_reference.permanent_escape_survival =
            ensemble.fr_permanent;

        ensemble_result.future_to_reference.full_survival =
            ensemble.fr_full;

        ensemble_result.longitudinal_survival =
            ensemble.longitudinal;

        CheckHistoryInvariants(
            ensemble_result
        );

        const TubeSurvivalFunction tube_rf_transverse =
            ComputeTubeSurvivalFunction(
                ensemble.rf_transverse);

        const TubeSurvivalFunction tube_rf_permanent =
            ComputeTubeSurvivalFunction(
                ensemble.rf_permanent);

        const TubeSurvivalFunction tube_rf_full =
            ComputeTubeSurvivalFunction(
                ensemble.rf_full);

        const TubeSurvivalFunction tube_fr_transverse =
            ComputeTubeSurvivalFunction(
                ensemble.fr_transverse);

        const TubeSurvivalFunction tube_fr_permanent =
            ComputeTubeSurvivalFunction(
                ensemble.fr_permanent);

        const TubeSurvivalFunction tube_fr_full =
            ComputeTubeSurvivalFunction(
                ensemble.fr_full);

        const TubeSurvivalFunction tube_longitudinal =
            ComputeTubeSurvivalFunction(
                ensemble.longitudinal);

        CheckLagZero(
            tube_rf_transverse,
            "RF transverse");

        CheckLagZero(
            tube_rf_permanent,
            "RF permanent");

        CheckLagZero(
            tube_rf_full,
            "RF full");

        CheckLagZero(
            tube_fr_transverse,
            "FR transverse");

        CheckLagZero(
            tube_fr_permanent,
            "FR permanent");

        CheckLagZero(
            tube_fr_full,
            "FR full");

        CheckLagZero(
            tube_longitudinal,
            "longitudinal");

        WriteSegmentSurvivalFunction(
            ensemble.rf_transverse,
            Out(
                output_directory,
                "directional_history_segment_rf_transverse.csv"
            )
        );

        WriteSegmentSurvivalFunction(
            ensemble.rf_permanent,
            Out(
                output_directory,
                "directional_history_segment_rf_permanent.csv"
            )
        );

        WriteSegmentSurvivalFunction(
            ensemble.rf_full,
            Out(
                output_directory,
                "directional_history_segment_rf_full.csv"
            )
        );

        WriteSegmentSurvivalFunction(
            ensemble.fr_transverse,
            Out(
                output_directory,
                "directional_history_segment_fr_transverse.csv"
            )
        );

        WriteSegmentSurvivalFunction(
            ensemble.fr_permanent,
            Out(
                output_directory,
                "directional_history_segment_fr_permanent.csv"
            )
        );

        WriteSegmentSurvivalFunction(
            ensemble.fr_full,
            Out(
                output_directory,
                "directional_history_segment_fr_full.csv"
            )
        );

        WriteSegmentSurvivalFunction(
            ensemble.longitudinal,
            Out(
                output_directory,
                "directional_history_segment_longitudinal.csv"
            )
        );

        WriteTubeSurvivalFunction(
            tube_rf_transverse,
            Out(
                output_directory,
                "directional_history_tube_rf_transverse.csv"
            )
        );

        WriteTubeSurvivalFunction(
            tube_rf_permanent,
            Out(
                output_directory,
                "directional_history_tube_rf_permanent.csv"
            )
        );

        WriteTubeSurvivalFunction(
            tube_rf_full,
            Out(
                output_directory,
                "directional_history_tube_rf_full.csv"
            )
        );

        WriteTubeSurvivalFunction(
            tube_fr_transverse,
            Out(
                output_directory,
                "directional_history_tube_fr_transverse.csv"
            )
        );

        WriteTubeSurvivalFunction(
            tube_fr_permanent,
            Out(
                output_directory,
                "directional_history_tube_fr_permanent.csv"
            )
        );

        WriteTubeSurvivalFunction(
            tube_fr_full,
            Out(
                output_directory,
                "directional_history_tube_fr_full.csv"
            )
        );

        WriteTubeSurvivalFunction(
            tube_longitudinal,
            Out(
                output_directory,
                "directional_history_tube_longitudinal.csv"
            )
        );

        WriteTubeComparisonSummary(
            output_directory,
            tube_rf_full,
            tube_fr_full,
            tube_rf_transverse,
            tube_fr_transverse,
            tube_rf_permanent,
            tube_fr_permanent,
            tube_longitudinal
        );

        const double segment_full_difference =
            MaxAbsDifference(
                ensemble.rf_full,
                ensemble.fr_full
            );

        const double segment_transverse_difference =
            MaxAbsDifference(
                ensemble.rf_transverse,
                ensemble.fr_transverse
            );

        const double segment_permanent_difference =
            MaxAbsDifference(
                ensemble.rf_permanent,
                ensemble.fr_permanent
            );

        const double tube_full_difference =
            MaxAbsDifference(
                tube_rf_full,
                tube_fr_full
            );

        const double tube_transverse_difference =
            MaxAbsDifference(
                tube_rf_transverse,
                tube_fr_transverse
            );

        const double tube_permanent_difference =
            MaxAbsDifference(
                tube_rf_permanent,
                tube_fr_permanent
            );

        const auto work_stop =
            std::chrono::steady_clock::now();

        const std::chrono::duration<double>
            parse_elapsed =
                parse_stop - parse_start;

        const std::chrono::duration<double>
            work_elapsed =
                work_stop - work_start;

        std::cout
            << "\n=============================================\n"
            << "RF vs FR DIFFERENCE SUMMARY\n"
            << "=============================================\n"
            << std::setprecision(8)
            << "max |delta phi|, transverse : "
            << segment_transverse_difference
            << '\n'
            << "max |delta phi|, permanent  : "
            << segment_permanent_difference
            << '\n'
            << "max |delta phi|, full       : "
            << segment_full_difference
            << '\n'
            << "max |delta mu|, transverse  : "
            << tube_transverse_difference
            << '\n'
            << "max |delta mu|, permanent   : "
            << tube_permanent_difference
            << '\n'
            << "max |delta mu|, full        : "
            << tube_full_difference
            << '\n'
            << "\nParsing time                : "
            << parse_elapsed.count()
            << " s\n"
            << "Calculation time            : "
            << work_elapsed.count()
            << " s\n"
            << "Total lag-zero origins      : "
            << tube_rf_full.sample_counts.front()
            << "\n\n"
            << "REAL-DATA DIRECTIONAL HISTORY COMPARISON "
            << "COMPLETED SUCCESSFULLY\n";

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
