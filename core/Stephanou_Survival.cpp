#include "Geometry_Utils.hpp"
#include "Parse_Z1_File.hpp"
#include "Survival_IO.hpp"
#include "Trajectory_Time.hpp"
#include "Tube_Survival.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


// ------------------------------------------------------------
// Parse "3.0,5.0,7.0" -> vector<double>
// ------------------------------------------------------------
std::vector<double> ParseDoubleList(const std::string& input)
{
    std::vector<double> values;

    std::stringstream stream(input);
    std::string token;

    while (std::getline(stream, token, ','))
    {
        values.push_back(std::stod(token));
    }

    return values;
}


// ------------------------------------------------------------
// Parse "20.6,20.6,39.2" -> Vec3
// ------------------------------------------------------------
Vec3 ParseVec3(const std::string& input)
{
    const std::vector<double> values =
        ParseDoubleList(input);

    if (values.size() != 3)
    {
        throw std::invalid_argument(
            "Box center must contain exactly three values: x,y,z"
        );
    }

    return Vec3{
        values[0],
        values[1],
        values[2]
    };
}


// ------------------------------------------------------------
// Reconstruct each box origin from a fixed physical box center.
//
// center = origin + H * (1/2, 1/2, 1/2)
// ------------------------------------------------------------
void SetFixedBoxCenter(
    std::vector<Box>& frame_boxes,
    const Vec3& center)
{
    const Vec3 fractional_center{
        0.5,
        0.5,
        0.5
    };

    for (Box& box : frame_boxes)
    {
        const Vec3 half_box =
            geometry::Multiply(
                box.matrix,
                fractional_center
            );

        box.origin = Vec3{
            center.x - half_box.x,
            center.y - half_box.y,
            center.z - half_box.z
        };
    }
}


// ------------------------------------------------------------
// Initialize an ensemble accumulator from one chain's result.
// ------------------------------------------------------------
void InitializeAccumulator(
    SegmentSurvivalFunction& accumulator,
    const SegmentSurvivalFunction& result)
{
    accumulator = result;

    std::fill(
        accumulator.survival.begin(),
        accumulator.survival.end(),
        0.0
    );

    std::fill(
        accumulator.sample_counts.begin(),
        accumulator.sample_counts.end(),
        0
    );
}


// ------------------------------------------------------------
// Add one chain to the ensemble.
//
// Weight by the number of valid time origins for each lag.
// ------------------------------------------------------------
void Accumulate(
    SegmentSurvivalFunction& accumulator,
    const SegmentSurvivalFunction& result)
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

        accumulator.sample_counts[lag_index] += count;

        for (std::size_t d = 0;
            d < num_diameters;
            ++d)
        {
            for (std::size_t s = 0;
                s < NUM_SAMPLE_POINTS;
                ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS
                    + lag_index * NUM_SAMPLE_POINTS
                    + s;

                accumulator.survival[index] +=
                    result.survival[index]
                    * static_cast<double>(count);
            }
        }
    }
}


// ------------------------------------------------------------
// Convert accumulated counts back into ensemble probabilities.
// ------------------------------------------------------------
void FinalizeAccumulator(
    SegmentSurvivalFunction& accumulator)
{
    const std::size_t num_lags =
        accumulator.lag_times.size();

    const std::size_t num_diameters =
        accumulator.tube_diameters.size();

    for (std::size_t lag_index = 0;
        lag_index < num_lags;
        ++lag_index)
    {
        const std::uint64_t count =
            accumulator.sample_counts[lag_index];

        if (count == 0) {
            continue;
        }

        for (std::size_t d = 0;
            d < num_diameters;
            ++d)
        {
            for (std::size_t s = 0;
                s < NUM_SAMPLE_POINTS;
                ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS
                    + lag_index * NUM_SAMPLE_POINTS
                    + s;

                accumulator.survival[index] /=
                    static_cast<double>(count);
            }
        }
    }
}


int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cerr
            << "Usage:\n"
            << "  Parallel_Segment_Survival "
            << "<Z1+SP.dat> "
            << "<tube_diameters> "
            << "<box_center_x,y,z> "
            << "<output_prefix>\n\n"

            << "Example:\n"
            << "  Parallel_Segment_Survival "
            << "Z1+SP.dat "
            << "\"3.0,5.0,7.0,9.0\" "
            << "\"20.6116,20.6103,39.226586\" "
            << "\"AMP_7.6_FREQ_1750_TRIAL_1\"\n";

        return 1;
    }

    try
    {
        // --------------------------------------------------------
        // Command-line arguments
        // --------------------------------------------------------

        const std::string z1_file =
            argv[1];

        const std::vector<double> tube_diameters =
            ParseDoubleList(argv[2]);

        const Vec3 box_center =
            ParseVec3(argv[3]);

        const std::string output_prefix =
            argv[4];


        // --------------------------------------------------------
        // Read Z1+ trajectory
        // --------------------------------------------------------

        PrimitivePathTrajectory parsed_trajectory =
            parse_z1_file(z1_file);

        auto& chain_trajectories =
            parsed_trajectory.chains;

        auto& frame_boxes =
            parsed_trajectory.frame_boxes;

        const std::size_t num_frames =
            frame_boxes.size();


        // --------------------------------------------------------
        // Restore the physical origin of every deforming box
        // --------------------------------------------------------

        SetFixedBoxCenter(
            frame_boxes,
            box_center
        );


        // --------------------------------------------------------
        // Z1+ currently does not contain physical timestamps.
        //
        // Assign frame numbers:
        //
        //     0, 1, 2, 3, ...
        //
        // so lag_time in the output currently means frame lag.
        // --------------------------------------------------------

        AssignUniformTimesteps(
            chain_trajectories,
            0,
            1
        );


        // --------------------------------------------------------
        // Lag schedule
        // --------------------------------------------------------

        std::vector<std::size_t> lag_frames{
            0, 1, 2, 4, 8, 16
        };

        for (std::size_t lag = 32;
            lag < num_frames;
            lag += 16)
        {
            lag_frames.push_back(lag);
        }


        // --------------------------------------------------------
        // Enable affine correction
        // --------------------------------------------------------

        const SegmentSurvivalOptions options{
            .apply_affine_correction = true
        };


        // --------------------------------------------------------
        // Ensemble accumulators
        // --------------------------------------------------------

        SegmentSurvivalFunction ensemble_raw;
        SegmentSurvivalFunction ensemble_affine;
        SegmentSurvivalFunction ensemble_history;
        SegmentSurvivalFunction ensemble_permanent;

        bool initialized = false;


        // --------------------------------------------------------
        // Analyze every chain
        // --------------------------------------------------------

        for (const ChainTrajectory& chain :
            chain_trajectories)
        {
            // Raw instantaneous survival
            const SegmentSurvivalFunction raw =
                ComputeSegmentSurvivalFunction(
                    chain,
                    lag_frames,
                    tube_diameters
                );


            // Affine-corrected instantaneous survival
            const SegmentSurvivalFunction affine =
                ComputeSegmentSurvivalFunction(
                    chain,
                    frame_boxes,
                    lag_frames,
                    tube_diameters,
                    options
                );


            // Affine-corrected history-dependent survival
            const HistoryDependentSurvivalResult history =
                ComputeHistoryDependentSurvivalFunction(
                    chain,
                    frame_boxes,
                    lag_frames,
                    tube_diameters,
                    options,
                    false
                );


            if (!initialized)
            {
                InitializeAccumulator(
                    ensemble_raw,
                    raw
                );

                InitializeAccumulator(
                    ensemble_affine,
                    affine
                );

                InitializeAccumulator(
                    ensemble_history,
                    history.stephanou_survival
                );

                InitializeAccumulator(
                    ensemble_permanent,
                    history.permanent_escape_survival
                );

                initialized = true;
            }


            Accumulate(
                ensemble_raw,
                raw
            );

            Accumulate(
                ensemble_affine,
                affine
            );

            Accumulate(
                ensemble_history,
                history.stephanou_survival
            );

            Accumulate(
                ensemble_permanent,
                history.permanent_escape_survival
            );
        }


        // --------------------------------------------------------
        // Finish ensemble averaging
        // --------------------------------------------------------

        FinalizeAccumulator(ensemble_raw);
        FinalizeAccumulator(ensemble_affine);
        FinalizeAccumulator(ensemble_history);
        FinalizeAccumulator(ensemble_permanent);


        // --------------------------------------------------------
        // Convert contour-resolved survival into tube survival
        // --------------------------------------------------------

        const TubeSurvivalFunction tube_raw =
            ComputeTubeSurvivalFunction(
                ensemble_raw
            );

        const TubeSurvivalFunction tube_affine =
            ComputeTubeSurvivalFunction(
                ensemble_affine
            );

        const TubeSurvivalFunction tube_history =
            ComputeTubeSurvivalFunction(
                ensemble_history
            );

        const TubeSurvivalFunction tube_permanent =
            ComputeTubeSurvivalFunction(
                ensemble_permanent
            );


        // --------------------------------------------------------
        // Write output
        // --------------------------------------------------------

        WriteSegmentSurvivalFunction(
            ensemble_raw,
            output_prefix + "_segment_raw.csv"
        );

        WriteSegmentSurvivalFunction(
            ensemble_affine,
            output_prefix + "_segment_affine.csv"
        );

        WriteSegmentSurvivalFunction(
            ensemble_history,
            output_prefix + "_segment_history.csv"
        );

        WriteSegmentSurvivalFunction(
            ensemble_permanent,
            output_prefix + "_segment_permanent_escape.csv"
        );


        WriteTubeSurvivalFunction(
            tube_raw,
            output_prefix + "_tube_raw.csv"
        );

        WriteTubeSurvivalFunction(
            tube_affine,
            output_prefix + "_tube_affine.csv"
        );

        WriteTubeSurvivalFunction(
            tube_history,
            output_prefix + "_tube_history.csv"
        );

        WriteTubeSurvivalFunction(
            tube_permanent,
            output_prefix + "_tube_permanent_escape.csv"
        );


        std::cout
            << "Analyzed "
            << chain_trajectories.size()
            << " chains.\n"
            << "Wrote results with prefix: "
            << output_prefix
            << '\n';

        return 0;
    }

    catch (const std::exception& error)
    {
        std::cerr
            << "ERROR: "
            << error.what()
            << '\n';

        return 1;
    }
}