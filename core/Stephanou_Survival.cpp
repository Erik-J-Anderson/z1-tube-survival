#include "Input_Parsers.hpp"
#include "Parse_Z1_File.hpp"
#include "Survival_Accumulator.hpp"
#include "Survival_IO.hpp"
#include "Trajectory_Time.hpp"
#include "Tube_Survival.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


int main(int argc, char* argv[])
{
    try
    {
        // ------------------------------------------------------------
        // Command-line arguments
        //
        // argv[1] : Z1+SP.dat
        // argv[2] : comma-separated tube diameters
        // argv[3] : box center x,y,z
        // argv[4] : output prefix
        // ------------------------------------------------------------
        if (argc != 5)
        {
            std::cerr
                << "Usage:\n"
                << "  " << argv[0]
                << " <Z1+SP.dat>"
                << " <tube_diameters>"
                << " <box_center_x,y,z>"
                << " <output_prefix>\n\n"
                << "Example:\n"
                << "  " << argv[0]
                << " Z1+SP.dat"
                << " \"5.0,7.0,9.0,9.2,9.4\""
                << " \"20.6116,20.6103,39.226586\""
                << " AMP_7.6_FREQ_1750_TRIAL_1\n";

            return 1;
        }


        const std::string z1_file =
            argv[1];

        const std::vector<double> tube_diameters =
            ParseDoubleList(argv[2]);

        const Vec3 box_center =
            ParseVec3(argv[3]);

        const std::string output_prefix =
            argv[4];


        if (tube_diameters.empty())
        {
            throw std::invalid_argument(
                "At least one tube diameter must be specified."
            );
        }


        // ------------------------------------------------------------
        // Parse Z1+ trajectory.
        // ------------------------------------------------------------
        std::cout
            << "Reading Z1+ trajectory...\n";

        PrimitivePathTrajectory parsed_trajectory =
            parse_z1_file(z1_file);


        std::vector<Box>& frame_boxes =
            parsed_trajectory.frame_boxes;

        std::vector<ChainTrajectory>& chain_trajectories =
            parsed_trajectory.chains;


        if (chain_trajectories.empty())
        {
            throw std::runtime_error(
                "No chain trajectories were parsed."
            );
        }

        if (frame_boxes.empty())
        {
            throw std::runtime_error(
                "No frame boxes were parsed."
            );
        }


        // ------------------------------------------------------------
        // Reconstruct the physical box origin using the fixed
        // simulation-box center supplied on the command line.
        // ------------------------------------------------------------
        SetFixedBoxCenter(
            frame_boxes,
            box_center
        );


        // ------------------------------------------------------------
        // For the current production interface, lag time is expressed
        // in saved-frame units.
        // ------------------------------------------------------------
        AssignUniformTimesteps(
            chain_trajectories,
            0,
            1
        );


        const std::size_t num_frames =
            frame_boxes.size();


        // ------------------------------------------------------------
        // Lag schedule.
        //
        // Dense early sampling, followed by every 16 frames.
        // ------------------------------------------------------------
        std::vector<std::size_t> lag_frames;

        const std::vector<std::size_t> early_lags{
            0, 1, 2, 4, 8, 16
        };

        for (const std::size_t lag : early_lags)
        {
            if (lag < num_frames) {
                lag_frames.push_back(lag);
            }
        }

        for (std::size_t lag = 32;
            lag < num_frames;
            lag += 16)
        {
            lag_frames.push_back(lag);
        }


        // ------------------------------------------------------------
        // Affine-corrected calculations.
        // ------------------------------------------------------------
        SegmentSurvivalOptions options{};

        options.apply_affine_correction = true;


        // ------------------------------------------------------------
        // Ensemble accumulators.
        // ------------------------------------------------------------
        SegmentSurvivalFunction ensemble_raw;
        SegmentSurvivalFunction ensemble_affine;
        SegmentSurvivalFunction ensemble_history;
        SegmentSurvivalFunction ensemble_permanent;


        bool raw_initialized = false;
        bool affine_initialized = false;
        bool history_initialized = false;
        bool permanent_initialized = false;


        // ------------------------------------------------------------
        // Process each chain.
        // ------------------------------------------------------------
        std::size_t chain_count = 0;

        for (const ChainTrajectory& chain : chain_trajectories)
        {
            // --------------------------------------------------------
            // Instantaneous raw survival.
            // --------------------------------------------------------
            const SegmentSurvivalFunction raw =
                ComputeSegmentSurvivalFunction(
                    chain,
                    lag_frames,
                    tube_diameters
                );


            // --------------------------------------------------------
            // Instantaneous affine-corrected survival.
            // --------------------------------------------------------
            const SegmentSurvivalFunction affine =
                ComputeSegmentSurvivalFunction(
                    chain,
                    frame_boxes,
                    lag_frames,
                    tube_diameters,
                    options
                );


            // --------------------------------------------------------
            // History-dependent affine-corrected survival.
            //
            // false = use the valid time-origin cohort associated
            // with each individual lag.
            // --------------------------------------------------------
            const HistoryDependentSurvivalResult history =
                ComputeHistoryDependentSurvivalFunction(
                    chain,
                    frame_boxes,
                    lag_frames,
                    tube_diameters,
                    options,
                    false
                );


            // --------------------------------------------------------
            // Initialize accumulators from the first chain.
            // --------------------------------------------------------
            if (!raw_initialized)
            {
                InitializeAccumulator(
                    ensemble_raw,
                    raw
                );

                raw_initialized = true;
            }

            if (!affine_initialized)
            {
                InitializeAccumulator(
                    ensemble_affine,
                    affine
                );

                affine_initialized = true;
            }

            if (!history_initialized)
            {
                InitializeAccumulator(
                    ensemble_history,
                    history.reference_to_future.transverse_survival
                );

                history_initialized = true;
            }

            if (!permanent_initialized)
            {
                InitializeAccumulator(
                    ensemble_permanent,
                    history.reference_to_future.permanent_escape_survival
                );

                permanent_initialized = true;
            }


            // --------------------------------------------------------
            // Accumulate this chain into the ensemble.
            // --------------------------------------------------------
            AccumulateSegmentSurvival(
                ensemble_raw,
                raw
            );

            AccumulateSegmentSurvival(
                ensemble_affine,
                affine
            );

            AccumulateSegmentSurvival(
                ensemble_history,
                history.reference_to_future.transverse_survival
            );

            AccumulateSegmentSurvival(
                ensemble_permanent,
                history.reference_to_future.permanent_escape_survival
            );


            // --------------------------------------------------------
            // Progress output.
            // --------------------------------------------------------
            ++chain_count;

            if (chain_count == 1 ||
                chain_count % 25 == 0 ||
                chain_count == chain_trajectories.size())
            {
                std::cout
                    << "Finished Chain "
                    << chain_count
                    << " / "
                    << chain_trajectories.size()
                    << std::endl;
            }
        }


        // ------------------------------------------------------------
        // Convert weighted accumulated counts into ensemble averages.
        // ------------------------------------------------------------
        FinalizeAccumulator(
            ensemble_raw
        );

        FinalizeAccumulator(
            ensemble_affine
        );

        FinalizeAccumulator(
            ensemble_history
        );

        FinalizeAccumulator(
            ensemble_permanent
        );


        // ------------------------------------------------------------
        // Integrate segment survival along the contour to obtain
        // tube survival.
        // ------------------------------------------------------------
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


        // ------------------------------------------------------------
        // Write segment-survival output.
        // ------------------------------------------------------------
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


        // ------------------------------------------------------------
        // Write tube-survival output.
        // ------------------------------------------------------------
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
            << "Survival analysis complete.\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "Error: "
            << error.what()
            << '\n';

        return 1;
    }
}
