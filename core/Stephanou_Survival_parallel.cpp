#include "Input_Parsers.hpp"
#include "MPI_Comms.hpp"
#include "Parse_Z1_File.hpp"
#include "Survival_Accumulator.hpp"
#include "Survival_IO.hpp"
#include "Trajectory_Time.hpp"
#include "Tube_Survival.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
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

    // ============================================================
    // Production lag schedule
    //
    // Keep the short-time resolution used previously:
    //
    //     0, 1, 2, 4, 8, 16
    //
    // then sample every 16 saved frames:
    //
    //     32, 48, 64, ...
    //
    // All lags are in saved-frame units because the driver assigns
    // timesteps 0,1,2,... to the parsed primitive-path frames.
    // ============================================================

    std::vector<std::size_t> BuildProductionLagFrames(
        std::size_t num_frames)
    {
        if (num_frames == 0)
        {
            throw std::runtime_error(
                "Cannot construct lag schedule for zero frames."
            );
        }

        const std::vector<std::size_t> initial_lags{
            0, 1, 2, 4, 8, 16
        };

        std::vector<std::size_t> lag_frames;

        for (const std::size_t lag : initial_lags)
        {
            if (lag < num_frames) {
                lag_frames.push_back(lag);
            }
        }

        for (std::size_t lag = 16;
            lag < num_frames;
            lag += 8)
        {
            lag_frames.push_back(lag);
        }

        if (lag_frames.empty())
        {
            throw std::runtime_error(
                "Production lag schedule is empty."
            );
        }

        return lag_frames;
    }


    // ============================================================
    // Verify that all four estimators used the same time-origin
    // cohorts at every lag.
    // ============================================================

    void CheckMatchingCohorts(
        const SegmentSurvivalFunction& raw,
        const SegmentSurvivalFunction& affine,
        const SegmentSurvivalFunction& history,
        const SegmentSurvivalFunction& permanent)
    {
        if (raw.sample_counts != affine.sample_counts ||
            affine.sample_counts != history.sample_counts ||
            history.sample_counts != permanent.sample_counts)
        {
            throw std::runtime_error(
                "Survival estimators used different lag-wise "
                "time-origin cohorts."
            );
        }
    }


    // ============================================================
    // History-dependent survival must never exceed the affine
    // instantaneous result.
    //
    // Permanent-escape survival must contain the Stephanou-like
    // history survival:
    //
    //     history <= permanent
    //
    // because permanent survival only asks whether the outer
    // absorbing boundary has ever been crossed.
    // ============================================================

    void CheckHistoryRelations(
        const SegmentSurvivalFunction& affine,
        const SegmentSurvivalFunction& history,
        const SegmentSurvivalFunction& permanent)
    {
        if (affine.survival.size() != history.survival.size() ||
            history.survival.size() != permanent.survival.size())
        {
            throw std::runtime_error(
                "History relation check received incompatible arrays."
            );
        }

        constexpr double tolerance = 1.0e-12;

        for (std::size_t i = 0;
            i < history.survival.size();
            ++i)
        {
            if (history.survival[i] >
                affine.survival[i] + tolerance)
            {
                throw std::runtime_error(
                    "History-dependent survival exceeded affine "
                    "instantaneous survival."
                );
            }

            if (history.survival[i] >
                permanent.survival[i] + tolerance)
            {
                throw std::runtime_error(
                    "History-dependent survival exceeded "
                    "permanent-escape survival."
                );
            }
        }
    }


    // ============================================================
    // Every survival definition should equal one at lag zero.
    // ============================================================

    void CheckLagZero(
        const TubeSurvivalFunction& result,
        const std::string& label)
    {
        if (result.lag_times.empty()) {
            throw std::runtime_error(
                label + ": empty lag array."
            );
        }

        if (std::abs(result.lag_times.front()) > 1.0e-12)
        {
            throw std::runtime_error(
                label + ": first lag is not zero."
            );
        }

        const std::size_t num_lags =
            result.lag_times.size();

        constexpr double tolerance = 1.0e-12;

        for (std::size_t diameter_index = 0;
            diameter_index < result.tube_diameters.size();
            ++diameter_index)
        {
            const std::size_t index =
                diameter_index * num_lags;

            if (std::abs(result.survival[index] - 1.0) >
                tolerance)
            {
                throw std::runtime_error(
                    label +
                    ": survival at lag zero is not unity."
                );
            }
        }
    }


    // ============================================================
    // Output all segment-level and tube-integrated functions.
    // ============================================================

    void WriteProductionOutputs(
        const std::string& prefix,
        const SegmentSurvivalFunction& raw,
        const SegmentSurvivalFunction& affine,
        const SegmentSurvivalFunction& history,
        const SegmentSurvivalFunction& permanent)
    {
        const TubeSurvivalFunction tube_raw =
            ComputeTubeSurvivalFunction(raw);

        const TubeSurvivalFunction tube_affine =
            ComputeTubeSurvivalFunction(affine);

        const TubeSurvivalFunction tube_history =
            ComputeTubeSurvivalFunction(history);

        const TubeSurvivalFunction tube_permanent =
            ComputeTubeSurvivalFunction(permanent);


        CheckLagZero(
            tube_raw,
            "Raw instantaneous tube survival"
        );

        CheckLagZero(
            tube_affine,
            "Affine instantaneous tube survival"
        );

        CheckLagZero(
            tube_history,
            "History-dependent tube survival"
        );

        CheckLagZero(
            tube_permanent,
            "Permanent-escape tube survival"
        );


        WriteSegmentSurvivalFunction(
            raw,
            prefix + "_segment_raw.csv"
        );

        WriteSegmentSurvivalFunction(
            affine,
            prefix + "_segment_affine.csv"
        );

        WriteSegmentSurvivalFunction(
            history,
            prefix + "_segment_history.csv"
        );

        WriteSegmentSurvivalFunction(
            permanent,
            prefix + "_segment_permanent_escape.csv"
        );


        WriteTubeSurvivalFunction(
            tube_raw,
            prefix + "_tube_raw.csv"
        );

        WriteTubeSurvivalFunction(
            tube_affine,
            prefix + "_tube_affine.csv"
        );

        WriteTubeSurvivalFunction(
            tube_history,
            prefix + "_tube_history.csv"
        );

        WriteTubeSurvivalFunction(
            tube_permanent,
            prefix + "_tube_permanent_escape.csv"
        );
    }

} // namespace



int main(
    int argc,
    char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int num_ranks = 1;

    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &rank
    );

    MPI_Comm_size(
        MPI_COMM_WORLD,
        &num_ranks
    );


    // ========================================================
    // Command line
    // ========================================================

    if (argc != 5)
    {
        if (rank == 0)
        {
            std::cerr
                << "Usage:\n"
                << "  Parallel_Segment_Survival "
                << "<Z1+SP.dat> "
                << "<tube_diameters> "
                << "<box_center_x,y,z> "
                << "<output_prefix>\n\n"
                << "Example:\n"
                << "  mpirun -np 16 "
                << "Parallel_Segment_Survival "
                << "Z1+SP.dat "
                << "\"5.0,7.0,9.0,9.2,9.4\" "
                << "\"20.6116,20.6103,39.226586\" "
                << "AMP_7.6_FREQ_1750_TRIAL_1\n";
        }

        MPI_Finalize();
        return 1;
    }


    try
    {
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
            throw std::runtime_error(
                "No tube diameters were supplied."
            );
        }

        for (const double diameter : tube_diameters)
        {
            if (diameter <= 0.0)
            {
                throw std::runtime_error(
                    "Tube diameters must be positive."
                );
            }
        }


        // ====================================================
        // Rank 0 parses the trajectory exactly once.
        // ====================================================

        std::vector<Box> frame_boxes;

        std::vector<ChainTrajectory> all_chains;

        std::uint64_t total_chain_count = 0;

        double parse_seconds = 0.0;


        if (rank == 0)
        {
            std::cout
                << "\n"
                << "========================================\n"
                << "MPI Z1+ Tube Survival Production Run\n"
                << "========================================\n"
                << "MPI ranks      : "
                << num_ranks
                << '\n'
                << "Input file     : "
                << z1_file
                << '\n'
                << "Output prefix  : "
                << output_prefix
                << '\n'
                << "Tube diameters : ";

            for (const double diameter :
            tube_diameters)
            {
                std::cout
                    << diameter
                    << ' ';
            }

            std::cout
                << "\n\n"
                << "Parsing trajectory on rank 0...\n";


            const double parse_start =
                MPI_Wtime();

            PrimitivePathTrajectory trajectory =
                parse_z1_file(z1_file);

            frame_boxes =
                std::move(
                    trajectory.frame_boxes
                );

            all_chains =
                std::move(
                    trajectory.chains
                );


            if (all_chains.empty())
            {
                throw std::runtime_error(
                    "Parsed trajectory contains zero chains."
                );
            }

            if (frame_boxes.empty())
            {
                throw std::runtime_error(
                    "Parsed trajectory contains zero frame boxes."
                );
            }


            // ------------------------------------------------
            // Z1+ stores box lengths but the origin must be
            // reconstructed from the known fixed box center.
            // Do this before broadcasting boxes.
            // ------------------------------------------------

            SetFixedBoxCenter(
                frame_boxes,
                box_center
            );


            // ------------------------------------------------
            // Saved-frame index becomes the trajectory time.
            //
            // Therefore CSV lag_time is currently in saved
            // frame units:
            //
            //     0, 1, 2, ...
            // ------------------------------------------------

            AssignUniformTimesteps(
                all_chains,
                0,
                1
            );


            total_chain_count =
                static_cast<std::uint64_t>(
                    all_chains.size()
                    );

            parse_seconds =
                MPI_Wtime() -
                parse_start;


            std::cout
                << "Parsed "
                << total_chain_count
                << " chains and "
                << frame_boxes.size()
                << " frames.\n"
                << "Parse time: "
                << std::fixed
                << std::setprecision(3)
                << parse_seconds
                << " s\n";
        }


        // ====================================================
        // Tell every rank how many chains exist before scatter.
        // ====================================================

        MPI_Bcast(
            &total_chain_count,
            1,
            MPI_UINT64_T,
            0,
            MPI_COMM_WORLD
        );


        if (total_chain_count <
            static_cast<std::uint64_t>(num_ranks))
        {
            throw std::runtime_error(
                "More MPI ranks were requested than chains. "
                "Use no more than one rank per chain."
            );
        }


        // ====================================================
        // Broadcast simulation boxes.
        // ====================================================

        BroadcastFrameBoxes(
            frame_boxes,
            0,
            MPI_COMM_WORLD
        );


        const std::size_t num_frames =
            frame_boxes.size();

        const std::vector<std::size_t> lag_frames =
            BuildProductionLagFrames(
                num_frames
            );


        // ====================================================
        // Scatter chains.
        //
        // Rank 0 keeps its own contiguous block and sends the
        // remaining blocks to the other ranks.
        // ====================================================

        std::vector<ChainTrajectory> local_chains =
            ScatterChainTrajectories(
                all_chains,
                0,
                MPI_COMM_WORLD
            );


        // Rank 0 no longer needs the complete trajectory after
        // distribution. Release that memory before calculation.
        if (rank == 0)
        {
            std::vector<ChainTrajectory>()
                .swap(all_chains);
        }


        const std::uint64_t local_chain_count =
            static_cast<std::uint64_t>(
                local_chains.size()
                );


        if (local_chain_count == 0)
        {
            throw std::runtime_error(
                "An MPI rank received zero chains."
            );
        }


        // ====================================================
        // Verify distributed chain count.
        // ====================================================

        std::uint64_t distributed_chain_count = 0;

        MPI_Reduce(
            &local_chain_count,
            &distributed_chain_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            0,
            MPI_COMM_WORLD
        );


        // ====================================================
        // Gather ownership information for a clean root-only
        // diagnostic instead of interleaved rank output.
        // ====================================================

        const std::array<std::uint64_t, 3>
            local_ownership{
                local_chain_count,
                static_cast<std::uint64_t>(
                    local_chains.front().chain_id
                ),
                static_cast<std::uint64_t>(
                    local_chains.back().chain_id
                )
        };


        std::vector<std::uint64_t>
            ownership_table;

        if (rank == 0)
        {
            ownership_table.resize(
                static_cast<std::size_t>(
                    3 * num_ranks
                    )
            );
        }


        MPI_Gather(
            local_ownership.data(),
            3,
            MPI_UINT64_T,

            rank == 0
            ? ownership_table.data()
            : nullptr,

            3,
            MPI_UINT64_T,

            0,
            MPI_COMM_WORLD
        );


        if (rank == 0)
        {
            if (distributed_chain_count !=
                total_chain_count)
            {
                throw std::runtime_error(
                    "Distributed chain count does not "
                    "match parsed chain count."
                );
            }


            std::cout
                << "\nChain ownership:\n";

            for (int r = 0;
                r < num_ranks;
                ++r)
            {
                const std::size_t offset =
                    static_cast<std::size_t>(
                        3 * r
                        );

                std::cout
                    << "  Rank "
                    << r
                    << ": "
                    << ownership_table[offset]
                    << " chains ["
                    << ownership_table[offset + 1]
                    << " ... "
                    << ownership_table[offset + 2]
                    << "]\n";
            }


            std::cout
                << "\nProduction lags ("
                << lag_frames.size()
                << " total):\n  ";

            for (std::size_t i = 0;
                i < lag_frames.size();
                ++i)
            {
                std::cout
                    << lag_frames[i];

                if (i + 1 != lag_frames.size()) {
                    std::cout << ' ';
                }

                if ((i + 1) % 16 == 0 &&
                    i + 1 != lag_frames.size())
                {
                    std::cout
                        << "\n  ";
                }
            }

            std::cout
                << "\n\nStarting distributed production calculation...\n";
        }


        // ====================================================
        // Synchronize before timing the expensive portion.
        // ====================================================

        MPI_Barrier(
            MPI_COMM_WORLD
        );

        const double distributed_start =
            MPI_Wtime();


        // ====================================================
        // Rank-local accumulators.
        //
        // IMPORTANT:
        // These remain UNFINALIZED until after MPI_Reduce.
        // ====================================================

        SegmentSurvivalFunction local_raw;
        SegmentSurvivalFunction local_affine;
        SegmentSurvivalFunction local_history;
        SegmentSurvivalFunction local_permanent;

        bool initialized = false;


        const SegmentSurvivalOptions affine_options{
            .apply_affine_correction = true
        };


        // ====================================================
        // Main chain loop.
        // ====================================================

        for (std::size_t chain_index = 0;
            chain_index < local_chains.size();
            ++chain_index)
        {
            const ChainTrajectory& chain =
                local_chains[chain_index];


            // ------------------------------------------------
            // 1. Raw instantaneous survival.
            // ------------------------------------------------

            const SegmentSurvivalFunction raw =
                ComputeSegmentSurvivalFunction(
                    chain,
                    lag_frames,
                    tube_diameters
                );


            // ------------------------------------------------
            // 2. Affine-corrected instantaneous survival.
            // ------------------------------------------------

            const SegmentSurvivalFunction affine =
                ComputeSegmentSurvivalFunction(
                    chain,
                    frame_boxes,
                    lag_frames,
                    tube_diameters,
                    affine_options
                );


            // ------------------------------------------------
            // 3 + 4.
            //
            // One history scan produces:
            //
            //     Stephanou-like history survival
            //     pure permanent-escape survival
            //
            // false -> use the normal varying lag-wise origin
            // cohort rather than one common origin cohort.
            // ------------------------------------------------

            const HistoryDependentSurvivalResult
                history_result =
                ComputeHistoryDependentSurvivalFunction(
                    chain,
                    frame_boxes,
                    lag_frames,
                    tube_diameters,
                    affine_options,
                    false
                );


            const SegmentSurvivalFunction& history =
                history_result.stephanou_survival;

            const SegmentSurvivalFunction& permanent =
                history_result
                .permanent_escape_survival;


            CheckMatchingCohorts(
                raw,
                affine,
                history,
                permanent
            );


            // ------------------------------------------------
            // Initialize empty weighted accumulators using the
            // first real result as the shape/template.
            // ------------------------------------------------

            if (!initialized)
            {
                InitializeAccumulator(
                    local_raw,
                    raw
                );

                InitializeAccumulator(
                    local_affine,
                    affine
                );

                InitializeAccumulator(
                    local_history,
                    history
                );

                InitializeAccumulator(
                    local_permanent,
                    permanent
                );

                initialized = true;
            }


            // ------------------------------------------------
            // Add normalized chain results back as weighted
            // origin counts.
            //
            // DO NOT FinalizeAccumulator here.
            // ------------------------------------------------

            AccumulateSegmentSurvival(
                local_raw,
                raw
            );

            AccumulateSegmentSurvival(
                local_affine,
                affine
            );

            AccumulateSegmentSurvival(
                local_history,
                history
            );

            AccumulateSegmentSurvival(
                local_permanent,
                permanent
            );


            // Keep production output relatively quiet.
            // Rank 0 reports its own progress as a proxy for
            // the approximately even chain distribution.
            if (rank == 0)
            {
                const std::size_t done =
                    chain_index + 1;

                if (done == 1 ||
                    done % 10 == 0 ||
                    done == local_chains.size())
                {
                    const double elapsed =
                        MPI_Wtime() -
                        distributed_start;

                    std::cout
                        << "Rank 0: finished "
                        << done
                        << " / "
                        << local_chains.size()
                        << " local chains"
                        << "  (elapsed "
                        << std::fixed
                        << std::setprecision(1)
                        << elapsed
                        << " s)\n";
                }
            }
        }


        if (!initialized)
        {
            throw std::runtime_error(
                "Local survival accumulator was never initialized."
            );
        }


        const double local_calculation_seconds =
            MPI_Wtime() -
            distributed_start;


        // ====================================================
        // Gather rank-local compute times for diagnostics.
        // ====================================================

        std::vector<double>
            rank_calculation_seconds;

        if (rank == 0)
        {
            rank_calculation_seconds.resize(
                static_cast<std::size_t>(
                    num_ranks
                    )
            );
        }


        MPI_Gather(
            &local_calculation_seconds,
            1,
            MPI_DOUBLE,

            rank == 0
            ? rank_calculation_seconds.data()
            : nullptr,

            1,
            MPI_DOUBLE,

            0,
            MPI_COMM_WORLD
        );


        // ====================================================
        // Global MPI reduction.
        //
        // These four calls sum:
        //
        //     weighted survival numerators
        //     sample_counts
        //
        // onto rank 0.
        //
        // Local accumulators are deliberately NOT finalized.
        // ====================================================

        SegmentSurvivalFunction global_raw;
        SegmentSurvivalFunction global_affine;
        SegmentSurvivalFunction global_history;
        SegmentSurvivalFunction global_permanent;


        ReduceSurvivalAccumulator(
            local_raw,
            global_raw,
            0,
            MPI_COMM_WORLD
        );

        ReduceSurvivalAccumulator(
            local_affine,
            global_affine,
            0,
            MPI_COMM_WORLD
        );

        ReduceSurvivalAccumulator(
            local_history,
            global_history,
            0,
            MPI_COMM_WORLD
        );

        ReduceSurvivalAccumulator(
            local_permanent,
            global_permanent,
            0,
            MPI_COMM_WORLD
        );


        // ====================================================
        // Rank 0 normalizes exactly once and writes output.
        // ====================================================

        if (rank == 0)
        {
            const double distributed_seconds =
                MPI_Wtime() -
                distributed_start;


            FinalizeAccumulator(
                global_raw
            );

            FinalizeAccumulator(
                global_affine
            );

            FinalizeAccumulator(
                global_history
            );

            FinalizeAccumulator(
                global_permanent
            );


            CheckMatchingCohorts(
                global_raw,
                global_affine,
                global_history,
                global_permanent
            );


            CheckHistoryRelations(
                global_affine,
                global_history,
                global_permanent
            );


            std::cout
                << "\nDistributed calculation complete.\n"
                << "Writing production CSV files...\n";


            WriteProductionOutputs(
                output_prefix,
                global_raw,
                global_affine,
                global_history,
                global_permanent
            );


            const double min_rank_time =
                *std::min_element(
                    rank_calculation_seconds.begin(),
                    rank_calculation_seconds.end()
                );

            const double max_rank_time =
                *std::max_element(
                    rank_calculation_seconds.begin(),
                    rank_calculation_seconds.end()
                );

            double mean_rank_time = 0.0;

            for (const double value :
            rank_calculation_seconds)
            {
                mean_rank_time += value;
            }

            mean_rank_time /=
                static_cast<double>(
                    rank_calculation_seconds.size()
                    );


            std::cout
                << "\n========================================\n"
                << "MPI production calculation complete\n"
                << "========================================\n"
                << "Chains analyzed        : "
                << total_chain_count
                << '\n'
                << "Frames                 : "
                << num_frames
                << '\n'
                << "Requested lags         : "
                << lag_frames.size()
                << '\n'
                << "MPI ranks              : "
                << num_ranks
                << '\n'
                << "Parse time             : "
                << std::fixed
                << std::setprecision(3)
                << parse_seconds
                << " s\n"
                << "Distributed wall time  : "
                << distributed_seconds
                << " s\n"
                << "Fastest rank compute   : "
                << min_rank_time
                << " s\n"
                << "Mean rank compute      : "
                << mean_rank_time
                << " s\n"
                << "Slowest rank compute   : "
                << max_rank_time
                << " s\n"
                << "Lag-0 total origins    : "
                << global_raw.sample_counts.front()
                << '\n'
                << "\nOutputs:\n"
                << "  "
                << output_prefix
                << "_segment_raw.csv\n"
                << "  "
                << output_prefix
                << "_segment_affine.csv\n"
                << "  "
                << output_prefix
                << "_segment_history.csv\n"
                << "  "
                << output_prefix
                << "_segment_permanent_escape.csv\n"
                << "  "
                << output_prefix
                << "_tube_raw.csv\n"
                << "  "
                << output_prefix
                << "_tube_affine.csv\n"
                << "  "
                << output_prefix
                << "_tube_history.csv\n"
                << "  "
                << output_prefix
                << "_tube_permanent_escape.csv\n"
                << "========================================\n";
        }


        MPI_Finalize();
        return 0;
    }


    // ========================================================
    // If any rank encounters an exception, abort the whole MPI
    // job rather than leaving the other ranks blocked in a
    // collective operation.
    // ========================================================

    catch (const std::exception& error)
    {
        std::cerr
            << "\nERROR on MPI rank "
            << rank
            << ": "
            << error.what()
            << '\n';

        MPI_Abort(
            MPI_COMM_WORLD,
            1
        );

        return 1;
    }
}