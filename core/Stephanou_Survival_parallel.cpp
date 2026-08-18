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
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{

// ------------------------------------------------------------
// Short lag schedule for the MPI accumulator validation.
//
// We deliberately keep this small so that the purpose of this
// executable is testing MPI correctness rather than doing the
// full production calculation.
// ------------------------------------------------------------
std::vector<std::size_t> BuildTestLagFrames(
    std::size_t num_frames)
{
    const std::vector<std::size_t> requested{
        0,
        1,
        2,
        4,
        8,
        16
    };

    std::vector<std::size_t> result;

    for (const std::size_t lag : requested)
    {
        if (lag < num_frames) {
            result.push_back(lag);
        }
    }

    return result;
}


// ------------------------------------------------------------
// Compare two finalized survival functions.
//
// Returns the largest absolute difference in phi.
// ------------------------------------------------------------
double CompareSurvivalFunctions(
    const SegmentSurvivalFunction& a,
    const SegmentSurvivalFunction& b)
{
    if (a.lag_times != b.lag_times)
    {
        throw std::runtime_error(
            "MPI/serial lag-time arrays do not match."
        );
    }

    if (a.tube_diameters != b.tube_diameters)
    {
        throw std::runtime_error(
            "MPI/serial tube-diameter arrays do not match."
        );
    }

    if (a.sample_counts != b.sample_counts)
    {
        throw std::runtime_error(
            "MPI/serial sample-count arrays do not match."
        );
    }

    if (a.survival.size() != b.survival.size())
    {
        throw std::runtime_error(
            "MPI/serial survival-array sizes do not match."
        );
    }


    double max_abs_difference = 0.0;

    for (std::size_t i = 0;
         i < a.survival.size();
         ++i)
    {
        max_abs_difference =
            std::max(
                max_abs_difference,
                std::abs(
                    a.survival[i] -
                    b.survival[i]
                )
            );
    }

    return max_abs_difference;
}


// ------------------------------------------------------------
// Check that lag-zero survival is unity.
// ------------------------------------------------------------
void CheckLagZero(
    const SegmentSurvivalFunction& result)
{
    if (result.lag_times.empty()) {
        throw std::runtime_error(
            "Lag-time array is empty."
        );
    }

    if (result.lag_times.front() != 0.0) {
        throw std::runtime_error(
            "First lag is not zero."
        );
    }


    const std::size_t num_lags =
        result.lag_times.size();


    for (std::size_t diameter_index = 0;
         diameter_index < result.tube_diameters.size();
         ++diameter_index)
    {
        for (std::size_t segment_index = 0;
             segment_index < NUM_SAMPLE_POINTS;
             ++segment_index)
        {
            const std::size_t flat_index =
                diameter_index *
                    num_lags *
                    NUM_SAMPLE_POINTS
                + segment_index;


            if (std::abs(
                    result.survival[flat_index] -
                    1.0) > 1.0e-12)
            {
                throw std::runtime_error(
                    "Lag-zero segment survival is not unity."
                );
            }
        }
    }
}

} // namespace


int main(int argc, char* argv[])
{
    MPI_Init(
        &argc,
        &argv
    );


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


    try
    {
        // ------------------------------------------------------------
        // Same command-line interface as the serial production driver.
        // ------------------------------------------------------------
        if (argc != 5)
        {
            if (rank == 0)
            {
                std::cerr
                    << "Usage:\n"
                    << "  " << argv[0]
                    << " <Z1+SP.dat>"
                    << " <tube_diameters>"
                    << " <box_center_x,y,z>"
                    << " <output_prefix>\n\n"

                    << "Example:\n"
                    << "  mpirun -np 4 "
                    << argv[0]
                    << " Z1+SP.dat"
                    << " \"5.0,7.0,9.0,9.2,9.4\""
                    << " \"20.6116,20.6103,39.226586\""
                    << " MPI_ACCUM_TEST\n";
            }

            MPI_Finalize();
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
                "At least one tube diameter is required."
            );
        }


        // ------------------------------------------------------------
        // Global data exists only on rank 0 before distribution.
        //
        // IMPORTANT:
        // For THIS TEST we intentionally keep all_chains alive on
        // rank 0 after scattering because rank 0 later computes an
        // independent serial reference.
        // ------------------------------------------------------------
        std::vector<Box> frame_boxes;
        std::vector<ChainTrajectory> all_chains;

        std::uint64_t total_chain_count = 0;


        if (rank == 0)
        {
            std::cout
                << "\nMPI accumulator validation test\n"
                << "--------------------------------\n"
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

            for (const double diameter : tube_diameters) {
                std::cout << diameter << ' ';
            }

            std::cout
                << "\n\nParsing trajectory on rank 0...\n";


            PrimitivePathTrajectory parsed =
                parse_z1_file(z1_file);


            frame_boxes =
                std::move(parsed.frame_boxes);

            all_chains =
                std::move(parsed.chains);


            if (frame_boxes.empty())
            {
                throw std::runtime_error(
                    "Parser returned no frame boxes."
                );
            }

            if (all_chains.empty())
            {
                throw std::runtime_error(
                    "Parser returned no chains."
                );
            }


            SetFixedBoxCenter(
                frame_boxes,
                box_center
            );


            AssignUniformTimesteps(
                all_chains,
                0,
                1
            );


            total_chain_count =
                static_cast<std::uint64_t>(
                    all_chains.size()
                );


            std::cout
                << "Parsed "
                << all_chains.size()
                << " chains and "
                << frame_boxes.size()
                << " frames.\n";
        }


        // ------------------------------------------------------------
        // Broadcast global chain count.
        // ------------------------------------------------------------
        MPI_Bcast(
            &total_chain_count,
            1,
            MPI_UINT64_T,
            0,
            MPI_COMM_WORLD
        );


        // ------------------------------------------------------------
        // Broadcast box geometry.
        // ------------------------------------------------------------
        BroadcastFrameBoxes(
            frame_boxes,
            0,
            MPI_COMM_WORLD
        );


        // ------------------------------------------------------------
        // Build identical test lag schedule on every rank.
        // ------------------------------------------------------------
        const std::vector<std::size_t> lag_frames =
            BuildTestLagFrames(
                frame_boxes.size()
            );


        if (lag_frames.empty())
        {
            throw std::runtime_error(
                "No valid test lags were generated."
            );
        }


        // ------------------------------------------------------------
        // Scatter chains.
        // ------------------------------------------------------------
        std::vector<ChainTrajectory> local_chains =
            ScatterChainTrajectories(
                all_chains,
                0,
                MPI_COMM_WORLD
            );


        if (local_chains.empty())
        {
            throw std::runtime_error(
                "MPI accumulator test requires at least "
                "one chain on every rank."
            );
        }


        // ------------------------------------------------------------
        // Report ownership cleanly through rank 0.
        //
        // Each rank sends:
        //
        //     [number of local chains,
        //      first chain id,
        //      last chain id]
        // ------------------------------------------------------------
        const std::array<std::uint64_t, 3> local_info{
            static_cast<std::uint64_t>(
                local_chains.size()
            ),
            static_cast<std::uint64_t>(
                local_chains.front().chain_id
            ),
            static_cast<std::uint64_t>(
                local_chains.back().chain_id
            )
        };


        std::vector<std::uint64_t> gathered_info;

        if (rank == 0)
        {
            gathered_info.resize(
                static_cast<std::size_t>(num_ranks) * 3
            );
        }


        MPI_Gather(
            local_info.data(),
            3,
            MPI_UINT64_T,

            rank == 0
                ? gathered_info.data()
                : nullptr,

            3,
            MPI_UINT64_T,
            0,
            MPI_COMM_WORLD
        );


        if (rank == 0)
        {
            std::cout
                << "\nChain ownership:\n";

            for (int r = 0;
                 r < num_ranks;
                 ++r)
            {
                const std::size_t base =
                    static_cast<std::size_t>(r) * 3;

                std::cout
                    << "  Rank "
                    << r
                    << ": "
                    << gathered_info[base + 0]
                    << " chains ["
                    << gathered_info[base + 1]
                    << " ... "
                    << gathered_info[base + 2]
                    << "]\n";
            }


            std::cout
                << "\nTest lags: ";

            for (const std::size_t lag : lag_frames) {
                std::cout << lag << ' ';
            }

            std::cout
                << "\n\nStarting distributed calculation...\n";
        }


        MPI_Barrier(
            MPI_COMM_WORLD
        );


        const auto parallel_start =
            std::chrono::steady_clock::now();


        // ============================================================
        // MPI PART
        //
        // Each rank independently builds an UNFINALIZED weighted
        // accumulator from its own chains.
        // ============================================================

        SegmentSurvivalFunction local_accumulator;

        bool local_initialized = false;


        for (std::size_t local_index = 0;
             local_index < local_chains.size();
             ++local_index)
        {
            const ChainTrajectory& chain =
                local_chains[local_index];


            const SegmentSurvivalFunction chain_result =
                ComputeSegmentSurvivalFunction(
                    chain,
                    lag_frames,
                    tube_diameters
                );


            if (!local_initialized)
            {
                InitializeAccumulator(
                    local_accumulator,
                    chain_result
                );

                local_initialized = true;
            }


            AccumulateSegmentSurvival(
                local_accumulator,
                chain_result
            );


            const std::size_t done =
                local_index + 1;


            if (done == 1 ||
                done % 50 == 0 ||
                done == local_chains.size())
            {
                std::cout
                    << "Rank "
                    << rank
                    << ": finished "
                    << done
                    << " / "
                    << local_chains.size()
                    << " local chains\n"
                    << std::flush;
            }
        }


        if (!local_initialized)
        {
            throw std::runtime_error(
                "Local accumulator was never initialized."
            );
        }


        // ------------------------------------------------------------
        // CRITICAL:
        //
        // DO NOT call FinalizeAccumulator(local_accumulator) here.
        //
        // local_accumulator currently contains:
        //
        //     sum(phi_chain * chain_sample_count)
        //
        // and:
        //
        //     sum(chain_sample_count)
        //
        // Those are exactly the quantities MPI_Reduce must add.
        // ------------------------------------------------------------


        SegmentSurvivalFunction mpi_accumulator;


        ReduceSurvivalAccumulator(
            local_accumulator,
            mpi_accumulator,
            0,
            MPI_COMM_WORLD
        );


        // ------------------------------------------------------------
        // Only after reduction do we normalize the global result.
        // ------------------------------------------------------------
        if (rank == 0)
        {
            FinalizeAccumulator(
                mpi_accumulator
            );

            CheckLagZero(
                mpi_accumulator
            );
        }


        MPI_Barrier(
            MPI_COMM_WORLD
        );


        const auto parallel_stop =
            std::chrono::steady_clock::now();


        // ============================================================
        // SERIAL REFERENCE
        //
        // Rank 0 now calculates the exact same observable using all
        // chains in its original trajectory.
        //
        // This gives us a direct end-to-end validation of:
        //
        //     scatter
        //       +
        //     rank-local accumulation
        //       +
        //     MPI reduction
        //       +
        //     final normalization
        // ============================================================

        int test_pass = 1;


        if (rank == 0)
        {
            std::cout
                << "\nDistributed calculation complete.\n"
                << "Computing independent serial reference...\n";


            SegmentSurvivalFunction serial_accumulator;

            bool serial_initialized = false;


            for (std::size_t chain_index = 0;
                 chain_index < all_chains.size();
                 ++chain_index)
            {
                const SegmentSurvivalFunction chain_result =
                    ComputeSegmentSurvivalFunction(
                        all_chains[chain_index],
                        lag_frames,
                        tube_diameters
                    );


                if (!serial_initialized)
                {
                    InitializeAccumulator(
                        serial_accumulator,
                        chain_result
                    );

                    serial_initialized = true;
                }


                AccumulateSegmentSurvival(
                    serial_accumulator,
                    chain_result
                );


                const std::size_t done =
                    chain_index + 1;


                if (done == 1 ||
                    done % 100 == 0 ||
                    done == all_chains.size())
                {
                    std::cout
                        << "Serial reference: finished "
                        << done
                        << " / "
                        << all_chains.size()
                        << " chains\n";
                }
            }


            if (!serial_initialized)
            {
                throw std::runtime_error(
                    "Serial reference accumulator "
                    "was never initialized."
                );
            }


            FinalizeAccumulator(
                serial_accumulator
            );


            CheckLagZero(
                serial_accumulator
            );


            // --------------------------------------------------------
            // Compare MPI result to serial result.
            // --------------------------------------------------------
            const double max_abs_difference =
                CompareSurvivalFunctions(
                    mpi_accumulator,
                    serial_accumulator
                );


            constexpr double tolerance =
                1.0e-12;


            if (max_abs_difference > tolerance) {
                test_pass = 0;
            }


            // --------------------------------------------------------
            // Also compare integrated tube-survival results.
            // --------------------------------------------------------
            const TubeSurvivalFunction mpi_tube =
                ComputeTubeSurvivalFunction(
                    mpi_accumulator
                );


            const TubeSurvivalFunction serial_tube =
                ComputeTubeSurvivalFunction(
                    serial_accumulator
                );


            double max_tube_difference = 0.0;


            if (mpi_tube.survival.size() !=
                serial_tube.survival.size())
            {
                throw std::runtime_error(
                    "MPI/serial tube-survival dimensions differ."
                );
            }


            for (std::size_t i = 0;
                 i < mpi_tube.survival.size();
                 ++i)
            {
                max_tube_difference =
                    std::max(
                        max_tube_difference,
                        std::abs(
                            mpi_tube.survival[i] -
                            serial_tube.survival[i]
                        )
                    );
            }


            if (max_tube_difference > tolerance) {
                test_pass = 0;
            }


            // --------------------------------------------------------
            // Write both results so they can also be inspected
            // manually if desired.
            // --------------------------------------------------------
            WriteSegmentSurvivalFunction(
                mpi_accumulator,
                output_prefix +
                    "_mpi_accumulator_segment.csv"
            );


            WriteSegmentSurvivalFunction(
                serial_accumulator,
                output_prefix +
                    "_serial_reference_segment.csv"
            );


            WriteTubeSurvivalFunction(
                mpi_tube,
                output_prefix +
                    "_mpi_accumulator_tube.csv"
            );


            WriteTubeSurvivalFunction(
                serial_tube,
                output_prefix +
                    "_serial_reference_tube.csv"
            );


            const std::chrono::duration<double>
                parallel_elapsed =
                    parallel_stop -
                    parallel_start;


            std::cout
                << std::setprecision(17)
                << "\nMaximum segment-survival difference : "
                << max_abs_difference
                << '\n'
                << "Maximum tube-survival difference    : "
                << max_tube_difference
                << '\n'
                << "MPI calculation time                : "
                << parallel_elapsed.count()
                << " s\n";


            if (test_pass)
            {
                std::cout
                    << "\nMPI accumulator validation: PASS\n";
            }
            else
            {
                std::cout
                    << "\nMPI accumulator validation: FAIL\n";
            }
        }


        // ------------------------------------------------------------
        // Tell every rank whether the validation passed.
        // ------------------------------------------------------------
        MPI_Bcast(
            &test_pass,
            1,
            MPI_INT,
            0,
            MPI_COMM_WORLD
        );


        MPI_Finalize();

        return test_pass ? 0 : 1;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "Rank "
            << rank
            << " error: "
            << error.what()
            << '\n'
            << std::flush;


        MPI_Abort(
            MPI_COMM_WORLD,
            1
        );

        return 1;
    }
}
