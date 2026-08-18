#include "MPI_Comms.hpp"

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <vector>


bool NearlyEqual(double a, double b, double tol = 1.0e-12)
{
    return std::abs(a - b) < tol;
}


ChainTrajectory MakeTestChain(std::size_t chain_id)
{
    ChainTrajectory chain;

    chain.chain_id = chain_id;

    chain.timesteps = {
        0,
        1
    };

    chain.frame_offsets = {
        0,
        2,
        4
    };

    const double x =
        static_cast<double>(chain_id);

    chain.nodes = {
        PPNode{x + 0.0, 1.0, 2.0, true},
        PPNode{x + 0.1, 1.1, 2.1, false},
        PPNode{x + 0.2, 1.2, 2.2, false},
        PPNode{x + 0.3, 1.3, 2.3, true}
    };

    return chain;
}


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    int num_ranks;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);


    if (num_ranks != 2)
    {
        if (rank == 0) {
            std::cerr
                << "Smoke test must be run with exactly 2 ranks.\n";
        }

        MPI_Finalize();
        return 1;
    }


    bool passed = true;


    // ========================================================
    // Test 1: BroadcastFrameBoxes
    // ========================================================

    std::vector<Box> frame_boxes;

    if (rank == 0)
    {
        frame_boxes.resize(2);

        frame_boxes[0].origin =
            Vec3{ 1.0, 2.0, 3.0 };

        frame_boxes[1].origin =
            Vec3{ 4.0, 5.0, 6.0 };


        for (std::size_t row = 0; row < 3; ++row)
        {
            for (std::size_t col = 0; col < 3; ++col)
            {
                frame_boxes[0].matrix.value[row][col] =
                    10.0 * row + col;

                frame_boxes[1].matrix.value[row][col] =
                    100.0 + 10.0 * row + col;
            }
        }

        frame_boxes[0].affine_geometry_available = true;
        frame_boxes[1].affine_geometry_available = false;
    }


    BroadcastFrameBoxes(
        frame_boxes,
        0,
        MPI_COMM_WORLD
    );


    if (frame_boxes.size() != 2) {
        passed = false;
    }

    if (!NearlyEqual(frame_boxes[0].origin.x, 1.0) ||
        !NearlyEqual(frame_boxes[0].origin.y, 2.0) ||
        !NearlyEqual(frame_boxes[0].origin.z, 3.0))
    {
        passed = false;
    }

    if (!NearlyEqual(
        frame_boxes[1].matrix.value[2][1],
        121.0))
    {
        passed = false;
    }

    if (!frame_boxes[0].affine_geometry_available ||
        frame_boxes[1].affine_geometry_available)
    {
        passed = false;
    }


    // ========================================================
    // Test 2: ScatterChainTrajectories
    // ========================================================

    std::vector<ChainTrajectory> all_chains;

    if (rank == 0)
    {
        all_chains.push_back(MakeTestChain(100));
        all_chains.push_back(MakeTestChain(101));
        all_chains.push_back(MakeTestChain(102));
        all_chains.push_back(MakeTestChain(103));
    }


    const std::vector<ChainTrajectory> local_chains =
        ScatterChainTrajectories(
            all_chains,
            0,
            MPI_COMM_WORLD
        );


    if (local_chains.size() != 2) {
        passed = false;
    }

    const std::size_t expected_first_id =
        rank == 0 ? 100 : 102;


    if (local_chains.size() == 2)
    {
        if (local_chains[0].chain_id != expected_first_id ||
            local_chains[1].chain_id != expected_first_id + 1)
        {
            passed = false;
        }

        if (local_chains[0].timesteps.size() != 2 ||
            local_chains[0].frame_offsets.size() != 3 ||
            local_chains[0].nodes.size() != 4)
        {
            passed = false;
        }

        const double expected_x =
            static_cast<double>(expected_first_id) + 0.2;

        if (!NearlyEqual(
            local_chains[0].nodes[2].x,
            expected_x))
        {
            passed = false;
        }
    }


    // ========================================================
    // Test 3: ReduceSurvivalAccumulator
    // ========================================================

    SegmentSurvivalFunction local;
    SegmentSurvivalFunction global;

    local.lag_times = { 0 };
    local.tube_diameters = { 9.0 };

    local.survival.resize(NUM_SAMPLE_POINTS);
    local.sample_counts = { 1 };


    const double local_value =
        static_cast<double>(rank + 1);

    for (double& value : local.survival) {
        value = local_value;
    }


    ReduceSurvivalAccumulator(
        local,
        global,
        0,
        MPI_COMM_WORLD
    );


    if (rank == 0)
    {
        // Rank 0 contributed 1.0
        // Rank 1 contributed 2.0
        // MPI_SUM should therefore give 3.0.

        if (global.sample_counts[0] != 2) {
            passed = false;
        }

        for (double value : global.survival)
        {
            if (!NearlyEqual(value, 3.0)) {
                passed = false;
            }
        }
    }


    // ========================================================
    // Combine test status from both ranks
    // ========================================================

    int local_pass =
        passed ? 1 : 0;

    int global_pass = 0;

    MPI_Allreduce(
        &local_pass,
        &global_pass,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD
    );


    if (rank == 0)
    {
        if (global_pass)
        {
            std::cout
                << "MPI communications smoke test: PASS\n";
        }
        else
        {
            std::cout
                << "MPI communications smoke test: FAIL\n";
        }
    }


    MPI_Finalize();

    return global_pass ? 0 : 1;
}