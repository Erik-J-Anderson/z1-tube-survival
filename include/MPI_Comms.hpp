#pragma once

#include <mpi.h>

#include "Primitive_Path_Trajectory.hpp"
#include "Tube_Survival.hpp"

#include <vector>


void BroadcastFrameBoxes(
    std::vector<Box>& frame_boxes,
    int root_rank,
    MPI_Comm comm
);

std::vector<ChainTrajectory> ScatterChainTrajectories(
    const std::vector<ChainTrajectory>& all_chains,
    int root_rank,
    MPI_Comm comm
);

void SendChainTrajectory(
    const ChainTrajectory& chain,
    int dest_rank,
    MPI_Comm comm
);

ChainTrajectory ReceiveChainTrajectory(
    int source_rank,
    MPI_Comm comm
);

void ReduceSurvivalAccumulator(
    const SegmentSurvivalFunction& local,
    SegmentSurvivalFunction& global,
    int root_rank,
    MPI_Comm comm
);