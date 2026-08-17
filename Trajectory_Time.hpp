#pragma once

#include "Primitive_Path_Trajectory.hpp"
#include <vector>


void AssignUniformTimesteps(
    std::vector<ChainTrajectory>& trajectories,
    long initial_timestep,
    long dump_interval
);


void AssignTimesteps(
    std::vector<ChainTrajectory>& trajectories,
    const std::vector<long>& frame_timesteps
);
