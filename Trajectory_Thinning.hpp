#pragma once

#include <cstddef>

#include "Primitive_Path_Trajectory.hpp"


ChainTrajectory MakeStridedTrajectory(
    const ChainTrajectory& trajectory,
    std::size_t frame_stride
);

