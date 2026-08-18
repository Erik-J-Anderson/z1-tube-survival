#pragma once

#include <array>
#include <cstddef>
#include <vector>


struct PPNode
{
    double x;
    double y;
    double z;
    bool end;
};


struct Segment
{
    double x;
    double y;
    double z;
};


struct Vec3
{
    double x{};
    double y{};
    double z{};
};


struct Mat3
{
    double value[3][3]{};
};


struct Box
{
    Vec3 origin{};
    Mat3 matrix{};

    // False for legacy Z1+ records containing only Lx Ly Lz.
    // Such records cannot describe shear tilt and must not be used
    // for affine correction of a sheared trajectory.
    bool affine_geometry_available{false};
};


struct ChainTrajectory
{
    std::size_t chain_id{};
    std::vector<long> timesteps;
    std::vector<std::size_t> frame_offsets;
    std::vector<PPNode> nodes;
};


struct PrimitivePathTrajectory
{
    std::vector<Box> frame_boxes;
    std::vector<ChainTrajectory> chains;
};


inline constexpr std::size_t NUM_SAMPLE_POINTS = 101;

using SampledPath = std::array<Segment, NUM_SAMPLE_POINTS>;


// Sample NUM_SAMPLE_POINTS positions uniformly in contour length
// along one primitive-path frame, including both chain ends.
SampledPath SamplePrimitivePath(
    const ChainTrajectory& trajectory,
    std::size_t frame_index
);
