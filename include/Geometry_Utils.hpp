#pragma once

#include <span>

#include "Primitive_Path_Trajectory.hpp"


namespace geometry
{

    Vec3 Multiply(
        const Mat3& matrix,
        const Vec3& vector);

    Mat3 Inverse(
        const Mat3& matrix);

    Vec3 CartesianToFractional(
        const Vec3& position,
        const Box& box);

    Vec3 FractionalToCartesian(
        const Vec3& fractional,
        const Box& box);

    Vec3 MapPositionBetweenBoxes(
        const Vec3& position,
        const Box& source_box,
        const Box& destination_box);

    void MapPPNodesBetweenBoxes(
        std::span<const PPNode> source_nodes,
        const Box& source_box,
        const Box& destination_box,
        std::span<PPNode> mapped_nodes);

}  // namespace geometry