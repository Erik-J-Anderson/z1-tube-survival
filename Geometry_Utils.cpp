#include "Geometry_Utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>


namespace geometry
{

    Vec3 Multiply(
        const Mat3& matrix,
        const Vec3& vector)
    {
        return Vec3{
            matrix.value[0][0] * vector.x
                + matrix.value[0][1] * vector.y
                + matrix.value[0][2] * vector.z,

            matrix.value[1][0] * vector.x
                + matrix.value[1][1] * vector.y
                + matrix.value[1][2] * vector.z,

            matrix.value[2][0] * vector.x
                + matrix.value[2][1] * vector.y
                + matrix.value[2][2] * vector.z
        };
    }


    Mat3 Inverse(const Mat3& matrix)
    {
        const double a = matrix.value[0][0];
        const double b = matrix.value[0][1];
        const double c = matrix.value[0][2];

        const double d = matrix.value[1][0];
        const double e = matrix.value[1][1];
        const double f = matrix.value[1][2];

        const double g = matrix.value[2][0];
        const double h = matrix.value[2][1];
        const double i = matrix.value[2][2];

        const double determinant =
            a * (e * i - f * h)
            - b * (d * i - f * g)
            + c * (d * h - e * g);

        const double matrix_scale = std::max({
            std::abs(a), std::abs(b), std::abs(c),
            std::abs(d), std::abs(e), std::abs(f),
            std::abs(g), std::abs(h), std::abs(i)
            });

        constexpr double relative_tolerance = 1.0e-12;

        const double determinant_tolerance =
            relative_tolerance
            * matrix_scale
            * matrix_scale
            * matrix_scale;

        if (matrix_scale == 0.0 ||
            std::abs(determinant) <= determinant_tolerance) {

            throw std::runtime_error(
                "Cannot invert a singular or nearly singular matrix.");
        }

        const double inverse_determinant = 1.0 / determinant;

        return Mat3{ {
            {
                (e * i - f * h) * inverse_determinant,
                (c * h - b * i) * inverse_determinant,
                (b * f - c * e) * inverse_determinant
            },
            {
                (f * g - d * i) * inverse_determinant,
                (a * i - c * g) * inverse_determinant,
                (c * d - a * f) * inverse_determinant
            },
            {
                (d * h - e * g) * inverse_determinant,
                (b * g - a * h) * inverse_determinant,
                (a * e - b * d) * inverse_determinant
            }
        } };
    }


    Vec3 CartesianToFractional(
        const Vec3& position,
        const Box& box)
    {
        /*
         * Fractional coordinates:
         *
         *     s = H^{-1}(r - origin)
         */
        const Vec3 relative_position{
            position.x - box.origin.x,
            position.y - box.origin.y,
            position.z - box.origin.z
        };

        const Mat3 inverse_box_matrix = Inverse(box.matrix);

        return Multiply(
            inverse_box_matrix,
            relative_position);
    }


    Vec3 FractionalToCartesian(
        const Vec3& fractional,
        const Box& box)
    {
        /*
         * Cartesian coordinates:
         *
         *     r = origin + Hs
         */
        const Vec3 relative_position =
            Multiply(box.matrix, fractional);

        return Vec3{
            box.origin.x + relative_position.x,
            box.origin.y + relative_position.y,
            box.origin.z + relative_position.z
        };
    }


    Vec3 MapPositionBetweenBoxes(
        const Vec3& position,
        const Box& source_box,
        const Box& destination_box)
    {
        /*
         * Preserve the position's fractional coordinates while
         * expressing it in the destination box:
         *
         *     r_destination =
         *         origin_destination
         *       + H_destination H_source^{-1}
         *         (r_source - origin_source)
         */
        const Vec3 fractional =
            CartesianToFractional(
                position,
                source_box);

        return FractionalToCartesian(
            fractional,
            destination_box);
    }


    void MapPPNodesBetweenBoxes(
        std::span<const PPNode> source_nodes,
        const Box& source_box,
        const Box& destination_box,
        std::span<PPNode> mapped_nodes)
    {
        if (source_nodes.size() != mapped_nodes.size()) {
            throw std::invalid_argument(
                "Source and mapped PPNode spans must have equal sizes.");
        }

        /*
         * Calculate the inverse once for the entire path. Calling
         * MapPositionBetweenBoxes() inside the loop would invert the
         * same source-box matrix for every node.
         */
        const Mat3 inverse_source_matrix =
            Inverse(source_box.matrix);

        for (std::size_t node_index = 0;
            node_index < source_nodes.size();
            ++node_index) {

            const PPNode& source_node =
                source_nodes[node_index];

            PPNode& mapped_node =
                mapped_nodes[node_index];

            // Preserve node IDs and any other stored metadata.
            mapped_node = source_node;

            const Vec3 relative_position{
                source_node.x - source_box.origin.x,
                source_node.y - source_box.origin.y,
                source_node.z - source_box.origin.z
            };

            const Vec3 fractional =
                Multiply(
                    inverse_source_matrix,
                    relative_position);

            const Vec3 mapped_relative_position =
                Multiply(
                    destination_box.matrix,
                    fractional);

            mapped_node.x =
                destination_box.origin.x + mapped_relative_position.x;

            mapped_node.y =
                destination_box.origin.y + mapped_relative_position.y;

            mapped_node.z =
                destination_box.origin.z + mapped_relative_position.z;
        }
    }

}  // namespace geometry
