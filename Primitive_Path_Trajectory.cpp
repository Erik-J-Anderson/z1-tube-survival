#include "Primitive_Path_Trajectory.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>


SampledPath SamplePrimitivePath(
    const ChainTrajectory& trajectory,
    std::size_t frame_index
)
{
    if (trajectory.frame_offsets.empty()) {
        throw std::invalid_argument(
            "SamplePrimitivePath: trajectory contains no frame offsets."
        );
    }

    const std::size_t num_frames =
        trajectory.frame_offsets.size() - 1;

    if (frame_index >= num_frames) {
        throw std::out_of_range(
            "SamplePrimitivePath: frame index is out of range."
        );
    }

    const std::size_t begin =
        trajectory.frame_offsets[frame_index];

    const std::size_t end =
        trajectory.frame_offsets[frame_index + 1];

    if (end > trajectory.nodes.size() || begin > end) {
        throw std::runtime_error(
            "SamplePrimitivePath: invalid frame offsets."
        );
    }

    const std::size_t num_nodes = end - begin;

    if (num_nodes < 2) {
        throw std::runtime_error(
            "SamplePrimitivePath: primitive path requires at least two nodes."
        );
    }

    const std::size_t num_edges = num_nodes - 1;

    std::vector<double> edge_lengths(num_edges, 0.0);

    double total_length = 0.0;

    for (std::size_t edge = 0; edge < num_edges; ++edge)
    {
        const PPNode& a = trajectory.nodes[begin + edge];
        const PPNode& b = trajectory.nodes[begin + edge + 1];

        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double dz = b.z - a.z;

        const double length =
            std::sqrt(dx * dx + dy * dy + dz * dz);

        edge_lengths[edge] = length;
        total_length += length;
    }

    if (total_length <= 0.0) {
        throw std::runtime_error(
            "SamplePrimitivePath: primitive path has zero contour length."
        );
    }

    SampledPath sampled_path{};

    // Endpoints are exact and avoid roundoff at s/L = 0 and 1.
    const PPNode& first_node = trajectory.nodes[begin];
    const PPNode& last_node = trajectory.nodes[end - 1];

    sampled_path.front() = Segment{
        first_node.x,
        first_node.y,
        first_node.z
    };

    sampled_path.back() = Segment{
        last_node.x,
        last_node.y,
        last_node.z
    };

    std::size_t edge = 0;
    double edge_start_distance = 0.0;

    for (std::size_t sample = 1;
         sample + 1 < NUM_SAMPLE_POINTS;
         ++sample)
    {
        const double fraction =
            static_cast<double>(sample) /
            static_cast<double>(NUM_SAMPLE_POINTS - 1);

        const double target_distance =
            fraction * total_length;

        while (edge + 1 < num_edges &&
               edge_start_distance + edge_lengths[edge] < target_distance)
        {
            edge_start_distance += edge_lengths[edge];
            ++edge;
        }

        // Skip any zero-length edges at the current location.
        while (edge + 1 < num_edges && edge_lengths[edge] <= 0.0)
        {
            edge_start_distance += edge_lengths[edge];
            ++edge;
        }

        const PPNode& a = trajectory.nodes[begin + edge];
        const PPNode& b = trajectory.nodes[begin + edge + 1];

        const double length = edge_lengths[edge];

        double t = 0.0;

        if (length > 0.0) {
            t =
                (target_distance - edge_start_distance) /
                length;
        }

        if (t < 0.0) {
            t = 0.0;
        }
        else if (t > 1.0) {
            t = 1.0;
        }

        sampled_path[sample] = Segment{
            a.x + t * (b.x - a.x),
            a.y + t * (b.y - a.y),
            a.z + t * (b.z - a.z)
        };
    }

    return sampled_path;
}
