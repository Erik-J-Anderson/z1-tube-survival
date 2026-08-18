#include "Tube_Survival.hpp"
#include "Geometry_Utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>


double ComputeDistanceToPolyline(
    const Segment& segment,
    std::span<const PPNode> polyline
)
{
    if (polyline.size() < 2) {
        throw std::invalid_argument(
            "ComputeDistanceToPolyline: polyline requires at least two nodes."
        );
    }

    double min_distance_squared =
        std::numeric_limits<double>::max();

    for (std::size_t i = 0; i + 1 < polyline.size(); ++i)
    {
        const PPNode& a = polyline[i];
        const PPNode& b = polyline[i + 1];

        const double ab_x = b.x - a.x;
        const double ab_y = b.y - a.y;
        const double ab_z = b.z - a.z;

        const double ap_x = segment.x - a.x;
        const double ap_y = segment.y - a.y;
        const double ap_z = segment.z - a.z;

        const double ab_squared =
            ab_x * ab_x +
            ab_y * ab_y +
            ab_z * ab_z;

        double t = 0.0;

        if (ab_squared > 0.0)
        {
            t =
                (
                    ap_x * ab_x +
                    ap_y * ab_y +
                    ap_z * ab_z
                ) /
                ab_squared;

            t = std::clamp(t, 0.0, 1.0);
        }

        const double closest_x = a.x + t * ab_x;
        const double closest_y = a.y + t * ab_y;
        const double closest_z = a.z + t * ab_z;

        const double dx = segment.x - closest_x;
        const double dy = segment.y - closest_y;
        const double dz = segment.z - closest_z;

        const double distance_squared =
            dx * dx +
            dy * dy +
            dz * dz;

        if (distance_squared < min_distance_squared) {
            min_distance_squared = distance_squared;
        }
    }

    return std::sqrt(min_distance_squared);
}


SegmentSurvivalFunction ComputeSegmentSurvivalFunction(
    const ChainTrajectory& chain_trajectory,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters
)
{
    return ComputeSegmentSurvivalFunction(
        chain_trajectory,
        std::span<const Box>{},
        lag_frames,
        tube_diameters,
        SegmentSurvivalOptions{}
    );
}


SegmentSurvivalFunction ComputeSegmentSurvivalFunction(
    const ChainTrajectory& chain_trajectory,
    std::span<const Box> frame_boxes,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters,
    const SegmentSurvivalOptions& options
)
{
    const std::size_t num_frames =
        chain_trajectory.timesteps.size();

    const std::size_t num_lags =
        lag_frames.size();

    const std::size_t num_diameters =
        tube_diameters.size();

    if (num_frames == 0) {
        throw std::invalid_argument(
            "ComputeSegmentSurvivalFunction: trajectory contains no timesteps."
        );
    }

    if (chain_trajectory.frame_offsets.size() != num_frames + 1) {
        throw std::invalid_argument(
            "ComputeSegmentSurvivalFunction: invalid frame_offsets array."
        );
    }

    if (lag_frames.empty()) {
        throw std::invalid_argument(
            "ComputeSegmentSurvivalFunction: no lag frames supplied."
        );
    }

    if (tube_diameters.empty()) {
        throw std::invalid_argument(
            "ComputeSegmentSurvivalFunction: no tube diameters supplied."
        );
    }

    for (double diameter : tube_diameters) {
        if (diameter <= 0.0) {
            throw std::invalid_argument(
                "ComputeSegmentSurvivalFunction: tube diameters must be positive."
            );
        }
    }

    if (options.apply_affine_correction)
    {
        if (frame_boxes.size() != num_frames) {
            throw std::invalid_argument(
                "ComputeSegmentSurvivalFunction: affine correction requires "
                "exactly one box per trajectory frame."
            );
        }

        /*
         * Legacy three-value Z1+ records (Lx Ly Lz) have
         * affine_geometry_available == false because they cannot encode
         * arbitrary shear tilt or a changing box origin.  They are still
         * sufficient for a trajectory that is KNOWN to remain orthorhombic,
         * such as the current z-only SLLOD/deform runs.  Therefore this
         * routine does not reject them merely because the flag is false.
         *
         * For a genuinely sheared trajectory, use box records that retain
         * the required tilt/origin information.
         */
    }

    SegmentSurvivalFunction result;

    result.options = options;
    result.tube_diameters = tube_diameters;
    result.lag_times.resize(num_lags);
    result.sample_counts.assign(num_lags, 0);

    result.survival.assign(
        num_diameters * num_lags * NUM_SAMPLE_POINTS,
        0.0
    );

    for (std::size_t lag_index = 0;
         lag_index < num_lags;
         ++lag_index)
    {
        const std::size_t lag = lag_frames[lag_index];

        if (lag >= num_frames) {
            throw std::invalid_argument(
                "ComputeSegmentSurvivalFunction: lag exceeds trajectory length."
            );
        }

        result.lag_times[lag_index] =
            static_cast<double>(
                chain_trajectory.timesteps[lag] -
                chain_trajectory.timesteps[0]
            );
    }

    // Reused for every frame pair. In threaded code, each thread must own
    // its own buffer.
    std::vector<PPNode> mapped_future_nodes;

    for (std::size_t reference_frame = 0;
         reference_frame < num_frames;
         ++reference_frame)
    {
        const SampledPath sampled_path =
            SamplePrimitivePath(
                chain_trajectory,
                reference_frame
            );

        for (std::size_t lag_index = 0;
             lag_index < num_lags;
             ++lag_index)
        {
            const std::size_t lag =
                lag_frames[lag_index];

            const std::size_t future_frame =
                reference_frame + lag;

            if (future_frame >= num_frames) {
                continue;
            }

            ++result.sample_counts[lag_index];

            const std::size_t begin =
                chain_trajectory.frame_offsets[future_frame];

            const std::size_t end =
                chain_trajectory.frame_offsets[future_frame + 1];

            if (end <= begin + 1 || end > chain_trajectory.nodes.size()) {
                throw std::runtime_error(
                    "ComputeSegmentSurvivalFunction: invalid primitive-path frame."
                );
            }

            const std::span<const PPNode> future_polyline(
                chain_trajectory.nodes.data() + begin,
                end - begin
            );

            std::span<const PPNode> comparison_polyline =
                future_polyline;

            if (options.apply_affine_correction)
            {
                mapped_future_nodes.resize(future_polyline.size());

                geometry::MapPPNodesBetweenBoxes(
                    future_polyline,
                    frame_boxes[future_frame],
                    frame_boxes[reference_frame],
                    mapped_future_nodes
                );

                comparison_polyline =
                    std::span<const PPNode>(mapped_future_nodes);
            }

            for (std::size_t segment_index = 0;
                 segment_index < NUM_SAMPLE_POINTS;
                 ++segment_index)
            {
                const double distance =
                    ComputeDistanceToPolyline(
                        sampled_path[segment_index],
                        comparison_polyline
                    );

                for (std::size_t diameter_index = 0;
                     diameter_index < num_diameters;
                     ++diameter_index)
                {
                    const double tube_radius =
                        0.5 * tube_diameters[diameter_index];

                    if (distance <= tube_radius)
                    {
                        const std::size_t index =
                            diameter_index *
                                num_lags *
                                NUM_SAMPLE_POINTS +
                            lag_index *
                                NUM_SAMPLE_POINTS +
                            segment_index;

                        result.survival[index] += 1.0;
                    }
                }
            }
        }
    }

    for (std::size_t diameter_index = 0;
         diameter_index < num_diameters;
         ++diameter_index)
    {
        for (std::size_t lag_index = 0;
             lag_index < num_lags;
             ++lag_index)
        {
            const std::uint64_t count =
                result.sample_counts[lag_index];

            if (count == 0) {
                continue;
            }

            for (std::size_t segment_index = 0;
                 segment_index < NUM_SAMPLE_POINTS;
                 ++segment_index)
            {
                const std::size_t index =
                    diameter_index *
                        num_lags *
                        NUM_SAMPLE_POINTS +
                    lag_index *
                        NUM_SAMPLE_POINTS +
                    segment_index;

                result.survival[index] /=
                    static_cast<double>(count);
            }
        }
    }

    return result;
}


HistoryDependentSurvivalResult ComputeHistoryDependentSurvivalFunction(
    const ChainTrajectory& chain_trajectory,
    std::span<const Box> frame_boxes,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters,
    const SegmentSurvivalOptions& options,
    bool use_common_origin_cohort
)
{
    const std::size_t num_frames =
        chain_trajectory.timesteps.size();

    const std::size_t num_lags =
        lag_frames.size();

    const std::size_t num_diameters =
        tube_diameters.size();

    if (num_frames == 0) {
        throw std::invalid_argument(
            "ComputeHistoryDependentSurvivalFunction: "
            "trajectory contains no timesteps."
        );
    }

    if (chain_trajectory.frame_offsets.size() != num_frames + 1) {
        throw std::invalid_argument(
            "ComputeHistoryDependentSurvivalFunction: "
            "invalid frame_offsets array."
        );
    }

    if (lag_frames.empty()) {
        throw std::invalid_argument(
            "ComputeHistoryDependentSurvivalFunction: "
            "no lag frames supplied."
        );
    }

    if (tube_diameters.empty()) {
        throw std::invalid_argument(
            "ComputeHistoryDependentSurvivalFunction: "
            "no tube diameters supplied."
        );
    }

    for (const double diameter : tube_diameters)
    {
        if (diameter <= 0.0) {
            throw std::invalid_argument(
                "ComputeHistoryDependentSurvivalFunction: "
                "tube diameters must be positive."
            );
        }
    }

    std::size_t max_lag = 0;

    for (const std::size_t lag : lag_frames)
    {
        if (lag >= num_frames) {
            throw std::invalid_argument(
                "ComputeHistoryDependentSurvivalFunction: "
                "lag exceeds trajectory length."
            );
        }

        max_lag = std::max(max_lag, lag);
    }

    if (options.apply_affine_correction &&
        frame_boxes.size() != num_frames)
    {
        throw std::invalid_argument(
            "ComputeHistoryDependentSurvivalFunction: "
            "affine correction requires exactly one box per "
            "trajectory frame."
        );
    }

    HistoryDependentSurvivalResult result;

    auto initialize_field =
        [&](SegmentSurvivalFunction& field)
        {
            field.options = options;
            field.tube_diameters = tube_diameters;
            field.lag_times.resize(num_lags);
            field.sample_counts.assign(num_lags, 0);

            field.survival.assign(
                num_diameters *
                num_lags *
                NUM_SAMPLE_POINTS,
                0.0
            );

            for (std::size_t lag_index = 0;
                 lag_index < num_lags;
                 ++lag_index)
            {
                const std::size_t lag =
                    lag_frames[lag_index];

                field.lag_times[lag_index] =
                    static_cast<double>(
                        chain_trajectory.timesteps[lag] -
                        chain_trajectory.timesteps[0]
                    );
            }
        };

    initialize_field(result.stephanou_survival);
    initialize_field(result.permanent_escape_survival);

    // Multiple requested output entries may refer to the same lag.
    // Store all output indices associated with each scanned delta.
    std::vector<std::vector<std::size_t>>
        output_indices_by_delta(max_lag + 1);

    for (std::size_t lag_index = 0;
         lag_index < num_lags;
         ++lag_index)
    {
        output_indices_by_delta[
            lag_frames[lag_index]
        ].push_back(lag_index);
    }

    std::vector<PPNode> mapped_future_nodes;

    // One byte per [diameter][contour sample].  This state is reset
    // for every time origin.
    std::vector<unsigned char> permanently_escaped(
        num_diameters * NUM_SAMPLE_POINTS,
        0
    );

    // Distance from every reference contour sample to the current
    // comparison primitive path.  Geometry is independent of diameter,
    // so calculate it only once per scanned frame.
    std::array<double, NUM_SAMPLE_POINTS> distances{};

    std::size_t reference_frame_end = num_frames;

    if (use_common_origin_cohort) {
        reference_frame_end =
            num_frames - max_lag;
    }

    for (std::size_t reference_frame = 0;
         reference_frame < reference_frame_end;
         ++reference_frame)
    {
        const SampledPath sampled_path =
            SamplePrimitivePath(
                chain_trajectory,
                reference_frame
            );

        std::fill(
            permanently_escaped.begin(),
            permanently_escaped.end(),
            static_cast<unsigned char>(0)
        );

        const std::size_t available_lag =
            num_frames - 1 - reference_frame;

        const std::size_t scan_lag =
            std::min(max_lag, available_lag);

        // IMPORTANT:
        // Scan every saved intermediate frame.  A permanent escape that
        // occurs between requested output lags must not be missed.
        for (std::size_t delta = 0;
             delta <= scan_lag;
             ++delta)
        {
            const std::size_t future_frame =
                reference_frame + delta;

            const std::size_t begin =
                chain_trajectory.frame_offsets[
                    future_frame
                ];

            const std::size_t end =
                chain_trajectory.frame_offsets[
                    future_frame + 1
                ];

            if (end <= begin + 1 ||
                end > chain_trajectory.nodes.size())
            {
                throw std::runtime_error(
                    "ComputeHistoryDependentSurvivalFunction: "
                    "invalid primitive-path frame."
                );
            }

            const std::span<const PPNode> future_polyline(
                chain_trajectory.nodes.data() + begin,
                end - begin
            );

            std::span<const PPNode> comparison_polyline =
                future_polyline;

            if (options.apply_affine_correction)
            {
                mapped_future_nodes.resize(
                    future_polyline.size()
                );

                geometry::MapPPNodesBetweenBoxes(
                    future_polyline,
                    frame_boxes[future_frame],
                    frame_boxes[reference_frame],
                    mapped_future_nodes
                );

                comparison_polyline =
                    std::span<const PPNode>(
                        mapped_future_nodes
                    );
            }

            for (std::size_t segment_index = 0;
                 segment_index < NUM_SAMPLE_POINTS;
                 ++segment_index)
            {
                distances[segment_index] =
                    ComputeDistanceToPolyline(
                        sampled_path[segment_index],
                        comparison_polyline
                    );
            }

            // Update the absorbing outer-boundary state BEFORE
            // accumulating this delta.  Therefore distance == a is
            // permanently escaped at the current frame.
            for (std::size_t diameter_index = 0;
                 diameter_index < num_diameters;
                 ++diameter_index)
            {
                const double escape_distance =
                    tube_diameters[diameter_index];

                for (std::size_t segment_index = 0;
                     segment_index < NUM_SAMPLE_POINTS;
                     ++segment_index)
                {
                    const std::size_t state_index =
                        diameter_index *
                            NUM_SAMPLE_POINTS +
                        segment_index;

                    if (distances[segment_index] >=
                        escape_distance)
                    {
                        permanently_escaped[
                            state_index
                        ] = 1;
                    }
                }
            }

            const auto& output_indices =
                output_indices_by_delta[delta];

            if (output_indices.empty()) {
                continue;
            }

            for (const std::size_t lag_index :
                 output_indices)
            {
                ++result.stephanou_survival
                    .sample_counts[lag_index];

                ++result.permanent_escape_survival
                    .sample_counts[lag_index];

                for (std::size_t diameter_index = 0;
                     diameter_index < num_diameters;
                     ++diameter_index)
                {
                    const double tube_radius =
                        0.5 *
                        tube_diameters[diameter_index];

                    for (std::size_t segment_index = 0;
                         segment_index < NUM_SAMPLE_POINTS;
                         ++segment_index)
                    {
                        const std::size_t state_index =
                            diameter_index *
                                NUM_SAMPLE_POINTS +
                            segment_index;

                        const bool alive =
                            permanently_escaped[
                                state_index
                            ] == 0;

                        const std::size_t output_index =
                            diameter_index *
                                num_lags *
                                NUM_SAMPLE_POINTS +
                            lag_index *
                                NUM_SAMPLE_POINTS +
                            segment_index;

                        // Pure absorbing outer-boundary survival.
                        if (alive) {
                            result
                                .permanent_escape_survival
                                .survival[output_index]
                                += 1.0;
                        }

                        // Stephanou-like history-dependent tube
                        // occupancy.  Temporary excursions into
                        // a/2 < d < a may recover, but an outer
                        // escape can never recover.
                        if (alive &&
                            distances[segment_index] <=
                                tube_radius)
                        {
                            result
                                .stephanou_survival
                                .survival[output_index]
                                += 1.0;
                        }
                    }
                }
            }
        }
    }

    auto normalize_field =
        [&](SegmentSurvivalFunction& field)
        {
            for (std::size_t lag_index = 0;
                 lag_index < num_lags;
                 ++lag_index)
            {
                const std::uint64_t count =
                    field.sample_counts[lag_index];

                if (count == 0) {
                    continue;
                }

                const double denominator =
                    static_cast<double>(count);

                for (std::size_t diameter_index = 0;
                     diameter_index < num_diameters;
                     ++diameter_index)
                {
                    for (std::size_t segment_index = 0;
                         segment_index <
                            NUM_SAMPLE_POINTS;
                         ++segment_index)
                    {
                        const std::size_t index =
                            diameter_index *
                                num_lags *
                                NUM_SAMPLE_POINTS +
                            lag_index *
                                NUM_SAMPLE_POINTS +
                            segment_index;

                        field.survival[index] /=
                            denominator;
                    }
                }
            }
        };

    normalize_field(result.stephanou_survival);
    normalize_field(result.permanent_escape_survival);

    return result;
}


TubeSurvivalFunction ComputeTubeSurvivalFunction(
    const SegmentSurvivalFunction& segment_survival
)
{
    const std::size_t num_lags =
        segment_survival.lag_times.size();

    const std::size_t num_diameters =
        segment_survival.tube_diameters.size();

    const std::size_t expected_size =
        num_diameters *
        num_lags *
        NUM_SAMPLE_POINTS;

    if (segment_survival.survival.size() != expected_size) {
        throw std::invalid_argument(
            "ComputeTubeSurvivalFunction: segment-survival array has incorrect dimensions."
        );
    }

    if (segment_survival.sample_counts.size() != num_lags) {
        throw std::invalid_argument(
            "ComputeTubeSurvivalFunction: sample_counts has incorrect size."
        );
    }

    TubeSurvivalFunction result;

    result.lag_times = segment_survival.lag_times;
    result.tube_diameters = segment_survival.tube_diameters;
    result.sample_counts = segment_survival.sample_counts;

    result.survival.assign(
        num_diameters * num_lags,
        0.0
    );

    for (std::size_t diameter_index = 0;
         diameter_index < num_diameters;
         ++diameter_index)
    {
        for (std::size_t lag_index = 0;
             lag_index < num_lags;
             ++lag_index)
        {
            double contour_sum = 0.0;

            for (std::size_t segment_index = 0;
                 segment_index < NUM_SAMPLE_POINTS;
                 ++segment_index)
            {
                const std::size_t flat_index =
                    diameter_index *
                        num_lags *
                        NUM_SAMPLE_POINTS +
                    lag_index *
                        NUM_SAMPLE_POINTS +
                    segment_index;

                double weight = 1.0;

                if (segment_index == 0 ||
                    segment_index == NUM_SAMPLE_POINTS - 1)
                {
                    weight = 0.5;
                }

                contour_sum +=
                    weight *
                    segment_survival.survival[flat_index];
            }

            const double tube_survival =
                contour_sum /
                static_cast<double>(NUM_SAMPLE_POINTS - 1);

            const std::size_t output_index =
                diameter_index * num_lags +
                lag_index;

            result.survival[output_index] =
                tube_survival;
        }
    }

    return result;
}

namespace
{

    struct PathProjection
    {
        double distance;
        double s_fraction;
    };


    // ------------------------------------------------------------
    // Return a non-owning view of one primitive-path frame.
    // ------------------------------------------------------------
    std::span<const PPNode> GetFrameNodes(
        const ChainTrajectory& trajectory,
        std::size_t frame_index
    )
    {
        const std::size_t num_frames =
            trajectory.frame_offsets.size() - 1;

        if (frame_index >= num_frames)
        {
            throw std::out_of_range(
                "GetFrameNodes: frame_index out of range."
            );
        }

        const std::size_t begin =
            trajectory.frame_offsets[frame_index];

        const std::size_t end =
            trajectory.frame_offsets[frame_index + 1];

        if (end <= begin)
        {
            throw std::runtime_error(
                "GetFrameNodes: primitive-path frame is empty."
            );
        }

        return std::span<const PPNode>(
            trajectory.nodes.data() + begin,
            end - begin
        );
    }


    // ------------------------------------------------------------
    // Project a point onto a primitive-path polyline.
    //
    // Returns:
    //
    //     distance
    //         Minimum Euclidean distance to the reference path.
    //
    //     s_fraction
    //         Contour position of the closest point,
    //         normalized from 0 to 1.
    //
    // This uses primitive-path contour length, NOT node index.
    // ------------------------------------------------------------
    PathProjection ProjectPointOntoPrimitivePath(
        const Segment& point,
        std::span<const PPNode> polyline
    )
    {
        if (polyline.size() < 2)
        {
            throw std::runtime_error(
                "ProjectPointOntoPrimitivePath: "
                "polyline must contain at least two nodes."
            );
        }


        // --------------------------------------------------------
        // First calculate total primitive-path contour length.
        // --------------------------------------------------------

        double total_length = 0.0;

        for (std::size_t i = 0; i + 1 < polyline.size(); ++i)
        {
            const double dx =
                polyline[i + 1].x - polyline[i].x;

            const double dy =
                polyline[i + 1].y - polyline[i].y;

            const double dz =
                polyline[i + 1].z - polyline[i].z;

            total_length += std::sqrt(
                dx * dx +
                dy * dy +
                dz * dz
            );
        }

        if (total_length <= 0.0)
        {
            throw std::runtime_error(
                "ProjectPointOntoPrimitivePath: "
                "reference path has zero contour length."
            );
        }


        // --------------------------------------------------------
        // Find the closest projection onto any finite line segment.
        // --------------------------------------------------------

        double best_distance_squared =
            std::numeric_limits<double>::infinity();

        double best_contour_position = 0.0;

        double contour_length_before_segment = 0.0;


        for (std::size_t i = 0; i + 1 < polyline.size(); ++i)
        {
            const PPNode& a = polyline[i];
            const PPNode& b = polyline[i + 1];


            const double abx = b.x - a.x;
            const double aby = b.y - a.y;
            const double abz = b.z - a.z;

            const double apx = point.x - a.x;
            const double apy = point.y - a.y;
            const double apz = point.z - a.z;


            const double segment_length_squared =
                abx * abx +
                aby * aby +
                abz * abz;

            const double segment_length =
                std::sqrt(segment_length_squared);


            if (segment_length_squared <= 0.0)
            {
                continue;
            }


            // Projection parameter onto the infinite line.
            double u =
                (
                    apx * abx +
                    apy * aby +
                    apz * abz
                    )
                / segment_length_squared;


            // Clamp to the finite line segment.
            u = std::clamp(u, 0.0, 1.0);


            const double closest_x =
                a.x + u * abx;

            const double closest_y =
                a.y + u * aby;

            const double closest_z =
                a.z + u * abz;


            const double dx =
                point.x - closest_x;

            const double dy =
                point.y - closest_y;

            const double dz =
                point.z - closest_z;


            const double distance_squared =
                dx * dx +
                dy * dy +
                dz * dz;


            if (distance_squared < best_distance_squared)
            {
                best_distance_squared = distance_squared;

                best_contour_position =
                    contour_length_before_segment
                    +
                    u * segment_length;
            }


            contour_length_before_segment += segment_length;
        }


        PathProjection result;

        result.distance =
            std::sqrt(best_distance_squared);

        result.s_fraction =
            best_contour_position / total_length;

        // Protect against tiny floating-point excursions.
        result.s_fraction =
            std::clamp(result.s_fraction, 0.0, 1.0);

        return result;
    }


    // ------------------------------------------------------------
    // Flat-array indices.
    // ------------------------------------------------------------

    std::size_t RetractionFieldIndex(
        std::size_t diameter_index,
        std::size_t lag_index,
        std::size_t s_index,
        std::size_t num_lags
    )
    {
        return
            (
                diameter_index * num_lags
                +
                lag_index
                )
            * NUM_SAMPLE_POINTS
            +
            s_index;
    }


    std::size_t RetractionSummaryIndex(
        std::size_t diameter_index,
        std::size_t lag_index,
        std::size_t num_lags
    )
    {
        return
            diameter_index * num_lags
            +
            lag_index;
    }

} // anonymous namespace



EndRetractionFunction ComputeEndRetractionFunction(
    const ChainTrajectory& chain_trajectory,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters
)
{
    if (lag_frames.empty())
    {
        throw std::invalid_argument(
            "ComputeEndRetractionFunction: "
            "lag_frames cannot be empty."
        );
    }

    if (tube_diameters.empty())
    {
        throw std::invalid_argument(
            "ComputeEndRetractionFunction: "
            "tube_diameters cannot be empty."
        );
    }

    if (chain_trajectory.frame_offsets.size() < 2)
    {
        throw std::invalid_argument(
            "ComputeEndRetractionFunction: "
            "trajectory contains fewer than one frame."
        );
    }


    const std::size_t num_frames =
        chain_trajectory.frame_offsets.size() - 1;

    const std::size_t num_lags =
        lag_frames.size();

    const std::size_t num_diameters =
        tube_diameters.size();


    // --------------------------------------------------------
    // Validate diameters.
    // --------------------------------------------------------

    for (double diameter : tube_diameters)
    {
        if (diameter <= 0.0)
        {
            throw std::invalid_argument(
                "ComputeEndRetractionFunction: "
                "tube diameters must be positive."
            );
        }
    }


    EndRetractionFunction result;

    result.lag_frames = lag_frames;
    result.tube_diameters = tube_diameters;


    // --------------------------------------------------------
    // Build the normalized contour grid.
    //
    // s = 0.00, 0.01, ..., 1.00
    // --------------------------------------------------------

    result.s_fraction.resize(NUM_SAMPLE_POINTS);

    for (std::size_t s = 0;
        s < NUM_SAMPLE_POINTS;
        ++s)
    {
        result.s_fraction[s] =
            static_cast<double>(s)
            /
            static_cast<double>(
                NUM_SAMPLE_POINTS - 1
                );
    }


    const std::size_t field_size =
        num_diameters
        *
        num_lags
        *
        NUM_SAMPLE_POINTS;

    const std::size_t summary_size =
        num_diameters
        *
        num_lags;


    result.end0_reached.assign(
        field_size,
        0.0
    );

    result.end1_reached.assign(
        field_size,
        0.0
    );

    result.combined_reached.assign(
        field_size,
        0.0
    );


    result.mean_max_depth_end0.assign(
        summary_size,
        0.0
    );

    result.mean_max_depth_end1.assign(
        summary_size,
        0.0
    );

    result.mean_max_depth_combined.assign(
        summary_size,
        0.0
    );


    result.sample_counts.assign(
        summary_size,
        0
    );


    // --------------------------------------------------------
    // We want to walk forward through the trajectory only once
    // per origin.
    //
    // Sort the requested lag values while remembering where
    // they belonged in the user's original lag array.
    // --------------------------------------------------------

    std::vector<std::pair<std::size_t, std::size_t>>
        ordered_lags;

    ordered_lags.reserve(num_lags);


    for (std::size_t lag_index = 0;
        lag_index < num_lags;
        ++lag_index)
    {
        ordered_lags.emplace_back(
            lag_frames[lag_index],
            lag_index
        );
    }


    std::sort(
        ordered_lags.begin(),
        ordered_lags.end(),
        [](const auto& a, const auto& b)
        {
            return a.first < b.first;
        }
    );


    // --------------------------------------------------------
    // Main time-origin loop.
    // --------------------------------------------------------

    for (std::size_t origin = 0;
        origin < num_frames;
        ++origin)
    {
        const std::span<const PPNode> reference_path =
            GetFrameNodes(
                chain_trajectory,
                origin
            );


        if (reference_path.size() < 2)
        {
            continue;
        }


        // One running maximum for each tube diameter.
        std::vector<double> max_depth_end0(
            num_diameters,
            0.0
        );

        std::vector<double> max_depth_end1(
            num_diameters,
            0.0
        );


        std::size_t previous_lag = 0;


        for (const auto& lag_entry : ordered_lags)
        {
            const std::size_t lag =
                lag_entry.first;

            const std::size_t lag_index =
                lag_entry.second;


            if (origin + lag >= num_frames)
            {
                // Since ordered_lags is sorted, all later
                // lags are also invalid for this origin.
                break;
            }


            // ------------------------------------------------
            // Examine every intermediate frame between the
            // previous requested lag and this requested lag.
            //
            // THIS is what makes the observable first-passage /
            // history dependent.
            // ------------------------------------------------

            for (std::size_t delta = previous_lag + 1;
                delta <= lag;
                ++delta)
            {
                const std::size_t future_frame =
                    origin + delta;


                const std::span<const PPNode> future_path =
                    GetFrameNodes(
                        chain_trajectory,
                        future_frame
                    );


                if (future_path.size() < 2)
                {
                    continue;
                }


                // Contour end 0 = first PP node.
                const Segment current_end0{
                    future_path.front().x,
                    future_path.front().y,
                    future_path.front().z
                };


                // Contour end 1 = last PP node.
                const Segment current_end1{
                    future_path.back().x,
                    future_path.back().y,
                    future_path.back().z
                };


                // Projection geometry itself does NOT depend
                // on tube diameter, so compute it only once.
                const PathProjection projection_end0 =
                    ProjectPointOntoPrimitivePath(
                        current_end0,
                        reference_path
                    );

                const PathProjection projection_end1 =
                    ProjectPointOntoPrimitivePath(
                        current_end1,
                        reference_path
                    );


                // ------------------------------------------------
                // Convert contour coordinate into penetration
                // depth measured inward from the corresponding end.
                //
                // end 0:
                //     depth = s
                //
                // end 1:
                //     depth = 1 - s
                // ------------------------------------------------

                const double depth_end0 =
                    projection_end0.s_fraction;

                const double depth_end1 =
                    1.0
                    -
                    projection_end1.s_fraction;


                for (std::size_t d = 0;
                    d < num_diameters;
                    ++d)
                {
                    const double radius =
                        0.5 * tube_diameters[d];


                    // Only count the projection as retraction
                    // if the current end is actually inside the
                    // reference tube neighborhood.
                    if (projection_end0.distance <= radius)
                    {
                        max_depth_end0[d] =
                            std::max(
                                max_depth_end0[d],
                                depth_end0
                            );
                    }


                    if (projection_end1.distance <= radius)
                    {
                        max_depth_end1[d] =
                            std::max(
                                max_depth_end1[d],
                                depth_end1
                            );
                    }
                }
            }


            previous_lag = lag;


            // ------------------------------------------------
            // Record the running maximum at this requested lag.
            // ------------------------------------------------

            for (std::size_t d = 0;
                d < num_diameters;
                ++d)
            {
                const std::size_t summary_index =
                    RetractionSummaryIndex(
                        d,
                        lag_index,
                        num_lags
                    );


                result.mean_max_depth_end0[
                    summary_index
                ] += max_depth_end0[d];

                result.mean_max_depth_end1[
                    summary_index
                ] += max_depth_end1[d];

                result.mean_max_depth_combined[
                    summary_index
                ] +=
                    0.5
                        *
                        (
                            max_depth_end0[d]
                            +
                            max_depth_end1[d]
                            );


                    ++result.sample_counts[
                        summary_index
                    ];


                    // --------------------------------------------
                    // First-passage / reached probability.
                    // --------------------------------------------

                    for (std::size_t s_index = 0;
                        s_index < NUM_SAMPLE_POINTS;
                        ++s_index)
                    {
                        const double s =
                            result.s_fraction[s_index];


                        const std::size_t field_index =
                            RetractionFieldIndex(
                                d,
                                lag_index,
                                s_index,
                                num_lags
                            );


                        if (max_depth_end0[d] >= s)
                        {
                            result.end0_reached[
                                field_index
                            ] += 1.0;
                        }


                        if (max_depth_end1[d] >= s)
                        {
                            result.end1_reached[
                                field_index
                            ] += 1.0;
                        }


                        // Pool the two ends.
                        result.combined_reached[
                            field_index
                        ] +=
                            (
                                (max_depth_end0[d] >= s)
                                ? 1.0
                                : 0.0
                                )
                                +
                                (
                                    (max_depth_end1[d] >= s)
                                    ? 1.0
                                    : 0.0
                                    );
                    }
            }
        }
    }


    // --------------------------------------------------------
    // Convert accumulated counts/sums into probabilities/means.
    // --------------------------------------------------------

    for (std::size_t d = 0;
        d < num_diameters;
        ++d)
    {
        for (std::size_t lag_index = 0;
            lag_index < num_lags;
            ++lag_index)
        {
            const std::size_t summary_index =
                RetractionSummaryIndex(
                    d,
                    lag_index,
                    num_lags
                );


            const std::uint64_t count =
                result.sample_counts[
                    summary_index
                ];


            if (count == 0)
            {
                continue;
            }


            const double count_double =
                static_cast<double>(count);


            result.mean_max_depth_end0[
                summary_index
            ] /= count_double;

            result.mean_max_depth_end1[
                summary_index
            ] /= count_double;

            result.mean_max_depth_combined[
                summary_index
            ] /= count_double;


            for (std::size_t s_index = 0;
                s_index < NUM_SAMPLE_POINTS;
                ++s_index)
            {
                const std::size_t field_index =
                    RetractionFieldIndex(
                        d,
                        lag_index,
                        s_index,
                        num_lags
                    );


                result.end0_reached[
                    field_index
                ] /= count_double;

                result.end1_reached[
                    field_index
                ] /= count_double;


                // Two ends per time origin.
                result.combined_reached[
                    field_index
                ] /=
                    (2.0 * count_double);
            }
        }
    }


    return result;
}
