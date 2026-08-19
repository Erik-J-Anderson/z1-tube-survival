#pragma once

#include "Primitive_Path_Trajectory.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>


struct SegmentSurvivalOptions
{
    // When true, each future primitive path is pulled back into the
    // box of its own time origin before distances are evaluated.
    bool apply_affine_correction{ false };
};


struct SegmentSurvivalFunction
{
    std::vector<double> lag_times;
    std::vector<double> tube_diameters;

    // Flattened storage: survival[a, lag, s]
    // Dimensions:
    //   num_diameters x num_lags x NUM_SAMPLE_POINTS
    std::vector<double> survival;

    // Number of valid time origins contributing to each lag.
    std::vector<std::uint64_t> sample_counts;

    SegmentSurvivalOptions options;
};


struct TubeSurvivalFunction
{
    std::vector<double> lag_times;
    std::vector<double> tube_diameters;

    // Flattened storage: survival[a, lag]
    std::vector<double> survival;

    // Carried through from the segment-survival calculation.
    std::vector<std::uint64_t> sample_counts;
};


// History-dependent survival for one directed distance definition.
//
// transverse_survival:
//   History-gated inner-tube occupancy:
//       distance <= a/2
//   provided the same contour sample has never previously reached
//       distance >= a.
//
// permanent_escape_survival:
//   Pure absorbing survival of the outer transverse boundary:
//       1 while max_history(distance) < a,
//       0 forever after distance >= a.
//
// full_survival:
//   transverse_survival with the additional longitudinal/end-retraction
//   first-passage gate applied.
//
// All fields use the standard [diameter][lag][s] layout.
struct DirectionalHistorySurvival
{
    SegmentSurvivalFunction transverse_survival;
    SegmentSurvivalFunction permanent_escape_survival;
    SegmentSurvivalFunction full_survival;
};


// Compare both directed geometric definitions under the same time origins,
// affine correction, transverse history gate, and longitudinal first-passage
// history:
//
// reference_to_future:
//     d_RF(s,t) = distance(r_reference(s), P_future)
//
// future_to_reference:
//     d_FR(s,t) = distance(r_future(s), P_reference)
//
// longitudinal_survival contains only the common end-retraction gate and is
// useful for diagnosing whether longitudinal escape dominates either result.
struct HistoryDependentSurvivalResult
{
    DirectionalHistorySurvival reference_to_future;
    DirectionalHistorySurvival future_to_reference;
    SegmentSurvivalFunction longitudinal_survival;
};


double ComputeDistanceToPolyline(
    const Segment& segment,
    std::span<const PPNode> polyline
);


// Legacy/uncorrected overload.  Kept so existing equilibrium smoke tests
// and callers continue to compile unchanged.
SegmentSurvivalFunction ComputeSegmentSurvivalFunction(
    const ChainTrajectory& chain_trajectory,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters
);


// Box-aware overload.  With apply_affine_correction=true, the future
// path at t0 + lag is mapped from frame_boxes[t0 + lag] into
// frame_boxes[t0].  Therefore the affine reference box changes with
// every time origin.
SegmentSurvivalFunction ComputeSegmentSurvivalFunction(
    const ChainTrajectory& chain_trajectory,
    std::span<const Box> frame_boxes,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters,
    const SegmentSurvivalOptions& options = {}
);


// History-dependent affine-aware comparison.
//
// Every saved intermediate frame is examined, not only requested output lags.
//
// Two directed distances are evaluated on the same normalized contour grid:
//
//   reference -> future:
//       distance from the reference sample r_reference(s) to the future PP.
//
//   future -> reference:
//       distance from the future sample r_future(s) to the reference PP.
//
// For either directed distance:
//
//   distance <= a/2
//       currently inside the inner tube.
//
//   a/2 < distance < a
//       temporarily outside the inner tube; recovery is allowed.
//
//   distance >= a
//       permanently escaped transversely for that time origin.
//
// Longitudinal escape is tracked independently by projecting each future
// chain end onto the reference primitive path and retaining the running
// maximum inward penetration depth from each end.  full_survival requires
// both the directed transverse criterion and the longitudinal gate.
//
// If use_common_origin_cohort is false (default), each requested lag uses
// every time origin valid for that lag.  If true, every lag uses only time
// origins valid through the largest requested lag.


HistoryDependentSurvivalResult ComputeHistoryDependentSurvivalFunction(
    const ChainTrajectory& chain_trajectory,
    std::span<const Box> frame_boxes,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters,
    const SegmentSurvivalOptions& options = {},
    bool use_common_origin_cohort = false,
    bool compute_reference_to_future = true
);


TubeSurvivalFunction ComputeTubeSurvivalFunction(
    const SegmentSurvivalFunction& segment_survival
);




struct EndRetractionFunction
{
    // Requested lag values, in saved-frame units.
    std::vector<std::size_t> lag_frames;

    // Tube diameters used to decide whether the end is still
    // inside the neighborhood of the reference primitive path.
    std::vector<double> tube_diameters;

    // Same normalized contour grid as the segment-survival function:
    // 0.00, 0.01, ..., 1.00
    std::vector<double> s_fraction;


    // ------------------------------------------------------------
    // First-passage / reached probabilities.
    //
    // Layout:
    //     [diameter][lag][s]
    //
    // end0_reached:
    //     P(max penetration depth from contour end 0 >= s)
    //
    // end1_reached:
    //     P(max penetration depth from contour end 1 >= s)
    //
    // combined_reached:
    //     Same probability, pooling both chain ends.
    // ------------------------------------------------------------

    std::vector<double> end0_reached;
    std::vector<double> end1_reached;
    std::vector<double> combined_reached;


    // ------------------------------------------------------------
    // Mean maximum penetration depth.
    //
    // Layout:
    //     [diameter][lag]
    // ------------------------------------------------------------

    std::vector<double> mean_max_depth_end0;
    std::vector<double> mean_max_depth_end1;
    std::vector<double> mean_max_depth_combined;


    // Number of valid time origins contributing to [diameter][lag].
    //
    // Note:
    // combined_reached has twice this many end samples because
    // each origin contributes two chain ends.
    std::vector<std::uint64_t> sample_counts;
};

EndRetractionFunction ComputeEndRetractionFunction(
    const ChainTrajectory& chain_trajectory,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters
);
