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
    bool apply_affine_correction{false};
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


// History-dependent tube survival.
//
// stephanou_survival:
//   A contour sample contributes 1 only when
//       distance <= a/2
//   AND it has never previously reached
//       distance >= a.
//   Motion through a/2 < distance < a is reversible: the sample may
//   return to the inner tube until the outer escape boundary is reached.
//
// permanent_escape_survival:
//   Pure absorbing survival of the outer escape boundary:
//       1 while max_history(distance) < a,
//       0 forever after distance >= a.
//
// Both fields use the same [diameter][lag][s] layout as
// SegmentSurvivalFunction.
struct HistoryDependentSurvivalResult
{
    SegmentSurvivalFunction stephanou_survival;
    SegmentSurvivalFunction permanent_escape_survival;
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


// History-dependent affine-aware survival.
//
// Every saved intermediate frame is examined, not only the requested
// output lags.  For each time origin, contour sample, and tube diameter:
//
//   distance <= a/2
//       inside the inner tube.
//
//   a/2 < distance < a
//       temporarily outside the inner tube; recovery is allowed.
//
//   distance >= a
//       permanently escaped for that time origin.  It can never
//       contribute to stephanou_survival again.
//
// If use_common_origin_cohort is false (the default), each requested lag
// uses every time origin valid for that lag, matching the ordinary
// instantaneous survival estimator.
//
// If true, only origins valid through the largest requested lag are used
// for every lag.  This is useful for strict absorbing-survival regression
// tests and first-passage analysis because the cohort is identical at
// all lags.
HistoryDependentSurvivalResult ComputeHistoryDependentSurvivalFunction(
    const ChainTrajectory& chain_trajectory,
    std::span<const Box> frame_boxes,
    const std::vector<std::size_t>& lag_frames,
    const std::vector<double>& tube_diameters,
    const SegmentSurvivalOptions& options = {},
    bool use_common_origin_cohort = false
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
