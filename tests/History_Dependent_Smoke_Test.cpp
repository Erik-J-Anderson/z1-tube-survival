#include "Tube_Survival.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{

    std::size_t idx(
        std::size_t diameter_index,
        std::size_t lag_index,
        std::size_t segment_index,
        std::size_t num_lags
    )
    {
        return diameter_index * num_lags * NUM_SAMPLE_POINTS
            + lag_index * NUM_SAMPLE_POINTS
            + segment_index;
    }


    void expect_near(
        double actual,
        double expected,
        const std::string& message
    )
    {
        constexpr double tolerance = 1.0e-12;

        if (std::abs(actual - expected) > tolerance)
        {
            throw std::runtime_error(
                message
                + " | expected = "
                + std::to_string(expected)
                + ", actual = "
                + std::to_string(actual)
            );
        }
    }


    /*
     * Add a straight two-node primitive path:
     *
     *      (x0, y, 0) -------- (x1, y, 0)
     *
     * Each call adds one saved frame.
     */
    void add_frame(
        ChainTrajectory& trajectory,
        double x0,
        double x1,
        double y
    )
    {
        trajectory.nodes.push_back(
            PPNode{ x0, y, 0.0, true }
        );

        trajectory.nodes.push_back(
            PPNode{ x1, y, 0.0, true }
        );

        trajectory.frame_offsets.push_back(
            trajectory.nodes.size()
        );
    }

} // namespace


int main()
{
    /*
     * Tube diameter:
     *
     *      a = 4
     *
     * Inner tube radius:
     *
     *      a/2 = 2
     *
     * Therefore:
     *
     *      d <= 2       -> inside
     *      2 < d < 4    -> temporarily outside, recovery allowed
     *      d >= 4       -> permanent transverse escape
     */

    const std::vector<double> tube_diameters{
        4.0
    };


    /*
     * Synthetic trajectory.
     *
     * Reference primitive path:
     *
     *      frame 0:
     *
     *          y = 0
     *
     *          x=0 ---------------- x=10
     *
     *
     * frame 1:
     *
     *      Entire PP translated transversely to y=3.
     *
     *      d = 3
     *
     *      This is OUTSIDE the inner radius (2),
     *      but INSIDE the absorbing boundary (4).
     *
     *      Therefore:
     *
     *          transverse occupancy = 0
     *          permanent survival   = 1
     *
     *
     * frame 2:
     *
     *      PP returns to y=1.
     *
     *      d = 1
     *
     *      It should RECOVER:
     *
     *          transverse occupancy = 1
     *          permanent survival   = 1
     *
     *
     * frame 3:
     *
     *      PP goes to y=5.
     *
     *      d = 5 >= 4
     *
     *      This crosses the absorbing boundary.
     *
     *      IMPORTANT:
     *
     *      frame 3 is deliberately NOT one of our requested
     *      output lags.
     *
     *      The history-dependent code must still see it.
     *
     *
     * frame 4:
     *
     *      PP returns to y=1.
     *
     *      Geometrically it is back inside the tube,
     *      but it must remain permanently dead because
     *      frame 3 crossed d >= 4.
     *
     *
     * frame 5:
     *
     *      PP becomes:
     *
     *          x=2 ---------------- x=12
     *
     *      This creates a longitudinal advance of the
     *      left end by s = 0.2 along the reference PP.
     *
     *      Therefore points near s=0.1 should be
     *      longitudinally dead, while s=0.5 remains alive.
     */


    ChainTrajectory trajectory;

    trajectory.chain_id = 0;

    trajectory.timesteps = {
        0,
        1,
        2,
        3,
        4,
        5
    };

    /*
     * frame_offsets always begins with zero.
     */
    trajectory.frame_offsets.push_back(0);


    // frame 0: reference
    add_frame(
        trajectory,
        0.0,
        10.0,
        0.0
    );


    // frame 1: temporary transverse excursion
    add_frame(
        trajectory,
        0.0,
        10.0,
        3.0
    );


    // frame 2: recovery into inner tube
    add_frame(
        trajectory,
        0.0,
        10.0,
        1.0
    );


    // frame 3: PERMANENT transverse escape
    //
    // Note that lag 3 will NOT be requested below.
    add_frame(
        trajectory,
        0.0,
        10.0,
        5.0
    );


    // frame 4: returns geometrically, but should remain dead
    add_frame(
        trajectory,
        0.0,
        10.0,
        1.0
    );


    // frame 5: longitudinal slide / end retraction
    add_frame(
        trajectory,
        2.0,
        12.0,
        0.0
    );


    /*
     * Deliberately skip lag 3.
     *
     * This is important:
     *
     * frame 3 contains the permanent escape event.
     *
     * If ComputeHistoryDependentSurvivalFunction only examined
     * requested output frames, it would miss the escape and
     * incorrectly report frame 4 as alive.
     */
    const std::vector<std::size_t> lag_frames{
        0,
        1,
        2,
        4,
        5
    };


    /*
     * No affine correction needed for this synthetic test.
     *
     * use_common_origin_cohort = true
     *
     * Since the maximum requested lag is 5 and there are
     * exactly 6 frames, only frame 0 can be used as a time
     * origin.
     *
     * That makes every expected value exactly 0 or 1.
     */
    const auto result =
        ComputeHistoryDependentSurvivalFunction(
            trajectory,
            std::span<const Box>{},
            lag_frames,
            tube_diameters,
            SegmentSurvivalOptions{},
            true
        );


    const std::size_t num_lags =
        lag_frames.size();


    /*
     * NUM_SAMPLE_POINTS = 101:
     *
     *      index 10 -> s = 0.10
     *      index 50 -> s = 0.50
     */
    constexpr std::size_t s10 = 10;
    constexpr std::size_t s50 = 50;


    /*
     * lag_frames:
     *
     *      lag index 0 -> frame 0
     *      lag index 1 -> frame 1
     *      lag index 2 -> frame 2
     *      lag index 3 -> frame 4
     *      lag index 4 -> frame 5
     */
    constexpr std::size_t lag0 = 0;
    constexpr std::size_t lag1 = 1;
    constexpr std::size_t lag2 = 2;
    constexpr std::size_t lag4 = 3;
    constexpr std::size_t lag5 = 4;


    /*
     * Convenience references.
     */
    const auto& rf =
        result.reference_to_future;

    const auto& fr =
        result.future_to_reference;

    const auto& longitudinal =
        result.longitudinal_survival;


    // ============================================================
    // 1. Initial state
    // ============================================================

    expect_near(
        rf.transverse_survival.survival[
            idx(0, lag0, s50, num_lags)
        ],
        1.0,
        "frame 0 RF transverse survival should be 1"
    );

    expect_near(
        fr.transverse_survival.survival[
            idx(0, lag0, s50, num_lags)
        ],
        1.0,
        "frame 0 FR transverse survival should be 1"
    );


    // ============================================================
    // 2. Temporary excursion: d = 3
    //
    //    radius = 2
    //    absorbing boundary = 4
    //
    //    Currently outside, but NOT permanently escaped.
    // ============================================================

    expect_near(
        rf.transverse_survival.survival[
            idx(0, lag1, s50, num_lags)
        ],
        0.0,
        "frame 1 RF should be temporarily outside inner tube"
    );

    expect_near(
        rf.permanent_escape_survival.survival[
            idx(0, lag1, s50, num_lags)
        ],
        1.0,
        "frame 1 RF should not yet be permanently escaped"
    );


    expect_near(
        fr.transverse_survival.survival[
            idx(0, lag1, s50, num_lags)
        ],
        0.0,
        "frame 1 FR should be temporarily outside inner tube"
    );

    expect_near(
        fr.permanent_escape_survival.survival[
            idx(0, lag1, s50, num_lags)
        ],
        1.0,
        "frame 1 FR should not yet be permanently escaped"
    );


    // ============================================================
    // 3. Recovery: d = 1
    //
    //    The temporary excursion should NOT kill the segment.
    //
    //    Expected history:
    //
    //        1 -> 0 -> 1
    // ============================================================

    expect_near(
        rf.transverse_survival.survival[
            idx(0, lag2, s50, num_lags)
        ],
        1.0,
        "frame 2 RF should recover after temporary excursion"
    );

    expect_near(
        rf.permanent_escape_survival.survival[
            idx(0, lag2, s50, num_lags)
        ],
        1.0,
        "frame 2 RF should still be permanently alive"
    );


    expect_near(
        fr.transverse_survival.survival[
            idx(0, lag2, s50, num_lags)
        ],
        1.0,
        "frame 2 FR should recover after temporary excursion"
    );

    expect_near(
        fr.permanent_escape_survival.survival[
            idx(0, lag2, s50, num_lags)
        ],
        1.0,
        "frame 2 FR should still be permanently alive"
    );


    // ============================================================
    // 4. Hidden permanent escape at frame 3
    //
    //    frame 3 itself was NOT requested.
    //
    //    At frame 4 the geometry is back at d = 1.
    //
    //    Nevertheless, the segment must remain DEAD because
    //    the history scanner should have observed d = 5 at
    //    intermediate frame 3.
    //
    //    Expected history:
    //
    //        1 -> 0 -> 1 -> [escape] -> 0
    // ============================================================

    expect_near(
        rf.transverse_survival.survival[
            idx(0, lag4, s50, num_lags)
        ],
        0.0,
        "frame 4 RF must remain dead after hidden frame-3 escape"
    );

    expect_near(
        rf.permanent_escape_survival.survival[
            idx(0, lag4, s50, num_lags)
        ],
        0.0,
        "frame 4 RF permanent survival must remain zero"
    );


    expect_near(
        fr.transverse_survival.survival[
            idx(0, lag4, s50, num_lags)
        ],
        0.0,
        "frame 4 FR must remain dead after hidden frame-3 escape"
    );

    expect_near(
        fr.permanent_escape_survival.survival[
            idx(0, lag4, s50, num_lags)
        ],
        0.0,
        "frame 4 FR permanent survival must remain zero"
    );


    /*
     * Longitudinal motion has not occurred yet at frame 4.
     *
     * Therefore the midpoint should still be longitudinally alive.
     *
     * This also confirms that full survival at this point is
     * being killed by TRANSVERSE history, not by the end gate.
     */
    expect_near(
        longitudinal.survival[
            idx(0, lag4, s50, num_lags)
        ],
        1.0,
        "frame 4 midpoint should still be longitudinally alive"
    );

    expect_near(
        rf.full_survival.survival[
            idx(0, lag4, s50, num_lags)
        ],
        0.0,
        "frame 4 RF full survival should be killed by transverse history"
    );


    // ============================================================
    // 5. Longitudinal/end-retraction gate
    //
    //    frame 5:
    //
    //        reference: [0,10]
    //        future:    [2,12]
    //
    //    Left end projects to:
    //
    //        s = 2/10 = 0.2
    //
    //    Therefore:
    //
    //        s = 0.10 -> longitudinally dead
    //        s = 0.50 -> longitudinally alive
    // ============================================================

    expect_near(
        longitudinal.survival[
            idx(0, lag5, s10, num_lags)
        ],
        0.0,
        "frame 5 s=0.10 should be killed by left-end retraction"
    );


    expect_near(
        longitudinal.survival[
            idx(0, lag5, s50, num_lags)
        ],
        1.0,
        "frame 5 s=0.50 should remain longitudinally alive"
    );


    // ============================================================
    // Final report
    // ============================================================

    std::cout
        << "History-dependent survival smoke test passed.\n"
        << "\n"
        << "Verified:\n"
        << "  [1] temporary transverse escape gives survival 0\n"
        << "      without permanent death\n"
        << "  [2] return inside the inner tube restores survival\n"
        << "  [3] crossing the outer boundary causes permanent death\n"
        << "  [4] permanent escape is detected at an unrequested\n"
        << "      intermediate frame\n"
        << "  [5] return after permanent escape does not restore survival\n"
        << "  [6] longitudinal end retraction independently kills\n"
        << "      the expected reference-tube region\n"
        << "  [7] both RF and FR transverse histories obey the\n"
        << "      same absorbing-history logic in this symmetric case\n";

    return 0;
}