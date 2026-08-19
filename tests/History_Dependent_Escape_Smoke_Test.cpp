#include "Tube_Survival.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

Box MakeUnitBox()
{
    Box box;
    box.origin = Vec3{0.0, 0.0, 0.0};
    box.matrix.value[0][0] = 20.0;
    box.matrix.value[1][1] = 20.0;
    box.matrix.value[2][2] = 20.0;
    box.affine_geometry_available = true;
    return box;
}

void RequireNear(
    double actual,
    double expected,
    const char* message)
{
    if (std::abs(actual - expected) > 1.0e-12)
    {
        throw std::runtime_error(
            std::string(message) +
            ": expected " +
            std::to_string(expected) +
            ", got " +
            std::to_string(actual)
        );
    }
}

} // namespace


int main()
{
    try
    {
        // One straight primitive path per frame.
        //
        // Tube diameter a = 10, so:
        //     inner tube radius = a/2 = 5
        //     permanent escape threshold = a = 10
        //
        // Relative shifts from the origin-frame path:
        //
        //     frame 0:  0  -> inside
        //     frame 1:  7  -> temporary excursion
        //     frame 2:  3  -> recovers
        //     frame 3: 11  -> permanent escape
        //     frame 4:  2  -> geometrically back inside, but must stay dead
        //
        // Expected Stephanou-like history:
        //     1, 0, 1, 0, 0
        //
        // Expected pure permanent-escape survival:
        //     1, 1, 1, 0, 0

        const std::vector<double> shifts{
            0.0, 7.0, 3.0, 11.0, 2.0
        };

        ChainTrajectory chain;
        chain.chain_id = 1;
        chain.timesteps = {0, 1, 2, 3, 4};
        chain.frame_offsets.push_back(0);

        for (const double shift : shifts)
        {
            chain.nodes.push_back(
                PPNode{0.0, shift, 0.0, true}
            );
            chain.nodes.push_back(
                PPNode{10.0, shift, 0.0, true}
            );

            chain.frame_offsets.push_back(
                chain.nodes.size()
            );
        }

        std::vector<Box> boxes(
            shifts.size(),
            MakeUnitBox()
        );

        const std::vector<std::size_t> lags{
            0, 1, 2, 3, 4
        };

        const std::vector<double> diameters{
            10.0
        };

        const SegmentSurvivalOptions options{
            .apply_affine_correction = true
        };

        const HistoryDependentSurvivalResult history =
            ComputeHistoryDependentSurvivalFunction(
                chain,
                boxes,
                lags,
                diameters,
                options,
                true
            );

        const std::vector<double> expected_stephanou{
            1.0, 0.0, 1.0, 0.0, 0.0
        };

        const std::vector<double> expected_permanent{
            1.0, 1.0, 1.0, 0.0, 0.0
        };

        for (std::size_t lag_index = 0;
             lag_index < lags.size();
             ++lag_index)
        {
            if (history.reference_to_future.transverse_survival
                    .sample_counts[lag_index] != 1 ||
                history.reference_to_future.permanent_escape_survival
                    .sample_counts[lag_index] != 1)
            {
                throw std::runtime_error(
                    "Common-origin smoke test did not use exactly one origin."
                );
            }

            for (std::size_t s = 0;
                 s < NUM_SAMPLE_POINTS;
                 ++s)
            {
                const std::size_t index =
                    lag_index *
                        NUM_SAMPLE_POINTS +
                    s;

                RequireNear(
                    history.reference_to_future.transverse_survival
                        .survival[index],
                    expected_stephanou[lag_index],
                    "Stephanou history mismatch"
                );

                RequireNear(
                    history.reference_to_future.permanent_escape_survival
                        .survival[index],
                    expected_permanent[lag_index],
                    "Permanent-escape history mismatch"
                );
            }
        }

        const TubeSurvivalFunction stephanou_tube =
            ComputeTubeSurvivalFunction(
                history.reference_to_future.transverse_survival
            );

        const TubeSurvivalFunction permanent_tube =
            ComputeTubeSurvivalFunction(
                history.reference_to_future.permanent_escape_survival
            );

        for (std::size_t lag_index = 0;
             lag_index < lags.size();
             ++lag_index)
        {
            RequireNear(
                stephanou_tube.survival[lag_index],
                expected_stephanou[lag_index],
                "Stephanou tube mismatch"
            );

            RequireNear(
                permanent_tube.survival[lag_index],
                expected_permanent[lag_index],
                "Permanent tube mismatch"
            );
        }

        // Pure absorbing outer-boundary survival must be nonincreasing
        // for a fixed/common cohort.
        for (std::size_t lag_index = 1;
             lag_index < lags.size();
             ++lag_index)
        {
            if (permanent_tube.survival[lag_index] >
                permanent_tube.survival[lag_index - 1] +
                    1.0e-12)
            {
                throw std::runtime_error(
                    "Permanent-escape survival increased with time."
                );
            }
        }

        std::cout
            << "History sequence         : 1 0 1 0 0\n"
            << "Permanent escape sequence: 1 1 1 0 0\n"
            << "ALL HISTORY-DEPENDENT ESCAPE CHECKS PASSED\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "FAILED: "
            << error.what()
            << '\n';

        return 1;
    }
}
