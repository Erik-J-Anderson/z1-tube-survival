#include "Primitive_Path_Trajectory.hpp"
#include "Tube_Survival.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>


namespace
{

Box MakeOrthorhombicBox(double lz)
{
    Box box{};

    box.origin = Vec3{0.0, 0.0, 0.0};
    box.matrix.value[0][0] = 10.0;
    box.matrix.value[1][1] = 10.0;
    box.matrix.value[2][2] = lz;

    // Deliberately false: this reproduces the legacy three-value
    // Z1+ box record used by the current z-only deforming trajectory.
    box.affine_geometry_available = false;

    return box;
}


std::size_t TubeIndex(
    std::size_t diameter_index,
    std::size_t lag_index,
    std::size_t num_lags)
{
    return diameter_index * num_lags + lag_index;
}

} // namespace


int main()
{
    try
    {
        // ------------------------------------------------------------
        // Construct a three-frame primitive path with NO non-affine
        // motion.  The path has fixed fractional z = 0.25, while Lz
        // changes strongly between frames.
        //
        //   frame 0: Lz = 10 -> z = 2.50
        //   frame 1: Lz = 20 -> z = 5.00
        //   frame 2: Lz =  5 -> z = 1.25
        //
        // For lag 1 there are TWO different time origins:
        //
        //   origin 0 -> frame 1
        //   origin 1 -> frame 2
        //
        // Correct affine handling must use the box at EACH origin.
        // Mapping every frame to frame 0 would fail the second pair.
        // ------------------------------------------------------------

        ChainTrajectory chain;
        chain.chain_id = 1;
        chain.timesteps = {0, 100, 200};
        chain.frame_offsets = {0, 2, 4, 6};

        chain.nodes = {
            PPNode{2.0, 0.0, 2.50, true},
            PPNode{8.0, 0.0, 2.50, true},

            PPNode{2.0, 0.0, 5.00, true},
            PPNode{8.0, 0.0, 5.00, true},

            PPNode{2.0, 0.0, 1.25, true},
            PPNode{8.0, 0.0, 1.25, true}
        };

        const std::vector<Box> boxes{
            MakeOrthorhombicBox(10.0),
            MakeOrthorhombicBox(20.0),
            MakeOrthorhombicBox(5.0)
        };

        const std::vector<std::size_t> lag_frames{0, 1};
        const std::vector<double> tube_diameters{0.2};

        const SegmentSurvivalFunction uncorrected_segments =
            ComputeSegmentSurvivalFunction(
                chain,
                boxes,
                lag_frames,
                tube_diameters,
                SegmentSurvivalOptions{
                    .apply_affine_correction = false
                }
            );

        const SegmentSurvivalFunction corrected_segments =
            ComputeSegmentSurvivalFunction(
                chain,
                boxes,
                lag_frames,
                tube_diameters,
                SegmentSurvivalOptions{
                    .apply_affine_correction = true
                }
            );

        const TubeSurvivalFunction uncorrected_tube =
            ComputeTubeSurvivalFunction(uncorrected_segments);

        const TubeSurvivalFunction corrected_tube =
            ComputeTubeSurvivalFunction(corrected_segments);

        constexpr std::size_t diameter_index = 0;
        constexpr std::size_t lag_index = 1;
        constexpr std::size_t num_lags = 2;
        constexpr double tolerance = 1.0e-12;

        const std::size_t index =
            TubeIndex(diameter_index, lag_index, num_lags);

        const double uncorrected =
            uncorrected_tube.survival[index];

        const double corrected =
            corrected_tube.survival[index];

        std::cout
            << "uncorrected lag-1 tube survival = "
            << uncorrected << '\n'
            << "corrected lag-1 tube survival   = "
            << corrected << '\n'
            << "lag-1 time origins              = "
            << corrected_tube.sample_counts[lag_index]
            << '\n';

        if (std::abs(uncorrected) > tolerance)
        {
            throw std::runtime_error(
                "Expected uncorrected lag-1 survival to be zero."
            );
        }

        if (std::abs(corrected - 1.0) > tolerance)
        {
            throw std::runtime_error(
                "Affine correction did not recover unit survival."
            );
        }

        if (corrected_tube.sample_counts[lag_index] != 2)
        {
            throw std::runtime_error(
                "Expected two distinct time origins at lag 1."
            );
        }

        std::cout
            << "ALL TIME-ORIGIN AFFINE TESTS PASSED\n";
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "TIME-ORIGIN AFFINE TEST FAILED:\n"
            << error.what()
            << '\n';

        return 1;
    }

    return 0;
}
