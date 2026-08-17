#include "Primitive_Path_Trajectory.hpp"
#include "Tube_Survival.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>


namespace
{

Box MakeBox(double xy)
{
    return Box{
        .origin = Vec3{0.0, 0.0, 0.0},
        .matrix = Mat3{{
            {10.0, xy,   0.0},
            {0.0,  10.0, 0.0},
            {0.0,  0.0,  10.0}
        }},
        .affine_geometry_available = true
    };
}


std::size_t SurvivalIndex(
    std::size_t diameter_index,
    std::size_t lag_index,
    std::size_t segment_index,
    std::size_t num_lags)
{
    return
        (diameter_index * num_lags + lag_index)
        * NUM_SAMPLE_POINTS
        + segment_index;
}

} // namespace


int main()
{
    /*
     * Frame 0 is a vertical primitive path:
     *
     *     r(s) = (0, 10s, 0).
     *
     * Frame 1 is exactly the same path after an affine xy shear with
     * gamma = xy/Ly = 5/10 = 0.5:
     *
     *     r'(s) = (5s, 10s, 0).
     *
     * There is no non-affine motion. Therefore affine-corrected segment
     * survival must equal one at every sampled contour position.
     */
    ChainTrajectory chain;
    chain.chain_id = 1;
    chain.timesteps = {0, 100};
    chain.frame_offsets = {0, 2, 4};
    chain.nodes = {
        PPNode{0.0, 0.0, 0.0, true},
        PPNode{0.0, 10.0, 0.0, true},
        PPNode{0.0, 0.0, 0.0, true},
        PPNode{5.0, 10.0, 0.0, true}
    };

    const std::vector<Box> boxes{
        MakeBox(0.0),
        MakeBox(5.0)
    };

    const std::vector<std::size_t> lag_frames{0, 1};
    const std::vector<double> tube_diameters{0.2};

    const SegmentSurvivalFunction uncorrected =
        ComputeSegmentSurvivalFunction(
            chain,
            boxes,
            lag_frames,
            tube_diameters,
            SegmentSurvivalOptions{
                .apply_affine_correction = false
            }
        );

    const SegmentSurvivalFunction corrected =
        ComputeSegmentSurvivalFunction(
            chain,
            boxes,
            lag_frames,
            tube_diameters,
            SegmentSurvivalOptions{
                .apply_affine_correction = true
            }
        );

    constexpr std::size_t diameter_index = 0;
    constexpr std::size_t lag_index = 1;
    constexpr std::size_t num_lags = 2;
    constexpr double tolerance = 1.0e-12;

    double minimum_corrected_survival = 1.0;
    double mean_uncorrected_survival = 0.0;
    double mean_corrected_survival = 0.0;

    std::ofstream output("affine_smoke_test.csv");

    if (!output) {
        throw std::runtime_error(
            "Could not open affine_smoke_test.csv for writing."
        );
    }

    output << "s_fraction,uncorrected,corrected\n";

    for (std::size_t segment_index = 0;
         segment_index < NUM_SAMPLE_POINTS;
         ++segment_index)
    {
        const std::size_t index =
            SurvivalIndex(
                diameter_index,
                lag_index,
                segment_index,
                num_lags
            );

        const double uncorrected_value =
            uncorrected.survival[index];

        const double corrected_value =
            corrected.survival[index];

        minimum_corrected_survival =
            std::min(
                minimum_corrected_survival,
                corrected_value
            );

        mean_uncorrected_survival += uncorrected_value;
        mean_corrected_survival += corrected_value;

        const double s_fraction =
            static_cast<double>(segment_index)
            / static_cast<double>(NUM_SAMPLE_POINTS - 1);

        output
            << s_fraction << ','
            << uncorrected_value << ','
            << corrected_value << '\n';
    }

    mean_uncorrected_survival /=
        static_cast<double>(NUM_SAMPLE_POINTS);

    mean_corrected_survival /=
        static_cast<double>(NUM_SAMPLE_POINTS);

    if (std::abs(minimum_corrected_survival - 1.0) > tolerance)
    {
        std::cerr
            << "FAIL: corrected survival was not one everywhere.\n"
            << "Minimum corrected survival: "
            << minimum_corrected_survival << '\n';

        return 1;
    }

    if (mean_uncorrected_survival >= 0.5)
    {
        std::cerr
            << "FAIL: uncorrected result did not detect the imposed shear.\n"
            << "Mean uncorrected survival: "
            << mean_uncorrected_survival << '\n';

        return 1;
    }

    std::cout
        << "PASS: pure affine shear was removed correctly.\n"
        << "Mean uncorrected survival: "
        << mean_uncorrected_survival << '\n'
        << "Mean corrected survival:   "
        << mean_corrected_survival << '\n'
        << "Wrote affine_smoke_test.csv\n";

    return 0;
}
