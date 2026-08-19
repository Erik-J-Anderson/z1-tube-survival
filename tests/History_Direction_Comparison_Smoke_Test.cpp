#include "Tube_Survival.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

static std::size_t idx(std::size_t d, std::size_t lag, std::size_t s,
    std::size_t nlags)
{
    return d * nlags * NUM_SAMPLE_POINTS + lag * NUM_SAMPLE_POINTS + s;
}

int main()
{
    ChainTrajectory tr;
    tr.timesteps = { 0, 1 };
    tr.frame_offsets = { 0, 2, 4 };

    // Reference PP: [0,10].
    tr.nodes.push_back(PPNode{ 0,0,0,true });
    tr.nodes.push_back(PPNode{ 10,0,0,true });

    // Future PP: [4,14], a pure longitudinal slide.
    tr.nodes.push_back(PPNode{ 4,0,0,true });
    tr.nodes.push_back(PPNode{ 14,0,0,true });

    const std::vector<std::size_t> lags{ 0,1 };
    const std::vector<double> diameters{ 4.0 }; // radius=2, outer escape=4

    auto r = ComputeHistoryDependentSurvivalFunction(
        tr,
        std::span<const Box>{},
        lags,
        diameters,
        SegmentSurvivalOptions{},
        true
    );

    const std::size_t s90 = 90;
    const std::size_t s10 = 10;
    const std::size_t i90 = idx(0, 1, s90, lags.size());
    const std::size_t i10 = idx(0, 1, s10, lags.size());

    const double rf90 =
        r.reference_to_future.transverse_survival.survival[i90];
    const double fr90 =
        r.future_to_reference.transverse_survival.survival[i90];
    const double long90 =
        r.longitudinal_survival.survival[i90];
    const double long10 =
        r.longitudinal_survival.survival[i10];
    const double full_rf90 =
        r.reference_to_future.full_survival.survival[i90];
    const double full_fr90 =
        r.future_to_reference.full_survival.survival[i90];

    if (std::abs(rf90 - 1.0) > 1e-12)
        throw std::runtime_error("reference->future s=.9 should survive transversely");

    if (std::abs(fr90 - 0.0) > 1e-12)
        throw std::runtime_error("future->reference s=.9 should be outside inner tube");

    if (std::abs(long90 - 1.0) > 1e-12)
        throw std::runtime_error("s=.9 should survive longitudinal gate");

    if (std::abs(long10 - 0.0) > 1e-12)
        throw std::runtime_error("s=.1 should be eaten by left-end retraction front");

    if (std::abs(full_rf90 - 1.0) > 1e-12 ||
        std::abs(full_fr90 - 0.0) > 1e-12)
        throw std::runtime_error("full directional comparison failed");

    std::cout << "Directional history comparison smoke test passed.\n";
    return 0;
}
