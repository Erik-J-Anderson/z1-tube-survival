#include "Survival_Accumulator.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>

// ------------------------------------------------------------
// Initialize an ensemble accumulator from one chain's result.
// ------------------------------------------------------------
void InitializeAccumulator(
    SegmentSurvivalFunction& accumulator,
    const SegmentSurvivalFunction& result)
{
    accumulator = result;

    std::fill(
        accumulator.survival.begin(),
        accumulator.survival.end(),
        0.0
    );

    std::fill(
        accumulator.sample_counts.begin(),
        accumulator.sample_counts.end(),
        0
    );
}


// ------------------------------------------------------------
// Add one chain to the ensemble.
//
// Weight by the number of valid time origins for each lag.
// ------------------------------------------------------------
void AccumulateSegmentSurvival(SegmentSurvivalFunction& accumulator,
    const SegmentSurvivalFunction& result)
{
    const std::size_t num_lags =
        result.lag_times.size();

    const std::size_t num_diameters =
        result.tube_diameters.size();

    for (std::size_t lag_index = 0;
        lag_index < num_lags;
        ++lag_index)
    {
        const std::uint64_t count =
            result.sample_counts[lag_index];

        accumulator.sample_counts[lag_index] += count;

        for (std::size_t d = 0;
            d < num_diameters;
            ++d)
        {
            for (std::size_t s = 0;
                s < NUM_SAMPLE_POINTS;
                ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS
                    + lag_index * NUM_SAMPLE_POINTS
                    + s;

                accumulator.survival[index] +=
                    result.survival[index]
                    * static_cast<double>(count);
            }
        }
    }
}


// ------------------------------------------------------------
// Convert accumulated counts back into ensemble probabilities.
// ------------------------------------------------------------
void FinalizeAccumulator(
    SegmentSurvivalFunction& accumulator)
{
    const std::size_t num_lags =
        accumulator.lag_times.size();

    const std::size_t num_diameters =
        accumulator.tube_diameters.size();

    for (std::size_t lag_index = 0;
        lag_index < num_lags;
        ++lag_index)
    {
        const std::uint64_t count =
            accumulator.sample_counts[lag_index];

        if (count == 0) {
            continue;
        }

        for (std::size_t d = 0;
            d < num_diameters;
            ++d)
        {
            for (std::size_t s = 0;
                s < NUM_SAMPLE_POINTS;
                ++s)
            {
                const std::size_t index =
                    d * num_lags * NUM_SAMPLE_POINTS
                    + lag_index * NUM_SAMPLE_POINTS
                    + s;

                accumulator.survival[index] /=
                    static_cast<double>(count);
            }
        }
    }
}
