#include "Survival_IO.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <stdexcept>


void WriteSegmentSurvivalFunction(
    const SegmentSurvivalFunction& result,
    const std::string& filename
)
{
    const std::size_t num_lags =
        result.lag_times.size();

    const std::size_t num_diameters =
        result.tube_diameters.size();

    const std::size_t expected_size =
        num_diameters *
        num_lags *
        NUM_SAMPLE_POINTS;

    if (result.survival.size() != expected_size) {
        throw std::invalid_argument(
            "WriteSegmentSurvivalFunction: survival array has incorrect size."
        );
    }

    if (result.sample_counts.size() != num_lags) {
        throw std::invalid_argument(
            "WriteSegmentSurvivalFunction: sample_counts has incorrect size."
        );
    }

    std::ofstream file(filename);

    if (!file) {
        throw std::runtime_error(
            "WriteSegmentSurvivalFunction: could not open output file: " + filename
        );
    }

    file << std::setprecision(17);

    file
        << "tube_diameter,"
        << "lag_time,"
        << "segment_index,"
        << "s_fraction,"
        << "phi,"
        << "sample_count\n";

    for (std::size_t diameter_index = 0;
         diameter_index < num_diameters;
         ++diameter_index)
    {
        for (std::size_t lag_index = 0;
             lag_index < num_lags;
             ++lag_index)
        {
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

                const double s_fraction =
                    static_cast<double>(segment_index) /
                    static_cast<double>(NUM_SAMPLE_POINTS - 1);

                file
                    << result.tube_diameters[diameter_index]
                    << ','
                    << result.lag_times[lag_index]
                    << ','
                    << segment_index
                    << ','
                    << s_fraction
                    << ','
                    << result.survival[index]
                    << ','
                    << result.sample_counts[lag_index]
                    << '\n';
            }
        }
    }
}


void WriteTubeSurvivalFunction(
    const TubeSurvivalFunction& result,
    const std::string& filename
)
{
    const std::size_t num_lags =
        result.lag_times.size();

    const std::size_t num_diameters =
        result.tube_diameters.size();

    const std::size_t expected_size =
        num_diameters * num_lags;

    if (result.survival.size() != expected_size) {
        throw std::invalid_argument(
            "WriteTubeSurvivalFunction: survival array has incorrect size."
        );
    }

    if (result.sample_counts.size() != num_lags) {
        throw std::invalid_argument(
            "WriteTubeSurvivalFunction: sample_counts has incorrect size."
        );
    }

    std::ofstream file(filename);

    if (!file) {
        throw std::runtime_error(
            "WriteTubeSurvivalFunction: could not open output file: " + filename
        );
    }

    file << std::setprecision(17);

    file
        << "tube_diameter,"
        << "lag_time,"
        << "mu,"
        << "sample_count\n";

    for (std::size_t diameter_index = 0;
         diameter_index < num_diameters;
         ++diameter_index)
    {
        for (std::size_t lag_index = 0;
             lag_index < num_lags;
             ++lag_index)
        {
            const std::size_t index =
                diameter_index * num_lags +
                lag_index;

            file
                << result.tube_diameters[diameter_index]
                << ','
                << result.lag_times[lag_index]
                << ','
                << result.survival[index]
                << ','
                << result.sample_counts[lag_index]
                << '\n';
        }
    }
}

void WriteEndRetractionFunction(
    const EndRetractionFunction& result,
    const std::string& filename
)
{
    const std::size_t num_diameters =
        result.tube_diameters.size();

    const std::size_t num_lags =
        result.lag_frames.size();

    const std::size_t num_s =
        result.s_fraction.size();


    const std::size_t expected_field_size =
        num_diameters
        *
        num_lags
        *
        num_s;

    const std::size_t expected_summary_size =
        num_diameters
        *
        num_lags;


    if (
        result.end0_reached.size()
        != expected_field_size
        ||
        result.end1_reached.size()
        != expected_field_size
        ||
        result.combined_reached.size()
        != expected_field_size
        )
    {
        throw std::runtime_error(
            "WriteEndRetractionFunction: "
            "retraction field dimensions are inconsistent."
        );
    }


    if (
        result.mean_max_depth_end0.size()
        != expected_summary_size
        ||
        result.mean_max_depth_end1.size()
        != expected_summary_size
        ||
        result.mean_max_depth_combined.size()
        != expected_summary_size
        ||
        result.sample_counts.size()
        != expected_summary_size
        )
    {
        throw std::runtime_error(
            "WriteEndRetractionFunction: "
            "retraction summary dimensions are inconsistent."
        );
    }


    std::ofstream output(filename);

    if (!output)
    {
        throw std::runtime_error(
            "WriteEndRetractionFunction: "
            "could not open output file."
        );
    }


    output << std::setprecision(12);


    output
        << "tube_diameter,"
        << "lag_frame,"
        << "s_index,"
        << "s_fraction,"
        << "end0_reached_probability,"
        << "end1_reached_probability,"
        << "combined_reached_probability,"
        << "mean_max_depth_end0,"
        << "mean_max_depth_end1,"
        << "mean_max_depth_combined,"
        << "sample_count\n";


    for (std::size_t d = 0;
        d < num_diameters;
        ++d)
    {
        for (std::size_t lag_index = 0;
            lag_index < num_lags;
            ++lag_index)
        {
            const std::size_t summary_index =
                d * num_lags
                +
                lag_index;


            for (std::size_t s_index = 0;
                s_index < num_s;
                ++s_index)
            {
                const std::size_t field_index =
                    (
                        d * num_lags
                        +
                        lag_index
                        )
                    * num_s
                    +
                    s_index;


                output
                    << result.tube_diameters[d]
                    << ","
                    << result.lag_frames[lag_index]
                    << ","
                    << s_index
                    << ","
                    << result.s_fraction[s_index]
                    << ","
                    << result.end0_reached[
                        field_index
                    ]
                    << ","
                    << result.end1_reached[
                        field_index
                    ]
                    << ","
                    << result.combined_reached[
                        field_index
                    ]
                    << ","
                    << result.mean_max_depth_end0[
                        summary_index
                    ]
                    << ","
                    << result.mean_max_depth_end1[
                        summary_index
                    ]
                    << ","
                    << result.mean_max_depth_combined[
                        summary_index
                    ]
                    << ","
                    << result.sample_counts[
                        summary_index
                    ]
                    << "\n";
            }
        }
    }
}