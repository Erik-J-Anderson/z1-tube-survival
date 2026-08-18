#include "Parse_Z1_File.hpp"
#include "Trajectory_Time.hpp"
#include "Tube_Survival.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

double BoxVectorLength(const Box& box, std::size_t column)
{
    const double x = box.matrix.value[0][column];
    const double y = box.matrix.value[1][column];
    const double z = box.matrix.value[2][column];
    return std::sqrt(x * x + y * y + z * z);
}

std::vector<std::size_t> MakeSmokeLags(std::size_t num_frames)
{
    const std::vector<std::size_t> candidates{
        0, 1, 2, 5, 10, 20, 40, 80, 160, 320
    };

    std::vector<std::size_t> lags;
    for (const std::size_t lag : candidates) {
        if (lag < num_frames) {
            lags.push_back(lag);
        }
    }
    return lags;
}

void CheckProbability(double value, const char* label)
{
    constexpr double eps = 1.0e-12;
    if (!std::isfinite(value) || value < -eps || value > 1.0 + eps) {
        throw std::runtime_error(
            std::string(label) + " is outside [0,1] or non-finite.");
    }
}

} // namespace

int main(int argc, char* argv[])
{
    try
    {
        if (argc < 2 || argc > 4)
        {
            std::cerr
                << "Usage: " << argv[0]
                << " <Z1+SP.dat> [chain_index=0] [tube_diameter=5.0]\n";
            return 2;
        }

        const std::string filename = argv[1];
        const std::size_t chain_index =
            (argc >= 3) ? static_cast<std::size_t>(std::stoull(argv[2])) : 0;
        const double tube_diameter =
            (argc >= 4) ? std::stod(argv[3]) : 5.0;

        if (tube_diameter <= 0.0) {
            throw std::invalid_argument("tube_diameter must be positive.");
        }

        std::cout
            << "========================================\n"
            << "REAL DEFORMING-TRAJECTORY AFFINE TEST\n"
            << "========================================\n\n"
            << "Parsing: " << filename << "\n";

        PrimitivePathTrajectory trajectory = parse_z1_file(filename);

        if (trajectory.chains.empty()) {
            throw std::runtime_error("No chains were parsed.");
        }
        if (trajectory.frame_boxes.empty()) {
            throw std::runtime_error("No frame boxes were parsed.");
        }
        if (chain_index >= trajectory.chains.size()) {
            throw std::out_of_range("chain_index exceeds parsed chain count.");
        }

        const std::size_t num_frames = trajectory.frame_boxes.size();

        for (const ChainTrajectory& chain : trajectory.chains) {
            if (chain.frame_offsets.size() != num_frames + 1) {
                throw std::runtime_error(
                    "Chain frame count does not match frame_boxes count.");
            }
        }

        // For this smoke test, physical time is not needed.  Assign one
        // arbitrary timestep unit per stored frame so the survival API has
        // a valid monotonically increasing time axis.
        AssignUniformTimesteps(trajectory.chains, 0, 1);

        double min_lz = BoxVectorLength(trajectory.frame_boxes.front(), 2);
        double max_lz = min_lz;
        for (const Box& box : trajectory.frame_boxes) {
            const double lz = BoxVectorLength(box, 2);
            min_lz = std::min(min_lz, lz);
            max_lz = std::max(max_lz, lz);
        }

        const ChainTrajectory& chain = trajectory.chains[chain_index];
        const std::vector<std::size_t> lag_frames = MakeSmokeLags(num_frames);
        const std::vector<double> tube_diameters{tube_diameter};

        std::cout
            << "Parsed chains : " << trajectory.chains.size() << "\n"
            << "Parsed frames : " << num_frames << "\n"
            << "Chain index   : " << chain_index << "\n"
            << "Chain ID      : " << chain.chain_id << "\n"
            << "Tube diameter : " << tube_diameter << "\n"
            << "Lz min/max    : " << min_lz << " / " << max_lz << "\n"
            << "Delta Lz      : " << (max_lz - min_lz) << "\n\n";

        if (max_lz - min_lz <= 0.0) {
            throw std::runtime_error(
                "This input does not appear to contain a deforming z box.");
        }

        std::cout << "Computing UNCORRECTED survival...\n";
        const SegmentSurvivalFunction uncorrected_segments =
            ComputeSegmentSurvivalFunction(
                chain,
                lag_frames,
                tube_diameters);

        std::cout << "Computing AFFINE-CORRECTED survival...\n";
        const SegmentSurvivalOptions affine_options{
            .apply_affine_correction = true
        };

        const SegmentSurvivalFunction corrected_segments =
            ComputeSegmentSurvivalFunction(
                chain,
                trajectory.frame_boxes,
                lag_frames,
                tube_diameters,
                affine_options);

        const TubeSurvivalFunction uncorrected_tube =
            ComputeTubeSurvivalFunction(uncorrected_segments);

        const TubeSurvivalFunction corrected_tube =
            ComputeTubeSurvivalFunction(corrected_segments);

        if (uncorrected_tube.sample_counts != corrected_tube.sample_counts) {
            throw std::runtime_error(
                "Corrected and uncorrected calculations used different sample counts.");
        }

        constexpr std::size_t center_s_index = NUM_SAMPLE_POINTS / 2;
        double max_abs_tube_delta = 0.0;

        std::ofstream csv("real_affine_comparison.csv");
        if (!csv) {
            throw std::runtime_error("Could not open real_affine_comparison.csv.");
        }

        csv
            << "lag_frame,sample_count,uncorrected_tube_survival,"
            << "corrected_tube_survival,delta_tube_survival,"
            << "uncorrected_center_survival,corrected_center_survival,"
            << "delta_center_survival\n";

        std::cout
            << "\n"
            << std::setw(10) << "lag"
            << std::setw(12) << "origins"
            << std::setw(16) << "uncorr tube"
            << std::setw(16) << "corr tube"
            << std::setw(14) << "delta"
            << std::setw(16) << "uncorr s=.5"
            << std::setw(16) << "corr s=.5"
            << "\n";

        std::cout << std::string(100, '-') << "\n";

        for (std::size_t lag_index = 0;
             lag_index < lag_frames.size();
             ++lag_index)
        {
            const double u_tube = uncorrected_tube.survival[lag_index];
            const double c_tube = corrected_tube.survival[lag_index];
            const double d_tube = c_tube - u_tube;

            const std::size_t center_flat_index =
                lag_index * NUM_SAMPLE_POINTS + center_s_index;

            const double u_center =
                uncorrected_segments.survival[center_flat_index];
            const double c_center =
                corrected_segments.survival[center_flat_index];
            const double d_center = c_center - u_center;

            CheckProbability(u_tube, "uncorrected tube survival");
            CheckProbability(c_tube, "corrected tube survival");
            CheckProbability(u_center, "uncorrected center survival");
            CheckProbability(c_center, "corrected center survival");

            max_abs_tube_delta =
                std::max(max_abs_tube_delta, std::abs(d_tube));

            const std::uint64_t count =
                corrected_tube.sample_counts[lag_index];

            std::cout
                << std::setw(10) << lag_frames[lag_index]
                << std::setw(12) << count
                << std::setw(16) << std::setprecision(7) << u_tube
                << std::setw(16) << std::setprecision(7) << c_tube
                << std::setw(14) << std::setprecision(7) << d_tube
                << std::setw(16) << std::setprecision(7) << u_center
                << std::setw(16) << std::setprecision(7) << c_center
                << "\n";

            csv
                << lag_frames[lag_index] << ','
                << count << ','
                << std::setprecision(17) << u_tube << ','
                << c_tube << ','
                << d_tube << ','
                << u_center << ','
                << c_center << ','
                << d_center << '\n';
        }

        // Zero lag must be the identity transformation and the path must
        // survive inside its own tube.
        if (!lag_frames.empty() && lag_frames.front() == 0)
        {
            constexpr double zero_lag_tolerance = 1.0e-12;
            const double u0 = uncorrected_tube.survival.front();
            const double c0 = corrected_tube.survival.front();

            if (std::abs(u0 - 1.0) > zero_lag_tolerance ||
                std::abs(c0 - 1.0) > zero_lag_tolerance ||
                std::abs(u0 - c0) > zero_lag_tolerance)
            {
                throw std::runtime_error(
                    "Zero-lag identity check failed.");
            }
        }

        std::cout
            << "\nMaximum |corrected - uncorrected| tube survival = "
            << max_abs_tube_delta << "\n";

        if (max_abs_tube_delta < 1.0e-6) {
            std::cout
                << "WARNING: affine correction made almost no difference for "
                << "this chain/diameter/lag set. This is not automatically a failure.\n";
        }

        std::cout
            << "Wrote real_affine_comparison.csv\n\n"
            << "ALL REAL-AFFINE SMOKE CHECKS PASSED\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "\nREAL-AFFINE SMOKE TEST FAILED:\n"
            << error.what() << '\n';
        return 1;
    }
}
