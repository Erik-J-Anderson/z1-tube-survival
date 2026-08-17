#include "Parse_Z1_File.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


namespace
{

double BoxVectorLength(
    const Box& box,
    std::size_t column)
{
    const double x = box.matrix.value[0][column];
    const double y = box.matrix.value[1][column];
    const double z = box.matrix.value[2][column];

    return std::sqrt(x * x + y * y + z * z);
}


double MaxOffDiagonalMagnitude(const Box& box)
{
    double max_value = 0.0;

    for (std::size_t row = 0; row < 3; ++row)
    {
        for (std::size_t col = 0; col < 3; ++col)
        {
            if (row == col) {
                continue;
            }

            max_value = std::max(
                max_value,
                std::abs(box.matrix.value[row][col])
            );
        }
    }

    return max_value;
}

} // namespace


int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cerr
            << "Usage:\n"
            << "  Serial_Smoke_Test_Deforming_Box "
            << "<Z1+SP.dat> "
            << "<expected_Lz_center> "
            << "<expected_amplitude> "
            << "<tolerance>\n\n"
            << "Example:\n"
            << "  ./Serial_Smoke_Test_Deforming_Box "
            << "Z1+SP.dat 76.1024 7.6 0.2\n";

        return 1;
    }

    try
    {
        const std::string filename = argv[1];
        const double expected_Lz_center = std::stod(argv[2]);
        const double expected_amplitude = std::abs(std::stod(argv[3]));
        const double tolerance = std::stod(argv[4]);

        if (expected_Lz_center <= 0.0) {
            throw std::invalid_argument(
                "expected_Lz_center must be positive."
            );
        }

        if (tolerance < 0.0) {
            throw std::invalid_argument(
                "tolerance must be nonnegative."
            );
        }

        const PrimitivePathTrajectory trajectory =
            parse_z1_file(filename);

        if (trajectory.chains.empty()) {
            throw std::runtime_error(
                "Parser returned zero chains."
            );
        }

        if (trajectory.chains.front().frame_offsets.empty()) {
            throw std::runtime_error(
                "First chain contains no frame offsets."
            );
        }

        const std::size_t num_frames =
            trajectory.chains.front().frame_offsets.size() - 1;

        if (num_frames == 0) {
            throw std::runtime_error(
                "Parser returned zero frames."
            );
        }

        if (trajectory.frame_boxes.size() != num_frames)
        {
            throw std::runtime_error(
                "Stored frame-box count does not match parsed frame count."
            );
        }

        for (const ChainTrajectory& chain : trajectory.chains)
        {
            if (chain.frame_offsets.size() != num_frames + 1)
            {
                throw std::runtime_error(
                    "Chains contain inconsistent frame counts."
                );
            }
        }

        double min_Lx = std::numeric_limits<double>::infinity();
        double max_Lx = -std::numeric_limits<double>::infinity();
        double min_Ly = std::numeric_limits<double>::infinity();
        double max_Ly = -std::numeric_limits<double>::infinity();
        double min_Lz = std::numeric_limits<double>::infinity();
        double max_Lz = -std::numeric_limits<double>::infinity();

        std::size_t min_Lz_frame = 0;
        std::size_t max_Lz_frame = 0;
        double largest_off_diagonal = 0.0;
        std::size_t affine_geometry_true_count = 0;

        for (std::size_t frame = 0;
             frame < trajectory.frame_boxes.size();
             ++frame)
        {
            const Box& box = trajectory.frame_boxes[frame];

            const double Lx = BoxVectorLength(box, 0);
            const double Ly = BoxVectorLength(box, 1);
            const double Lz = BoxVectorLength(box, 2);

            if (!std::isfinite(Lx) ||
                !std::isfinite(Ly) ||
                !std::isfinite(Lz))
            {
                throw std::runtime_error(
                    "Non-finite box length found."
                );
            }

            if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0)
            {
                throw std::runtime_error(
                    "Non-positive box length found."
                );
            }

            min_Lx = std::min(min_Lx, Lx);
            max_Lx = std::max(max_Lx, Lx);
            min_Ly = std::min(min_Ly, Ly);
            max_Ly = std::max(max_Ly, Ly);

            if (Lz < min_Lz)
            {
                min_Lz = Lz;
                min_Lz_frame = frame;
            }

            if (Lz > max_Lz)
            {
                max_Lz = Lz;
                max_Lz_frame = frame;
            }

            largest_off_diagonal =
                std::max(
                    largest_off_diagonal,
                    MaxOffDiagonalMagnitude(box)
                );

            if (box.affine_geometry_available) {
                ++affine_geometry_true_count;
            }
        }

        const double recovered_center =
            0.5 * (max_Lz + min_Lz);

        const double recovered_amplitude =
            0.5 * (max_Lz - min_Lz);

        const double delta_Lx = max_Lx - min_Lx;
        const double delta_Ly = max_Ly - min_Ly;

        const double center_error =
            std::abs(recovered_center - expected_Lz_center);

        const double amplitude_error =
            std::abs(recovered_amplitude - expected_amplitude);

        std::cout
            << "========================================\n"
            << "DEFORMING-BOX PARSER SMOKE TEST\n"
            << "========================================\n\n"
            << "File: " << filename << "\n\n"
            << "Parsed chains : " << trajectory.chains.size() << '\n'
            << "Parsed frames : " << num_frames << '\n'
            << "Stored boxes  : " << trajectory.frame_boxes.size() << "\n\n"
            << "Transverse dimensions:\n"
            << "  Lx min/max = " << min_Lx << " / " << max_Lx << '\n'
            << "  Ly min/max = " << min_Ly << " / " << max_Ly << '\n'
            << "  Delta Lx   = " << delta_Lx << '\n'
            << "  Delta Ly   = " << delta_Ly << "\n\n"
            << "Deforming z dimension:\n"
            << "  Lz min     = " << min_Lz
            << "  (frame " << min_Lz_frame << ")\n"
            << "  Lz max     = " << max_Lz
            << "  (frame " << max_Lz_frame << ")\n"
            << "  center     = " << recovered_center << '\n'
            << "  amplitude  = " << recovered_amplitude << "\n\n"
            << "Expected:\n"
            << "  Lz center  = " << expected_Lz_center << '\n'
            << "  amplitude  = " << expected_amplitude << '\n'
            << "  tolerance  = " << tolerance << "\n\n"
            << "Errors:\n"
            << "  center error    = " << center_error << '\n'
            << "  amplitude error = " << amplitude_error << "\n\n"
            << "Box-matrix diagnostics:\n"
            << "  largest off-diagonal element = "
            << largest_off_diagonal << '\n'
            << "  affine_geometry_available = true for "
            << affine_geometry_true_count << " / "
            << trajectory.frame_boxes.size() << " frames\n\n";

        if (delta_Lx > tolerance)
        {
            throw std::runtime_error(
                "Lx changed more than expected for z-only deformation."
            );
        }

        if (delta_Ly > tolerance)
        {
            throw std::runtime_error(
                "Ly changed more than expected for z-only deformation."
            );
        }

        if (recovered_amplitude <= 0.0)
        {
            throw std::runtime_error(
                "No z-box deformation was detected."
            );
        }

        if (center_error > tolerance)
        {
            throw std::runtime_error(
                "Recovered z-box center does not match expected value."
            );
        }

        if (amplitude_error > tolerance)
        {
            throw std::runtime_error(
                "Recovered z deformation amplitude does not match expected amplitude."
            );
        }

        std::cout
            << "ALL DEFORMING-BOX PARSER TESTS PASSED\n";
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "\nDEFORMING-BOX SMOKE TEST FAILED:\n"
            << error.what() << '\n';

        return 1;
    }

    return 0;
}
