#include "Input_Parsers.hpp"
#include "Geometry_Utils.hpp"

#include <sstream>
#include <stdexcept>


// ------------------------------------------------------------
// Parse a comma-separated list of doubles.
// Example: "5.0,7.0,9.0".
// ------------------------------------------------------------
std::vector<double> ParseDoubleList(
    const std::string& input)
{
    std::vector<double> values;

    std::stringstream stream(input);
    std::string token;

    while (std::getline(stream, token, ','))
    {
        if (token.empty()) {
            continue;
        }

        values.push_back(
            std::stod(token)
        );
    }

    return values;
}


// ------------------------------------------------------------
// Parse three comma-separated values into a Vec3.
// ------------------------------------------------------------
Vec3 ParseVec3(
    const std::string& input)
{
    const std::vector<double> values =
        ParseDoubleList(input);

    if (values.size() != 3)
    {
        throw std::invalid_argument(
            "Box center must contain exactly three values: x,y,z"
        );
    }

    return Vec3{
        values[0],
        values[1],
        values[2]
    };
}


// ------------------------------------------------------------
// Reconstruct each box origin from a fixed physical box center.
// ------------------------------------------------------------
void SetFixedBoxCenter(
    std::vector<Box>& frame_boxes,
    const Vec3& center)
{
    const Vec3 fractional_center{
        0.5,
        0.5,
        0.5
    };

    for (Box& box : frame_boxes)
    {
        const Vec3 half_box =
            geometry::Multiply(
                box.matrix,
                fractional_center
            );

        box.origin = Vec3{
            center.x - half_box.x,
            center.y - half_box.y,
            center.z - half_box.z
        };
    }
}
