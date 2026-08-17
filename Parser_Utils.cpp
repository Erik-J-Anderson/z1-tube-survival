#include "Parser_Utils.hpp"

#include <charconv>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <system_error>


namespace
{

bool is_whitespace(char c)
{
    return
        c == ' '  ||
        c == '\t' ||
        c == '\n' ||
        c == '\r' ||
        c == '\f' ||
        c == '\v';
}


std::string_view trim_left(
    std::string_view text
)
{
    while (!text.empty() && is_whitespace(text.front())) {
        text.remove_prefix(1);
    }

    return text;
}


std::string_view trim(
    std::string_view text
)
{
    text = trim_left(text);

    while (!text.empty() && is_whitespace(text.back())) {
        text.remove_suffix(1);
    }

    return text;
}


bool consume_double(
    std::string_view& text,
    double& value
)
{
    text = trim_left(text);

    if (text.empty()) {
        return false;
    }

    const char* begin = text.data();
    const char* end = text.data() + text.size();

    const auto result =
        std::from_chars(
            begin,
            end,
            value,
            std::chars_format::general
        );

    if (result.ec != std::errc{}) {
        return false;
    }

    if (!std::isfinite(value)) {
        return false;
    }

    const std::size_t consumed =
        static_cast<std::size_t>(result.ptr - begin);

    text.remove_prefix(consumed);

    return true;
}

} // namespace


bool is_blank_line(
    std::string_view line
)
{
    return trim(line).empty();
}


bool parse_count_line(
    std::string_view line,
    std::size_t& count
)
{
    line = trim(line);

    if (line.empty()) {
        return false;
    }

    const char* begin = line.data();
    const char* end = line.data() + line.size();

    const auto result =
        std::from_chars(
            begin,
            end,
            count
        );

    if (result.ec != std::errc{}) {
        return false;
    }

    std::string_view remainder(
        result.ptr,
        static_cast<std::size_t>(end - result.ptr)
    );

    return trim(remainder).empty();
}


bool parse_box_line(
    std::string_view line,
    Box& box
)
{
    double values[9]{};
    std::size_t count = 0;

    while (count < 9 && consume_double(line, values[count])) {
        ++count;
    }

    if (!trim(line).empty()) {
        return false;
    }

    if (count != 3 && count != 6 && count != 9) {
        return false;
    }

    double xlo = 0.0;
    double ylo = 0.0;
    double zlo = 0.0;
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;

    if (count == 3) {
        lx = values[0];
        ly = values[1];
        lz = values[2];
    }
    else if (count == 6) {
        lx = values[0];
        ly = values[1];
        lz = values[2];
        xy = values[3];
        xz = values[4];
        yz = values[5];
    }
    else {
        xlo = values[0];
        ylo = values[1];
        zlo = values[2];
        lx = values[3];
        ly = values[4];
        lz = values[5];
        xy = values[6];
        xz = values[7];
        yz = values[8];
    }

    if (lx <= 0.0 || ly <= 0.0 || lz <= 0.0) {
        return false;
    }

    box.origin = Vec3{xlo, ylo, zlo};
    box.matrix = Mat3{{
        {lx,  xy,  xz},
        {0.0, ly,  yz},
        {0.0, 0.0, lz}
    }};
    box.affine_geometry_available = (count == 6 || count == 9);

    return true;
}


bool parse_node(
    std::string_view line,
    double& x,
    double& y,
    double& z
)
{
    if (!consume_double(line, x)) {
        return false;
    }

    if (!consume_double(line, y)) {
        return false;
    }

    if (!consume_double(line, z)) {
        return false;
    }

    // Additional Z1+ node columns, if present, are intentionally ignored.
    return true;
}
