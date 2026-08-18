#pragma once

#include "Primitive_Path_Trajectory.hpp"

#include <cstddef>
#include <string_view>
#include <vector>

enum class ParseState
{
    ChainCount,
    Box,
    ChainNodeCount,
    Nodes
};


bool is_blank_line(
    std::string_view line
);


bool parse_count_line(
    std::string_view line,
    std::size_t& count
);


bool parse_box_line(
    std::string_view line,
    Box& box
);


bool parse_node(
    std::string_view line,
    double& x,
    double& y,
    double& z
);
