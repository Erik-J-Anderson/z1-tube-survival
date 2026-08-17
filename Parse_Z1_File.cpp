#include "Parse_Z1_File.hpp"
#include "Parser_Utils.hpp"

#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{

[[noreturn]]
void throw_parse_error(
    const std::string& filename,
    std::size_t line_number,
    const std::string& message
)
{
    std::ostringstream error;

    error
        << "Error parsing "
        << filename
        << " at line "
        << line_number
        << ": "
        << message;

    throw std::runtime_error(error.str());
}

} // namespace


//std::vector<ChainTrajectory>
PrimitivePathTrajectory
parse_z1_file(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error(
            "parse_z1_file: could not open file: " + filename
        );
    }

    PrimitivePathTrajectory parsed_trajectory;

    std::vector<ChainTrajectory>& trajectories =
        parsed_trajectory.chains;

    ParseState state = ParseState::ChainCount;

    std::string line;
    std::size_t line_number = 0;
    std::size_t frame_index = 0;
    std::size_t chains_in_frame = 0;
    std::size_t current_chain = 0;
    std::size_t nodes_expected = 0;
    std::size_t nodes_read = 0;

    while (std::getline(file, line))
    {
        ++line_number;

        if (is_blank_line(line)) {
            continue;
        }

        switch (state)
        {
        case ParseState::ChainCount:
        {
            std::size_t parsed_chain_count = 0;

            if (!parse_count_line(line, parsed_chain_count)) {
                throw_parse_error(
                    filename,
                    line_number,
                    "expected number of chains."
                );
            }

            if (parsed_chain_count == 0) {
                throw_parse_error(
                    filename,
                    line_number,
                    "frame contains zero chains."
                );
            }

            if (frame_index == 0)
            {
                chains_in_frame = parsed_chain_count;
                trajectories.resize(chains_in_frame);

                for (std::size_t chain = 0;
                     chain < chains_in_frame;
                     ++chain)
                {
                    trajectories[chain].chain_id = chain + 1;
                    trajectories[chain].frame_offsets.push_back(0);
                }
            }
            else if (parsed_chain_count != chains_in_frame)
            {
                throw_parse_error(
                    filename,
                    line_number,
                    "number of chains changed between frames."
                );
            }

            current_chain = 0;
            state = ParseState::Box;
            break;
        }

        case ParseState::Box:
        {
            Box box;

            if (!parse_box_line(line, box)) {
                throw_parse_error(
                    filename,
                    line_number,
                    "expected 3, 6, or 9 box values."
                );
            }

            parsed_trajectory.frame_boxes.push_back(box);

            state = ParseState::ChainNodeCount;
            break;
        }

        case ParseState::ChainNodeCount:
        {
            if (!parse_count_line(line, nodes_expected)) {
                throw_parse_error(
                    filename,
                    line_number,
                    "expected number of primitive-path nodes for current chain."
                );
            }

            if (nodes_expected < 2) {
                throw_parse_error(
                    filename,
                    line_number,
                    "primitive path contains fewer than two nodes."
                );
            }

            nodes_read = 0;
            state = ParseState::Nodes;
            break;
        }

        case ParseState::Nodes:
        {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;

            if (!parse_node(line, x, y, z)) {
                throw_parse_error(
                    filename,
                    line_number,
                    "invalid primitive-path node record."
                );
            }

            const bool is_end =
                nodes_read == 0 ||
                nodes_read + 1 == nodes_expected;

            ChainTrajectory& trajectory =
                trajectories[current_chain];

            trajectory.nodes.push_back(
                PPNode{x, y, z, is_end}
            );

            ++nodes_read;

            if (nodes_read == nodes_expected)
            {
                trajectory.frame_offsets.push_back(
                    trajectory.nodes.size()
                );

                ++current_chain;

                if (current_chain == chains_in_frame)
                {
                    ++frame_index;
                    state = ParseState::ChainCount;
                }
                else
                {
                    state = ParseState::ChainNodeCount;
                }
            }

            break;
        }
        }
    }

    if (state != ParseState::ChainCount) {
        throw std::runtime_error(
            "parse_z1_file: unexpected end of file while parsing a frame."
        );
    }

    if (frame_index == 0) {
        throw std::runtime_error(
            "parse_z1_file: file contained no complete frames."
        );
    }

    for (const ChainTrajectory& trajectory : trajectories)
    {
        if (trajectory.frame_offsets.size() != frame_index + 1) {
            throw std::runtime_error(
                "parse_z1_file: inconsistent frame_offsets array across chains."
            );
        }

        if (!trajectory.timesteps.empty()) {
            throw std::runtime_error(
                "parse_z1_file: parser unexpectedly modified trajectory timesteps."
            );
        }
    }

    return parsed_trajectory;
}
