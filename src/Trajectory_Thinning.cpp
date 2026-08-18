#include "Trajectory_Thinning.hpp"

#include <cstddef>
#include <stdexcept>


ChainTrajectory MakeStridedTrajectory(
    const ChainTrajectory& trajectory,
    std::size_t frame_stride
)
{
    if (frame_stride == 0)
    {
        throw std::invalid_argument(
            "MakeStridedTrajectory: frame_stride must be >= 1."
        );
    }


    // --------------------------------------------------------
    // Validate trajectory structure
    // --------------------------------------------------------

    if (trajectory.frame_offsets.empty())
    {
        throw std::invalid_argument(
            "MakeStridedTrajectory: trajectory has no frame offsets."
        );
    }


    const std::size_t num_frames =
        trajectory.frame_offsets.size() - 1;


    if (trajectory.timesteps.size() != num_frames)
    {
        throw std::invalid_argument(
            "MakeStridedTrajectory: timesteps.size() does not match "
            "the number of frames."
        );
    }


    if (trajectory.frame_offsets.back() != trajectory.nodes.size())
    {
        throw std::invalid_argument(
            "MakeStridedTrajectory: final frame offset does not match "
            "nodes.size()."
        );
    }


    // --------------------------------------------------------
    // Construct output trajectory
    // --------------------------------------------------------

    ChainTrajectory result;

    result.chain_id = trajectory.chain_id;


    // Number of frames that will survive:
    //
    // physical frames:
    //
    //     0, stride, 2*stride, ...
    //
    const std::size_t num_output_frames =
        (num_frames + frame_stride - 1)
        / frame_stride;


    result.timesteps.reserve(
        num_output_frames
    );

    result.frame_offsets.reserve(
        num_output_frames + 1
    );


    // The first frame always begins at node zero.
    result.frame_offsets.push_back(0);


    // --------------------------------------------------------
    // Copy selected frames
    // --------------------------------------------------------

    for (
        std::size_t frame = 0;
        frame < num_frames;
        frame += frame_stride
        )
    {
        const std::size_t begin =
            trajectory.frame_offsets[frame];

        const std::size_t end =
            trajectory.frame_offsets[frame + 1];


        if (begin > end)
        {
            throw std::runtime_error(
                "MakeStridedTrajectory: invalid frame offsets."
            );
        }


        if (end > trajectory.nodes.size())
        {
            throw std::runtime_error(
                "MakeStridedTrajectory: frame offset exceeds nodes.size()."
            );
        }


        // Preserve the actual physical timestep.
        result.timesteps.push_back(
            trajectory.timesteps[frame]
        );


        // Copy all primitive-path nodes belonging to this frame.
        result.nodes.insert(
            result.nodes.end(),
            trajectory.nodes.begin()
            + static_cast<std::ptrdiff_t>(begin),
            trajectory.nodes.begin()
            + static_cast<std::ptrdiff_t>(end)
        );


        // Record where the NEXT frame will begin.
        result.frame_offsets.push_back(
            result.nodes.size()
        );
    }


    return result;
}