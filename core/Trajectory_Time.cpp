#include "Trajectory_Time.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>


namespace
{

std::size_t GetNumberOfFrames(
    const std::vector<ChainTrajectory>& trajectories
)
{
    if (trajectories.empty()) {
        return 0;
    }

    const ChainTrajectory& first_chain = trajectories.front();

    if (first_chain.frame_offsets.empty()) {
        throw std::runtime_error(
            "GetNumberOfFrames: trajectory contains no frame offsets."
        );
    }

    return first_chain.frame_offsets.size() - 1;
}


void ValidateFrameCounts(
    const std::vector<ChainTrajectory>& trajectories,
    std::size_t expected_frames
)
{
    for (const ChainTrajectory& trajectory : trajectories)
    {
        if (trajectory.frame_offsets.empty()) {
            throw std::runtime_error(
                "ValidateFrameCounts: trajectory contains no frame offsets."
            );
        }

        const std::size_t num_frames =
            trajectory.frame_offsets.size() - 1;

        if (num_frames != expected_frames) {
            throw std::runtime_error(
                "ValidateFrameCounts: chains contain different numbers of frames."
            );
        }
    }
}

} // namespace


void AssignUniformTimesteps(
    std::vector<ChainTrajectory>& trajectories,
    long initial_timestep,
    long dump_interval
)
{
    if (trajectories.empty()) {
        return;
    }

    if (dump_interval <= 0) {
        throw std::invalid_argument(
            "AssignUniformTimesteps: dump_interval must be positive."
        );
    }

    const std::size_t num_frames =
        GetNumberOfFrames(trajectories);

    ValidateFrameCounts(trajectories, num_frames);

    for (ChainTrajectory& trajectory : trajectories)
    {
        trajectory.timesteps.resize(num_frames);

        for (std::size_t frame = 0;
             frame < num_frames;
             ++frame)
        {
            trajectory.timesteps[frame] =
                initial_timestep +
                static_cast<long>(frame) * dump_interval;
        }
    }
}


void AssignTimesteps(
    std::vector<ChainTrajectory>& trajectories,
    const std::vector<long>& frame_timesteps
)
{
    if (trajectories.empty()) {
        return;
    }

    const std::size_t num_frames =
        GetNumberOfFrames(trajectories);

    ValidateFrameCounts(trajectories, num_frames);

    if (frame_timesteps.size() != num_frames) {
        throw std::invalid_argument(
            "AssignTimesteps: number of supplied timesteps does not match number of parsed frames."
        );
    }

    for (std::size_t i = 1;
         i < frame_timesteps.size();
         ++i)
    {
        if (frame_timesteps[i] <= frame_timesteps[i - 1]) {
            throw std::invalid_argument(
                "AssignTimesteps: timesteps must be strictly increasing."
            );
        }
    }

    for (ChainTrajectory& trajectory : trajectories) {
        trajectory.timesteps = frame_timesteps;
    }
}
