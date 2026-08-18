#include "MPI_Comms.hpp"

#include <array>
#include <cstdint>
#include <vector>
#include <algorithm>


namespace
{
    constexpr int CHAIN_COUNT_TAG = 100;
    constexpr int CHAIN_META_TAG = 101;
    constexpr int TIMESTEP_TAG = 102;
    constexpr int OFFSET_TAG = 103;
    constexpr int NODE_XYZ_TAG = 104;
    constexpr int NODE_END_TAG = 105;
}


// ------------------------------------------------------------
// Broadcast frame boxes to every rank.
// ------------------------------------------------------------
void BroadcastFrameBoxes(
    std::vector<Box>& frame_boxes,
    int root_rank,
    MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::uint64_t num_boxes =
        static_cast<std::uint64_t>(frame_boxes.size());

    MPI_Bcast(
        &num_boxes,
        1,
        MPI_UINT64_T,
        root_rank,
        comm
    );

    if (rank != root_rank) {
        frame_boxes.resize(num_boxes);
    }

    // 3 origin values + 9 matrix values per box.
    std::vector<double> box_data(12 * num_boxes);

    std::vector<unsigned char> affine_flags(num_boxes);

    if (rank == root_rank)
    {
        for (std::size_t i = 0; i < num_boxes; ++i)
        {
            const Box& box = frame_boxes[i];

            const std::size_t base = 12 * i;

            box_data[base + 0] = box.origin.x;
            box_data[base + 1] = box.origin.y;
            box_data[base + 2] = box.origin.z;

            std::size_t k = 3;

            for (std::size_t row = 0; row < 3; ++row)
            {
                for (std::size_t col = 0; col < 3; ++col)
                {
                    box_data[base + k] =
                        box.matrix.value[row][col];

                    ++k;
                }
            }

            affine_flags[i] =
                box.affine_geometry_available ? 1 : 0;
        }
    }

    MPI_Bcast(
        box_data.data(),
        static_cast<int>(box_data.size()),
        MPI_DOUBLE,
        root_rank,
        comm
    );

    MPI_Bcast(
        affine_flags.data(),
        static_cast<int>(affine_flags.size()),
        MPI_UNSIGNED_CHAR,
        root_rank,
        comm
    );

    if (rank != root_rank)
    {
        for (std::size_t i = 0; i < num_boxes; ++i)
        {
            Box& box = frame_boxes[i];

            const std::size_t base = 12 * i;

            box.origin.x = box_data[base + 0];
            box.origin.y = box_data[base + 1];
            box.origin.z = box_data[base + 2];

            std::size_t k = 3;

            for (std::size_t row = 0; row < 3; ++row)
            {
                for (std::size_t col = 0; col < 3; ++col)
                {
                    box.matrix.value[row][col] =
                        box_data[base + k];

                    ++k;
                }
            }

            box.affine_geometry_available =
                affine_flags[i] != 0;
        }
    }
}


// ------------------------------------------------------------
// Send one ChainTrajectory.
// ------------------------------------------------------------
void SendChainTrajectory(
    const ChainTrajectory& chain,
    int dest_rank,
    MPI_Comm comm)
{
    std::array<std::uint64_t, 4> metadata{
        static_cast<std::uint64_t>(chain.chain_id),
        static_cast<std::uint64_t>(chain.timesteps.size()),
        static_cast<std::uint64_t>(chain.frame_offsets.size()),
        static_cast<std::uint64_t>(chain.nodes.size())
    };

    MPI_Send(
        metadata.data(),
        metadata.size(),
        MPI_UINT64_T,
        dest_rank,
        CHAIN_META_TAG,
        comm
    );

    MPI_Send(
        chain.timesteps.data(),
        chain.timesteps.size(),
        MPI_LONG,
        dest_rank,
        TIMESTEP_TAG,
        comm
    );


    // size_t is not guaranteed to have a portable MPI datatype,
    // so convert offsets to uint64_t.
    std::vector<std::uint64_t> offsets(
        chain.frame_offsets.begin(),
        chain.frame_offsets.end()
    );

    MPI_Send(
        offsets.data(),
        offsets.size(),
        MPI_UINT64_T,
        dest_rank,
        OFFSET_TAG,
        comm
    );


    // Pack PPNode fields explicitly.
    std::vector<double> node_xyz(
        3 * chain.nodes.size()
    );

    std::vector<unsigned char> node_end(
        chain.nodes.size()
    );

    for (std::size_t i = 0; i < chain.nodes.size(); ++i)
    {
        node_xyz[3 * i + 0] = chain.nodes[i].x;
        node_xyz[3 * i + 1] = chain.nodes[i].y;
        node_xyz[3 * i + 2] = chain.nodes[i].z;

        node_end[i] =
            chain.nodes[i].end ? 1 : 0;
    }

    MPI_Send(
        node_xyz.data(),
        node_xyz.size(),
        MPI_DOUBLE,
        dest_rank,
        NODE_XYZ_TAG,
        comm
    );

    MPI_Send(
        node_end.data(),
        node_end.size(),
        MPI_UNSIGNED_CHAR,
        dest_rank,
        NODE_END_TAG,
        comm
    );
}


// ------------------------------------------------------------
// Receive one ChainTrajectory.
// ------------------------------------------------------------
ChainTrajectory ReceiveChainTrajectory(
    int source_rank,
    MPI_Comm comm)
{
    std::array<std::uint64_t, 4> metadata{};

    MPI_Recv(
        metadata.data(),
        metadata.size(),
        MPI_UINT64_T,
        source_rank,
        CHAIN_META_TAG,
        comm,
        MPI_STATUS_IGNORE
    );

    ChainTrajectory chain;

    chain.chain_id =
        static_cast<std::size_t>(metadata[0]);

    const std::size_t num_timesteps =
        static_cast<std::size_t>(metadata[1]);

    const std::size_t num_offsets =
        static_cast<std::size_t>(metadata[2]);

    const std::size_t num_nodes =
        static_cast<std::size_t>(metadata[3]);


    chain.timesteps.resize(num_timesteps);
    chain.frame_offsets.resize(num_offsets);
    chain.nodes.resize(num_nodes);


    MPI_Recv(
        chain.timesteps.data(),
        chain.timesteps.size(),
        MPI_LONG,
        source_rank,
        TIMESTEP_TAG,
        comm,
        MPI_STATUS_IGNORE
    );


    std::vector<std::uint64_t> offsets(num_offsets);

    MPI_Recv(
        offsets.data(),
        offsets.size(),
        MPI_UINT64_T,
        source_rank,
        OFFSET_TAG,
        comm,
        MPI_STATUS_IGNORE
    );

    for (std::size_t i = 0; i < num_offsets; ++i) {
        chain.frame_offsets[i] =
            static_cast<std::size_t>(offsets[i]);
    }


    std::vector<double> node_xyz(
        3 * num_nodes
    );

    std::vector<unsigned char> node_end(
        num_nodes
    );

    MPI_Recv(
        node_xyz.data(),
        node_xyz.size(),
        MPI_DOUBLE,
        source_rank,
        NODE_XYZ_TAG,
        comm,
        MPI_STATUS_IGNORE
    );

    MPI_Recv(
        node_end.data(),
        node_end.size(),
        MPI_UNSIGNED_CHAR,
        source_rank,
        NODE_END_TAG,
        comm,
        MPI_STATUS_IGNORE
    );


    for (std::size_t i = 0; i < num_nodes; ++i)
    {
        chain.nodes[i] = PPNode{
            node_xyz[3 * i + 0],
            node_xyz[3 * i + 1],
            node_xyz[3 * i + 2],
            node_end[i] != 0
        };
    }

    return chain;
}


// ------------------------------------------------------------
// Divide chains approximately evenly among ranks.
// ------------------------------------------------------------
std::vector<ChainTrajectory> ScatterChainTrajectories(
    const std::vector<ChainTrajectory>& all_chains,
    int root_rank,
    MPI_Comm comm)
{
    int rank;
    int num_ranks;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_ranks);

    std::vector<ChainTrajectory> local_chains;


    if (rank == root_rank)
    {
        const std::size_t total_chains =
            all_chains.size();

        const std::size_t base =
            total_chains / num_ranks;

        const std::size_t remainder =
            total_chains % num_ranks;

        std::size_t begin = 0;

        for (int target = 0;
            target < num_ranks;
            ++target)
        {
            const std::size_t count =
                base +
                (
                    static_cast<std::size_t>(target) < remainder
                    ? 1
                    : 0
                    );

            if (target == root_rank)
            {
                local_chains.insert(
                    local_chains.end(),
                    all_chains.begin() + begin,
                    all_chains.begin() + begin + count
                );
            }
            else
            {
                const std::uint64_t send_count =
                    static_cast<std::uint64_t>(count);

                MPI_Send(
                    &send_count,
                    1,
                    MPI_UINT64_T,
                    target,
                    CHAIN_COUNT_TAG,
                    comm
                );

                for (std::size_t i = begin;
                    i < begin + count;
                    ++i)
                {
                    SendChainTrajectory(
                        all_chains[i],
                        target,
                        comm
                    );
                }
            }

            begin += count;
        }
    }
    else
    {
        std::uint64_t recv_count = 0;

        MPI_Recv(
            &recv_count,
            1,
            MPI_UINT64_T,
            root_rank,
            CHAIN_COUNT_TAG,
            comm,
            MPI_STATUS_IGNORE
        );

        local_chains.reserve(recv_count);

        for (std::size_t i = 0;
            i < recv_count;
            ++i)
        {
            local_chains.push_back(
                ReceiveChainTrajectory(
                    root_rank,
                    comm
                )
            );
        }
    }

    return local_chains;
}


// ------------------------------------------------------------
// Sum rank-local accumulators onto root.
//
// IMPORTANT:
// local should NOT be finalized before this call.
// ------------------------------------------------------------
void ReduceSurvivalAccumulator(
    const SegmentSurvivalFunction& local,
    SegmentSurvivalFunction& global,
    int root_rank,
    MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == root_rank)
    {
        global = local;

        std::fill(
            global.survival.begin(),
            global.survival.end(),
            0.0
        );

        std::fill(
            global.sample_counts.begin(),
            global.sample_counts.end(),
            0
        );
    }

    MPI_Reduce(
        local.survival.data(),
        rank == root_rank ? global.survival.data() : nullptr,
        static_cast<int>(local.survival.size()),
        MPI_DOUBLE,
        MPI_SUM,
        root_rank,
        comm
    );

    MPI_Reduce(
        local.sample_counts.data(),
        rank == root_rank ? global.sample_counts.data() : nullptr,
        static_cast<int>(local.sample_counts.size()),
        MPI_UINT64_T,
        MPI_SUM,
        root_rank,
        comm
    );
}