cmake_minimum_required(VERSION 3.20)

project(
    z1_tube_survival
    LANGUAGES CXX
)

# ------------------------------------------------------------
# C++ standard
# ------------------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


# ------------------------------------------------------------
# MPI
# ------------------------------------------------------------

find_package(MPI REQUIRED)


# ============================================================
# Main serial survival executable
# ============================================================

add_executable(
    Parallel_Segment_Survival

    core/Stephanou_Survival.cpp

    core/Input_Parsers.cpp
    core/Survival_Accumulator.cpp

    core/Geometry_Utils.cpp
    core/Parse_Z1_File.cpp
    core/Parser_Utils.cpp
    core/Primitive_Path_Trajectory.cpp
    core/Survival_IO.cpp
    core/Trajectory_Thinning.cpp
    core/Trajectory_Time.cpp
    core/Tube_Survival.cpp
)

target_include_directories(
    Parallel_Segment_Survival
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)





# ============================================================
# MPI communications smoke test
# ============================================================

add_executable(
    MPI_Comms_Smoke_Test

    tests/MPI_Comms_Smoke_Test.cpp
    core/MPI_Comms.cpp
)

target_include_directories(
    MPI_Comms_Smoke_Test
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    MPI_Comms_Smoke_Test
    PRIVATE
        MPI::MPI_CXX
)


add_executable(
    Parallel_Segment_Survival_MPI

    core/Stephanou_Survival_parallel.cpp

    core/MPI_Comms.cpp
    core/Input_Parsers.cpp
    core/Survival_Accumulator.cpp

    core/Geometry_Utils.cpp
    core/Parse_Z1_File.cpp
    core/Parser_Utils.cpp
    core/Primitive_Path_Trajectory.cpp
    core/Survival_IO.cpp
    core/Trajectory_Time.cpp
    core/Tube_Survival.cpp
)

target_include_directories(
    Parallel_Segment_Survival_MPI
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    Parallel_Segment_Survival_MPI
    PRIVATE
        MPI::MPI_CXX
)

# ------------------------------------------------------------
# Compiler warnings
# ------------------------------------------------------------




if(MSVC)

    target_compile_options(
        Parallel_Segment_Survival
        PRIVATE
            /W4
    )

    target_compile_options(
        MPI_Comms_Smoke_Test
        PRIVATE
            /W4
    )

else()

    target_compile_options(
        Parallel_Segment_Survival
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )

    target_compile_options(
        MPI_Comms_Smoke_Test
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )

endif()
