cmake_minimum_required(VERSION 3.20)

project(
    z1_tube_survival
    VERSION 0.1.0
    LANGUAGES CXX
)

# ------------------------------------------------------------
# C++ language settings
# ------------------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Put all executables in build/bin.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Useful for IDEs and code-analysis tools.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# ------------------------------------------------------------
# Core tube-survival library
# ------------------------------------------------------------

add_library(tube_survival_core STATIC
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
    tube_survival_core
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

if(MSVC)
    target_compile_options(
        tube_survival_core
        PRIVATE
            /W4
    )
else()
    target_compile_options(
        tube_survival_core
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()


# ------------------------------------------------------------
# Main executable
# ------------------------------------------------------------

add_executable(
    Parallel_Segment_Survival
    core/Parallel_Segment_Survival.cpp
)

target_link_libraries(
    Parallel_Segment_Survival
    PRIVATE
        tube_survival_core
)


# ------------------------------------------------------------
# Regression / smoke tests
# ------------------------------------------------------------

include(CTest)

if(BUILD_TESTING)

    add_executable(
        affine_correction_smoke_test
        tests/Affine_Correction_Smoke_Test.cpp
    )
    target_link_libraries(
        affine_correction_smoke_test
        PRIVATE
            tube_survival_core
    )
    add_test(
        NAME affine_correction_smoke_test
        COMMAND affine_correction_smoke_test
    )

    add_executable(
        affine_time_origin_smoke_test
        tests/Affine_Time_Origin_Smoke_Test.cpp
    )
    target_link_libraries(
        affine_time_origin_smoke_test
        PRIVATE
            tube_survival_core
    )
    add_test(
        NAME affine_time_origin_smoke_test
        COMMAND affine_time_origin_smoke_test
    )

    add_executable(
        history_dependent_escape_smoke_test
        tests/History_Dependent_Escape_Smoke_Test.cpp
    )
    target_link_libraries(
        history_dependent_escape_smoke_test
        PRIVATE
            tube_survival_core
    )
    add_test(
        NAME history_dependent_escape_smoke_test
        COMMAND history_dependent_escape_smoke_test
    )

endif()

