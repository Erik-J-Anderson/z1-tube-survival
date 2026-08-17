cmake_minimum_required(VERSION 3.20)

project(
    z1_tube_survival
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Put the executable in build/bin.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_library(tube_survival_core STATIC
    Geometry_Utils.cpp
    Parse_Z1_File.cpp
    Parser_Utils.cpp
    Primitive_Path_Trajectory.cpp
    Survival_IO.cpp
    Trajectory_Thinning.cpp
    Trajectory_Time.cpp
    Tube_Survival.cpp
)

target_include_directories(
    tube_survival_core
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
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

# Main program.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Parallel_Segment_Survival.cpp")
    add_executable(
        Parallel_Segment_Survival
        Parallel_Segment_Survival.cpp
    )

    target_link_libraries(
        Parallel_Segment_Survival
        PRIVATE
            tube_survival_core
    )
endif()


# ------------------------------------------------------------
# Regression / smoke tests
# ------------------------------------------------------------

include(CTest)

if(BUILD_TESTING)

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/History_Dependent_Escape_Smoke_Test.cpp")
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

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/Affine_Correction_Smoke_Test.cpp")
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
    endif()

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/Affine_Time_Origin_Smoke_Test.cpp")
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
    endif()

endif()
