#!/usr/bin/env bash
set -euo pipefail

cxx="${CXX:-g++}"

"$cxx" -std=c++20 -O2 -Wall -Wextra -pedantic -I. \
    tests/Affine_Time_Origin_Smoke_Test.cpp \
    Tube_Survival.cpp \
    Primitive_Path_Trajectory.cpp \
    Geometry_Utils.cpp \
    -o Affine_Time_Origin_Smoke_Test

./Affine_Time_Origin_Smoke_Test
