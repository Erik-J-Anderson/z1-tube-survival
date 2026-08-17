#!/usr/bin/env bash
set -euo pipefail

CXX=g++

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage: $0 <Z1+SP.dat> [chain_index=0] [tube_diameter=5.0]"
    exit 2
fi

INPUT=$1
CHAIN_INDEX=${2:-0}
TUBE_DIAMETER=${3:-5.0}

echo "Using compiler:"
"$CXX" --version | head -n 1

echo
echo "Compiling real deforming-trajectory affine smoke test..."
"$CXX" \
    -std=c++20 \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -I. \
    tests/Serial_Smoke_Test_Real_Affine.cpp \
    Parse_Z1_File.cpp \
    Parser_Utils.cpp \
    Primitive_Path_Trajectory.cpp \
    Trajectory_Time.cpp \
    Geometry_Utils.cpp \
    Tube_Survival.cpp \
    -o Serial_Smoke_Test_Real_Affine

echo
echo "Running on real trajectory..."
./Serial_Smoke_Test_Real_Affine \
    "$INPUT" \
    "$CHAIN_INDEX" \
    "$TUBE_DIAMETER"
