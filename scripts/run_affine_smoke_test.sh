#!/usr/bin/env bash

set -euo pipefail

# Force the modern system GCC rather than the old Conda compiler.
CXX=g++

echo "Using compiler:"
"$CXX" --version | head -n 1
echo

export MPLCONFIGDIR="${TMPDIR:-/tmp}/z1-matplotlib-${USER:-user}"
mkdir -p "${MPLCONFIGDIR}"

"$CXX" \
    -std=c++20 \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -I. \
    tests/Affine_Correction_Smoke_Test.cpp \
    Geometry_Utils.cpp \
    Primitive_Path_Trajectory.cpp \
    Tube_Survival.cpp \
    -o affine_correction_smoke_test

echo
echo "Compilation succeeded."
echo "Running affine correction smoke test..."
echo

./affine_correction_smoke_test

echo
echo "Generating plot..."
python3 tests/Plot_Affine_Smoke_Test.py

echo
echo "Smoke test and plot completed successfully."
