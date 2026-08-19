#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-g++}"

INPUT_FILE="${1:-Z1+SP.dat}"
FRAME_TIME="${2:-20}"
TUBE_DIAMETER="${3:-9.3}"

BOX_CENTER_Z="${BOX_CENTER_Z:-}"
MAX_CHAINS="${MAX_CHAINS:-0}"
LAGS="${LAGS:-0,1,2,4,8,16,32,64,128,256}"
OUTPUT_DIR="${OUTPUT_DIR:-directional_history_output}"

if [[ -z "$BOX_CENTER_Z" ]]; then
    echo "ERROR: BOX_CENTER_Z is required." >&2
    echo >&2
    echo "Example:" >&2
    echo "  BOX_CENTER_Z=39.22658609485464 \\" >&2
    echo "  bash scripts/run_directional_history_comparison.sh \\" >&2
    echo "      Z1+SP.dat 20 9.3" >&2
    exit 2
fi

required=(
    include/Geometry_Utils.hpp
    include/Parse_Z1_File.hpp
    include/Primitive_Path_Trajectory.hpp
    include/Survival_IO.hpp
    include/Trajectory_Time.hpp
    include/Tube_Survival.hpp
    core/Geometry_Utils.cpp
    core/Parse_Z1_File.cpp
    core/Parser_Utils.cpp
    core/Primitive_Path_Trajectory.cpp
    core/Survival_IO.cpp
    core/Trajectory_Time.cpp
    core/Tube_Survival.cpp
    tests/Directional_History_Real_Data_Comparison.cpp
    scripts/plot_directional_history_comparison.py
)

for path in "${required[@]}"; do
    if [[ ! -f "$path" ]]; then
        echo "ERROR: run this from the project root. Missing: $path" >&2
        exit 1
    fi
done

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: input file not found: $INPUT_FILE" >&2
    exit 1
fi

mkdir -p build/bin
mkdir -p "$OUTPUT_DIR"

EXE="build/bin/directional_history_real_data_comparison"

echo "Compiler:"
"$CXX" --version | head -n 1
echo

echo "Input file      : $INPUT_FILE"
echo "Frame time      : $FRAME_TIME"
echo "Tube diameter   : $TUBE_DIAMETER"
echo "Fixed z center  : $BOX_CENTER_Z"
echo "Max chains      : $MAX_CHAINS (0 = all)"
echo "Lags            : $LAGS"
echo "Output directory: $OUTPUT_DIR"
echo

echo "Compiling real-data directional comparison..."

"$CXX" \
    -std=c++20 \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -Iinclude \
    tests/Directional_History_Real_Data_Comparison.cpp \
    core/Geometry_Utils.cpp \
    core/Parse_Z1_File.cpp \
    core/Parser_Utils.cpp \
    core/Primitive_Path_Trajectory.cpp \
    core/Survival_IO.cpp \
    core/Trajectory_Time.cpp \
    core/Tube_Survival.cpp \
    -o "$EXE"

echo
echo "Running directional comparison..."

"$EXE" \
    "$INPUT_FILE" \
    "$FRAME_TIME" \
    "$TUBE_DIAMETER" \
    "$BOX_CENTER_Z" \
    "$MAX_CHAINS" \
    "$LAGS" \
    "$OUTPUT_DIR"

echo
echo "Plotting..."

python3 \
    scripts/plot_directional_history_comparison.py \
    "$OUTPUT_DIR"

echo
echo "Done."
echo
echo "Most useful outputs:"
echo "  $OUTPUT_DIR/directional_history_tube_comparison.csv"
echo "  $OUTPUT_DIR/directional_history_comparison_dashboard.png"
echo "  $OUTPUT_DIR/directional_history_delta_phi_heatmap.png"
echo "  $OUTPUT_DIR/directional_history_largest_lag_profile.png"
