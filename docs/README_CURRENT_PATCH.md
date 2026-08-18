# Current affine-correction patch

This bundle contains the complete current C++ source tree plus the parser and affine smoke tests used in the August 17 update.

## What changed

1. `Primitive_Path_Trajectory.hpp`
   - `PrimitivePathTrajectory` stores one global `Box` per frame plus all chain trajectories.
   - `Box` uses `origin + matrix` so the same representation can support orthorhombic and triclinic boxes.

2. `Parse_Z1_File.cpp` / `Parser_Utils.cpp`
   - One box is retained for every parsed Z1+ frame.
   - Legacy `Lx Ly Lz` records become diagonal box matrices.

3. `Geometry_Utils.hpp/.cpp`
   - Provides source-box -> destination-box mapping.
   - `MapPPNodesBetweenBoxes()` works with the current `PPNode{x,y,z,end}` layout.

4. `Tube_Survival.hpp/.cpp`
   - The old three-argument `ComputeSegmentSurvivalFunction()` remains available and retains uncorrected behavior.
   - The box-aware overload optionally removes homogeneous affine deformation.
   - For every pair `(t0, t0 + lag)`, the future PP is mapped into the **box at that time origin**:

     r_mapped = o(t0) + H(t0) H(t0+lag)^(-1) [r_future - o(t0+lag)]

   - Legacy three-value boxes are accepted for affine correction when the trajectory is known to remain orthorhombic (the current z-only deformation case). They do not contain enough information for an unknown sheared/tilted box.

5. Existing end-retraction code is otherwise unchanged in this patch.

## Synthetic affine time-origin smoke test

```bash
g++ -std=c++20 -O2 -Wall -Wextra -pedantic -I. \
    tests/Affine_Time_Origin_Smoke_Test.cpp \
    Tube_Survival.cpp \
    Primitive_Path_Trajectory.cpp \
    Geometry_Utils.cpp \
    -o Affine_Time_Origin_Smoke_Test

./Affine_Time_Origin_Smoke_Test
```

Expected:

```text
uncorrected lag-1 tube survival = 0
corrected lag-1 tube survival   = 1
lag-1 time origins              = 2
ALL TIME-ORIGIN AFFINE TESTS PASSED
```

The test deliberately uses two different lag-1 time origins, so a wrong implementation that maps every frame to frame 0 will fail.

## Real deforming-box parser smoke test

```bash
g++ -std=c++20 -O2 -Wall -Wextra -pedantic \
    Serial_Smoke_Test_Deforming_Box.cpp \
    Parse_Z1_File.cpp \
    Parser_Utils.cpp \
    -o Serial_Smoke_Test_Deforming_Box

./Serial_Smoke_Test_Deforming_Box \
    Z1+SP.dat \
    76.1024 \
    7.6 \
    0.2
```

The user's A=7.6 trajectory passed with 770 chains, 2000 frames, fixed Lx/Ly, recovered Lz amplitude 7.59959, and center 76.1024.
