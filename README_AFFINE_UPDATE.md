# Z1+ affine-correction update

This source bundle compiles with C++20 and adds optional, time-origin-dependent affine remapping to the segment-survival calculation.

## Box-line formats

The parser accepts whitespace-separated box records containing 3, 6, or 9 values.

### 3 values

`Lx Ly Lz`

Legacy orthorhombic format. The parser stores a diagonal box matrix. This is sufficient for affine correction when the trajectory is known to remain orthorhombic, including the current z-only oscillatory deformation.

It is **not** sufficient to reconstruct unknown shear tilt or a changing box origin.

### 6 values

`Lx Ly Lz xy xz yz`

Restricted-triclinic box with origin `(0,0,0)`.

### 9 values

`xlo ylo zlo Lx Ly Lz xy xz yz`

Restricted-triclinic box with explicit origin.

The box matrix is

```text
H = [ Lx   xy   xz ]
    [  0   Ly   yz ]
    [  0    0   Lz ]
```

## Time-origin-dependent correction

For each survival time origin `t0` and future frame `t1 = t0 + lag`, the future primitive path is pulled back into the box at `t0`:

```text
r_mapped = o0 + H0 H1^{-1} (r1 - o1)
```

This means the reference box changes with every time origin; frame 0 is not used as a universal reference.

## Box-aware segment-survival call

```cpp
SegmentSurvivalOptions options{
    .apply_affine_correction = true
};

const SegmentSurvivalFunction result =
    ComputeSegmentSurvivalFunction(
        chain_trajectory,
        primitive_path_trajectory.frame_boxes,
        lag_frames,
        tube_diameters,
        options
    );
```

The original uncorrected call is also retained:

```cpp
const SegmentSurvivalFunction result =
    ComputeSegmentSurvivalFunction(
        chain_trajectory,
        lag_frames,
        tube_diameters
    );
```

## Important first-passage note

The current segment/tube survival remains the existing instantaneous geometric observable. The next history-dependent/MFPT patch should be implemented separately so the old observable remains available as a regression baseline. See `README_FIRST_PASSAGE_NEXT.md`.
