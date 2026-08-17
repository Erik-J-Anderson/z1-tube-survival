# Z1+ Primitive-Path Tube Survival

C++ tools for analyzing tube survival from Z1+ primitive-path trajectories.

The project is currently under active research development. The present codebase includes
primitive-path parsing and sampling, contour-resolved and tube-averaged survival calculations,
affine correction for deforming simulation boxes, history-dependent permanent tube escape,
and end-retraction / first-passage diagnostics.

## Current analysis definitions

For tube diameter \(a\), a sampled point on the reference primitive path is compared with the
future primitive path using its minimum Euclidean distance \(d\).

### Instantaneous tube survival

A sample is inside the tube when

\[
d \le a/2.
\]

This criterion is instantaneous, so a sample may leave and later re-enter.

### Affine correction

For a deforming box, future primitive-path coordinates can be pulled back into the box of the
time origin before distances are evaluated:

\[
\mathbf r_{t\rightarrow t_0}
=
\mathbf o(t_0)
+
H(t_0)H(t)^{-1}
\left[\mathbf r(t)-\mathbf o(t)\right].
\]

Here \(\mathbf o(t)\) is the box origin and \(H(t)\) is the box matrix.

For the current centered, orthorhombic \(z\)-deformation smoke tests, the per-frame origin is
reconstructed from the fixed physical box center and the Z1+ box length. This manual
`BOX_CENTER_Z` path is a validation mechanism, not the intended long-term production interface.

### History-dependent survival

The current history-dependent rule uses:

- `d <= a/2`: inside the inner tube;
- `a/2 < d < a`: temporarily outside, with recovery allowed;
- `d >= a`: permanently escaped for that time origin.

Every saved intermediate frame is scanned so an outer-boundary crossing between requested
output lags is not missed.

The code also computes the pure absorbing survival

\[
S_{\mathrm{perm}}(t)
=
P\!\left[\max_{\tau\le t} d(\tau) < a\right],
\]

which is the natural survival function for first-passage analysis.

## Repository layout

The repository is intentionally kept close to the current working project layout while the
scientific implementation is still evolving.

```text
.
├── Geometry_Utils.cpp/.hpp
├── Parse_Z1_File.cpp/.hpp
├── Parser_Utils.cpp/.hpp
├── Primitive_Path_Trajectory.cpp/.hpp
├── Survival_IO.cpp/.hpp
├── Trajectory_Thinning.cpp/.hpp
├── Trajectory_Time.cpp/.hpp
├── Tube_Survival.cpp/.hpp
├── Parallel_Segment_Survival.cpp
├── tests/
└── run_*.sh
```

Large simulation trajectories and generated plots/CSVs should not be committed. See
`.gitignore`.

## Requirements

- C++20 compiler
- CMake 3.20 or newer for the CMake build
- Python 3 + Matplotlib for plotting scripts

GCC 12+ or Clang 15+ is recommended. The code uses C++20 facilities including `std::span`.

## Build with CMake

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The main executable, when present, will be written to:

```text
build/bin/Parallel_Segment_Survival
```

To compile and run the synthetic regression tests:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Current history-dependent all-chain smoke workflow

The current validation script can be run from the project root with a known fixed box center:

```bash
BOX_CENTER_Z=<z-center> \
bash run_all_chain_history_escape_comparison.sh Z1+SP.dat <frame-time> <tube-diameter>
```

Example output includes tube-survival CSVs and contour/time-domain plots comparing:

1. raw instantaneous survival;
2. affine-corrected instantaneous survival;
3. affine-corrected history-dependent survival;
4. pure absorbing permanent-escape survival.

## Data policy

Do **not** commit production Z1+, LAMMPS trajectory, or large analysis-output files to this
repository. The `.gitignore` excludes common large inputs and generated results.

If a regression test needs input data, add a deliberately small fixture under `tests/data/`.
The ignore rules explicitly permit small test fixtures there.

## Development status

This is research software under active development. Numerical definitions and interfaces may
change while the affine correction, permanent-escape model, first-passage analysis, and
parallel implementation are validated.

A future production cleanup should remove the manual fixed-center input by reading complete
per-frame box-origin information from the corresponding LAMMPS trajectory or another
authoritative box-history source.

## License

No open-source license is included yet. If the repository is private, that is fine. Choose a
license deliberately before making the repository public.
