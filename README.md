# Z1+ Primitive-Path Tube Survival

C++ tools for analyzing tube survival from Z1+ primitive-path trajectories.

The project is under active research development. The current codebase includes primitive-path parsing and sampling, contour-resolved and tube-averaged survival calculations, affine correction for deforming simulation boxes, history-dependent permanent tube escape, and end-retraction / first-passage diagnostics.

## Current analysis definitions

For tube diameter $a$, a sampled point on the reference primitive path is compared with the future primitive path using its minimum Euclidean distance $d$.

### Instantaneous tube survival

A sample is inside the tube when

```math
d \le \frac{a}{2}.
```

This criterion is instantaneous, so a sample may leave and later re-enter.

### Affine correction

For a deforming simulation box, future primitive-path coordinates can be pulled back into the box of the time origin before distances are evaluated:

```math
\mathbf{r}_{t\rightarrow t_0}
=
\mathbf{o}(t_0)
+
H(t_0)H(t)^{-1}
\left[
\mathbf{r}(t)-\mathbf{o}(t)
\right].
```

Here $\mathbf{o}(t)$ is the box origin and $H(t)$ is the box matrix.

For the current centered, orthorhombic $z$-deformation validation tests, the per-frame origin is reconstructed from the fixed physical box center and the Z1+ box length. The manual `BOX_CENTER_Z` path is a validation mechanism, not the intended long-term production interface.

### History-dependent survival

The current history-dependent rule uses:

- $d \le a/2$: inside the inner tube;
- $a/2 < d < a$: temporarily outside, with recovery allowed;
- $d \ge a$: permanently escaped for that time origin.

Every saved intermediate frame is scanned so that an outer-boundary crossing between requested output lags is not missed.

The code also computes the pure absorbing survival

```math
S_{\mathrm{perm}}(t)
=
P\!\left[
\max_{\tau \le t} d(\tau) < a
\right],
```

which is the natural survival function for first-passage analysis.

## Repository layout

The project is currently organized as:

```text
.
├── CMakeLists.txt
├── README.md
├── .gitignore
├── .gitattributes
│
├── core/
│   ├── Geometry_Utils.cpp
│   ├── Parse_Z1_File.cpp
│   ├── Parser_Utils.cpp
│   ├── Primitive_Path_Trajectory.cpp
│   ├── Survival_IO.cpp
│   ├── Trajectory_Thinning.cpp
│   ├── Trajectory_Time.cpp
│   ├── Tube_Survival.cpp
│   ├── Parallel_Segment_Survival.cpp
│   └── additional diagnostic / smoke-test sources
│
├── include/
│   ├── Geometry_Utils.hpp
│   ├── Parse_Z1_File.hpp
│   ├── Parser_Utils.hpp
│   ├── Primitive_Path_Trajectory.hpp
│   ├── Survival_IO.hpp
│   ├── Trajectory_Thinning.hpp
│   ├── Trajectory_Time.hpp
│   └── Tube_Survival.hpp
│
├── tests/
│   ├── Affine_Correction_Smoke_Test.cpp
│   ├── Affine_Time_Origin_Smoke_Test.cpp
│   ├── History_Dependent_Escape_Smoke_Test.cpp
│   └── plotting / real-data smoke-test helpers
│
├── scripts/
│   └── run_*.sh
│
└── docs/
    └── development notes and patch documentation
```

Large simulation trajectories and generated plots/CSVs should not be committed. See `.gitignore`.

## Requirements

- C++20 compiler
- CMake 3.20 or newer
- Python 3 with Matplotlib for plotting scripts

The code uses C++20 facilities including `std::span`.

## Build with CMake

From the repository root:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON

cmake --build build -j
```

The main executable, when enabled by the current `CMakeLists.txt`, is written under:

```text
build/bin/
```

## Run regression tests

After building:

```bash
ctest --test-dir build --output-on-failure
```

The current synthetic regression suite includes tests for:

- affine correction;
- time-origin affine mapping;
- history-dependent permanent escape.

## Current all-chain history-dependent validation workflow

The current validation script can be run from the repository root with a known fixed box center:

```bash
BOX_CENTER_Z=<z-center> \
bash scripts/run_all_chain_history_escape_comparison.sh \
    Z1+SP.dat <frame-time> <tube-diameter>
```

The comparison includes:

1. raw instantaneous survival;
2. affine-corrected instantaneous survival;
3. affine-corrected history-dependent survival;
4. pure absorbing permanent-escape survival.

The manual `BOX_CENTER_Z` input is currently used for validation of centered deformation. A production implementation should obtain complete per-frame box origins from an authoritative simulation-box history rather than require the center to be entered manually.

## Data policy

Do **not** commit production Z1+, LAMMPS trajectory, or large generated analysis files to this repository.

Examples of files that should normally remain outside version control include:

```text
Z1+SP.dat
*.dump
*.data
*.log
generated *.csv
generated *.png
compiled executables
build/
```

If a regression test needs input data, add a deliberately small fixture under `tests/data/`.

## Development status

This is research software under active development. Numerical definitions and interfaces may change while the affine correction, permanent-escape model, first-passage analysis, and parallel implementation are validated.

## License

No open-source license is currently included. Choose a license deliberately before making the repository public.
