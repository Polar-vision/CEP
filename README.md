# Building From Source

**Collinearity Equation Parameterization:** A Systematic Benchmarking Study for Bundle Adjustment

This branch is configured for a CPU build of CEP with Ceres Solver and SuiteSparse. CUDA/cuDSS is not required for the default build.

## Tested Platform

- Operating system: Windows 11
- Compiler: Visual Studio 2022, x64 MSVC toolchain
- Build system: CMake + Ninja
- Package manager: vcpkg, `x64-windows`

## Dependencies

The current local configuration has been tested with:

- Ceres Solver: 2.3.0-dev, built from source
- SuiteSparse: 7.12.2, from vcpkg
- BLAS/LAPACK: OpenBLAS + LAPACK, from vcpkg
- Eigen: 5.0.1, from vcpkg
- Sophus: 1.24.6, from vcpkg

See [Build Ceres Solver 2.3.0-dev on Windows](docs/CERES_WINDOWS_BUILD.md) for the Ceres source build used by this branch.

Install the vcpkg packages used by CEP and the Ceres build:

```bat
cd C:\zuo\Projects\repos\vcpkg
vcpkg install eigen3:x64-windows abseil:x64-windows glog:x64-windows gflags:x64-windows openblas:x64-windows lapack:x64-windows suitesparse:x64-windows suitesparse-spqr:x64-windows sophus:x64-windows
```

Ceres should be built and installed before building CEP. The CMake files default to:

```text
C:\zuo\Projects\repos\ceres-solver-install\lib\cmake\Ceres
```

If Ceres is installed somewhere else, pass `-DCeres_DIR=<ceres_install>\lib\cmake\Ceres` when configuring `ba_v2` and `example_v2`.

## Build

Open an x64 Visual Studio developer environment before running CMake. In a normal `cmd.exe` terminal:

```bat
call "C:\zuo\Program Files\VS\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=x64
```

If another Ninja appears first on `PATH`, pass `-DCMAKE_MAKE_PROGRAM` as shown below.

### Build `ba_v2`

```bat
cmake -S C:\zuo\Projects\repos\CEP\ba_v2 ^
  -B C:\zuo\Projects\repos\CEP\ba_v2\build ^
  -G Ninja ^
  -DCMAKE_MAKE_PROGRAM="C:\zuo\Program Files\VS\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe" ^
  -DCMAKE_TOOLCHAIN_FILE=C:\zuo\Projects\repos\vcpkg\scripts\buildsystems\vcpkg.cmake ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_FLAGS=/utf-8 ^
  -DCMAKE_CXX_FLAGS=/utf-8

cmake --build C:\zuo\Projects\repos\CEP\ba_v2\build -j
```

The build produces `ba.dll`, `ba.lib`, and `ba_static.lib` in `ba_v2\build`.

### Build `example_v2`

```bat
cmake -S C:\zuo\Projects\repos\CEP\example_v2 ^
  -B C:\zuo\Projects\repos\CEP\example_v2\build ^
  -G Ninja ^
  -DCMAKE_MAKE_PROGRAM="C:\zuo\Program Files\VS\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe" ^
  -DCMAKE_TOOLCHAIN_FILE=C:\zuo\Projects\repos\vcpkg\scripts\buildsystems\vcpkg.cmake ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_FLAGS=/utf-8 ^
  -DCMAKE_CXX_FLAGS=/utf-8

cmake --build C:\zuo\Projects\repos\CEP\example_v2\build -j
```

The example build copies `ba.dll` and the required runtime DLLs next to `example.exe`.

Running `example.exe` evaluates the nine object-point parameterizations configured in `example_v2/example.cpp` and writes comparison tables next to the selected dataset:

- `BA-comparison.csv`: full numeric summary for post-processing
- `BA-comparison.md`: compact Markdown table for quick inspection

## Solver Backend

`ba_v2/src/PBAImp_v2.cpp` uses:

```cpp
options.linear_solver_type = ceres::SPARSE_SCHUR;
options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
```

This matches the Ceres build with SuiteSparse and avoids the previous CUDA Sparse/cuDSS requirement.

## Synthetic Tests

See [Synthetic Parallax Short-Baseline Test](docs/SYNTHETIC_PARALLAX_TEST.md) for a fixed-camera synthetic check of the CEP parallax parameterization under weak-baseline geometry.

## Datasets

This repository benchmarks bundle adjustment datasets under the CEP framework. Download links:

| Dataset | Description | Download Link |
| --- | --- | --- |
| Close-Range (CR) | Small-scale indoor/outdoor scenes | [Download CR](https://drive.google.com/drive/folders/1mvsQEFGBvZ-VcfxIJ3hqXiBV3tVTPMSe) |
| UAV (UAV) | UAV imagery for oblique photogrammetry | [Download UAV](https://drive.google.com/drive/folders/1-VA-JrVe03PVZnswjuLAvx_EC7vqzSXN) |
| Vehicle (KD) | Vehicle-mounted multi-camera data | [Download KD](https://drive.google.com/drive/folders/1_GID2a5O5CSfUn5QfoWFhhnA_NwI5kmi) |
| Multi-camera Oblique (LM) | Multi-camera oblique aerial datasets | [Download LM](https://drive.google.com/drive/folders/1NDOMrSZocyTG7JEdLQ7KdujHUKDikmQa) |

For the current Windows workspace, place datasets here:

```text
C:\zuo\Projects\BA Datasets
```

`3rdparty/dataPath.h` uses paths relative to the directory containing `example.exe`. For the default build this is `example_v2\build`, so the example can be launched from any working directory. For example:

```cpp
static const char* dt = "../../../../BA Datasets/Close-Range/CR1-problem-11-9611/Initial Value/cal.txt";
```

To switch datasets, edit `3rdparty/dataPath.h` and enable one `dt` line.

## Data Format

Each dataset folder contains:

- `cal.txt`: camera intrinsics, including `fx`, `fy`, `cx`, and `cy`
- `Cam.txt`: camera poses, including Euler angles `ey`, `ex`, `ez`, camera center `(Xc, Yc, Zc)`, and camera ID
- `XYZ.txt`: object point coordinates `X`, `Y`, `Z`
- `Feature.txt`: feature tracks, image indices, and observed image coordinates `(u, v)`
