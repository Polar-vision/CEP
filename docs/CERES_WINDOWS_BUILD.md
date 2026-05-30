# Build Ceres Solver 2.3.0-dev on Windows

This page records the Ceres build used by the CEP SuiteSparse branch. It builds Ceres from source with CPU SuiteSparse, OpenBLAS/LAPACK, and Eigen from vcpkg. CUDA/cuDSS is disabled.

## Tested Configuration

- Windows 11
- Visual Studio 2022, x64 MSVC toolchain
- CMake + Ninja
- vcpkg triplet: `x64-windows`
- Ceres Solver: 2.3.0-dev
- Ceres source revision used locally: `0ba987ac` on `master`
- Eigen: 5.0.1
- SuiteSparse: 7.12.2
- BLAS/LAPACK: OpenBLAS + LAPACK

## Directory Layout

The commands below assume this layout:

```text
C:\zuo\Projects\repos
  ceres-solver
  ceres-solver-build
  ceres-solver-install
  vcpkg
```

Adjust the paths if your workspace is different.

## 1. Prepare an x64 MSVC Environment

Run the build from a normal `cmd.exe` terminal after initializing the Visual Studio toolchain:

```bat
call "C:\zuo\Program Files\VS\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=x64
```

Check that the Windows resource tools are visible:

```bat
where cl
where rc
where mt
```

If `rc` or `mt` is missing, the terminal is not using a complete Visual Studio developer environment.

## 2. Install vcpkg Dependencies

From the vcpkg root:

```bat
cd C:\zuo\Projects\repos\vcpkg
vcpkg install eigen3:x64-windows ^
  abseil:x64-windows ^
  gtest:x64-windows ^
  glog:x64-windows ^
  gflags:x64-windows ^
  openblas:x64-windows ^
  lapack:x64-windows ^
  suitesparse:x64-windows ^
  suitesparse-cholmod[matrixops,modify,partition,supernodal]:x64-windows ^
  suitesparse-spqr:x64-windows ^
  sophus:x64-windows
```

For CEP, `sophus` is needed by the BA project. For Ceres itself, the important packages are Eigen, Abseil, GTest, OpenBLAS/LAPACK, and SuiteSparse.

This build used SuiteSparse with CHOLMOD/SPQR available. If Ceres reports that SuiteSparse is disabled because SPQR or LAPACK is missing, reinstall the missing vcpkg packages and configure Ceres again from a clean build directory.

If the Ceres compile step fails with missing identifiers such as `cholmod_scale` or `cholmod_sdmult`, CHOLMOD was built without the MatrixOps module. Install `suitesparse-cholmod[matrixops,modify,partition,supernodal]:x64-windows`, delete the Ceres build directory, and configure Ceres again.

## 3. Clone Ceres

Clone the upstream repository:

```bat
cd C:\zuo\Projects\repos
git clone https://github.com/ceres-solver/ceres-solver.git
```

If direct GitHub access is unreliable, download the source archive from GitHub or use a working mirror. Avoid committing mirror URLs into project files.

## 4. Configure Ceres

Use the Visual Studio bundled Ninja explicitly if another Ninja appears first on `PATH`.

```bat
cmake -S C:\zuo\Projects\repos\ceres-solver ^
  -B C:\zuo\Projects\repos\ceres-solver-build ^
  -G Ninja ^
  -DCMAKE_MAKE_PROGRAM="C:\zuo\Program Files\VS\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe" ^
  -DCMAKE_TOOLCHAIN_FILE=C:\zuo\Projects\repos\vcpkg\scripts\buildsystems\vcpkg.cmake ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_FLAGS=/utf-8 ^
  -DCMAKE_CXX_FLAGS=/utf-8 ^
  -DBUILD_TESTING=ON ^
  -DBUILD_EXAMPLES=ON ^
  -DUSE_CUDA=OFF ^
  -DCMAKE_INSTALL_PREFIX=C:\zuo\Projects\repos\ceres-solver-install
```

Expected CMake summary items include:

```text
Found Ceres version: 2.3.0
Found Eigen version 5.0.1
Found BLAS: ...\openblas.lib
Found LAPACK library: ...\lapack.lib;...\openblas.lib
Found SuiteSparse 7.12.2, building with SuiteSparse.
Building without CUDA.
```

The following missing components are acceptable for this CEP CPU build:

- METIS missing: disables Eigen METIS ordering only
- Intel TBB missing: SuiteSparseQR is assumed to be built without TBB
- SuiteSparse Partition missing: Ceres disables CHOLMOD partition support
- Google benchmark missing: Ceres benchmarks are not built

They are not required for `SPARSE_SCHUR` with `SUITE_SPARSE`.

## 5. Build and Install Ceres

```bat
cmake --build C:\zuo\Projects\repos\ceres-solver-build -j
cmake --install C:\zuo\Projects\repos\ceres-solver-build
```

After installation, CEP expects:

```text
C:\zuo\Projects\repos\ceres-solver-install\lib\cmake\Ceres\CeresConfig.cmake
```

## 6. Optional Test Data Fix

Some Ceres tests look for the test data relative to the build layout. If bundle adjustment tests fail because data files are missing, create a junction from the repo-level `data` path to the Ceres source data directory:

```bat
mklink /J C:\zuo\Projects\repos\data C:\zuo\Projects\repos\ceres-solver\data
```

Then rerun the tests:

```bat
ctest --test-dir C:\zuo\Projects\repos\ceres-solver-build -C Release --output-on-failure
```

This step is only for Ceres tests. CEP itself does not need this junction.

## 7. Use This Ceres Build in CEP

CEP's CMake files default to this Ceres install location:

```text
C:\zuo\Projects\repos\ceres-solver-install\lib\cmake\Ceres
```

If your Ceres install path is different, pass it when configuring CEP:

```bat
-DCeres_DIR=<your_ceres_install>\lib\cmake\Ceres
```

CEP is configured to use:

```cpp
options.linear_solver_type = ceres::SPARSE_SCHUR;
options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
```

So the Ceres configure step must find SuiteSparse and LAPACK successfully.

## Notes

- Use `/utf-8` for MSVC to avoid code page warnings from third-party headers.
- Do not configure from a terminal where `where rc` or `where mt` fails.
- If CMake accidentally uses `C:\Strawberry\c\bin\ninja.exe`, pass `-DCMAKE_MAKE_PROGRAM` to the Visual Studio Ninja path shown above.
- When changing vcpkg packages or Ceres options, delete `C:\zuo\Projects\repos\ceres-solver-build` before configuring again.
