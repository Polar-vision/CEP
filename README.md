# CEP
**Collinearity Equation Parameterization:** A Systematic Benchmarking Study for Bundle Adjustment

---

## Dependencies
- **Ceres-solver**: 2.2.0  
- **cuDSS**  
- **Sophus**  
- **Eigen**: 3.3.5  

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/Polar-vision/CEP.git
cd CEP
```

### 2. Build the BA library (ba_v2)
```bash
cd ba_v2
mkdir build
cmake -B build -S . ^
  -DCMAKE_TOOLCHAIN_FILE="<vcpkg_root>\scripts\buildsystems\vcpkg.cmake" ^
  -DEIGEN3_INCLUDE_DIR="<path_to_third_party>/eigen3" ^
  -DSOPHUS_INCLUDE_DIR="<path_to_third_party>/Sophus/install/include" ^
  -DCeres_DIR="<path_to_Ceres>/lib/cmake/Ceres" ^
  -Dcudss_DIR="<path_to_cudss>/lib/12/cmake/cudss" ^
  -G "Visual Studio 17 2022" -A x64

cmake --build build --config Release
```

### 3. Build the demo (example_v2)
```bash
cd ../example_v2
mkdir build
cmake -B build -S . ^
  -DCMAKE_TOOLCHAIN_FILE="<vcpkg_root>\scripts\buildsystems\vcpkg.cmake" ^
  -DCeres_DIR="<path_to_Ceres>/lib/cmake/Ceres" ^
  -Dcudss_DIR="<path_to_cudss>/lib/12/cmake/cudss" ^
  -G "Visual Studio 17 2022" -A x64

cmake --build build --config Release
```
