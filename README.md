# CEP
Collinearity Equation Parameterization: A Systematic Benchmarking Study for Bundle Adjustment

# Dependencies
Ceres-solver: 2.2.0
Cudss
Sophus
Eigen: 3.3.5

# Installation
```bash
git clone https://github.com/Polar-vision/CEP.git  
cd CEP  
cd ba_v2  
mkdir build  
cmake -B build -S .  
-DCMAKE_TOOLCHAIN_FILE="<VCPKG_ROOT>\scripts\buildsystems\vcpkg.cmake"  
-DEIGEN3_INCLUDE_DIR="<PATH_TO_THIRD_PARTY>/eigen3"  
-DSOPHUS_INCLUDE_DIR="<PATH_TO_THIRD_PARTY>/Sophus/install/include"  
-DCeres_DIR="<PATH_TO_CERES>/lib/cmake/Ceres"  
-Dcudss_DIR="<PATH_TO_CUDSS>/lib/12/cmake/cudss"  
-G "Visual Studio 17 2022" -A x64  

cmake --build build --config Release -j8  

cd ..  
cd example_v2  
mkdir build  
cmake -B build -S .  
-DCMAKE_TOOLCHAIN_FILE="<VCPKG_ROOT>\scripts\buildsystems\vcpkg.cmake"  
-DCeres_DIR="<PATH_TO_CERES>/lib/cmake/Ceres"  
-Dcudss_DIR="<PATH_TO_CUDSS>/lib/12/cmake/cudss"  
-G "Visual Studio 17 2022" -A x64  

cmake --build build --config Release -j8
