# CEP
Collinearity Equation Parameterization: A Systematic Benchmarking Study for Bundle Adjustment

# Dependencies
Ceres, Cudss, Sophus, Eigen

# Installation
git clone https://github.com/Polar-vision/CEP.git \\
cd CEP \\
cd ba_v2 \\
mkdir build \\
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 \\
cmake --build build --config Release -j8 \\

cd .. \\
cd example_v2 \\
mkdir build \\
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 \\
cmake --build build --config Release -j8
