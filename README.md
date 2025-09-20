# 🚀 CEP
**Collinearity Equation Parameterization:** A Systematic Benchmarking Study for Bundle Adjustment
**Parameterization Overview:**

🌳 **CEP**
├─ 🖼️ Image Point Parameterization
│ ├─ 🔹 UV-coordinate
│ └─ 🔹 Light cone
├─ 📦 Object Point Parameterization
│ ├─ 🔸 Cartesian coordinate
│ ├─ 🔸 Spherical coordinate
│ ├─ 🔸 Inverse distance
│ ├─ 🔸 Inverse depth (0 anchor)
│ ├─ 🔹 Cartesian coordinate (1 anchor)
│ ├─ 🔹 Spherical coordinate (1 anchor)
│ ├─ 🔹 Inverse distance (1 anchor)
│ ├─ 🔹 Inverse depth (1 anchor)
│ └─ 🔹 Parallax angle (2 anchors)
└─ 🔄 3D Rotation
├─ ⚙️ Euler angle
├─ ⚙️ Axis angle
├─ ⚙️ Quaternion
└─ ⚙️ Rotation matrix
---

## 🖥 Tested Platforms
This project has been tested on:

- **Operating System:** Windows 11  
- **Compiler / IDE:** Visual Studio 2022 (x64)  
- **CMake version:** ≥ 3.10  

> Note: Other platforms may work, but they have not been verified.


## 📦 Dependencies
- **Ceres-solver**: 2.2.0  
- **cuDSS**: 0.5.0  
- **Sophus**: 1.24.6  
- **Eigen**: 3.3.5  

---

## ⚙️ Installation

### 1️⃣ Clone the repository
```bash
git clone https://github.com/Polar-vision/CEP.git
cd CEP
```

### 2️⃣ Build the BA library (ba_v2)
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

### 3️⃣ Build the demo (example_v2)
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
