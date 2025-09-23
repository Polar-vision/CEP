# 🚀 Building from Source
**Collinearity Equation Parameterization:** A Systematic Benchmarking Study for Bundle Adjustment  

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
# 🚀 Datasets

This repository provides a set of datasets used for benchmarking **Bundle Adjustment** under the Collinearity Equation Parameterization (CEP) framework.  

---

## 📥 Download

You can download the datasets from the following links:

| Dataset | Description | Download Link |
|---------|-------------|---------------|
| Close-Range (CR) | Small-scale indoor/outdoor scenes | [Download CR](https://drive.google.com/drive/folders/1mvsQEFGBvZ-VcfxIJ3hqXiBV3tVTPMSe) |
| UAV (UAV) | UAV imagery for oblique photogrammetry | [Download UAV](https://drive.google.com/drive/folders/1-VA-JrVe03PVZnswjuLAvx_EC7vqzSXN) |
| Vehicle (KD) | Vehicle-mounted multi-camera data | [Download KD](https://drive.google.com/drive/folders/1_GID2a5O5CSfUn5QfoWFhhnA_NwI5kmi) |
| Multi-camera Oblique (LM) | Multi-camera oblique aerial datasets | [Download LM](https://drive.google.com/drive/folders/1NDOMrSZocyTG7JEdLQ7KdujHUKDikmQa) |
---
## 🗂 Data Format

Each dataset is organized in the following structure:

- **Intrinsics (`cal.txt`)**  
  Contains the camera intrinsic parameters: `fx`, `fy`, `cx`, `cy` for each camera.

- **Extrinsics (`Cam.txt`)**  
  Contains the camera poses, including Euler angles (`ey`, `ex`, `ez`), the perspective center `(Xc, Yc, Zc)`, and the camera ID.  
  - `ey` : rotation around the y-axis  
  - `ex` : rotation around the x-axis  
  - `ez` : rotation around the z-axis

- **3D Points (`XYZ.txt`)**  
  Lists the 3D coordinates of the object points: `X`, `Y`, `Z`.

- **Feature Tracks (`Feature.txt`)**  
  Each line represents a feature track, including the number of views, the corresponding image indices, and the (u, v) coordinates in each image.


---
