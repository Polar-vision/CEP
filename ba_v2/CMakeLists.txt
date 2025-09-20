cmake_minimum_required(VERSION 3.10)
project(ba)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Please specify EIGEN3_INCLUDE_DIR, e.g., -DEIGEN3_INCLUDE_DIR=/path/to/eigen3")
endif()


if(NOT SOPHUS_INCLUDE_DIR)
    message(FATAL_ERROR "Please specify SOPHUS_INCLUDE_DIR, e.g., -DSOPHUS_INCLUDE_DIR=/path/to/sophus/include")
endif()


if (VCPKG_INSTALLED_DIR AND VCPKG_TARGET_TRIPLET)
    list(INSERT CMAKE_PREFIX_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
    list(INSERT CMAKE_MODULE_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
endif()


find_package(Ceres REQUIRED)


set(SOURCES
    src/BAExporter_v2.cpp
    src/IBA_v2.cpp
    src/PBAImp_v2.cpp
)


set(BA_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")


add_library(ba STATIC ${SOURCES})
target_compile_definitions(ba PRIVATE BA_STATIC)
target_include_directories(ba PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${BA_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${SOPHUS_INCLUDE_DIR}
)
target_link_libraries(ba PRIVATE ${CERES_LIBRARIES})


add_library(ba_shared SHARED ${SOURCES})
target_compile_definitions(ba_shared PRIVATE BA_EXPORTS)
target_include_directories(ba_shared PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${BA_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${SOPHUS_INCLUDE_DIR}
)
target_link_libraries(ba_shared PRIVATE ${CERES_LIBRARIES})
set_target_properties(ba_shared PROPERTIES OUTPUT_NAME "ba")
add_custom_command(TARGET ba_shared POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_RUNTIME_DLLS:ba_shared> $<TARGET_FILE_DIR:ba_shared>
    COMMAND_EXPAND_LISTS
)


# cmake_minimum_required(VERSION 3.10)
# project(ba)
# set(3RDPARTY_DIR "E:/zuo/projects/backup/3rdparty")
# set(EIGEN3_INCLUDE_DIR "${3RDPARTY_DIR}/eigen3")
# set(SOPHUS_INCLUDE_DIR "${3RDPARTY_DIR}/Sophus/install/include")
# if (MSVC)
#     message(STATUS "Using MSVC compiler")
#     string(FIND "${CMAKE_CXX_FLAGS_RELEASE}" "/O2" O2_FOUND)
#     if(O2_FOUND EQUAL -1)
#         set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
#     endif()
#     message(STATUS "CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")

#     string(FIND "${CMAKE_C_FLAGS_RELEASE}" "/O2" O2_FOUND)
#     if(O2_FOUND EQUAL -1)
#         set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2")
#     endif()
#     message(STATUS "CMAKE_C_FLAGS_RELEASE = ${CMAKE_C_FLAGS_RELEASE}")
    
#     # 确保 Vcpkg 的路径优先于系统路径
#     set(VCPKG_INSTALLED_DIR "E:/zuo/vcpkg/installed") 
#     set(VCPKG_TARGET_TRIPLET "x64-windows")
#     list(INSERT CMAKE_PREFIX_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
#     list(INSERT CMAKE_MODULE_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
#     # 查找 Ceres
#     set(cudss_DIR "E:/zuo/Program Files/NVIDIA cuDSS/v0.5/lib/12/cmake/cudss")
#     set(Ceres_DIR "E:/zuo/projects/ceres-solver/install/lib/cmake/Ceres")
#     find_package(Ceres REQUIRED)
# endif()


# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # 源文件
# set(SOURCES
#     "src/BAExporter_v2.cpp"
#     "src/IBA_v2.cpp"
#     "src/PBAImp_v2.cpp"
# )

# # 头文件路径
# set(BA_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
# message(STATUS "BA_INCLUDE_DIR路径: ${BA_INCLUDE_DIR}")
# # ----------------------------
# # 1. 生成静态库（.lib）
# # ----------------------------
# add_library(ba STATIC ${SOURCES})
# target_compile_definitions(ba PRIVATE BA_STATIC)  # 静态库不需要导出符号
# target_include_directories(ba PRIVATE
#     ${CMAKE_SOURCE_DIR}
#     ${BA_INCLUDE_DIR}
#     ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include  # 添加 Vcpkg 包含路径
#     ${EIGEN3_INCLUDE_DIR}
#     ${SOPHUS_INCLUDE_DIR}
# )
# target_link_libraries(ba PRIVATE ${CERES_LIBRARIES})


# # ----------------------------
# # 2. 生成动态库（.dll + .lib 导入库）
# # ----------------------------
# add_library(ba_shared SHARED ${SOURCES})
# set_target_properties(ba_shared PROPERTIES LINK_FLAGS "/SAFESEH:NO")
# target_compile_definitions(ba_shared PRIVATE BA_EXPORTS)  # 动态库需要导出符号
# target_include_directories(ba_shared PRIVATE
#     ${CMAKE_SOURCE_DIR}
#     ${BA_INCLUDE_DIR}
#     ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include  # 添加 Vcpkg 包含路径
#     ${EIGEN3_INCLUDE_DIR}
#     ${SOPHUS_INCLUDE_DIR}
# )

# target_link_libraries(ba_shared PRIVATE ${CERES_LIBRARIES})
# # 动态库的导出文件名（可选）
# set_target_properties(ba_shared PROPERTIES OUTPUT_NAME "ba")
# add_custom_command(TARGET ba_shared POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different
#     $<TARGET_RUNTIME_DLLS:ba_shared> $<TARGET_FILE_DIR:ba_shared>
#     COMMAND_EXPAND_LISTS
# )