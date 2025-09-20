cmake_minimum_required(VERSION 3.10)
project(example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:CONSOLE")

find_package(Ceres REQUIRED)

set(BA_PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ba_v2")
set(BA_INCLUDE_DIR "${BA_PROJECT_DIR}/include")
set(BA_BUILD_DIR   "${BA_PROJECT_DIR}/build/$<CONFIG>")  # Debug/Release 自动切换

set(SOURCES
    example.cpp
)

add_executable(example ${SOURCES})
set_target_properties(example PROPERTIES LINK_FLAGS "/SAFESEH:NO")

target_include_directories(example PRIVATE
    ${BA_INCLUDE_DIR}
    ${BA_PROJECT_DIR}/../3rdparty   # 假设 3rdparty 在 workspace 下
)

target_link_libraries(example PRIVATE
    "${BA_BUILD_DIR}/ba.lib"
    ${CERES_LIBRARIES}
)

add_custom_command(TARGET example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${BA_BUILD_DIR}/ba.dll"
    $<TARGET_FILE_DIR:example>
)

add_custom_command(TARGET example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_RUNTIME_DLLS:example> $<TARGET_FILE_DIR:example>
    COMMAND_EXPAND_LISTS
)



# cmake_minimum_required(VERSION 3.10)
# project(example)
# set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:CONSOLE")
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(3RDPARTY_DIR "E:/zuo/projects/backup/3rdparty")
# set(BA_LIBS "E:/zuo/projects/backup/ba_v2/build/Release")

# # 确保 Vcpkg 的路径优先于系统路径
# set(VCPKG_INSTALLED_DIR "E:/zuo/vcpkg/installed") 
# set(VCPKG_TARGET_TRIPLET "x64-windows")
# list(INSERT CMAKE_PREFIX_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
# list(INSERT CMAKE_MODULE_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
# # 查找 Ceres
# set(cudss_DIR "E:/zuo/Program Files/NVIDIA cuDSS/v0.5/lib/12/cmake/cudss")
# set(Ceres_DIR "E:/zuo/projects/ceres-solver/install/lib/cmake/Ceres")
# find_package(Ceres REQUIRED)


# # 源文件
# set(SOURCES
#     example.cpp
# )

# # 创建可执行程序
# add_executable(example ${SOURCES})
# set_target_properties(example PROPERTIES LINK_FLAGS "/SAFESEH:NO")
# set(BA_INCLUDE_DIR "E:/zuo/projects/backup/ba_v2/include")

# # 头文件路径
# target_include_directories(example PRIVATE
#     ${BA_INCLUDE_DIR}
#     ${3RDPARTY_DIR}
#     ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include  # 添加 Vcpkg 包含路径
# )

# target_link_libraries(example PRIVATE 
#     "${BA_LIBS}/ba.lib"
#     ${CERES_LIBRARIES}
# )

# add_custom_command(TARGET example POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different
#     "${BA_LIBS}/ba.dll"
#     $<TARGET_FILE_DIR:example>)

# add_custom_command(TARGET example POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different
#     $<TARGET_RUNTIME_DLLS:example> $<TARGET_FILE_DIR:example>
#     COMMAND_EXPAND_LISTS
# )