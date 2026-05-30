cmake_minimum_required(VERSION 3.10)
project(example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:CONSOLE")

if (VCPKG_INSTALLED_DIR AND VCPKG_TARGET_TRIPLET)
    list(INSERT CMAKE_PREFIX_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
    list(INSERT CMAKE_MODULE_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
endif()

set(DEFAULT_CERES_DIR "C:/zuo/Projects/repos/ceres-solver-install/lib/cmake/Ceres")
if(NOT Ceres_DIR AND EXISTS "${DEFAULT_CERES_DIR}/CeresConfig.cmake")
    set(Ceres_DIR "${DEFAULT_CERES_DIR}" CACHE PATH "Path to CeresConfig.cmake" FORCE)
endif()

find_package(Ceres REQUIRED)

set(BA_PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ba_v2")
set(BA_INCLUDE_DIR "${BA_PROJECT_DIR}/include")
if(CMAKE_CONFIGURATION_TYPES)
    set(BA_LIB_DIR "${BA_PROJECT_DIR}/build/$<CONFIG>")
else()
    set(BA_LIB_DIR "${BA_PROJECT_DIR}/build")
endif()

set(COPY_BA_RUNTIME_DLLS_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/copy_ba_runtime_dlls.cmake")
file(WRITE "${COPY_BA_RUNTIME_DLLS_SCRIPT}" [=[
if(NOT EXISTS "${src_dir}")
    message(FATAL_ERROR "BA runtime directory does not exist: ${src_dir}")
endif()

file(GLOB runtime_dlls "${src_dir}/*.dll")
foreach(runtime_dll IN LISTS runtime_dlls)
    file(COPY "${runtime_dll}" DESTINATION "${dst_dir}")
endforeach()
]=])

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
    "${BA_LIB_DIR}/ba.lib"
    Ceres::ceres
)

add_custom_command(TARGET example POST_BUILD
    COMMAND ${CMAKE_COMMAND}
    "-Dsrc_dir=${BA_LIB_DIR}"
    "-Ddst_dir=$<TARGET_FILE_DIR:example>"
    -P "${COPY_BA_RUNTIME_DLLS_SCRIPT}"
    VERBATIM
)

add_custom_command(TARGET example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_RUNTIME_DLLS:example> $<TARGET_FILE_DIR:example>
    COMMAND_EXPAND_LISTS
)

add_executable(synthetic_parallax_test synthetic_parallax_test.cpp)
set_target_properties(synthetic_parallax_test PROPERTIES LINK_FLAGS "/SAFESEH:NO")
if(MSVC)
    target_compile_options(synthetic_parallax_test PRIVATE /EHsc)
endif()
target_include_directories(synthetic_parallax_test PRIVATE
    ${BA_INCLUDE_DIR}
    ${BA_PROJECT_DIR}/src
    ${BA_PROJECT_DIR}/../3rdparty
)
target_link_libraries(synthetic_parallax_test PRIVATE
    "${BA_LIB_DIR}/ba.lib"
    Ceres::ceres
)
add_custom_command(TARGET synthetic_parallax_test POST_BUILD
    COMMAND ${CMAKE_COMMAND}
    "-Dsrc_dir=${BA_LIB_DIR}"
    "-Ddst_dir=$<TARGET_FILE_DIR:synthetic_parallax_test>"
    -P "${COPY_BA_RUNTIME_DLLS_SCRIPT}"
    VERBATIM
)
add_custom_command(TARGET synthetic_parallax_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_RUNTIME_DLLS:synthetic_parallax_test> $<TARGET_FILE_DIR:synthetic_parallax_test>
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
