# # CMakeList.txt: Blunt_Solver 的 CMake 项目，在此处包括源代码并定义
# # 项目特定的逻辑。

cmake_minimum_required(VERSION 3.8)

project("Blunt_Solver")

# 添加OpenMP支持
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "Found OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-D_OPENMPH)
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()

set(CMAKE_CXX_COMPILER g++)

# 包含目录
include_directories("./include")

set(CMAKE_CXX_STANDARD 11)

# 将源代码添加到此项目的可执行文件
add_executable(Blunt_Solver 
    src/Blunt_Solver.cpp
    src/QlQr_Solver.cpp
    src/Flux_Solver.cpp
    src/Global.cpp
    src/Geometry.cpp
    src/Spatial_Derivative.cpp
    src/Time_Integral.cpp
    src/Post_Process.cpp
    src/Time_Step.cpp
    src/Compute_Boundary.cpp
    src/Init_Flow.cpp
)

# 添加OpenMP库链接
target_link_libraries(Blunt_Solver OpenMP::OpenMP_CXX)

# 设置OpenMP线程数
# set(ENV{OMP_NUM_THREADS} 1)