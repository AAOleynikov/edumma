cmake_minimum_required(VERSION 3.5)
project(HeatSolverUpdated)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(heat_solver main.cpp)
add_executable(mesh_gen mesh_gen.cpp)

target_compile_options(heat_solver PRIVATE -O3)
target_compile_options(mesh_gen PRIVATE -O3)
