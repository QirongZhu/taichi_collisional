cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS "-O3")

project(mTaichi)

add_executable(Taichi src/main.cpp src/hilbert.h src/allvars.cpp src/allvars.h src/tree.h)

find_package(MPI REQUIRED)

target_link_libraries(Taichi PRIVATE MPI::MPI_CXX)