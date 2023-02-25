cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)

project(mTaichi)

add_executable(Taichi src/main.cpp src/hilbert.h src/allvars.h)

find_package(MPI REQUIRED)

target_link_libraries(Taichi PUBLIC MPI::MPI_CXX)