cmake_minimum_required(VERSION 3.12)

project(trees VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(tree tree.cxx tree.h)


