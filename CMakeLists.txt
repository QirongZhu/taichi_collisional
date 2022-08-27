cmake_minimum_required(VERSION 3.12)

project(trees VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


add_definitions(-D_GLIBCXX_PARALLEL)

find_package(OpenMP)

if(OPENMP_FOUND)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

add_executable(tree tree.cxx tree.h)

target_link_libraries(tree PRIVATE OpenMP::OpenMP_CXX)

target_compile_options(tree PRIVATE "${OpenMP_CXX_FLAGS}")


