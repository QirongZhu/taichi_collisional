cmake_minimum_required(VERSION 3.12)

project(trees VERSION 1.0)

# choose to use octree instead of the default binary tree
#add_compile_definitions(USE_OCTREE)

# set the bucket size to construct leaf cells  
add_compile_definitions(LEAFSIZE=96)

# multipole EXPANSION order as a macro 
add_compile_definitions(EXPANSION=10)

#add_compile_definitions(DEBUG)
add_compile_definitions(SIMD_P2P)
#add_compile_definitions(DOUBLE_P2P)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# CGAL and Boost directories
set(CGAL_DIR "${CMAKE_SOURCE_DIR}/include")
include_directories(${CGAL_DIR})

set(BOOST_DIR "/usr/local/include")
include_directories(${BOOST_DIR})

set(TBB_DIR "/usr/local/include/tbb/../oneapi/tbb/")
include_directories(${TBB_DIR})

#parallel STL and OpenMP stuff
add_definitions(-D_GLIBCXX_PARALLEL)
find_package(OpenMP)
if(OPENMP_FOUND)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()


find_package(TBB REQUIRED)
if(TBB_FOUND)
   message("TBB found")
endif()

#build the target
add_executable(tree main.cxx tree.cxx tree.h common.h kernel.cxx radixtree.cxx radixtree.h)

target_link_libraries(tree PRIVATE OpenMP::OpenMP_CXX TBB::tbb TBB::tbbmalloc)
target_compile_options(tree PRIVATE "${OpenMP_CXX_FLAGS}")


