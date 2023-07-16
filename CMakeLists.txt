cmake_minimum_required(VERSION 3.13)
project(quatro)

### set compiler
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++17 -pthread -fexceptions -g -ggdb") #c++17, for parallel programming 
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
  message("ERROR: OpenMP could not be found.")
endif(OPENMP_FOUND)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra ${CMAKE_CXX_FLAGS}")

### get packages
find_package(PCL REQUIRED)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(catkin REQUIRED)
catkin_package( #this automatically adds/installs the results of build
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
)

### get packages' headers
include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  ${PCL_INCLUDE_DIRS}
)

include(3rdparty/find_dependencies.cmake)

###########
## Build ##
###########
### teaser
set(teaser_src
  src/graph.cc
  src/teaser_utils/feature_matcher.cc
  src/teaser_utils/fpfh.cc
)

target_link_libraries(run_example
  PUBLIC
  ${PCL_LIBRARY_DIRS}
  ${catkin_LIBRARIES}
  stdc++fs
  pmc::pmc
)