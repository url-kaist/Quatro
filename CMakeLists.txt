cmake_minimum_required(VERSION 3.10)
project(quatro)

set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

find_package(catkin REQUIRED COMPONENTS
        roscpp
        rospy
        std_msgs
        geometry_msgs
        message_generation
        sensor_msgs
        nav_msgs
        pcl_ros
        cv_bridge
)

find_package(OpenCV REQUIRED)
find_package(OpenCV REQUIRED QUIET)
find_package(OpenMP)
find_package(PCL 1.8 REQUIRED)
find_package(Boost 1.54 REQUIRED)
find_package(Eigen3 3.2 QUIET REQUIRED NO_MODULE)
find_package(catkin QUIET)

add_message_files(
        FILES
        cloud_info.msg
)

generate_messages(
        DEPENDENCIES
        geometry_msgs
        std_msgs
        nav_msgs
)

catkin_package(
        INCLUDE_DIRS
        LIBRARIES
        CATKIN_DEPENDS roscpp rospy std_msgs
)

include_directories(
        ${catkin_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
        include
)

include(3rdparty/find_dependencies.cmake)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()



####### Executable ####### 
set(teaser_src
        src/graph.cc
        src/teaser_utils/feature_matcher.cc
        src/teaser_utils/fpfh.cc
        )

add_executable(run_example examples/run_global_registration.cpp ${teaser_src})

target_link_libraries(run_example
        PUBLIC
        ${PCL_LIBRARY_DIRS}
        ${catkin_LIBRARIES}
        stdc++fs
        pmc::pmc
        )
