#ifndef _UTILITY_LIDAR_ODOMETRY_H_
#define _UTILITY_LIDAR_ODOMETRY_H_


#include <ros/ros.h>
#include <ros/package.h>

#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>

#include "quatro/cloud_info.h"

#include <opencv/cv.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/range_image/range_image.h>
#include <pcl/filters/filter.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/common.h>
#include <pcl/registration/icp.h>

#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>
 
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <deque>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include <iterator>
#include <sstream>
#include <string>
#include <limits>
#include <iomanip>
#include <array>
#include <thread>
#include <mutex>
#include <experimental/filesystem>

using namespace std;
using PointType = pcl::PointXYZ;

typedef pcl::PointXYZI  PointTypeIP;

//Velodyne 64 HDE
// extern const int N_SCAN = 64;
// extern const int Horizon_SCAN = 1800; //1028~4500
// extern const float ang_res_x = 360.0/float(Horizon_SCAN);
// extern const float ang_res_y = 26.9/float(N_SCAN-1);//28.0/float(N_SCAN-1);
// extern const float ang_bottom = 25.0;
// extern const int groundScanInd = 60;    // 60 ;

// VLP-16
// extern const int N_SCAN = 16;
// extern const int Horizon_SCAN = 1800;
// extern const float ang_res_x = 0.2;
// extern const float ang_res_y = 2.0;
// extern const float ang_bottom = 15.0+0.1;
// extern const int groundScanInd = 7;

// HDL-32E
// extern const int N_SCAN = 32;
// extern const int Horizon_SCAN = 1800;
// extern const float ang_res_x = 360.0/float(Horizon_SCAN);
// extern const float ang_res_y = 41.33/float(N_SCAN-1);
// extern const float ang_bottom = 30.67;
// extern const int groundScanInd = 20;

// Ouster users may need to uncomment line 159 in imageProjection.cpp
// Usage of Ouster imu data is not fully supported yet, please just publish point cloud data
// Ouster OS1-16
// extern const int N_SCAN = 16;
// extern const int Horizon_SCAN = 1024;
// extern const float ang_res_x = 360.0/float(Horizon_SCAN);
// extern const float ang_res_y = 33.2/float(N_SCAN-1);
// extern const float ang_bottom = 16.6+0.1;
// extern const int groundScanInd = 7;

// Ouster OS1-64
//  extern const int N_SCAN = 64;
//  extern const int Horizon_SCAN = 1024;
//  extern const float ang_res_x = 360.0/float(Horizon_SCAN);
//  extern const float ang_res_y = 33.2/float(N_SCAN-1);
//  extern const float ang_bottom = 16.6+0.1;
//  extern const int groundScanInd = 15;

extern const bool loopClosureEnableFlag = true;
extern const double mappingProcessInterval = 0.3;

extern const float scanPeriod = 0.1;
extern const int systemDelay = 0;
extern const int imuQueLength = 200;

// extern const float sensorMountAngle = 0.0;

// //segmentation threshold
// extern const float segmentTheta = 60.0/180.0*M_PI; // decrese this value may improve accuracy       //60.0/180.0*M_PI

// //If number of segment is below than 30, check line number. this for minimum number of point for it
// extern const int segmentValidPointNum = 5;

// //if number of segment is small, number of line is checked, this is threshold for it.
// extern const int segmentValidLineNum = 3;

// extern const float segmentAlphaX = ang_res_x / 180.0 * M_PI;
// extern const float segmentAlphaY = ang_res_y / 180.0 * M_PI;


//extern const int edgeFeatureNum = 2;
//extern const int surfFeatureNum = 4;
//extern const int sectionsTotal = 6;
//extern const float edgeThreshold = 0.1;
//extern const float surfThreshold = 0.1;
//extern const float nearestFeatureSearchSqDist = 25;
//
//
//// Mapping Params
//extern const float surroundingKeyframeSearchRadius = 50.0; // key frame that is within n meters from current pose will be considerd for scan-to-map optimization (when loop closure disabled)
//extern const int   surroundingKeyframeSearchNum = 50; // submap size (when loop closure enabled)
//// history key frames (history submap for loop closure)
//extern const float historyKeyframeSearchRadius = 7.0; // key frame that is within n meters from current pose will be considerd for loop closure
//extern const int   historyKeyframeSearchNum = 25; // 2n+1 number of hostory key frames will be fused into a submap for loop closure
//extern const float historyKeyframeFitnessScore = 0.3; // the smaller the better alignment
//
//extern const float globalMapVisualizationSearchRadius = 500.0; // key frames with in n meters will be visualized


struct smoothness_t{ 
    float value;
    size_t ind;
};

struct by_value{ 
    bool operator()(smoothness_t const &left, smoothness_t const &right) { 
        return left.value < right.value;
    }
};

/*
    * A point cloud type that has 6D pose info ([x,y,z,roll,pitch,yaw] intensity is time stamp)
    */
struct PointXYZIRPYT
{
    PCL_ADD_POINT4D
    PCL_ADD_INTENSITY;
    float roll;
    float pitch;
    float yaw;
    double time;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} EIGEN_ALIGN16;

POINT_CLOUD_REGISTER_POINT_STRUCT (PointXYZIRPYT,
                                   (float, x, x) (float, y, y)
                                   (float, z, z) (float, intensity, intensity)
                                   (float, roll, roll) (float, pitch, pitch) (float, yaw, yaw)
                                   (double, time, time)
)

typedef PointXYZIRPYT  PointTypePose;

void setCorrespondenceMarker(const pcl::PointCloud<PointType>& src_matched, const pcl::PointCloud<PointType>& tgt_matched,
                             visualization_msgs::Marker& marker, float thickness=0.1, std::vector<float> rgb_color={0.0, 0.0, 0.0}, int id=0){
    if (!marker.points.empty()) marker.points.clear();
    marker.header.frame_id = "map";
    marker.header.stamp = ros::Time();
    marker.ns = "my_namespace";
    marker.id = id; // To avoid overlap
    marker.type = visualization_msgs::Marker::LINE_LIST;
    marker.action = visualization_msgs::Marker::ADD;

    marker.scale.x = thickness; // thickness
    marker.color.r = rgb_color[0];
    marker.color.g = rgb_color[1];
    marker.color.b = rgb_color[2];
    marker.color.a = 1.0; // Don't forget to set the alpha!

    geometry_msgs::Point srcP;
    geometry_msgs::Point tgtP;
    assert(src_matched.size() == tgt_matched.size());
    for (int idx = 0; idx < src_matched.size(); ++idx){
        PointType sP = src_matched[idx];
        PointType sT = tgt_matched[idx];
        srcP.x = sP.x;     srcP.y = sP.y;     srcP.z = sP.z;
        tgtP.x = sT.x;     tgtP.y = sT.y;     tgtP.z = sT.z;

        marker.points.emplace_back(srcP);
        marker.points.emplace_back(tgtP);
    }
}

#endif
