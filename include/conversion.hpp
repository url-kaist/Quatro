#ifndef CONVERSION_HPP
#define CONVERSION_HPP

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_types_conversion.h>
#include <pcl_conversions/pcl_conversions.h>
#include <teaser/geometry.h>

template<typename T>
void pcl2teaser(const pcl::PointCloud<T> &pcl_raw, teaser::PointCloud &cloud) {
    cloud.clear();
    for (const auto &pt: pcl_raw.points) {
        cloud.push_back({pt.x, pt.y, pt.z});
    }
}

template<typename T>
sensor_msgs::PointCloud2 cloud2msg(pcl::PointCloud<T> cloud, std::string frame_id = "map") {
    sensor_msgs::PointCloud2 cloud_ROS;
    pcl::toROSMsg(cloud, cloud_ROS);
    cloud_ROS.header.frame_id = frame_id;
    return cloud_ROS;
}

void xyzi2xyz(pcl::PointCloud<pcl::PointXYZI>::Ptr XYZI, pcl::PointCloud<pcl::PointXYZ>::Ptr XYZ) {
    (*XYZ).points.resize((*XYZI).size());
    for (size_t i = 0; i < (*XYZI).points.size(); i++) {
        (*XYZ).points[i].x = (*XYZI).points[i].x;
        (*XYZ).points[i].y = (*XYZI).points[i].y;
        (*XYZ).points[i].z = (*XYZI).points[i].z;
    }
}

template<typename T>
void pcl2eigen(const pcl::PointCloud<T> &pcl_raw, Eigen::Matrix<double, 3, Eigen::Dynamic> &cloud) {
    int N = pcl_raw.points.size();
    cloud.resize(3, N);
    for (int i = 0; i < N; ++i) {
        cloud.col(i) << pcl_raw.points[i].x, pcl_raw.points[i].y, pcl_raw.points[i].z;
    }
}

template<typename T>
void eigen2pcl(const Eigen::Matrix<double, 3, Eigen::Dynamic> &src, pcl::PointCloud<T> &cloud) {
    int num_pc = src.cols();
    T   pt_tmp;
    if (!cloud.empty()) cloud.clear();
    for (int i = 0; i < num_pc; ++i) {
        pt_tmp.x = src(0, i);
        pt_tmp.y = src(1, i);
        pt_tmp.z = src(2, i);
        cloud.points.emplace_back(pt_tmp);
    }
}

#endif // CONVERSION_HPP
