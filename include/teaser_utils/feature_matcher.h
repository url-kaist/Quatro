/**
 * Copyright 2020, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Jingnan Shi, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 */

#pragma once

#include <flann/flann.hpp>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include "teaser/geometry.h"
#include "./fpfh.h"

namespace teaser {

class Matcher {
public:
  typedef std::vector<Eigen::VectorXf> Feature;
  typedef flann::Index<flann::L2<float>> KDTree;

  // New methods
  // Public methods:
  // 1. calculateCorrespondences
  //    input: source point cloud, target point cloud
  //    output: correspondences
  Matcher() = default;

  /**
   * Calculate correspondences based on given features and point clouds.
   * @param source_points
   * @param target_points
   * @param use_absolute_scaleS
   * @param use_crosscheck
   * @param use_tuple_test
   * @return
   */

//  template <typename T>
  std::vector<std::pair<int, int>> calculateCorrespondences(teaser::PointCloud& source_points, teaser::PointCloud& target_points,
                           FPFHCloud & source_features, FPFHCloud& target_features,
                           bool use_absolute_scale = true, bool use_crosscheck = true,
                           bool use_tuple_test = true, float tuple_scale = 0){
      Feature cloud_features;
      pointcloud_.push_back(source_points);
      pointcloud_.push_back(target_points);

      // It compute the global_scale_ required to set correctly the search radius
      normalizePoints(use_absolute_scale);
      int size_descriptor = source_features[0].descriptorSize();

      for (auto& f : source_features) {
          Eigen::VectorXf fpfh(size_descriptor);
          for (int i = 0; i < size_descriptor; i++)
              fpfh(i) = f.histogram[i];
          cloud_features.push_back(fpfh);
      }
      features_.push_back(cloud_features);

      cloud_features.clear();
      for (auto& f : target_features) {
          Eigen::VectorXf fpfh(size_descriptor);
          for (int i = 0; i < size_descriptor; i++)
              fpfh(i) = f.histogram[i];
          cloud_features.push_back(fpfh);
      }
      features_.push_back(cloud_features);

      advancedMatching(use_crosscheck, use_tuple_test, tuple_scale);

      return corres_;
  }

private:
  template <typename T> void buildKDTree(const std::vector<T>& data, KDTree* tree);

  template <typename T>
  void searchKDTree(KDTree* tree, const T& input, std::vector<int>& indices,
                    std::vector<float>& dists, int nn);

  void advancedMatching(bool use_crosscheck, bool use_tuple_test, float tuple_scale);

  void normalizePoints(bool use_absolute_scale);

  std::vector<std::pair<int, int>> corres_;
  std::vector<teaser::PointCloud> pointcloud_;
  std::vector<Feature> features_;
  std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > means_; // for normalization
  float global_scale_;
};

} // namespace teaser
