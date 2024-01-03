/**
 * Copyright 2020, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Jingnan Shi, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 */
#pragma once

#include <chrono>
#include <numeric>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <flann/flann.hpp>
#include <teaser/geometry.h>
#include <quatro/fpfh.h>

#ifdef TBB_EN
#include <tbb/tbb.h>
#endif

namespace teaser {

class Matcher {
public:
  typedef std::vector<Eigen::VectorXf> Feature;
  typedef flann::Index<flann::L2<float>> KDTree;

  Matcher() = default;
  /**
   * Calculate correspondences based on given features and point clouds.
   * input: source point cloud, target point cloud
   * output: correspondences
   */
  std::vector<std::pair<int, int>>
  calculateCorrespondences(const teaser::PointCloud& source_points, const teaser::PointCloud& target_points,
                           const teaser::FPFHCloud& source_features, const teaser::FPFHCloud& target_features,
                           const bool& use_absolute_scale = false, const bool& use_crosscheck = true,
                           const bool& use_tuple_test = true, const float& tuple_scale = 0.95, 
                           const bool& use_optimized_matching = true, const float& thr_dist = 30.0, const int& num_max_corres = 600);

private:
  // variables
  std::vector<std::pair<int, int>> corres_;
  std::vector<teaser::PointCloud> pointcloud_;
  std::vector<Feature> features_;
  std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > means_; // for normalization
  float global_scale_;
  
  // methods
  template <typename T> void buildKDTree(const std::vector<T>& data, KDTree* tree);
#ifdef TBB_EN
  template <typename T> void buildKDTreeWithTBB(const std::vector<T>& data, KDTree* tree);
#endif
  template <typename T>
  void searchKDTree(KDTree* tree, const T& input, std::vector<int>& indices,
                    std::vector<float>& dists, int nn);
  template <typename T>
  void searchKDTreeAll(Matcher::KDTree* tree, const std::vector<T>& inputs,
                              std::vector<int>& indices, std::vector<float>& dists, int nn);
  void advancedMatching(const bool& use_crosscheck, const bool& use_tuple_test, const float& tuple_scale);
  void optimizedMatching(const float& thr_dist, const int& num_max_corres, const float& tuple_scale);
  void normalizePoints(const bool& use_absolute_scale);
};

} // namespace teaser
