#pragma once


///// Eigen
#include <Eigen/Core>
///// PCL
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
///// imported headers - teaser
#include <teaser/registration.h>
#include <teaser/geometry.h>
#include <teaser/utils.h>
///// this package header
#include <quatro/matcher.h>

using namespace std;

//////////////////////////////////////////////////////
template <typename PointType>
class quatro
{
	private:
		int m_rotation_max_iter = 100, m_num_max_corres = 200;
		double m_normal_radius = 0.02, m_fpfh_radius = 0.04, m_distance_threshold = 30.0;
		double m_noise_bound = 0.25, m_rotation_gnc_factor = 1.39, m_rotation_cost_thr = 0.0001;
		bool m_estimate_scale = false, m_use_optimized_matching = true;
		teaser::RobustRegistrationSolver::Params m_quatro_params;
	public:
		quatro(){};
		quatro(const double &fpfh_normal_radi, const double &fpfh_radi, const double noise_bound, const double &rot_gnc_fact,
				const double &rot_cost_thr, const int &rot_max_iter, const bool &estimat_scale,
				const bool& use_optimized_matching = true, const double& distance_threshold = 30.0, const int& num_max_corres = 200);
		Eigen::Matrix4d align(const pcl::PointCloud<PointType> &src, const pcl::PointCloud<PointType> &dst, bool &if_valid);
	private:
		void set_params();
		teaser::PointCloud pcl_to_teaser_pcl(const pcl::PointCloud<PointType> &cloud_in);
};