#ifndef QUATRO_MODULE_H
#define QUATRO_MODULE_H


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
class quatro
{
	private:
		double m_normal_radius = 0.02, m_fpfh_radius = 0.04;
		double m_noise_bound = 0.25, m_rotation_gnc_factor = 1.39, m_rotation_cost_thr = 0.0001;
		int m_rotation_max_iter = 100;
		bool m_estimate_scale = false;
		teaser::RobustRegistrationSolver::Params m_quatro_params;
	public:
		quatro(){};
		quatro(const double &fpfh_normal_radi, const double &fpfh_radi, const double noise_bound, const double &rot_gnc_fact, const double &rot_cost_thr, const int &rot_max_iter, const bool &estimat_scale);
		template <typename T>
		Eigen::Matrix4d align(const pcl::PointCloud<T> &src, const pcl::PointCloud<T> &dst, bool &if_valid);
	private:
		void set_params();
		template <typename T>
		teaser::PointCloud pcl_to_teaser_pcl(const pcl::PointCloud<T> &cloud_in);
};


#endif