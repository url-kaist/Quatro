quatro::quatro(const double &fpfh_normal_radi, const double &fpfh_radi, const double noise_bound, const double &rot_gnc_fact, const double &rot_cost_thr, const int &rot_max_iter, const bool &estimat_scale)
{
	m_normal_radius = fpfh_normal_radi;
	m_fpfh_radius = fpfh_radi;
	m_noise_bound = noise_bound;
	m_rotation_gnc_factor = rot_gnc_fact;
	m_rotation_cost_thr = rot_cost_thr;
	m_rotation_max_iter = rot_max_iter;
	m_estimate_scale = estimat_scale;
}

template <typename T>
teaser::PointCloud quatro::pcl_to_teaser_pcl(const pcl::PointCloud<T> &cloud_in)
{
	teaser::PointCloud t_pcl_out_;
	if (cloud_in.size() > 0 ) t_pcl_out_.reserve(cloud_in.size());
  for (size_t i = 0; i < cloud_in.size(); ++i)
  {
    tgt_cloud.push_back({cloud_in.points[i].x, cloud_in.points[i].y, cloud_in.points[i].z});
  }
	return t_pcl_out_;
}

void quatro::reset_params()
{
  m_quatro_params.noise_bound = m_noise_bound;
  m_quatro_params.estimate_scaling = m_estimate_scale;
  m_quatro_params.rotation_max_iterations = m_rotation_max_iter;
  m_quatro_params.rotation_gnc_factor = m_rotation_gnc_factor;
  m_quatro_params.rotation_cost_threshold = m_rotation_cost_thr;
  m_quatro_params.cbar2 = 1;
  m_quatro_params.rotation_estimation_algorithm = teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::QUATRO;
  m_quatro_params.inlier_selection_mode == teaser::RobustRegistrationSolver::INLIER_SELECTION_MODE::PMC_HEU;
	return;
}

template <typename T>
Eigen::Matrix4d quatro::align(const pcl::PointCloud<T> &src, const pcl::PointCloud<T> &dst)
{
	Eigen::Matrix4d out_tf_ = Eigen::Matrix4d::Identity();

	teaser::PointCloud src_cloud_ = pcl_to_teaser_pcl(src);
	teaser::PointCloud tgt_cloud_ = pcl_to_teaser_pcl(dst);
  teaser::FPFHEstimation fpfh_;
  teaser::FPFHCloudPtr obj_descriptors_ = fpfh_.computeFPFHFeatures(src_cloud_, m_normal_radius, m_fpfh_radius);
  teaser::FPFHCloudPtr scene_descriptors_ = fpfh_.computeFPFHFeatures(tgt_cloud_, m_normal_radius, m_fpfh_radius);

  teaser::Matcher matcher_;
  std::vector<std::pair<int, int>> correspondences_ = matcher_.calculateCorrespondences(src_cloud_, tgt_cloud_, *obj_descriptors_, *scene_descriptors_, false, true, false, 0.95);
  
  reset_params();
	teaser::RobustRegistrationSolver Quatro_(m_quatro_params);

	Quatro_.solve(src_cloud_, tgt_cloud_, correspondences_);
	auto solution_by_quatro = Quatro_.getSolution();
	
  out_tf_.block<3, 3>(0, 0) = solution_by_quatro.rotation;
  out_tf_.block<3, 1>(0, 3) = solution_by_quatro.translation;

	return out_tf_;
}