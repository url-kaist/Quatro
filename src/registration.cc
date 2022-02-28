/**
 * Copyright 2020, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Jingnan Shi, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 */

#include "teaser/registration.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <iterator>
#include <ctime>
#include "teaser/utils.h"
#include "teaser/graph.h"
#include "teaser/macros.h"
#include <fstream>
#include <iomanip>

double findMedian(std::vector<int> a, int n) {

    // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Applying nth_element
        // on (n-1)/2 th index
        nth_element(a.begin(), a.begin() + (n - 1) / 2, a.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double) (a[(n - 1) / 2] + a[n / 2]) / 2.0;
    }

        // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Value at index (N/2)th
        // is the median
        return (double) a[n / 2];
    }
}

// - - - - - - - - - - - - - - - - - - - -
bool        SAVE_RESULTS    = false;
int         COUNT_EXECUTION = 0;
std::string SAVE_ABS_PATH   = "/home/shapelim/viz_for_paper";
// - - - - - - - - - - - - - - - - - - - -

void teaser::ScalarTLSEstimator::estimate(const Eigen::RowVectorXd &X,
                                          const Eigen::RowVectorXd &ranges, double*estimate,
                                          Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers,
                                          bool using_median_selection) {
    // check input parameters
    bool dimension_inconsistent = (X.rows() != ranges.rows()) || (X.cols() != ranges.cols());
    if (inliers) {
        dimension_inconsistent |= ((inliers->rows() != 1) || (inliers->cols() != ranges.cols()));
    }
    bool only_one_element = (X.rows() == 1) && (X.cols() == 1);
    assert(!dimension_inconsistent);
    assert(!only_one_element); // TODO: admit a trivial solution

    int                                 N = X.cols();
    std::vector<std::pair<double, int>> h;
    for (size_t                         i = 0; i < N; ++i) {
        h.push_back(std::make_pair(X(i) - ranges(i), i + 1));
        h.push_back(std::make_pair(X(i) + ranges(i), -i - 1));
    }

    // ascending order
    std::sort(h.begin(), h.end(),
              [](std::pair<double, int> a, std::pair<double, int> b) { return a.first < b.first; });

    // calculate weights
    Eigen::RowVectorXd weights = ranges.array().square();
    weights = weights.array().inverse();
    int                nr_centers      = 2 * N;
    Eigen::RowVectorXd x_hat           = Eigen::MatrixXd::Zero(1, nr_centers);
    Eigen::RowVectorXd x_cost          = Eigen::MatrixXd::Zero(1, nr_centers);
    Eigen::RowVectorXi set_cardinality = Eigen::MatrixXi::Zero(1, nr_centers);

    double ranges_inverse_sum     = ranges.sum();
    double dot_X_weights          = 0;
    double dot_weights_consensus  = 0;
    int    consensus_set_cardinal = 0;
    double sum_xi                 = 0;
    double sum_xi_square          = 0;

    // To save outputs
    static int component_idx = 0;
    component_idx = component_idx % 3;

    std::string   h_path, estimate_path, idx_path;
    std::ofstream h_txt, estimate_txt, idx_txt;
    // if (SAVE_RESULTS) {
    //     h_path =
    //             (boost::format("%s/%04d_h_%d.txt") % SAVE_ABS_PATH % COUNT_EXECUTION % component_idx).str();

    //     estimate_path = (boost::format("%s/%04d_estimates_%d.txt") % SAVE_ABS_PATH % COUNT_EXECUTION %
    //                      component_idx)
    //             .str();

    //     idx_path = (boost::format("%s/%04d_final_idx_%d.txt") % SAVE_ABS_PATH % COUNT_EXECUTION %
    //                 component_idx)
    //             .str();
    //     estimate_txt.open(estimate_path);
    //     estimate_txt.close();

    //     idx_txt.open(estimate_path);
    //     idx_txt.close();

    //     h_txt.open(h_path);
    //     for (int ijk = 0; ijk < h.size(); ++ijk) {
    //         h_txt << h[ijk].first << " " << h[ijk].second << "\n";
    //     }
    //     h_txt.close();
    // }

    for (size_t i = 0; i < nr_centers; ++i) {

        int idx     = int(std::abs(h.at(i).second)) - 1; // Indices starting at 1
        int epsilon = (h.at(i).second > 0) ? 1 : -1;

        consensus_set_cardinal += epsilon;
        dot_weights_consensus += epsilon * weights(idx);
        dot_X_weights += epsilon * weights(idx) * X(idx); // X(idx): v_ij in my paper
        ranges_inverse_sum -= epsilon * ranges(idx);      // handling truncated loss!
        sum_xi += epsilon * X(idx);
        sum_xi_square += epsilon * X(idx) * X(idx);

        if (using_median_selection) {
            // To set median values
            set_cardinality(i) = consensus_set_cardinal;
        }

        x_hat(i) = dot_X_weights / dot_weights_consensus;
        // sum_xi_square: already includes consensus_set_cardinal
        double residual =
                       consensus_set_cardinal * x_hat(i) * x_hat(i) + sum_xi_square - 2 * sum_xi * x_hat(i);
        x_cost(i) = residual + ranges_inverse_sum;

        // if (SAVE_RESULTS) {
        //     estimate_txt.open(estimate_path, std::ios::app);
        //     estimate_txt << consensus_set_cardinal << " " << idx << " " << X(idx) << " " << x_hat(i)
        //                  << " " << x_cost(i) << "\n";
        //     estimate_txt.close();
        // }
    }

    size_t min_idx;
    x_cost.minCoeff(&min_idx);

    double median;
    if (using_median_selection) {
        int n_card = set_cardinality(min_idx);
        std::cout << "\033[1;32m[Median Selection] num: " << n_card << "\033[0m" << std::endl;
        if (n_card > 0) {
            std::vector<double> candidates;
            candidates.reserve(set_cardinality(min_idx));
            for (int j = 0; j < n_card; ++j) {
                int idx = int(std::abs(h.at(min_idx - j).second)) - 1;
                assert(idx >= 0);
                candidates.push_back(X(idx));
            }
            sort(candidates.begin(), candidates.end());
            median = static_cast<double>(candidates[candidates.size() / 2 - 1] +
                                         candidates[candidates.size() / 2]) /
                     2.0;
        }
    }

    // if (SAVE_RESULTS) {
    //     idx_txt.open(idx_path);
    //     idx_txt << min_idx << std::endl;
    //     idx_txt.close();
    // }
    double estimate_temp;
    if (using_median_selection) {
        estimate_temp = median;
    } else {
        estimate_temp = x_hat(min_idx);
    }

    if (estimate) {
        // update estimate output if it's not nullptr
        *estimate = estimate_temp;
    }
    if (inliers) {
        // update inlier output if it's not nullptr
        *inliers = (X.array() - estimate_temp).array().abs() <= ranges.array();
    }
    component_idx++;
}

void teaser::ScalarTLSEstimator::estimate_tiled(const Eigen::RowVectorXd &X,
                                                const Eigen::RowVectorXd &ranges, const int &s,
                                                double*estimate,
                                                Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers) {
    // check input parameters
    bool dimension_inconsistent = (X.rows() != ranges.rows()) || (X.cols() != ranges.cols());
    if (inliers) {
        dimension_inconsistent |= ((inliers->rows() != 1) || (inliers->cols() != ranges.cols()));
    }
    bool only_one_element = (X.rows() == 1) && (X.cols() == 1);
    assert(!dimension_inconsistent);
    assert(!only_one_element); // TODO: admit a trivial solution

    // Prepare variables for calculations
    int                N          = X.cols();
    Eigen::RowVectorXd h(N * 2);
    h << X - ranges, X + ranges;
    // ascending order
    std::sort(h.data(), h.data() + h.cols(), [](double a, double b) { return a < b; });
    // calculate interval centers
    Eigen::RowVectorXd h_centers  = (h.head(h.cols() - 1) + h.tail(h.cols() - 1)) / 2;
    auto               nr_centers = h_centers.cols();

    // calculate weights
    Eigen::RowVectorXd weights = ranges.array().square();
    weights = weights.array().inverse();

    Eigen::RowVectorXd x_hat  = Eigen::MatrixXd::Zero(1, nr_centers);
    Eigen::RowVectorXd x_cost = Eigen::MatrixXd::Zero(1, nr_centers);

    // loop tiling
    size_t ih_bound = ((nr_centers) & ~((s) - 1));
    size_t jh_bound = ((N) & ~((s) - 1));

    std::vector<double>              ranges_inverse_sum_vec(nr_centers, 0);
    std::vector<double>              dot_X_weights_vec(nr_centers, 0);
    std::vector<double>              dot_weights_consensus_vec(nr_centers, 0);
    std::vector<std::vector<double>> X_consensus_table(nr_centers, std::vector<double>());

    auto inner_loop_f = [&](const size_t &i, const size_t &jh, const size_t &jl_lower_bound,
                            const size_t &jl_upper_bound) {
        double              &ranges_inverse_sum    = ranges_inverse_sum_vec[i];
        double              &dot_X_weights         = dot_X_weights_vec[i];
        double              &dot_weights_consensus = dot_weights_consensus_vec[i];
        std::vector<double> &X_consensus_vec       = X_consensus_table[i];

        size_t      j  = 0;
        for (size_t jl = jl_lower_bound; jl < jl_upper_bound; ++jl) {
            j = jh + jl;
            bool consensus = std::abs(X(j) - h_centers(i)) <= ranges(j);
            if (consensus) {
                dot_X_weights += X(j) * weights(j);
                dot_weights_consensus += weights(j);
                X_consensus_vec.push_back(X(j));
            } else {
                ranges_inverse_sum += ranges(j);
            }
        }

        if (j == N - 1) {
            // x_hat(i) = dot(X(consensus), weights(consensus)) / dot(weights, consensus);
            x_hat(i)                             = dot_X_weights / dot_weights_consensus;

            // residual = X(consensus)-x_hat(i);
            Eigen::Map<Eigen::VectorXd> X_consensus(X_consensus_vec.data(), X_consensus_vec.size());
            Eigen::VectorXd             residual = X_consensus.array() - x_hat(i);

            // x_cost(i) = dot(residual,residual) + sum(ranges(~consensus));
            x_cost(i) = residual.squaredNorm() + ranges_inverse_sum;
        }
    };

#pragma omp parallel for default(none) shared(                                                     \
    jh_bound, ih_bound, ranges_inverse_sum_vec, dot_X_weights_vec, dot_weights_consensus_vec, \
    X_consensus_table, h_centers, weights, N, X, x_hat, x_cost, s, ranges, inner_loop_f)
    for (size_t ih = 0; ih < ih_bound; ih += s) {
        for (size_t jh = 0; jh < jh_bound; jh += s) {
            for (size_t il = 0; il < s; ++il) {
                size_t i = ih + il;
                inner_loop_f(i, jh, 0, s);
            }
        }
    }

    // finish the left over entries
    // 1. Finish the unfinished js
#pragma omp parallel for default(none)                                                             \
    shared(jh_bound, ih_bound, ranges_inverse_sum_vec, dot_X_weights_vec, \
           dot_weights_consensus_vec, X_consensus_table, h_centers, weights, N, X, x_hat, x_cost, \
           s, ranges, nr_centers, inner_loop_f)
    for (size_t i = 0; i < nr_centers; ++i) {
        inner_loop_f(i, 0, jh_bound, N);
    }

    // 2. Finish the unfinished is
#pragma omp parallel for default(none)                                                             \
    shared(jh_bound, ih_bound, ranges_inverse_sum_vec, dot_X_weights_vec, \
           dot_weights_consensus_vec, X_consensus_table, h_centers, weights, N, X, x_hat, x_cost, \
           s, ranges, nr_centers, inner_loop_f)
    for (size_t i = ih_bound; i < nr_centers; ++i) {
        inner_loop_f(i, 0, 0, N);
    }

    size_t min_idx;
    x_cost.minCoeff(&min_idx);
    double estimate_temp = x_hat(min_idx);
    if (estimate) {
        // update estimate output if it's not nullptr
        *estimate = estimate_temp;
    }
    if (inliers) {
        // update inlier output if it's not nullptr
        *inliers = (X.array() - estimate_temp).array().abs() <= ranges.array();
    }
}

void teaser::FastGlobalRegistrationSolver::solveForRotation(
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &src,
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &dst,
        Eigen::Matrix3d*rotation,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers) {
    assert(rotation);                 // make sure R is not a nullptr
    assert(src.cols() == dst.cols()); // check dimensions of input data
    assert(params_.gnc_factor > 1);   // make sure mu will decrease
    assert(params_.noise_bound != 0); // make sure noise bound is not zero
    if (inliers) {
        assert(inliers->cols() == src.cols());
    }

    // Prepare some intermediate variables
    double noise_bound_sq = std::pow(params_.noise_bound, 2);
    size_t match_size     = src.cols();
    cost_ = std::numeric_limits<double>::infinity();

    // Calculate the initial mu
    double src_diameter  = teaser::utils::calculateDiameter<double, 3>(src);
    double dest_diameter = teaser::utils::calculateDiameter<double, 3>(dst);
    double global_scale  = src_diameter > dest_diameter ? src_diameter : dest_diameter;
    global_scale /= noise_bound_sq;
    double mu = std::pow(global_scale, 2) / noise_bound_sq;

    // stopping condition for mu
    double min_mu = 1.0;
    *rotation = Eigen::Matrix3d::Identity(3, 3); // rotation matrix
    Eigen::Matrix<double, 1, Eigen::Dynamic> l_pq(1, match_size);
    l_pq.setOnes(1, match_size);

    // Assumptions of the two inputs:
    // they should be of the same scale,
    // outliers should be removed as much as possible
    // input vectors should contain TIM vectors (if only estimating rotation)
    for (size_t i = 0; i < params_.max_iterations; ++i) {
        double scaled_mu              = mu * noise_bound_sq;

        // 1. Optimize for line processes weights
        Eigen::Matrix<double, 3, 1> q, p, rpq;
        for (size_t                 j = 0; j < match_size; ++j) {
            // p = Rq
            q   = src.col(j);
            p   = dst.col(j);
            rpq = p - (*rotation) * q;
            l_pq(j) = std::pow(scaled_mu / (scaled_mu + rpq.squaredNorm()), 2);
        }

        // 2. Optimize for Rotation Matrix
        *rotation = teaser::utils::svdRot(src, dst, l_pq, i);

        // update cost
        Eigen::Matrix<double, 3, Eigen::Dynamic> diff = (dst - (*rotation) * src).array().square();
        cost_ = ((scaled_mu * diff.colwise().sum()).array() /
                 (scaled_mu + diff.colwise().sum().array()).array())
                .sum();

        // additional termination conditions
        if (cost_ < params_.cost_threshold || mu < min_mu) {
            TEASER_DEBUG_INFO_MSG("Convergence condition met.");
            TEASER_DEBUG_INFO_MSG("Iterations: " << i);
            TEASER_DEBUG_INFO_MSG("Mu: " << mu);
            TEASER_DEBUG_INFO_MSG("Cost: " << cost_);
            break;
        }

        // update mu
        mu /= params_.gnc_factor;
    }

    if (inliers) {
        *inliers = l_pq.cast<bool>();
    }
}

void teaser::FastGlobalRegistrationSolver::solveForRotation2D(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &src,
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &dst, Eigen::Matrix2d*rotation,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers) {
    throw std::invalid_argument("Not implemented!");
}

void teaser::TLSScaleSolver::solveForScale(const Eigen::Matrix<double, 3, Eigen::Dynamic> &src,
                                           const Eigen::Matrix<double, 3, Eigen::Dynamic> &dst,
                                           double*scale,
                                           Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers) {

    Eigen::Matrix<double, 1, Eigen::Dynamic> v1_dist =
                                                     src.array().square().colwise().sum().array().sqrt();
    Eigen::Matrix<double, 1, Eigen::Dynamic> v2_dist =
                                                     dst.array().square().colwise().sum().array().sqrt();

    Eigen::Matrix<double, 1, Eigen::Dynamic> raw_scales = v2_dist.array() / v1_dist.array();
    double                                   beta       = 2 * noise_bound_ * sqrt(cbar2_);
    Eigen::Matrix<double, 1, Eigen::Dynamic> alphas     = beta * v1_dist.cwiseInverse();

    tls_estimator_.estimate(raw_scales, alphas, scale, inliers);
}

void teaser::ScaleInliersSelector::solveForScale(
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &src,
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &dst,
        double*scale,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers) {
    // We assume no scale difference between the two vectors of points.
    *scale = 1;

    Eigen::Matrix<double, 1, Eigen::Dynamic> v1_dist =
                                                     src.array().square().colwise().sum().array().sqrt();
    Eigen::Matrix<double, 1, Eigen::Dynamic> v2_dist =
                                                     dst.array().square().colwise().sum().array().sqrt();
    double                                   beta    = 2 * noise_bound_ * sqrt(cbar2_);

    // A pair-wise correspondence is an inlier if it passes the following two tests:
    // 1. dst / src is within maximum allowed error
    // 2. src / dst is within maximum allowed error
    Eigen::Matrix<double, 1, Eigen::Dynamic> alphas_forward     = beta * v1_dist.cwiseInverse();
    Eigen::Matrix<double, 1, Eigen::Dynamic> raw_scales_forward = v2_dist.array() / v1_dist.array();
    Eigen::Matrix<bool, 1, Eigen::Dynamic>   inliers_forward    =
                                                     (raw_scales_forward.array() - *scale).array().abs() <=
                                                     alphas_forward.array();

    Eigen::Matrix<double, 1, Eigen::Dynamic> alphas_reverse     = beta * v2_dist.cwiseInverse();
    Eigen::Matrix<double, 1, Eigen::Dynamic> raw_scales_reverse = v1_dist.array() / v2_dist.array();
    Eigen::Matrix<bool, 1, Eigen::Dynamic>   inliers_reverse    =
                                                     (raw_scales_reverse.array() - *scale).array().abs() <=
                                                     alphas_reverse.array();

    // element-wise AND using component-wise product (Eigen 3.2 compatible)
    *inliers = inliers_forward.cwiseProduct(inliers_reverse);
}

void teaser::TLSTranslationSolver::solveForTranslation(
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &src,
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &dst,
        Eigen::Vector3d*translation,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers,
        bool using_median_selection) {
    assert(src.cols() == dst.cols());
    if (inliers) {
        assert(inliers->cols() == src.cols());
    }

    // Raw translation
    Eigen::Matrix<double, 3, Eigen::Dynamic> raw_translation = dst - src;

    // Error bounds for each measurements
    int                                      N      = src.cols();
    double                                   beta   = noise_bound_ * sqrt(cbar2_);      
    Eigen::Matrix<double, 1, Eigen::Dynamic> alphas = beta * Eigen::MatrixXd::Ones(1, N);

    // Estimate x, y, and z component of translation: perform TLS on each row
    *inliers = Eigen::Matrix<bool, 1, Eigen::Dynamic>::Ones(1, N);
    Eigen::Matrix<bool, 1, Eigen::Dynamic> inliers_temp(1, N);

    for (size_t i = 0; i < raw_translation.rows(); ++i) {
        tls_estimator_.estimate(raw_translation.row(i), alphas, &((*translation)(i)), &inliers_temp,
                                using_median_selection);
        // element-wise AND using component-wise product (Eigen 3.2 compatible)
        // a point is an inlier iff. x,y,z are all inliers
        *inliers = (*inliers).cwiseProduct(inliers_temp);
    }
}

teaser::RobustRegistrationSolver::RobustRegistrationSolver(
        const teaser::RobustRegistrationSolver::Params &params) {
    reset(params);
}

Eigen::Matrix<double, 3, Eigen::Dynamic>
teaser::RobustRegistrationSolver::computeTIMs(const Eigen::Matrix<double, 3, Eigen::Dynamic> &v,
                                              Eigen::Matrix<int, 2, Eigen::Dynamic>*map) {

    auto                                     N = v.cols();
    // std::cout<<"\033[1;32m"<< "N:"<<N<<"\033[0m"<<std::endl;
    // std::cout<<"\033[1;32m"<< "N2:"<<N * (N - 1) / 2<<"\033[0m"<<std::endl;
    Eigen::Matrix<double, 3, Eigen::Dynamic> vtilde(3, N * (N - 1) / 2);
    map->resize(2, N * (N - 1) / 2);

    // std::cout<<"\033[1;32m"<< "Debug 1-3"<<"\033[0m"<<std::endl;

#pragma omp parallel for default(none) shared(N, v, vtilde, map)
    for (size_t i = 0; i < N - 1; i++) {
        // Calculate some important indices
        // For each measurement, we compute the TIMs between itself and all the measurements after it.
        // For examples:
        // i=0: add N-1 TIMs
        // i=1: add N-2 TIMs
        // etc..
        // i=k: add N-1-k TIMs
        // And by arithmatic series, we can get the starting index of each segment be:
        // k*N - k*(k+1)/2
        size_t segment_start_idx = i * N - i * (i + 1) / 2;
        size_t segment_cols      = N - 1 - i;

        // calculate TIM
        Eigen::Matrix<double, 3, 1>              m    = v.col(i);
        Eigen::Matrix<double, 3, Eigen::Dynamic> temp = v - m * Eigen::MatrixXd::Ones(1, N);

        // concatenate to the end of the tilde vector
        vtilde.middleCols(segment_start_idx, segment_cols) = temp.rightCols(segment_cols);

        // populate the index map
        Eigen::Matrix<int, 2, Eigen::Dynamic> map_addition(2, N);
        for (size_t                           j            = 0; j < N; ++j) {
            map_addition(0, j) = i;
            map_addition(1, j) = j;
        }
        map->middleCols(segment_start_idx, segment_cols)   = map_addition.rightCols(segment_cols);
    }

    return vtilde;
}

teaser::RegistrationSolution
teaser::RobustRegistrationSolver::solve(const teaser::PointCloud &src_cloud,
                                        const teaser::PointCloud &dst_cloud,
                                        const std::vector<std::pair<int, int>> correspondences) {
    Eigen::Matrix<double, 3, Eigen::Dynamic> src(3, correspondences.size());
    Eigen::Matrix<double, 3, Eigen::Dynamic> dst(3, correspondences.size());
    for (size_t                              i = 0; i < correspondences.size(); ++i) {
        auto src_idx = std::get<0>(correspondences[i]);
        auto dst_idx = std::get<1>(correspondences[i]);
        src.col(i) << src_cloud[src_idx].x, src_cloud[src_idx].y, src_cloud[src_idx].z;
        dst.col(i) << dst_cloud[dst_idx].x, dst_cloud[dst_idx].y, dst_cloud[dst_idx].z;
    }
    return solve(src, dst);
}

teaser::RegistrationSolution
teaser::RobustRegistrationSolver::solve(const Eigen::Matrix<double, 3, Eigen::Dynamic> &src,
                                        const Eigen::Matrix<double, 3, Eigen::Dynamic> &dst) {
    assert(scale_solver_ && rotation_solver_ && translation_solver_);
                                               
    // Handle deprecated params
    // 나중에 주석풀기
    /*if (!params_.use_max_clique) {
        TEASER_DEBUG_INFO_MSG(
                "Using deprecated param field use_max_clique. Switch to inlier_selection_mode instead.");
        params_.inlier_selection_mode = INLIER_SELECTION_MODE::NONE;
    }
    if (!params_.max_clique_exact_solution) {
        TEASER_DEBUG_INFO_MSG("Using deprecated param field max_clique_exact_solution. Switch to "
                              "inlier_selection_mode instead.");
        params_.inlier_selection_mode = INLIER_SELECTION_MODE::PMC_HEU;
    }*/

    /**
     * Steps to estimate T/R/s
     *
     * Estimate Scale
     * -- compute TIMs
     *
     * Remove outliers
     * De-scale the TIMs
     *        v2tilde = v2tilde/s_est; % correct scale from v2 side, more stable
     *
     * Estimate rotation
     *
     * Estimate Translation
     */

    std::time_t start, middle, end;
    start     = clock();
    src_tims_ = computeTIMs(src, &src_tims_map_);
    dst_tims_ = computeTIMs(dst, &dst_tims_map_);
    end = clock();
    double tim_time = (double) (end - start) / CLOCKS_PER_SEC;
    TEASER_DEBUG_INFO_MSG("Starting scale solver.");
    solveForScale(src_tims_, dst_tims_);
    TEASER_DEBUG_INFO_MSG("Scale estimation complete.");

    // Calculate Maximum Clique
    // Note: the max_clique_ vector holds the indices of original measurements that are within the
    // max clique of the built inlier graph.
    double gen_inlier_time = 0;
    double max_clique_time = 0;
    double maxclique_time  = 0;
    double sort_time       = 0;
    if (params_.inlier_selection_mode != INLIER_SELECTION_MODE::NONE) {
        // Create inlier graph: A graph with (indices of) original measurements as vertices, and edges
        // only when the TIM between two measurements are inliers. Note: src_tims_map_ is the same as
        // dst_tim_map_
        start = clock();

        inlier_graph_.populateVertices(src.cols());
        int         disgrrrrrcp_idx = 5;

        for (size_t i        = 0; i < scale_inliers_mask_.cols(); ++i) {
            // if (i % static_cast<int>(scale_inliers_mask_.cols()/disp_idx) == 0){ std::cout<<i <<" / "<<
            // scale_inliers_mask_.cols()<<std::endl; }

            if (scale_inliers_mask_(0, i)) {
                inlier_graph_.addEdge(src_tims_map_(0, i), src_tims_map_(1, i));
            }

        }

        end                                     = clock();
        gen_inlier_time                         = (double) (end - start) / CLOCKS_PER_SEC;

        // std::cout<<"\033[1;32m"<< "Debug 2"<<"\033[0m"<<std::endl;
        teaser::MaxCliqueSolver::Params clique_params;

        if (params_.inlier_selection_mode == INLIER_SELECTION_MODE::PMC_EXACT) {
            clique_params.solver_mode = teaser::MaxCliqueSolver::CLIQUE_SOLVER_MODE::PMC_EXACT;
        } else if (params_.inlier_selection_mode == INLIER_SELECTION_MODE::PMC_HEU) {               // SetParams 에서
            clique_params.solver_mode = teaser::MaxCliqueSolver::CLIQUE_SOLVER_MODE::PMC_HEU;
        } else {
            clique_params.solver_mode = teaser::MaxCliqueSolver::CLIQUE_SOLVER_MODE::KCORE_HEU;
        }

        clique_params.time_limit                = params_.max_clique_time_limit;
        clique_params.kcore_heuristic_threshold = params_.kcore_heuristic_threshold;

        teaser::MaxCliqueSolver clique_solver(clique_params);
        start       = clock();
        max_clique_ = clique_solver.findMaxClique(inlier_graph_);
        middle      = clock();
        std::sort(max_clique_.begin(), max_clique_.end());
        end            = clock();
        maxclique_time = (double) (middle - start) / CLOCKS_PER_SEC;
        sort_time      = (double) (end - middle) / CLOCKS_PER_SEC;
        TEASER_DEBUG_INFO_MSG("Max Clique of scale estimation inliers: ");
#ifndef NDEBUG
        std::copy(max_clique_.begin(), max_clique_.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
#endif

        if (max_clique_.size() <= 1) {
            TEASER_DEBUG_INFO_MSG("Clique size too small. Abort.");
            solution_.valid = false;
            return solution_;
        }
    } 
    //나중에 주석풀기
    // else {
    //     // not using clique filtering is equivalent to saying all measurements are in the max clique
    //     max_clique_.reserve(src.cols());
    //     for (size_t i = 0; i < src.cols(); ++i) {
    //         max_clique_.push_back(i);
    //     }
    // }

    // Calculate new measurements & TIMs based on max clique inliers
    start = clock();
    if (params_.rotation_tim_graph == INLIER_GRAPH_FORMULATION::CHAIN) {        //default = CHAIN
        // ==============>
        // chain graph
        //    TEASER_DEBUG_INFO_MSG("Using chain graph for GNC rotation.");

        std::cout << "\033[1;32mNum. of maximum cliques: " << max_clique_.size() << "\033[0m"
                  << std::endl;
        num_maxclique_ = max_clique_.size();

        pruned_src_tims_.resize(3, max_clique_.size());
        pruned_dst_tims_.resize(3, max_clique_.size());
        src_tims_map_rotation_.resize(2, max_clique_.size());
        dst_tims_map_rotation_.resize(2, max_clique_.size());
        for (size_t i = 0; i < max_clique_.size(); ++i) {
            const auto &root = max_clique_[i];
            int leaf;
            if (i != max_clique_.size() - 1) {
                leaf = max_clique_[i + 1];
            } else {
                leaf = max_clique_[0];
            }
            pruned_src_tims_.col(i) = src.col(leaf) - src.col(root);
            pruned_dst_tims_.col(i) = dst.col(leaf) - dst.col(root);

            // populate the TIMs map
            dst_tims_map_rotation_(0, i) = leaf;
            dst_tims_map_rotation_(1, i) = root;
            src_tims_map_rotation_(0, i) = leaf;
            src_tims_map_rotation_(1, i) = root;
        }
    }
    // 나중에 주석풀기
    /*else {
        // complete graph
        //    TEASER_DEBUG_INFO_MSG("Using complete graph for GNC rotation.");
        // select the inlier measurements with max clique
        Eigen::Matrix<double, 3, Eigen::Dynamic> src_inliers(3, max_clique_.size());
        Eigen::Matrix<double, 3, Eigen::Dynamic> dst_inliers(3, max_clique_.size());
        for (size_t                              i = 0; i < max_clique_.size(); ++i) {
            src_inliers.col(i) = src.col(max_clique_[i]);
            dst_inliers.col(i) = dst.col(max_clique_[i]);
        }
        // construct the TIMs
        pruned_dst_tims_ = computeTIMs(dst_inliers, &dst_tims_map_rotation_);
        pruned_src_tims_ = computeTIMs(src_inliers, &src_tims_map_rotation_);
    }*/
    end   = clock();
    double prune_time = (double) (end - start) / CLOCKS_PER_SEC;

    // Remove scaling for rotation estimation
    pruned_dst_tims_ *= (1 / solution_.scale);

    // Update GNC rotation solver's noise bound with the new information
    // Note: this implicitly assumes that rotation_solver_'s noise bound
    // is set to the original noise bound of the measurements.
    auto params = rotation_solver_->getParams();
    params.noise_bound *= (2 / solution_.scale);
    rotation_solver_->setParams(params);

    // Solve for rotation
    start = clock();
    TEASER_DEBUG_INFO_MSG(
            (boost::format("Starting rotation solver w/ %d TIMs.") % pruned_src_tims_.cols()).str());
    solveForRotation(pruned_src_tims_, pruned_dst_tims_);


    /** HT Lim: 알고리즘 가시화를 위한 저장용 **/
    if (SAVE_RESULTS) {
        std::string   src_tim_path =
                              (boost::format("%s/%04d_tims_src.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        std::string   dst_tim_path =
                              (boost::format("%s/%04d_tims_dst.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        std::ofstream src_tims_txt(src_tim_path);
        std::ofstream tgt_tims_txt(dst_tim_path);
        for (int      l            = 0; l < pruned_src_tims_.cols(); ++l) {
            src_tims_txt << pruned_src_tims_(0, l) << " " << pruned_src_tims_(1, l) << " "
                         << pruned_src_tims_(2, l) << "\n";
            tgt_tims_txt << pruned_dst_tims_(0, l) << " " << pruned_dst_tims_(1, l) << " "
                         << pruned_dst_tims_(2, l) << "\n";
        }
        src_tims_txt.close();
        tgt_tims_txt.close();
    }

    rotation_inliers_.clear();
    // Save indices of inlier TIMs from GNC rotation estimation
    // By Hyungtae Lim
    // rotation_inliers_에 들어가는 건 max_clique index임!!
    for (size_t i = 0; i < rotation_inliers_mask_.cols(); ++i) {
        if (i == 0) {
            // Check (N-1)-th and 0-th
            if (rotation_inliers_mask_(0, rotation_inliers_mask_.cols() - 1) &&
                rotation_inliers_mask_(0, i)) {
                rotation_inliers_.emplace_back(i);
            }
        } else {
            // Check i-1th and i-th
            if (rotation_inliers_mask_(0, i - 1) && rotation_inliers_mask_(0, i)) {
                rotation_inliers_.emplace_back(i);
            }
        }
    }
    TEASER_DEBUG_INFO_MSG(
            (boost::format("Rotation estimation complete: %d inliers") % rotation_inliers_.size()).str());
    num_rot_inliers_ = rotation_inliers_.size();

    int                                      N_R = rotation_inliers_.size();
    Eigen::Matrix<double, 3, Eigen::Dynamic> rotation_pruned_src;
    Eigen::Matrix<double, 3, Eigen::Dynamic> rotation_pruned_dst;

    if (params_.using_rot_inliers_when_estimating_cote && (N_R > 0)) {
        /**
         * Rotation을 구한 후 해당 inliers만 가져옴
         * 그런데 Quatro를 통해 rotation을 했을 때에는 inliers들 중 quasi-inliers들이 섞여 있는데,
         * quasi-inliers이 가끔 consensus를 이뤄서 부정확해질 때가 있음
         **/
        for (size_t i = 0; i < N_R; ++i) {
            rotation_pruned_src.resize(3, N_R);
            rotation_pruned_dst.resize(3, N_R);
            rotation_pruned_src.col(i) = src.col(max_clique_[rotation_inliers_[i]]);
            rotation_pruned_dst.col(i) = dst.col(max_clique_[rotation_inliers_[i]]);
        }
    } else {
        int         N_MC = max_clique_.size();
        for (size_t i    = 0; i < N_MC; ++i) {
            rotation_pruned_src.resize(3, N_MC);
            rotation_pruned_dst.resize(3, N_MC);
            if (reg_name_ == "Quatro") {
                rotation_pruned_src.col(i) = estimated_RyRx_ * src.col(max_clique_[i]);
            } else if (reg_name_ == "TEASER") {
                rotation_pruned_src.col(i) = src.col(max_clique_[i]);
            }
            rotation_pruned_dst.col(i) = dst.col(max_clique_[i]);
        }
    }

    end = clock();
    double rot_time = (double) (end - start) / CLOCKS_PER_SEC;

    // Solve for translation
    TEASER_DEBUG_INFO_MSG(
            (boost::format("Starting trainslation solver w/ %d pruned_src.") % rotation_pruned_src.cols())
                    .str());
    /**
     * COTE 모드: "weighted_mean"과 "median"모드가 있음
     */
    if (params_.cote_mode == "weighted_mean") {
        solveForTranslation(solution_.scale * solution_.rotation * rotation_pruned_src,
                            rotation_pruned_dst);
    } else if (params_.cote_mode == "median") {
        solveForTranslation(solution_.scale * solution_.rotation * rotation_pruned_src,
                            rotation_pruned_dst, true);
    }

    TEASER_DEBUG_INFO_MSG("Translation estimation complete.");
    start = clock();
    double trans_time = (double) (start - end) / CLOCKS_PER_SEC;

    // Find the final inliers
    translation_inliers_ = utils::findNonzero<bool>(translation_inliers_mask_);

    final_inliers_.clear();
    final_inliers_.reserve(translation_inliers_.size());
    // final_inliers: 매칭쌍 상의 index 순서를 나타내야 함!
    // input src, dst의 indices를 리턴해줌
    if (params_.using_rot_inliers_when_estimating_cote && (N_R > 0)) {
        for (const auto &idx : translation_inliers_) {
            //    final_inliers_.push_back(max_clique_[idx]);
            final_inliers_.push_back(max_clique_[rotation_inliers_[idx]]);
        }
    } else {
        for (const auto &idx : translation_inliers_) {
            final_inliers_.push_back(max_clique_[idx]); // original!!
        }
    }

    // Update validity flag
    solution_.valid = true;

    // - - - - - - - - - - - - - - - - - - - -
    // Save clique results

    if (SAVE_RESULTS) {
        std::string final_inlier_path;
        final_inlier_path = (boost::format("%s/%04d_%s_final_inliers.txt") % SAVE_ABS_PATH %
                             params_.reg_name % COUNT_EXECUTION)
                .str();

        std::ofstream fi_txt(final_inlier_path);
        for (int      ijk = 0; ijk < translation_inliers_mask_.cols(); ++ijk) {
            if (translation_inliers_mask_[ijk]) {
                if (ijk != translation_inliers_mask_.size() - 1) {
                    fi_txt << ijk << " ";
                } else {
                    fi_txt << ijk << "\n";
                }
            }
        }
        fi_txt.close();
    }
    COUNT_EXECUTION++;

    return solution_;
}

double teaser::RobustRegistrationSolver::solveForScale(
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &v1,
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &v2) {
    scale_inliers_mask_.resize(1, v1.cols());
    scale_solver_->solveForScale(v1, v2, &(solution_.scale), &scale_inliers_mask_);
    std::cout << "Scale: " << solution_.scale << std::endl;
    return solution_.scale;
}

Eigen::Vector3d teaser::RobustRegistrationSolver::solveForTranslation(
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &v1,
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &v2, 
        bool using_median_selection)                            
{
    translation_inliers_mask_.resize(1, v1.cols());
    translation_solver_->solveForTranslation(v1, v2, &(solution_.translation), &translation_inliers_mask_, using_median_selection);

    return solution_.translation;
}

Eigen::Matrix3d teaser::RobustRegistrationSolver::solveForRotation(
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &v1,
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &v2) {
    rotation_inliers_mask_.resize(1, v1.cols());

    if (reg_name_ == "TEASER") {
        rotation_solver_->solveForRotation(v1, v2, &(solution_.rotation), &rotation_inliers_mask_);

    } 
    // 이 부분만 옮겨두됨
    else if (reg_name_ == "Quatro") {
        Eigen::Matrix2d                          rotation_2d;
        Eigen::Matrix<double, 2, Eigen::Dynamic> src_2d;
        Eigen::Matrix<double, 2, Eigen::Dynamic> dst_2d;
        // XY Coordinates for calculate yaw
        src_2d.resize(2, v1.cols());
        dst_2d.resize(2, v2.cols());
        src_2d = v1.topRows(2);
        dst_2d = v2.topRows(2);

        /**
         * XY에 대해 rotation을 구한 후, R3x3의 2x2 부분에 yaw rotation을 채워 줌
         */
        Eigen::Matrix3d rot_yaw = Eigen::Matrix3d::Identity();
        rotation_2d = Eigen::Matrix2d::Identity();
        rotation_solver_->solveForRotation2D(src_2d, dst_2d, &rotation_2d, &rotation_inliers_mask_);
        rot_yaw.block<2, 2>(0, 0) = rotation_2d;
        solution_.rotation = rot_yaw;

        /*
         * trial 1
         * Calculation of uncertainty of optimization
         * 하지만, optimization된 값의 covariance로 실패/성공을 구분해낼 수 없었음
         */
        //    double information = 0;
        //    int n_pairs = v1.cols();
        //    // c_h: cos hat, s_h: sine hat
        //    double c_h = rotation_2d(0, 0);
        //    double s_h = rotation_2d(1, 0);
        //    for (int i = 0; i < n_pairs; ++i) {
        //      double jacob1 = s_h * src_2d(0, i) + c_h * src_2d(1, i);
        //      double jacob2 = (-c_h) * src_2d(0, i) + s_h * src_2d(1, i);
        //      double jacob_sum = pow(jacob1, 2) + pow(jacob2, 2);
        //      information += jacob_sum;
        //    }
        //    rot_cov = 1 / (information + 0.000000001);
        //
        //    std::string costPath = (boost::format("%s/rot_cov.txt") %
        //    "/home/shapelim/dummy/00").str(); std::ofstream costTxt(costPath, std::ios::app); costTxt
        //    << rot_cov << std::endl; costTxt.close();

        /*
         * trial 2
         * yaw를 구한 후 ceres solver를 이용한 fine tuning?
         * 하지만 약간 성능이 개선되기는 하지만, 오히려 논문의 contribution을 떨어뜨려 배제됨
         * 경미한 수준이어서 IMU measurement 쓰는 게 더 나음
         */
        //    std::vector<int> inlier_indices;
        //    inlier_indices.reserve(v1.cols());
        //    for (int n = 0; n < rotation_inliers_mask_.size(); ++n) {
        //      if (rotation_inliers_mask_[n]) {
        //        inlier_indices.push_back(n);
        //      }
        //    }
        //    if (inlier_indices.size() > 3){
        //      ceres::Problem problem;
        //      double rot_angle_est = 0;
        //      for (int i = 0; i < inlier_indices.size(); ++i) {
        //        int l = inlier_indices[i] ;
        //        ceres::CostFunction *cost_function =
        //            new ceres::AutoDiffCostFunction<RotResidual, 1, 1>(new RotResidual(v1.col(l),
        //            v2.col(l)));
        //        problem.AddResidualBlock(cost_function, new ceres::CauchyLoss(0.5), &rot_angle_est);
        //      }
        //      // Run the solver!
        //      ceres::Solver::Options options;
        //      options.minimizer_progress_to_stdout = false;
        //      ceres::Solver::Summary summary;
        //      Solve(options, &problem, &summary);
        //
        //      std::cout << summary.BriefReport() << "\n";
        //      Eigen::Matrix3d rot_pitch = Eigen::Matrix3d::Identity();
        //      if (abs(rot_angle_est) < 7.0 * M_PI / 180.0 ){
        //        double c_ = cos(rot_angle_est);
        //        double s_ = sin(rot_angle_est);
        //        rot_pitch(0, 0) = c_;
        //        rot_pitch(2, 2) = c_;
        //        rot_pitch(0, 2) = s_;
        //        rot_pitch(2, 0) = -s_;
        //        solution_.rotation *= rot_pitch;
        //      }
        //    }

    } else
        throw std::invalid_argument(
                "[solveForRotation] The param is wrong! It should be 'TEASER' or 'Quatro'");

    /**
     * solveForRotation 함수로 rotation을 구한 후,
     * 미리 계산된 RyRx로 pose를 보정하는 과정
     * setPreEstimatedRyRx()함수를 호출하면 자동으로 using_pre_estimated_RyRX_가 true가 됨
     */
    if (using_pre_estimated_RyRx_) {
        std::cout << "\033[1;32m";
        std::cout << "IMU CORRECTION!!!" << reg_name_ << std::endl;
        std::cout << "IMU CORRECTION!!!\033[0m" << std::endl;
        if (reg_name_ == "TEASER") {
            std::cout << "TEASER++: Previous RyRx is discarded" << std::endl;
            double          yaw_est = atan2(solution_.rotation(1, 0), solution_.rotation(0, 0));
            Eigen::Matrix3d rot_yaw = Eigen::Matrix3d::Identity();
            double          c_y     = static_cast<double>(cos(yaw_est));
            double          s_y     = static_cast<double>(sin(yaw_est));
            rot_yaw(0, 0) = c_y;
            rot_yaw(0, 1) = -s_y;
            rot_yaw(1, 0) = s_y;
            rot_yaw(1, 1) = c_y;
            solution_.rotation = rot_yaw * estimated_RyRx_;
        } else if (reg_name_ == "Quatro") {
            std::cout << "Quatro: Pre-estimated RyRx is updated" << std::endl;
            solution_.rotation = solution_.rotation * estimated_RyRx_;
        } else {
            throw std::invalid_argument("Wrong reg type name is coming!");
        }
    }
    return solution_.rotation;
}

// Original
void teaser::GNCTLSRotationSolver::solveForRotation(    
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &src,
        const Eigen::Matrix<double, 3, Eigen::Dynamic> &dst,
        Eigen::Matrix3d*rotation,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers) {
    assert(rotation);                 // make sure R is not a nullptr
    assert(src.cols() == dst.cols()); // check dimensions of input data
    assert(params_.gnc_factor > 1);   // make sure mu will increase
    assert(params_.noise_bound != 0); // make sure noise sigma is not zero

    std::cout << "Conducting \033[1;32m3D\033[0m Rot solver" << std::endl;

    if (inliers) {
        assert(inliers->cols() == src.cols());
    }

    /**
     * Loop: terminate when:
     *    1. the change in cost in two consecutive runs is smaller than a user-defined threshold
     *    2. # iterations exceeds the maximum allowed
     *
     * Within each loop:
     * 1. fix weights and solve for R
     * 2. fix R and solve for weights
     */

    // Prepare some variables
    size_t match_size = src.cols(); // number of correspondences

    double mu = 1; // arbitrary starting mu

    double prev_cost = std::numeric_limits<double>::infinity();
    cost_ = std::numeric_limits<double>::infinity();
    double noise_bound_sq = std::pow(params_.noise_bound, 2);
    //  double noise_bound_sq = 1.25;
    if (noise_bound_sq < 1e-16) {
        noise_bound_sq = 1e-2;
    }
    TEASER_DEBUG_INFO_MSG("GNC rotation estimation noise bound:" << params_.noise_bound);
    TEASER_DEBUG_INFO_MSG("GNC rotation estimation noise bound squared:" << noise_bound_sq);

    std::string   costPath, estWeightPath, estRotPath;
    std::ofstream costTxt, weightTxt, rotTxt;
    if (SAVE_RESULTS) {
        costPath = (boost::format("%s/%04d_TEASER_cost.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        costTxt.open(costPath);
        estWeightPath =
                (boost::format("%s/%04d_TEASER_weights.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        weightTxt.open(estWeightPath);
        estRotPath = (boost::format("%s/%04d_TEASER_rot.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        rotTxt.open(estRotPath);
    }
    // - - - - - - - - - - - - - - - - - - - -

    Eigen::Matrix<double, 3, Eigen::Dynamic> diffs(3, match_size);
    Eigen::Matrix<double, 1, Eigen::Dynamic> weights(1, match_size);
    weights.setOnes(1, match_size);
    Eigen::Matrix<double, 1, Eigen::Dynamic> residuals_sq(1, match_size);

    // Loop for performing GNC-TLS
    for (size_t i = 0; i < params_.max_iterations; ++i) {

        // Fix weights and perform SVD rotation estimation
        *rotation = teaser::utils::svdRot(src, dst, weights, COUNT_EXECUTION);

        // Calculate residuals squared
        diffs        = (dst - (*rotation) * src).array().square();
        residuals_sq = diffs.colwise().sum();
        if (i == 0) {
            // Initialize rule for mu
            double max_residual = residuals_sq.maxCoeff();
            mu                  = 1 / (2 * max_residual / noise_bound_sq - 1);
            // Degenerate case: mu = -1 because max_residual is very small
            // i.e., little to none noise
            if (mu <= 0) {
                TEASER_DEBUG_INFO_MSG(
                        "GNC-TLS terminated because maximum residual at initialization is very small.");
                break;
            }
        }

        // Fix R and solve for weights in closed form
        double th1 = (mu + 1) / mu * noise_bound_sq;
        double th2 = mu / (mu + 1) * noise_bound_sq;
        cost_ = 0;
        for (size_t j = 0; j < match_size; ++j) {
            // Also calculate cost in this loop
            // Note: the cost calculated is using the previously solved weights
            cost_ += weights(j) * residuals_sq(j);

            if (residuals_sq(j) >= th1) {
                weights(j) = 0;
            } else if (residuals_sq(j) <= th2) {
                weights(j) = 1;
            } else {
                weights(j) = sqrt(noise_bound_sq * mu * (mu + 1) / residuals_sq(j)) - mu;
                assert(weights(j) >= 0 && weights(j) <= 1);
            }
            // - - - - - - - - - - - - - - - - - - - -
            if (SAVE_RESULTS) {
                if (j == match_size - 1) { // Last
                    weightTxt << weights(j) << "\n";
                } else {
                    weightTxt << weights(j) << " ";
                }
            }
            // - - - - - - - - - - - - - - - - - - - -
        }

        // Calculate cost
        double cost_diff = std::abs(cost_ - prev_cost);
        // - - - - - - - - - - - - - - - - - - - -
        if (SAVE_RESULTS) {
            for (int ii = 0; ii < 3; ++ii) {
                for (int jj = 0; jj < 3; ++jj) {
                    if ((ii == 2) && (jj == 2)) {
                        rotTxt << (*rotation)(ii, jj) << "\n";
                    } else {
                        rotTxt << (*rotation)(ii, jj) << " ";
                    }
                }
            }
            costTxt << cost_ << "\n";
        }
        // - - - - - - - - - - - - - - - - - - - -

        // Increase mu
        mu        = mu * params_.gnc_factor;
        prev_cost = cost_;

        if (cost_diff < params_.cost_threshold) {
            TEASER_DEBUG_INFO_MSG("GNC-TLS solver terminated due to cost convergence.");
            TEASER_DEBUG_INFO_MSG("Cost diff: " << cost_diff);
            TEASER_DEBUG_INFO_MSG("Iterations: " << i);
            break;
        }
    }

    // - - - - - - - - - - - - - - - - - - - -
    if (SAVE_RESULTS) {
        weightTxt.close();
        rotTxt.close();
        costTxt.close();
    }
    // - - - - - - - - - - - - - - - - - - - -

    if (inliers) {
        for (size_t i = 0; i < weights.cols(); ++i) {
            (*inliers)(0, i) = weights(0, i) >= 0.5;
        }
    }
}

// 2D VERSION
void teaser::GNCTLSRotationSolver::solveForRotation2D(                  // quatro는 이걸로
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &src,
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &dst, Eigen::Matrix2d*rotation,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>*inliers) {
    assert(rotation);                 // make sure R is not a nullptr
    assert(src.cols() == dst.cols()); // check dimensions of input data
    assert(params_.gnc_factor > 1);   // make sure mu will increase
    assert(params_.noise_bound != 0); // make sure noise sigma is not zero

    std::cout << "\033[1;32m=>Conducting 2D Rot solver w/ " << src.cols() << " \033[0m" << std::endl;

    if (inliers) {
        assert(inliers->cols() == src.cols());
    }
    /**
     * Loop: terminate when:
     *    1. the change in cost in two consecutive runs is smaller than a user-defined threshold
     *    2. # iterations exceeds the maximum allowed
     *
     * Within each loop:
     * 1. fix weights and solve for R
     * 2. fix R and solve for weights
     */

    std::string   costPath, estWeightPath, estRotPath;
    std::ofstream costTxt, weightTxt, rotTxt;
    if (SAVE_RESULTS) {
        costPath = (boost::format("%s/%04d_SONNY_cost.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        costTxt.open(costPath);
        estWeightPath =
                (boost::format("%s/%04d_SONNY_weights.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        weightTxt.open(estWeightPath);
        estRotPath = (boost::format("%s/%04d_SONNY_rot.txt") % SAVE_ABS_PATH % COUNT_EXECUTION).str();
        rotTxt.open(estRotPath);
    }
    // - - - - - - - - - - - - - - - - - - - -

    // Prepare some variables
    size_t match_size = src.cols(); // number of correspondences

    double mu = 1; // arbitrary starting mu

    double prev_cost = std::numeric_limits<double>::infinity();
    cost_ = std::numeric_limits<double>::infinity();
    //  double noise_bound_sq = std::pow(params_.noise_bound, 2);
    static double rot_noise_bound = params_.noise_bound;
    static double noise_bound_sq  = std::pow(rot_noise_bound, 2);
    if (noise_bound_sq < 1e-16) {
        noise_bound_sq = 1e-2;
    }
    TEASER_DEBUG_INFO_MSG("GNC rotation estimation noise bound:" << rot_noise_bound);
    //  TEASER_DEBUG_INFO_MSG("GNC rotation estimation noise bound:" << rot_noise_bound);
    TEASER_DEBUG_INFO_MSG("GNC rotation estimation noise bound squared:" << noise_bound_sq);


    Eigen::Matrix<double, 2, Eigen::Dynamic> diffs(2, match_size);
    Eigen::Matrix<double, 1, Eigen::Dynamic> weights(1, match_size);
    weights.setOnes(1, match_size);
    Eigen::Matrix<double, 1, Eigen::Dynamic> residuals_sq(1, match_size);

    // Loop for performing GNC-TLS
    for (size_t i = 0; i < params_.max_iterations; ++i) {

        // Fix weights and perform SVD rotation estimation
        *rotation = teaser::utils::svdRot2d(src, dst, weights);

        // Calculate residuals squared
        diffs        = (dst - (*rotation) * src).array().square();
        residuals_sq = diffs.colwise().sum();
        if (i == 0) {
            // Initialize rule for mu
            double max_residual = residuals_sq.maxCoeff();
            mu                  = 1 / (2 * max_residual / noise_bound_sq - 1);
            // Degenerate case: mu = -1 because max_residual is very small
            // i.e., little to none noise
            if (mu <= 0) {
                TEASER_DEBUG_INFO_MSG(
                        "GNC-TLS terminated because maximum residual at initialization is very small.");
                break;
            }
        }

        // Fix R and solve for weights in closed form
        double th1 = (mu + 1) / mu * noise_bound_sq;
        double th2 = mu / (mu + 1) * noise_bound_sq;
        cost_ = 0;
        for (size_t j = 0; j < match_size; ++j) {
            // Also calculate cost in this loop
            // Note: the cost calculated is using the previously solved weights
            cost_ += weights(j) * residuals_sq(j);

            if (residuals_sq(j) >= th1) {
                weights(j) = 0;
            } else if (residuals_sq(j) <= th2) {
                weights(j) = 1;
            } else {
                weights(j) = sqrt(noise_bound_sq * mu * (mu + 1) / residuals_sq(j)) - mu;
                assert(weights(j) >= 0 && weights(j) <= 1);
            }
            // - - - - - - - - - - - - - - - - - - - -
            if (SAVE_RESULTS) {
                if (j == match_size - 1) { // Last
                    weightTxt << weights(j) << "\n";
                } else {
                    weightTxt << weights(j) << " ";
                }
            }
            // - - - - - - - - - - - - - - - - - - - -
        }

        // Calculate cost
        double cost_diff = std::abs(cost_ - prev_cost);
        // - - - - - - - - - - - - - - - - - - - -
        if (SAVE_RESULTS) {
            for (int ii = 0; ii < 2; ++ii) {
                for (int jj = 0; jj < 2; ++jj) {
                    if ((ii == 1) && (jj == 1)) {
                        rotTxt << (*rotation)(ii, jj) << "\n";
                    } else {
                        rotTxt << (*rotation)(ii, jj) << " ";
                    }
                }
            }
            costTxt << cost_ << "\n";
        }
        // - - - - - - - - - - - - - - - - - - - -

        // Increase mu
        mu        = mu * params_.gnc_factor;
        prev_cost = cost_;

        if (cost_diff < params_.cost_threshold) {
            TEASER_DEBUG_INFO_MSG("GNC-TLS solver terminated due to cost convergence.");
            TEASER_DEBUG_INFO_MSG("Cost diff: " << cost_diff);
            TEASER_DEBUG_INFO_MSG("Iterations: " << i);
            break;
        }
    }

    // - - - - - - - - - - - - - - - - - - - -
    if (SAVE_RESULTS) {
        weightTxt.close();
        rotTxt.close();
        costTxt.close();
    }
    // - - - - - - - - - - - - - - - - - - - -

    if (inliers) {
        for (size_t i = 0; i < weights.cols(); ++i) {
            (*inliers)(0, i) = weights(0, i) >= 0.4;
        }
    }
}
