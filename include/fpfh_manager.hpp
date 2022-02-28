//
// Created by Hyungtae Lim on 1/24/22.
//

#ifndef FPFH_MANAGER_H
#define FPFH_MANAGER_H


#include <unistd.h>
#include <geometry_msgs/Pose.h>
#include <iostream>
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include "conversion.hpp"
#include "teaser_utils/fpfh.h"
#include "teaser_utils/feature_matcher.h"
#include "boost/format.hpp"
#include <chrono>
#include <stdexcept>
#include "utility.h"

using namespace std;

class FPFHManager
{
private:
    teaser::FPFHCloudPtr obj_descriptors;
    teaser::FPFHCloudPtr scene_descriptors;

public:
    teaser::PointCloud src_cloud, tgt_cloud;
    double normal_radius_;
    double fpfh_radius_;
    bool is_initial_ = true;
    bool is_odometry_test_ = false;
    int interval_;
    std::string savedir_;
    std::string loaddir_;
    std::vector<std::pair<int, int>> corr; // Correspondence

    Eigen::Matrix<double, 3, Eigen::Dynamic> src_matched;
    Eigen::Matrix<double, 3, Eigen::Dynamic> tgt_matched;

    Eigen::Matrix<double, 3, Eigen::Dynamic> src_normals; // Not in use
    Eigen::Matrix<double, 3, Eigen::Dynamic> tgt_normals;
    pcl::PointCloud<PointType> src_matched_pcl;
    pcl::PointCloud<PointType> tgt_matched_pcl;

    FPFHManager(double normal_radius, double fpfh_radius, int interval=1): normal_radius_(normal_radius), fpfh_radius_(fpfh_radius), interval_(interval) {
        obj_descriptors.reset(new pcl::PointCloud<pcl::FPFHSignature33>());
        scene_descriptors.reset(new pcl::PointCloud<pcl::FPFHSignature33>());
        is_initial_ = true;
    }

    FPFHManager(){
        obj_descriptors.reset(new pcl::PointCloud<pcl::FPFHSignature33>());
        scene_descriptors.reset(new pcl::PointCloud<pcl::FPFHSignature33>());
        is_initial_ = true;
    }

    ~FPFHManager(){
        normal_radius_ = 0.5;
        fpfh_radius_ = 0.6;
        interval_ = 1;
        obj_descriptors.reset(new pcl::PointCloud<pcl::FPFHSignature33>());
        scene_descriptors.reset(new pcl::PointCloud<pcl::FPFHSignature33>());
        is_initial_ = true;
    }
    void flushAllFeatures(){
        is_initial_ = true;
    }

    void swapTgt2Src(){
        src_cloud = tgt_cloud;
        *obj_descriptors = *scene_descriptors;
    }
    void setParams(float normal_radius, float fpfh_radius, int interval){
        normal_radius_ = normal_radius;
        fpfh_radius_ = fpfh_radius;
        interval_ = interval;
    }

    void clearInputs(){
        is_initial_ = true;
        src_cloud.clear();
        tgt_cloud.clear();
        obj_descriptors->clear();
        scene_descriptors->clear();
    }
    void setLoadDir(std::string loaddir){
        loaddir_ = loaddir;
    }
    void setSaveDir(std::string savedir){
        savedir_ = savedir;
    }

    void setFeaturePair(pcl::PointCloud<PointType>::Ptr src, pcl::PointCloud<PointType>::Ptr target){
        if (normal_radius_ > fpfh_radius_){
            std::cout<<normal_radius_<<" <-> "<< fpfh_radius_<<std::endl;
            throw invalid_argument("[FPFHManager]: Normal should be lower than fpfh_radius!!!!");
        }

        pcl::PointCloud<pcl::Normal>::Ptr src_normals_raw(new pcl::PointCloud<pcl::Normal>);
        pcl::PointCloud<pcl::Normal>::Ptr tgt_normals_raw(new pcl::PointCloud<pcl::Normal>);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        // Set FPFH
        static teaser::FPFHEstimation fpfh;
        if (is_initial_ && !is_odometry_test_) {
            pcl2teaser(*src, src_cloud);
            obj_descriptors = fpfh.computeFPFHFeatures(src_cloud, *src_normals_raw, normal_radius_, fpfh_radius_);
            is_initial_ = false;
        }else{
            // To reduce computational cost on odometry test
            swapTgt2Src();
        }
        pcl2teaser(*target, tgt_cloud);

        scene_descriptors = fpfh.computeFPFHFeatures(tgt_cloud, *tgt_normals_raw, normal_radius_, fpfh_radius_);

        std::chrono::system_clock::time_point end_extraction = std::chrono::system_clock::now();

        teaser::Matcher matcher;
        corr = matcher.calculateCorrespondences(
                src_cloud, tgt_cloud, *obj_descriptors, *scene_descriptors, true, true, true, 0.95);

        std::chrono::system_clock::time_point end_matching = std::chrono::system_clock::now();

        // ---------------------------
        src_matched.resize(3, corr.size());
        tgt_matched.resize(3, corr.size());
        tgt_normals.resize(3, corr.size());

        for (size_t i = 0; i < corr.size(); ++i) {
            auto src_idx = std::get<0>(corr[i]);
            auto dst_idx = std::get<1>(corr[i]);
            src_matched.col(i) << src_cloud[src_idx].x, src_cloud[src_idx].y, src_cloud[src_idx].z;
            tgt_matched.col(i) << tgt_cloud[dst_idx].x, tgt_cloud[dst_idx].y, tgt_cloud[dst_idx].z;

            tgt_normals.col(i) << double(tgt_normals_raw->points[dst_idx].normal_x), double(tgt_normals_raw->points[dst_idx].normal_y), double(tgt_normals_raw->points[dst_idx].normal_z);
        }

        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::chrono::duration<double> sec_extraction = end_extraction - start;
        std::chrono::duration<double> sec_matching = end_matching - end_extraction;
        std::cout << setprecision(4) << "FPFH takes: " << sec.count() << " sec. ";
        std::cout << "(Extraction: " << sec_extraction.count() << " sec. + Matching: " << sec_matching.count() << " sec.)" << std::endl;

        eigen2pcl(src_matched, src_matched_pcl);
        eigen2pcl(tgt_matched, tgt_matched_pcl);
    }
    Eigen::Matrix<double, 3, Eigen::Dynamic> getSrcMatched(){
        return src_matched;
    }
    Eigen::Matrix<double, 3, Eigen::Dynamic> getTgtMatched(){
        return tgt_matched;
    }

    Eigen::Matrix<double, 3, Eigen::Dynamic> getTgtNormals(){
        return tgt_normals;
    }

    pcl::PointCloud<pcl::FPFHSignature33> getObjDescriptor(){
        return *obj_descriptors;
    }
    pcl::PointCloud<pcl::FPFHSignature33> getSceneDescriptor(){
        return *scene_descriptors;
    }

    pcl::PointCloud<PointType> getSrcKps(){
        return src_matched_pcl;
    }
    pcl::PointCloud<PointType> getTgtKps(){
        return tgt_matched_pcl;
    }

    void saveFeaturePair(int src_idx, int tgt_idx, bool verbose=false){
        if (savedir_.empty()){
            throw invalid_argument("Save dir. is not set");
        }
        std::string pcdname = (boost::format("%s/%06d_to_%06d.pcd") % savedir_ % src_idx %tgt_idx).str();
//        std::string txtname = (boost::format("%s/%06d_to_%06d.txt") % savedir_ % src_idx %tgt_idx).str();
        if (verbose){
            std::cout<<"[SAVER]: "<<pcdname<<std::endl;
            std::cout<<src_matched_pcl.points.size()<< " + " << tgt_matched_pcl.points.size()<<std::endl;
        }
        /***
         * Source is the first!!!
         */
        pcl::PointCloud<PointType> merge = src_matched_pcl + tgt_matched_pcl;
        pcl::io::savePCDFile(pcdname, merge);

//        ofstream num_check;
//        num_check.open (txtname);
//        num_check << src_matched_pcl.points.size() << " " << tgt_matched_pcl.points.size()<<std::endl;
//        num_check.close();
    }

    void loadFeaturePair(int src_idx, int tgt_idx, bool verbose=false){
        time_t start = clock();
        if (loaddir_.empty()){
            throw invalid_argument("Load dir. is not set");
        }
        std::string pcdname = (boost::format("%s/%06d_to_%06d.pcd") % loaddir_ % src_idx %tgt_idx).str();
        pcl::PointCloud<PointType> merge;

        if (pcl::io::loadPCDFile(pcdname, merge) == -1){
            throw invalid_argument("[FPFHManager]: Load feature set failed.");
        }

        if (verbose){
            std::cout<<"[LOADER]: Loaded data from "<<pcdname<<"..."<<std::endl;
            std::cout<<merge.points.size();
        }

        src_matched_pcl.clear();
        tgt_matched_pcl.clear();
        for (int i = 0; i< merge.points.size();++i){
            if (i < merge.points.size()/2){
                src_matched_pcl.push_back(merge.points[i]);
            }else{
                tgt_matched_pcl.push_back(merge.points[i]);
            }
        }
        time_t end = clock();
        if (verbose){
            std::cout<<"=>"<< src_matched_pcl.points.size()<<" " <<tgt_matched_pcl.points.size()<<std::endl;
            std::cout<<"Takes "<<(double)(end-start)/CLOCKS_PER_SEC<<"s..."<<std::endl;
        }
    }

    std::vector<std::pair<int, int>> getCorrespondences(){
        return corr;
    }

};

#endif //FPFH_MANAGER_H
