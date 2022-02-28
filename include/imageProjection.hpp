
#include "utility.h"
#include "conversion.hpp"

#define NOT_ASSIGNED -1
#define USE_FOUR_NEIGHBORS_PIXELS 0
#define USE_EIGHT_NEIGHBORS_PIXELS 1
#define USE_FOUR_CROSS_NEIGHBORS_PIXELS 2

class ImageProjection {
private:
    pcl::PointCloud<PointTypeIP>::Ptr laserCloudIn;
    pcl::PointCloud<PointTypeIP>::Ptr fullCloud; // projected velodyne raw cloud, but saved in the form of 1-D matrix
    pcl::PointCloud<PointTypeIP>::Ptr fullInfoCloud; // same as fullCloud, but with intensity is replaced with range

    pcl::PointCloud<PointTypeIP>::Ptr groundCloud;
    pcl::PointCloud<PointTypeIP>::Ptr segmentedCloud;
    pcl::PointCloud<PointTypeIP>::Ptr segmentedCloudPure;
    pcl::PointCloud<PointTypeIP>::Ptr outlierCloud;

    PointTypeIP nanPoint; // fill in fullCloud at each iteration

    cv::Mat rangeMat; // range matrix for range image
    cv::Mat labelMat; // label matrix for segmentaiton marking
    cv::Mat groundMat; // ground matrix for ground cloud marking

    int labelCount;
    int numMinPtsForSubclustering = 30;

    float startOrientation;
    float endOrientation;

    quatro::cloud_info segMsg; // info of segmented cloud
    std_msgs::Header       cloudHeader;

    std::vector<std::pair<int8_t, int8_t> > neighborIterator; // neighbor iterator for segmentaiton process

    uint16_t *allPushedIndX; // array for tracking points of a segmented object
    uint16_t *allPushedIndY;

    uint16_t *queueIndX; // array for breadth-first search process of segmentation
    uint16_t *queueIndY;

public:
    int neighborPixelsSelectionMode = NOT_ASSIGNED;
    string groundSegMode;

    int N_SCAN;
    int Horizon_SCAN; //1028~4500
    float ang_res_x;
    float ang_res_y; //28.0/float(N_SCAN-1);
    float ang_bottom;
    int groundScanInd; 

    float segmentAlphaX;
    float segmentAlphaY;

    const float sensorMountAngle = 0.0;
    //segmentation threshold
    const float segmentTheta = 60.0/180.0*M_PI; // decrese this value may improve accuracy
    //If number of segment is below than 30, check line number. this for minimum number of point for it
    const int   segmentValidPointNum = 5;
    //if number of segment is small, number of line is checked, this is threshold for it.
    const int   segmentValidLineNum = 3;
    

    ImageProjection() {
        nanPoint.x         = std::numeric_limits<float>::quiet_NaN();
        nanPoint.y         = std::numeric_limits<float>::quiet_NaN();
        nanPoint.z         = std::numeric_limits<float>::quiet_NaN();
        nanPoint.intensity = -1;

        allocateMemory();
        resetParameters();
    }

    ImageProjection(
            string lidarType, string neighborSelectionMode, string groundSegmentationMode, int numSubclusteringCriteria = 30) {
        nanPoint.x         = std::numeric_limits<float>::quiet_NaN();
        nanPoint.y         = std::numeric_limits<float>::quiet_NaN();
        nanPoint.z         = std::numeric_limits<float>::quiet_NaN();
        nanPoint.intensity = -1;

        std::cout << "[Lidar Type]: " << lidarType << std::endl;
        if(lidarType == "Velodyne-64-HDE"){
            N_SCAN = 64;
            Horizon_SCAN = 1800; //1028~4500
            ang_res_x = 360.0/float(Horizon_SCAN);
            ang_res_y = 26.9/float(N_SCAN-1); //28.0/float(N_SCAN-1);
            ang_bottom = 25.0;
            groundScanInd = 60; 
        }
        else if(lidarType == "VLP-16"){
            N_SCAN = 16;
            Horizon_SCAN = 1800;
            ang_res_x = 0.2;
            ang_res_y = 2.0;
            ang_bottom = 15.0+0.1;
            groundScanInd = 7; 
        }
        else if(lidarType == "HDL-32E"){
            N_SCAN = 32;
            Horizon_SCAN = 1800;
            ang_res_x = 360.0/float(Horizon_SCAN);
            ang_res_y = 41.33/float(N_SCAN-1);
            ang_bottom = 30.67;
            groundScanInd = 20; 
        }
        else if(lidarType == "Ouster-OS1-16"){
            N_SCAN = 16;
            Horizon_SCAN = 1024;
            ang_res_x = 360.0/float(Horizon_SCAN);
            ang_res_y = 33.2/float(N_SCAN-1);
            ang_bottom = 16.6+0.1;
            groundScanInd = 7; 
        }
        else if(lidarType == "Ouster-OS1-64"){
            N_SCAN = 64;
            Horizon_SCAN = 1024;
            ang_res_x = 360.0/float(Horizon_SCAN);
            ang_res_y = 33.2/float(N_SCAN-1);
            ang_bottom = 16.6+0.1;
            groundScanInd = 15;
        }        
        else throw invalid_argument("[ImageProjection]:Check your paramter. Lidar Type is wrong!");
        segmentAlphaX = ang_res_x / 180.0 * M_PI;
        segmentAlphaY = ang_res_y / 180.0 * M_PI;

        std::cout << "[NEIGHBOR Mode]: " << neighborSelectionMode << std::endl;
        if (neighborSelectionMode == "4Neighbor") neighborPixelsSelectionMode = USE_FOUR_NEIGHBORS_PIXELS;
        else if (neighborSelectionMode == "8Neighbor") neighborPixelsSelectionMode = USE_EIGHT_NEIGHBORS_PIXELS;
        else if (neighborSelectionMode == "4CrossNeighbor") neighborPixelsSelectionMode = USE_FOUR_CROSS_NEIGHBORS_PIXELS;
        else throw invalid_argument("[ImageProjection]:Check your paramter. Neighbor selection mode is wrong!");

        groundSegMode = groundSegmentationMode;
        if (!(groundSegMode == "LeGO-LOAM") && !(groundSegMode == "Patchwork") ) {
            throw invalid_argument("[ImageProjection]: Check your paramter. Ground Segmentation mode is wrong!");
        }

        numMinPtsForSubclustering = numSubclusteringCriteria;


        allocateMemory();
        resetParameters();
    }

    void allocateMemory() {

        laserCloudIn.reset(new pcl::PointCloud<PointTypeIP>());

        fullCloud.reset(new pcl::PointCloud<PointTypeIP>());
        fullInfoCloud.reset(new pcl::PointCloud<PointTypeIP>());

        groundCloud.reset(new pcl::PointCloud<PointTypeIP>());
        segmentedCloud.reset(new pcl::PointCloud<PointTypeIP>());
        segmentedCloudPure.reset(new pcl::PointCloud<PointTypeIP>());
        outlierCloud.reset(new pcl::PointCloud<PointTypeIP>());

        fullCloud->points.resize(N_SCAN * Horizon_SCAN);
        fullInfoCloud->points.resize(N_SCAN * Horizon_SCAN);

        segMsg.startRingIndex.assign(N_SCAN, 0);
        segMsg.endRingIndex.assign(N_SCAN, 0);

        segMsg.segmentedCloudGroundFlag.assign(N_SCAN * Horizon_SCAN, false);
        segMsg.segmentedCloudColInd.assign(N_SCAN * Horizon_SCAN, 0);
        segMsg.segmentedCloudRange.assign(N_SCAN * Horizon_SCAN, 0);

        // 4 neigbor labeling
        cout << "=> neighborPixelsSelectionMode: ";
        if (neighborPixelsSelectionMode == USE_FOUR_NEIGHBORS_PIXELS) {
            cout << "USE_FOUR_NEIGHBORS_PIXELS" << endl;
            neighborIterator.push_back({-1, 0});
            neighborIterator.push_back({0, 1});
            neighborIterator.push_back({0, -1});
            neighborIterator.push_back({1, 0});

        } else if (neighborPixelsSelectionMode == USE_EIGHT_NEIGHBORS_PIXELS) {
            cout << "USE_EIGHT_NEIGHBORS_PIXELS" << endl;
            neighborIterator.push_back({-1, 0});
            neighborIterator.push_back({0, 1});
            neighborIterator.push_back({0, -1});
            neighborIterator.push_back({1, 0});
            neighborIterator.push_back({-1, -1});
            neighborIterator.push_back({-1, 1});
            neighborIterator.push_back({1, 1});
            neighborIterator.push_back({1, -1});

        } else if (neighborPixelsSelectionMode == USE_FOUR_CROSS_NEIGHBORS_PIXELS) {
            cout << "USE_FOUR_CROSS_NEIGHBORS_PIXELS" << endl;
            neighborIterator.push_back({-1, -1});
            neighborIterator.push_back({-1, 1});
            neighborIterator.push_back({1, 1});
            neighborIterator.push_back({1, -1});
        } else { cout << endl; }

        allPushedIndX = new uint16_t[N_SCAN * Horizon_SCAN];
        allPushedIndY = new uint16_t[N_SCAN * Horizon_SCAN];

        queueIndX = new uint16_t[N_SCAN * Horizon_SCAN];
        queueIndY = new uint16_t[N_SCAN * Horizon_SCAN];
    }

    void resetParameters() {
        laserCloudIn->clear();
        groundCloud->clear();
        segmentedCloud->clear();
        segmentedCloudPure->clear();
        outlierCloud->clear();

        rangeMat   = cv::Mat(N_SCAN, Horizon_SCAN, CV_32F, cv::Scalar::all(FLT_MAX));
        groundMat  = cv::Mat(N_SCAN, Horizon_SCAN, CV_8S, cv::Scalar::all(0));
        labelMat   = cv::Mat(N_SCAN, Horizon_SCAN, CV_32S, cv::Scalar::all(0));
        labelCount = 1;

        std::fill(fullCloud->points.begin(), fullCloud->points.end(), nanPoint);
        std::fill(fullInfoCloud->points.begin(), fullInfoCloud->points.end(), nanPoint);
    }

    void getValidSegments(pcl::PointCloud<PointXYZ>& output) {
        pcl::copyPointCloud(*segmentedCloudPure, output);
    }

    void getValidSegments(pcl::PointCloud<pcl::PointXYZI>& output) {
        output = *segmentedCloudPure;
    }

    void getGround(pcl::PointCloud<pcl::PointXYZ>& output) {
        pcl::copyPointCloud(*groundCloud, output);
    }

    void getGround(pcl::PointCloud<pcl::PointXYZI>& output) {
        output = *groundCloud;
    }

    void getOutliers(pcl::PointCloud<pcl::PointXYZ>& output) {
        pcl::copyPointCloud(*outlierCloud, output);
    }

    void getOutliers(pcl::PointCloud<pcl::PointXYZI>& output) {
        output = *outlierCloud;
    }

    cv::Mat getrangeMat() {
        return rangeMat;
    }


    ~ImageProjection() {}

    void setPointCloud(const pcl::PointCloud<PointType>::Ptr srcPtr) {
        (*laserCloudIn).points.resize((*srcPtr).size());
        for (size_t i = 0; i < (*srcPtr).size(); i++) {
            (*laserCloudIn).points[i].x = (*srcPtr).points[i].x;
            (*laserCloudIn).points[i].y = (*srcPtr).points[i].y;
            (*laserCloudIn).points[i].z = (*srcPtr).points[i].z;
        }
    }

    void copyPointCloud(const pcl::PointCloud<PointType>::Ptr srcPtr) {
        setPointCloud(srcPtr);

        // Remove Nan points
        std::vector<int> indices;
        pcl::removeNaNFromPointCloud(*laserCloudIn, *laserCloudIn, indices);
    }

    // (Ground Removal + Sub Clustering)
    /***
     * Ground Removal + Subclustering
     * @param pcPtr
     */
    void segmentCloud(const pcl::PointCloud<PointType>::Ptr pcPtr) {
        // 1. Copy pcl point cloud
        copyPointCloud(pcPtr);
        // 2. Start and end angle of a scan
        findStartEndAngle();
        // 3. Range image projection
        projectPointCloud();
        // 4. Mark ground points
        /***
         * If groundSegMode == "Patchwork",
         * then the ground points are already rejected by Patchwork
         * i.e. srcPtr does not contain ground points
         * Thus, only masking labelMat is performed
         */
        if (groundSegMode == "LeGO-LOAM") {
            groundRemoval();
        } else if (groundSegMode == "Patchwork") {
            maskGround();
        }
        // 5. Point cloud segmentation
        cloudSegmentation();
    }

    void findStartEndAngle() {
        // start and end orientation of this cloud
        segMsg.startOrientation = -atan2(laserCloudIn->points[0].y, laserCloudIn->points[0].x);
        segMsg.endOrientation   = -atan2(laserCloudIn->points[laserCloudIn->points.size() - 1].y,
                                         laserCloudIn->points[laserCloudIn->points.size() - 1].x) + 2 * M_PI;
        if (segMsg.endOrientation - segMsg.startOrientation > 3 * M_PI) {
            segMsg.endOrientation -= 2 * M_PI;
        } else if (segMsg.endOrientation - segMsg.startOrientation < M_PI)
            segMsg.endOrientation += 2 * M_PI;
        segMsg.orientationDiff  = segMsg.endOrientation - segMsg.startOrientation;
    }

    void projectPointCloud() {
        // range image projection
        float       verticalAngle, horizonAngle, range;
        size_t      rowIdn, columnIdn, index, cloudSize;
        PointTypeIP thisPoint;

        cloudSize = laserCloudIn->points.size();

        for (size_t i = 0; i < cloudSize; ++i) {

            thisPoint.x = laserCloudIn->points[i].x;
            thisPoint.y = laserCloudIn->points[i].y;
            thisPoint.z = laserCloudIn->points[i].z;

            // find the row and column index in the iamge for this point
            verticalAngle = atan2(thisPoint.z, sqrt(thisPoint.x * thisPoint.x + thisPoint.y * thisPoint.y)) * 180 / M_PI;
            rowIdn        = (verticalAngle + ang_bottom) / ang_res_y;

            if (rowIdn < 0 || rowIdn >= N_SCAN)
                continue;

            horizonAngle = atan2(thisPoint.x, thisPoint.y) * 180 / M_PI;
            columnIdn    = -round((horizonAngle - 90.0) / ang_res_x) + Horizon_SCAN / 2;
            if (columnIdn >= Horizon_SCAN)
                columnIdn -= Horizon_SCAN;

            if (columnIdn < 0 || columnIdn >= Horizon_SCAN)
                continue;

            range = sqrt(thisPoint.x * thisPoint.x + thisPoint.y * thisPoint.y + thisPoint.z * thisPoint.z);
            if (range < 0.1)
                continue;

            rangeMat.at<float>(rowIdn, columnIdn) = range;

            // ???
            thisPoint.intensity = (float) rowIdn + (float) columnIdn / 10000.0;

            index = columnIdn + rowIdn * Horizon_SCAN;
            fullCloud->points[index]     = thisPoint;
            fullInfoCloud->points[index] = thisPoint;
            fullInfoCloud->points[index].intensity = range; // the corresponding range of a point is saved as "intensity"
        }

    }

    void maskGround() {
        // The pixels where the points are not projected are masked (i.e. rangeMat.at<float>(i, j) == FLT_MAX)
        for (size_t i = 0; i < N_SCAN; ++i) {
            for (size_t j = 0; j < Horizon_SCAN; ++j) {
                if (groundMat.at<int8_t>(i, j) == 1 || rangeMat.at<float>(i, j) == FLT_MAX) {
                    labelMat.at<int>(i, j) = -1;
                }
            }
        }
    }

    void groundRemoval() {
        size_t      lowerInd, upperInd;
        float       diffX, diffY, diffZ, angle;
        // groundMat
        // -1, no valid info to check if ground of not
        //  0, initial value, after validation, means not ground
        //  1, ground
        for (size_t j = 0; j < Horizon_SCAN; ++j) {
            for (size_t i = 0; i < groundScanInd; ++i) {

                // lowerInd : 현재 채널에서, 수평한 방향으로 포인트 인덱스
                // upperInd : lowerInd에서 z방향 바로 위 포인트 인덱스
                lowerInd = j + (i) * Horizon_SCAN;
                upperInd = j + (i + 1) * Horizon_SCAN;

                if (fullCloud->points[lowerInd].intensity == -1 ||
                    fullCloud->points[upperInd].intensity == -1) {
                    // no info to check, invalid points
                    groundMat.at<int8_t>(i, j) = -1;
                    continue;
                }

                // lowerInd에서 upperInd까지 이루는 각도 : angle
                diffX = fullCloud->points[upperInd].x - fullCloud->points[lowerInd].x;
                diffY = fullCloud->points[upperInd].y - fullCloud->points[lowerInd].y;
                diffZ = fullCloud->points[upperInd].z - fullCloud->points[lowerInd].z;
                angle = atan2(diffZ, sqrt(diffX * diffX + diffY * diffY)) * 180 / M_PI;

                // angle이 theshold를 넘지 못하면 ground로 판단.
                if (abs(angle - sensorMountAngle) <= 10) {
                    groundMat.at<int8_t>(i, j)     = 1;
                    groundMat.at<int8_t>(i + 1, j) = 1;
                }
            }
        }

        // extract ground cloud (groundMat == 1)
        // mark entry that doesn't need to label (ground and invalid point) for segmentation
        // note that ground remove is from 0~N_SCAN-1, need rangeMat for mark label matrix for the 16th scan

        // 거리가 매우 먼 포인트와 ground 판단되는 포인트들은 invalid로 설정
        for (size_t i = 0; i < N_SCAN; ++i) {
            for (size_t j = 0; j < Horizon_SCAN; ++j) {
                if (groundMat.at<int8_t>(i, j) == 1 || rangeMat.at<float>(i, j) == FLT_MAX) {
                    labelMat.at<int>(i, j) = -1;
                }
            }
        }

        // if (pubGroundCloud.getNumSubscribers() != 0){
        for (size_t i = 0; i <= groundScanInd; ++i) {
            for (size_t j = 0; j < Horizon_SCAN; ++j) {
                if (groundMat.at<int8_t>(i, j) == 1)
                    groundCloud->push_back(fullCloud->points[j + i * Horizon_SCAN]);
            }
        }
        // }
    }

    void cloudSegmentation() {
        // segmentation process
        for (size_t i = 0; i < N_SCAN; ++i)
            for (size_t j = 0; j < Horizon_SCAN; ++j)
                if (labelMat.at<int>(i, j) == 0)
                    labelComponents(i, j);

        // label: 999999 -> invalid
        // label: over 0 -> valid

        int         sizeOfSegCloud = 0;
        // extract segmented cloud for lidar odometry
        for (size_t i              = 0; i < N_SCAN; ++i) {

            segMsg.startRingIndex[i] = sizeOfSegCloud - 1 + 5;

            for (size_t j          = 0; j < Horizon_SCAN; ++j) {
                if (labelMat.at<int>(i, j) > 0 || groundMat.at<int8_t>(i, j) == 1) {
                    // outliers that will not be used for optimization (always continue)
                    if (labelMat.at<int>(i, j) == 999999) {
                        if (j % 1 == 0) {    //  i > groundScanInd  && j % 5 == 0
                            outlierCloud->push_back(fullCloud->points[j + i * Horizon_SCAN]);
                            continue;
                        } else {
                            continue;
                        }
                    }
                    // majority of ground points are skipped
                    if (groundMat.at<int8_t>(i, j) == 1) {
                        if (j % 5 != 0 && j > 5 && j < Horizon_SCAN - 5)
                            continue;
                    }
                    // mark ground points so they will not be considered as edge features later
                    segMsg.segmentedCloudGroundFlag[sizeOfSegCloud] = (groundMat.at<int8_t>(i, j) == 1);
                    // mark the points' column index for marking occlusion later
                    segMsg.segmentedCloudColInd[sizeOfSegCloud]     = j;
                    // save range info
                    segMsg.segmentedCloudRange[sizeOfSegCloud]      = rangeMat.at<float>(i, j);
                    // save seg cloud
                    segmentedCloud->push_back(fullCloud->points[j + i * Horizon_SCAN]);
                    // size of seg cloud
                    ++sizeOfSegCloud;
                }
            }
            segMsg.endRingIndex[i] = sizeOfSegCloud - 1 - 5;
        }

        // extract segmented cloud for visualization
        // 해당토픽을 subscribe하지 않으면, 처리하지 않고 publish안함(아래코드)
        // if (pubSegmentedCloudPure.getNumSubscribers() != 0){
        for (size_t i = 0; i < N_SCAN; ++i) {
            for (size_t j = 0; j < Horizon_SCAN; ++j) {
                if (labelMat.at<int>(i, j) > 0 && labelMat.at<int>(i, j) != 999999) {        //
                    segmentedCloudPure->push_back(fullCloud->points[j + i * Horizon_SCAN]);
                    segmentedCloudPure->points.back().intensity = labelMat.at<int>(i, j);
                }
            }
        }
        // }
    }

    void labelComponents(int row, int col) {
        // use std::queue std::vector std::deque will slow the program down greatly
        float d1, d2, alpha, angle;
        int   fromIndX, fromIndY, thisIndX, thisIndY;
        bool  lineCountFlag[N_SCAN] = {false};

        queueIndX[0] = row;
        queueIndY[0] = col;
        int queueSize     = 1;
        int queueStartInd = 0;
        int queueEndInd   = 1;

        allPushedIndX[0] = row;
        allPushedIndY[0] = col;
        int allPushedIndSize = 1;

        while (queueSize > 0) {
            // Pop point
            fromIndX = queueIndX[queueStartInd];
            fromIndY = queueIndY[queueStartInd];
            --queueSize;
            ++queueStartInd;
            // Mark popped point
            labelMat.at<int>(fromIndX, fromIndY) = labelCount;
            // Loop through all the neighboring grids of popped grid
            for (auto iter = neighborIterator.begin(); iter != neighborIterator.end(); ++iter) {
                // new index
                thisIndX = fromIndX + (*iter).first;
                thisIndY = fromIndY + (*iter).second;
                // index should be within the boundary
                if (thisIndX < 0 || thisIndX >= N_SCAN)
                    continue;
                // at range image margin (left or right side)
                if (thisIndY < 0)
                    thisIndY = Horizon_SCAN - 1;
                if (thisIndY >= Horizon_SCAN)
                    thisIndY = 0;
                // prevent infinite loop (caused by put already examined point back)
                if (labelMat.at<int>(thisIndX, thisIndY) != 0)
                    continue;


                //두포인트의 range 값
                d1 = std::max(rangeMat.at<float>(fromIndX, fromIndY),
                              rangeMat.at<float>(thisIndX, thisIndY));
                d2 = std::min(rangeMat.at<float>(fromIndX, fromIndY),
                              rangeMat.at<float>(thisIndX, thisIndY));

                if ((*iter).first == 0)
                    alpha = segmentAlphaX;      // velodyne64 : segmentAlphaX = 0.2deg
                else
                    alpha = segmentAlphaY;

                // beta값
                angle = atan2(d2 * sin(alpha), (d1 - d2 * cos(alpha)));

                if (angle > segmentTheta) {

                    queueIndX[queueEndInd] = thisIndX;
                    queueIndY[queueEndInd] = thisIndY;
                    ++queueSize;
                    ++queueEndInd;

                    labelMat.at<int>(thisIndX, thisIndY) = labelCount;
                    lineCountFlag[thisIndX] = true;

                    allPushedIndX[allPushedIndSize] = thisIndX;
                    allPushedIndY[allPushedIndSize] = thisIndY;
                    ++allPushedIndSize;
                }
            }
        }
        // check if this segment is valid
        bool feasibleSegment = false;
        /***
         * once the # of allPushedIndSize is smaller than numMinPtsForSubclustering,
         * it is considered as outliers (rejected)
         */
        if (allPushedIndSize >= numMinPtsForSubclustering)
            feasibleSegment = true;
        else if (allPushedIndSize >= segmentValidPointNum) {
            int         lineCount = 0;
            for (size_t i         = 0; i < N_SCAN; ++i)
                if (lineCountFlag[i] == true)
                    ++lineCount;
            if (lineCount >= segmentValidLineNum)
                feasibleSegment = true;
        }
        // segment is valid, mark these points
        if (feasibleSegment == true) {
            ++labelCount;
        } else { // segment is invalid, mark these points
            for (size_t i = 0; i < allPushedIndSize; ++i) {
                labelMat.at<int>(allPushedIndX[i], allPushedIndY[i]) = 999999;
            }
        }
    }
};