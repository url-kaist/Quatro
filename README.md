# Quatro 
+ This branch is for being used as a third party module in other packages
+ `Quatro` is from [here, official repo](https://github.com/url-kaist/Quatro) and [here, author's module repository](https://github.com/LimHyungTae/quatro-cpp-fpfh/)
    + `Quatro` is for the robust and global registration of pointclouds to avoid degeneracy in urban environments


### Dependencies
+ `PCL` >= 1.8
+ `C++` >= 17
+ `Boost` >= 1.54
+ `Eigen` >= 3.2
+ `OpenMP` >= 4.5
+ `Teaser++` (tested with commit ver `e415c0d`, May 22, 2023)
    ```shell
    git clone https://github.com/MIT-SPARK/TEASER-plusplus.git
    cd TEASER-plusplus && mkdir build && cd build
    cmake .. -DENABLE_DIAGNOSTIC_PRINT=OFF
    sudo make install -j16
    sudo ldconfig
    ```
+ (Optional, but recommended) `tbb` for parallelization and ***x10 faster computation in matching***
    ```shell
    sudo apt install libtbb-dev
    ```

### How to build
+ Git clone and `catkin build` this repository
```shell
cd ~/your_workspace/src
git clone https://github.com/engcang/Quatro
cd ..
catkin build
```
+ (Optional, but recommended) with `tbb`
```shell
catkin build -DQUATRO_TBB=ON
```
+ (Optional) to enable `std::cout` for debugging,
```shell
catkin build -DQUATRO_DEBUG=ON
```

### Use case
+ Refer - [here](https://github.com/engcang/FAST-LIO-SAM-QN)
+ Or, refer the example as follows:
    1. Make sure that you have all dependencies and built the `quatro` properly with `catkin build` as above
    2. In the `CMakeLists.txt` of your wanted package, import `quatro` as a component of `catkin`
        ```CMake
        find_package(catkin REQUIRED COMPONENTS
            ...
            quatro #Include here
            ...
        )
        include_directories(
            ...
            ${catkin_INCLUDE_DIRS} #Header files are included in catkin_INCLUDE_DIRS
            ...
        )
        add_library(some_library src/some_src.cpp)
        target_link_libraries(some_library ${catkin_LIBRARIES}) #Libraries are included in catkin_LIBRARIES
        ```
    4. Use in the source file of your wanted package as:
        ```c++
        #include <quatro/quatro_module.h>
        using QuatroPointType = pcl::PointXYZI; //can be changed

        shared_ptr<quatro<QuatroPointType>> m_quatro_handler = nullptr;
        m_quatro_handler = std::make_shared<quatro<QuatroPointType>>(fpfh_normal_radius_, 
                            fpfh_radius_, noise_bound_, rot_gnc_factor_, rot_cost_diff_thr_,
                            quatro_max_iter_, estimat_scale_, optimize_matching, 
                            distance_threshold, max_correspondences); //refer https://github.com/engcang/FAST-LIO-SAM-QN/blob/master/fast_lio_sam_qn/config/config.yaml#L28

        ////// use
        pcl::PointCloud<QuatroPointType> src_; //this should be not empty but the real data
        pcl::PointCloud<QuatroPointType> dst_; //this should be not empty but the real data
        bool if_valid_;
        Eigen::Matrix4d output_tf_ = m_quatro_handler->align(src_, dst_, if_valid_);
        if (if_valid_)
        {
            //do something wiht output_tf_
        }
        ```

### License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>
- This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/)

### Acknowledgement
- This work was supported by the Industry Core Technology Development Project, 20005062, Development of Artificial Intelligence Robot Autonomous Navigation Technology for Agile Movement in Crowded Space, funded by the Ministry of Trade, Industry & Energy (MOTIE, Republic of Korea) and by the research project “Development of A.I. based recognition, judgement and control solution for autonomous vehicle corresponding to atypical driving environment,” which is financed from the Ministry of Science and ICT (Republic of Korea) Contract No. 2019-0-00399. The student was supported by the BK21 FOUR from the Ministry of Education (Republic of Korea).

### Copyright
- All codes on this page are copyrighted by KAIST published under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 License. You must attribute the work in the manner specified by the author. You may not use the work for commercial purposes, and you may only distribute the resulting work under the same license if you alter, transform, or create the work.

### Citation
- If our research has been helpful, please cite the below paper:
```tex
@article{lim2022quatro,
    title={A Single Correspondence Is Enough: Robust Global Registration to Avoid Degeneracy in Urban Environments},
    author={Lim, Hyungtae and Yeon, Suyong and Ryu, Suyong and Lee, Yonghan and Kim, Youngji and Yun, Jaeseong and Jung, Euigon and Lee, Donghwan and Myung, Hyun},
    booktitle={Proc. IEEE Int. Conf. Robot. Autom.},
    year={2022},
    pages={Accepted. To appear}
    }
```
```tex
@article{lim2023quatro-plusplus,
  author    = {Lim, Hyungtae and Kim, Beomsoo and Kim, Daebeom and Lee, Eungchang and Myung, Hyun},
  title     = {Quatro++: Robust Global Registration Exploiting Ground Segmentation for Loop Closing in LiDAR SLAM},
  booktitle = {International Journal of Robotics Research},
  year      = {2023},
}
```