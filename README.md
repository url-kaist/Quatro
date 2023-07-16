# Quatro 
+ This branch is for being used as a module in other packages
+ Quatro is from [here, official repo](https://github.com/url-kaist/Quatro) and for the robust and global registration of pointclouds to avoid degeneracy in urban environments


### Dependencies
+ PCL >= 1.8
+ C++ >= 17
+ Boost >= 1.54
+ Eigen >= 3.2
+ OpenMP >= 4.5

### Use case
+ Refer - [here](https://github.com/engcang/FAST-LIO-SAM-QN)
+ Or, refer the example as follows:
```c++
quatro.setInputSource(srcMatched);
quatro.setInputTarget(tgtMatched);
Eigen::Matrix4d output;
quatro.computeTransformation(output);
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