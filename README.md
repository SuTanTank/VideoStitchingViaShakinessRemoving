# Dynamic Video Stitching via Shakiness Removing

Demo code for the works below:
 - Nie, Yongwei, et al. "[Dynamic Video Stitching via Shakiness Removing](http:://ieeexplore.ieee.org/document/8003352/)." IEEE Transactions on Image Processing 27.1 (2018): 164-178.
 - Su, Tan, et al. "[Video stitching for handheld inputs via combined video stabilization](https://dl.acm.org/citation.cfm?id=3005383)." SIGGRAPH ASIA 2016 Technical Briefs. ACM, 2016.

[Video Demo](https://youtu.be/IDzFvqRb40Y) on YouTube:

[![Video Demo](https://img.youtube.com/vi/IDzFvqRb40Y/0.jpg)](https://youtu.be/IDzFvqRb40Y)

## Dataset
___Including video clips from previous works, see Readme.MD for more details)___

[OneDrive](https://1drv.ms/u/s!AuSZxO3eMeoIhvA9dUmE091Fw0QzrQ?e=UeaGU5)

[BaiduYun](https://pan.baidu.com/s/1bZd1k6#list/path=%2F)

***Notes on your own test data:** The proposed method is using traditional SIFT + Mesh-based warping to conduct stitching, which cannot handle the case where the parallax is too large, or the light condition is too different. To put it in a simple way, if a single pair of frames from the input video pair can't be stitched by Content-Preserve Warping, then this video stitching method will also fail.

## External libraries and code used:
 - Shankar Rao's Motion Segmentation Code: http://perception.csl.illinois.edu/coding/motion/#Software
 - CVX: http://cvxr.com/cvx/
 - vlfeat: http://www.vlfeat.org/
 - peter kovesi matlab toolbox: http://www.peterkovesi.com/matlabfns/
 - Liu Shuaicheng's As-similar-as-possible Warping code: http://www.liushuaicheng.org/SIGGRAPH2013/index.htm

## How to use this demo code:

1. In the folder `/case-cuhk_lib`, extract video frames of ***case17-l.mp4*** to folder `/left`, and extract video frames of ***case17-r.mp4*** to folder `/right`. After the frame extracation, each folder should contain 400 jpg files. The file names should be indexed properly. (e.g. ***001.jpg 002.jpg*** ...). A script `video2frames.m` is provided. Modify the filename

2. You may need to install `CVX` if you have not.

3. Set MATLAB path to `/Stitching-1.1.0`, run `RunStitching.m`. The generated output frames will be in the auto-created subfolder under `/case-cuhk_lib`. (`res_demo` if you didn't change the output path). **Warning:** the motion segmentation code from Shankar Rao is very slow, if your video doesn't have large foreground object, set `SKIP_BACKGROUND_SEGMENTATION = true` in `RunStitching.m` and treat all trajectories as background. For further details, see the code comments.

4. Build the OpenCV project in `/SeamCut` (You need to set OpenCV's include and library path manually) and copy the executable (e.g. `SeamCut.exe`) to the folder containing the output frames (***1.jpg***, ***A1.jpg***, ***B1.jpg***, ***2.jpg***, ***A2.jpg***, ***B2.jpg***, ...). This program finds the continuous optimal seam by GraphCut algorithm and use OpenCV's multi-band blending function.

5. Run `./SeamCut 5 1 400 1 0.2` to see the final result. Blended frames are saved as ***D1.jpg***, ***D2.jpg***, ...

For more details, please read the comments

## Please cite our papers if you use the code or our dataset.

@article{nie2018dynamic,
  title={Dynamic Video Stitching via Shakiness Removing},
  author={Nie, Yongwei and Su, Tan and Zhang, Zhensong and Sun, Hanqiu and Li, Guiqing},
  journal={IEEE Transactions on Image Processing},
  volume={27},
  number={1},
  pages={164--178},
  year={2018},
  publisher={IEEE}
}

@inproceedings{su2016video,
  title={Video stitching for handheld inputs via combined video stabilization},
  author={Su, Tan and Nie, Yongwei and Zhang, Zhensong and Sun, Hanqiu and Li, Guiqing},
  booktitle={SIGGRAPH ASIA 2016 Technical Briefs},
  pages={25},
  year={2016},
  organization={ACM}
}

Code by Tan Su, Zhensong Zhang and Yongwei Nie. For research purpose ONLY.
[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu)
