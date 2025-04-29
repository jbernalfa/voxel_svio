# Voxel-SVIO

**Voxel-SVIO** (Voxel-SVIO: Stereo Visual-Inertial Odometry based on Voxel Map) is a MSCKF based visual-inertial odometry (VIO) that utilizes voxel-based map management. The voxel-based map management enables VIO systems to efficiently retrieve the most suitable points for optimizer inclusion, thereby ensuring optimal allocation of computational resources to the variables most critical for optimization.

## Related Work

Voxel-SVIO: Stereo Visual-Inertial Odometry based on Voxel Map

Authors: [*Zikang Yuan*](https://scholar.google.com/citations?hl=zh-CN&user=acxdM9gAAAAJ), [*Fengtian Lang*](https://scholar.google.com/citations?hl=zh-CN&user=zwgGSkEAAAAJ&view_op=list_works&gmla=ABEO0Yrl4-YPuowyntSYyCW760yxM5-IWkF8FGV4t9bs9qz1oWrqnlHmPdbt7LMcMDc04kl2puqRR4FaZvaCUONsX7MQhuAC6a--VS2pTsuwj-CyKgWp3iWDP2TS0I__Zui5da4), *Jie Deng*, *Hongcheng Luo* and [*Xin Yang*](https://scholar.google.com/citations?user=lsz8OOYAAAAJ&hl=zh-CN)

## Installation

### 1. Requirements

> GCC >= 7.5.0
>
> Cmake >= 3.16.0
> 
> [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.3.4
>
> [OpenCV](https://github.com/opencv/opencv) == 4.2.0
>
> [PCL](https://pointclouds.org/downloads/) == 1.8 for Ubuntu 18.04, and == 1.10 for Ubuntu 20.04
>
> [Ceres](http://ceres-solver.org/installation.html) >= 1.14
>
> [ROS](http://wiki.ros.org/ROS/Installation)
>
> [Livox_ROS_Driver](https://github.com/Livox-SDK/livox_ros_driver)

##### Have Tested On:

| OS    | GCC  | Cmake | Eigen3 | OpenCV | Ceres |
|:-:|:-:|:-:|:-:|:-:|:-:|
| Ubuntu 20.04 | 9.4.0  | 3.16.3 | 3.3.7 | 4.2.0 | 1.14 |

### 2. Create ROS workspace

```bash
mkdir -p ~/Voxel-SVIO/src
cd Voxel-SVIO/src
```

### 3. Clone the directory and build

```bash
git clone https://github.com/ZikangYuan/sr_livo.git
cd ..
catkin_make
```
