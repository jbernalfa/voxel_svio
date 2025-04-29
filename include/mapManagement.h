#pragma once

// c++ include
#include <iostream>
#include <memory>
#include <mutex>
#include <vector>

// lib include
#include <Eigen/Core>

// function include
#include "mapPoint.h"

class mapManagement
{
public:

  static bool addPointToVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr, double voxel_size, int max_num_points_in_voxel, double min_distance_points);

  static void changeHostVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr, double voxel_size, int max_num_points_in_voxel, double min_distance_points);

  static void deleteFromVoxel(voxelHashMap &voxel_map, std::shared_ptr<mapPoint> map_point_ptr);

private:

  mapManagement();
};