#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>


std::vector<double> GetCoeffFromPoints(std::vector<Eigen::Vector3d> &points, bool is3D = true);

