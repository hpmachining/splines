#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen\Dense>


Eigen::Affine3d create_rotation_matrix(double ax, double ay, double az) {
  Eigen::Affine3d rx =
    Eigen::Affine3d(Eigen::AngleAxisd(ax, Eigen::Vector3d(1, 0, 0)));
  Eigen::Affine3d ry =
    Eigen::Affine3d(Eigen::AngleAxisd(ay, Eigen::Vector3d(0, 1, 0)));
  Eigen::Affine3d rz =
    Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
  return rz * ry * rx;
}

using namespace Eigen;

int main(void) {
  Vector4d xAxis(1.0, 0.0, 0.0, 0.0);
  Vector4d yAxis(0.0, 1.0, 0.0, 0.0);
  Vector4d zAxis(0.0, 0.0, 1.0, 0.0);
  Vector4d origin(0.0, 0.0, 0.0, 1.0);
  Matrix4d identity;
  identity.setIdentity();

  Affine3d r = create_rotation_matrix(0.0, 0.0, .25 * M_PI);
  Affine3d t(Translation3d(Vector3d(0.0, 0.0, 0.0)));
  Matrix4d transform = (t * r).matrix();

  std::cout << transform << '\n';

  return 0;
}