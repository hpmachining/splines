#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen\Dense>

//using namespace Eigen;

int main(void) {
  Eigen::Affine3d transform;
  transform.setIdentity();
  transform.rotate(Eigen::AngleAxisd(.25 * M_PI, Eigen::Vector3d::UnitZ()));
  
  Eigen::Affine3d offset(Eigen::Translation3d(Eigen::Vector3d(1.0, 0.5, -1.0)));
  transform = offset * transform;
  std::cout << "Transformation Matrix Rotated and Offset\n" << transform.matrix() << "\n\n";
  transform.rotate(Eigen::AngleAxisd(.25 * M_PI, Eigen::Vector3d::UnitZ()));
  std::cout << "Transformation Matrix Rotated Again\n" << transform.matrix() << "\n\n";

  Eigen::Vector4d p1(1, 0, 0, 0), p2;
  std::cout << "Point\n" << p1 << "\n\n";
  p2 = transform * p1;
  std::cout << "Point\n" << p2 << "\n\n";
  p1.w() = 1;
  p2 = transform * p1;
  std::cout << "Point\n" << p2 << "\n\n";

  return 0;
}