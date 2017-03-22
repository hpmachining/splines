#define _USE_MATH_DEFINES
#define EIGEN_MPL2_ONLY

#include <cmath>
#include <iostream>
#include <Eigen\Dense>

//using namespace Eigen;

int main(void) {
  Eigen::Affine3d transform;
  transform.setIdentity();
  std::cout << "Transformation Matrix Set to Identity\n" << transform.matrix() << "\n\n";
 
  Eigen::AngleAxisd rotation(.25 * M_PI, Eigen::Vector3d::UnitZ());
  Eigen::Vector3d offset(1.0, 0.5, -1.0);

  transform.rotate(rotation);
  transform.pretranslate(offset);
  std::cout << "Transformation Matrix Rotated and Offset\n" << transform.matrix() << "\n\n";

  Eigen::Vector4d p1(1, 0, 0, 1), p2;
  std::cout << "Point\n" << p1 << "\n\n";
  p2 = transform * p1;
  std::cout << "Rotated and Translated Point\n" << p2 << "\n\n";
  //Eigen::Translation3d test(transform.translation());
  //p2 = test.inverse() * p2;
  //std::cout << "Test\n" << p2 << "\n\n";
  p2 = transform.inverse() * p2;
  //transform.pretranslate(-offset);
  //std::cout << "Transformation Matrix Translated back\n" << transform.matrix() << "\n\n";
  //p2 -= (Eigen::Vector4d(offset, 0));
  std::cout << "Point Rotated and Translated back\n" << p2 << "\n\n";
  return 0;
}