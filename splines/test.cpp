#define EIGEN_MPL2_ONLY
#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <Eigen/Dense>

//using namespace Eigen;

int main(void) {
  Eigen::Affine3d transform;
  transform.setIdentity();
  std::cout << "Transformation Matrix Set to Identity\n" << transform.matrix() << "\n\n";
 
  Eigen::Vector3d p1(1, 0, 0);
  std::cout << "Point 1\n" << p1 << "\n\n";

  Eigen::AngleAxisd rotation(.25 * M_PI, Eigen::Vector3d::UnitZ());
  Eigen::Vector3d offset(0.75, 0.5, -1.0);

  transform.rotate(rotation).pretranslate(offset);
  std::cout << "Transformation Matrix: rotate & pretranslate\n" << transform.matrix() << "\n\n";
  Eigen::Vector3d p2 = transform * p1;
  std::cout << "Point 2: Point 1 transformed using transformation matrix\n" << p2 << "\n\n";

  p2 = transform.inverse() * p2;
  std::cout << "Point 2: Transformed using inverse transformation matrix\n" << p2 << "\n\n";

  transform.rotate(rotation.inverse()).pretranslate(-offset);
  std::cout << "Transformation Matrix: Inverse Rotate and Translate\n" << transform.matrix() << "\n\n";

  transform.translate(offset).rotate(rotation);
  std::cout << "Transformation Matrix: translate\n" << transform.matrix() << "\n\n";
  p2 = transform * p1;
  std::cout << "Point 2: Point 1 transformed using transformation matrix\n" << p2 << "\n\n";

  transform.prerotate(rotation);
  std::cout << "Transformation Matrix: translate\n" << transform.matrix() << "\n\n";
  p2 = transform * p1;
  std::cout << "Point 2: Point 1 transformed using transformation matrix\n" << p2 << "\n\n";

  return 0;
}