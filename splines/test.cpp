#define EIGEN_MPL2_ONLY
#define _USE_MATH_DEFINES

#include "splines.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/StdVector>

//using namespace Eigen;

//int main(void) {
//  Eigen::Affine3d transform;
//  transform.setIdentity();
//  std::cout << "Transformation Matrix Set to Identity\n" << transform.matrix() << "\n\n";
// 
//  Eigen::Eigen::Vector3d p1(1, 0, 0);
//  std::cout << "Point 1\n" << p1 << "\n\n";
//
//  Eigen::AngleAxisd rotation(.25 * M_PI, Eigen::Eigen::Vector3d::UnitZ());
//  Eigen::Eigen::Vector3d offset(0.75, 0.5, -1.0);
//
//  transform.rotate(rotation).pretranslate(offset);
//  std::cout << "Transformation Matrix: rotate & pretranslate\n" << transform.matrix() << "\n\n";
//  Eigen::Eigen::Vector3d p2 = transform * p1;
//  std::cout << "Point 2: Point 1 transformed using transformation matrix\n" << p2 << "\n\n";
//
//  p2 = transform.inverse() * p2;
//  std::cout << "Point 2: Transformed using inverse transformation matrix\n" << p2 << "\n\n";
//
//  transform.rotate(rotation.inverse()).pretranslate(-offset);
//  std::cout << "Transformation Matrix: Inverse Rotate and Translate\n" << transform.matrix() << "\n\n";
//
//  transform.translate(offset).rotate(rotation);
//  std::cout << "Transformation Matrix: translate\n" << transform.matrix() << "\n\n";
//  p2 = transform * p1;
//  std::cout << "Point 2: Point 1 transformed using transformation matrix\n" << p2 << "\n\n";
//
//  transform.prerotate(rotation);
//  std::cout << "Transformation Matrix: translate\n" << transform.matrix() << "\n\n";
//  p2 = transform * p1;
//  std::cout << "Point 2: Point 1 transformed using transformation matrix\n" << p2 << "\n\n";
//  return 0;
//}

int main(void) {
	Eigen::Vector3d controlPoint;
	std::vector<Eigen::Vector3d> controlPoints;
	std::ifstream dataFile("ctrl.dat");
	
	while (std::cin >> controlPoint.x() >> controlPoint.y() >> controlPoint.z()) {
	//while (dataFile >> controlPoint.x() >> controlPoint.y() >> controlPoint.z()) {
		controlPoints.push_back(controlPoint);
	}
	std::cout.precision(15);
	//controlPoints.pop_back();
	std::vector<double> coefficients =
		hpm::CalculateCoefficientsFromPoints<double, Eigen::Vector3d>(hpm::kCubic, controlPoints, false);
	for (auto i : coefficients) {
		std::cout << i << '\n';
	}

	Eigen::Vector3d coordinate = hpm::CalculateCoordinate(hpm::kCubic, 1, .5, controlPoints);
	
	auto uniqueEnd = std::unique(controlPoints.begin(), controlPoints.end(),
		[](const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs) { return lhs.isApprox(rhs); });
	controlPoints.erase(uniqueEnd, controlPoints.end());
	
	//std::cout << '\n';
	//for (auto i : controlPoints) {
	//	std::cout << i.x() << '\t' << i.y() << '\t' << i.z() << '\n';
	//}
	return 0;
}