#define EIGEN_MPL2_ONLY
#define _USE_MATH_DEFINES

#include "splines.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/StdVector>

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
	using bezier::kCubic;
	using bezier::kQuadratic;
	using bezier::k2d;	
	using bezier::k3d;

	//Eigen::Vector3d controlPoint;
	//std::vector<Eigen::Vector3d> controlPoints;
	//while (std::cin >> controlPoint.x() >> controlPoint.y() >> controlPoint.z()) {
	//	controlPoints.push_back(controlPoint);
	//}
	//std::vector<double> coefficients = bezier::CalculateCoefficients<double, k3d, kCubic>(controlPoints);
	//for (auto i : coefficients) {
	//	std::cout << i << '\n';
	//}

	//Eigen::Vector2d controlPoint;
	//std::vector<Eigen::Vector2d> controlPoints;
	//while (std::cin >> controlPoint.x() >> controlPoint.y()) {
	//	controlPoints.push_back(controlPoint);
	//}
	//std::vector<double> coefficients = bezier::CalculateCoefficients<double, k2d, kQuadratic>(controlPoints, 1);
	//for (auto i : coefficients) {
	//	std::cout << i << '\n';
	//}

	//// Routine for eliminating duplicate points to make b-spline control points
	//auto uniqueEnd = std::unique(controlPoints.begin(), controlPoints.end(),
	//	[](const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs) { return lhs.isApprox(rhs); });
	//controlPoints.erase(uniqueEnd, controlPoints.end());
	
	//Eigen::Vector2d controlPoint;
	//std::vector<Eigen::Vector2d> controlPoints;
	//while (std::cin >> controlPoint.x() >> controlPoint.y()) {
	//	controlPoints.push_back(controlPoint);
	//}
	//Eigen::Vector2d coordinate;
	//for (auto i = 1; i <= controlPoints.size() / 3; ++i) {
	//	coordinate = bezier::CalculateCoordinate<double, k2d, kQuadratic>(controlPoints, i, .25);
	//	std::cout.precision(5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::CalculateCoordinate<double, k2d, kQuadratic>(controlPoints, i, .5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::CalculateCoordinate<double, k2d, kQuadratic>(controlPoints, i, .75);
	//	std::cout << coordinate.transpose() << '\n';
	//}

	//Eigen::Vector3d controlPoint;
	//std::vector<Eigen::Vector3d> controlPoints;
	//while (std::cin >> controlPoint.x() >> controlPoint.y() >> controlPoint.z()) {
	//	controlPoints.push_back(controlPoint);
	//}
	//Eigen::Vector3d coordinate;
	//for (auto i = 1; i <= controlPoints.size() / 4; ++i) {
	//	coordinate = bezier::CalculateCoordinate<double, k3d, kCubic>(controlPoints, i, .25);
	//	std::cout.precision(5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::CalculateCoordinate<double, k3d, kCubic>(controlPoints, i, .5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::CalculateCoordinate<double, k3d, kCubic>(controlPoints, i, .75);
	//	std::cout << coordinate.transpose() << '\n';
	//}

	using Point = Eigen::Vector3d;
	Point control_point;
	std::vector<Point> control_points;
	Eigen::Matrix3Xd C;
	size_t col = 0;
	while (std::cin >> control_point.x() >> control_point.y() >> control_point.z()) {
		C.conservativeResize(Eigen::NoChange, col + 1);
		control_points.push_back(control_point);
		C.block(0, col, 3, 1) = control_point;
		++col;
	}
	size_t columns = C.cols();
	std::cout << std::fixed << std::setprecision(6);
	std::cout << C.transpose() << '\n';
	//for (auto i : control_points) {
	//	std::cout << i.transpose() << '\n';
	//}
	Eigen::Vector3d coordinate;
	Eigen::Vector3d tangent;
	std::cout << std::fixed << std::setprecision(14);
	for (auto i = 1; i <= control_points.size() / 4; ++i) {
		coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
		std::cout << coordinate.transpose() << '\n';
		//tangent = bezier::CalculateTangent<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
		//tangent.normalize();
		//tangent += coordinate;
		//std::cout << tangent.transpose() << '\n';
	}

	return 0;
}