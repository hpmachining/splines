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

template <typename Scalar, size_t rows, size_t cols>
void WriteFile(const std::string& file, const Eigen::Matrix<Scalar, rows, cols>& out) {
	std::ofstream out_file(file);
	file << out;
}

void TestMatrix() {
	typedef double Scalar;
	const size_t degree = 3;
	const size_t order = degree + 1;
	const size_t dimension = 3;
	double t = .75;
	double t1 = 1.0 - t;
	Eigen::Matrix4d lh;
	Eigen::Matrix4d rh;

	Eigen::Matrix<Scalar, dimension, 1> point;
	std::vector<Eigen::Matrix<Scalar, dimension, 1>> points;
	while (std::cin >> point.x() >> point.y() >> point.z()) {
		points.push_back(point);
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	Eigen::Matrix<Scalar, order, dimension> P;
	//for (auto i = 0; i < points.size() / order; ++i) {
	size_t i = 0;
	const size_t segment_index = (i)* order;
	for (size_t j = segment_index, p = 0; j < segment_index + order; ++j, ++p) {
		for (size_t k = 0; k < dimension; ++k) {
			P(p, k) = points[j][k];
		}
	}

	// Create and fill the power basis (t) matrix
	Eigen::Matrix<Scalar, order, order> Z;
	Z.setZero();
	Eigen::Matrix<Scalar, order, order> Z1;
	Z1.setZero();
	Eigen::Matrix<Scalar, 1, order> power;
	Eigen::Matrix<Scalar, 1, order> power1;
	for (size_t i = 0; i < order; ++i) {
		power(0, i) = std::pow(t, i);
		power1(0, i) = std::pow(t1, i);
		Z.diagonal()[i] = power(0, i);
		Z1.diagonal()[i] = power1(0, i);
	}

	Eigen::Matrix<double, order, order> M = bezier::GetPowerCoefficients<double, degree>();
	Eigen::Matrix<double, order, order> M1 = M.inverse();

	//std::cout <<
	//	"T\n" << power << "\n\n" <<
	//	"M1\n" << M1 << "\n\n" <<
	//	"Z\n" << Z << "\n\n" <<
	//	"M\n" << M << "\n\n" <<
	//	"P\n" << P << "\n\n";

	Eigen::Matrix<Scalar, order, dimension> result;
	Eigen::Matrix<Scalar, order, order> Q = M1 * Z * M;
	Eigen::Matrix<Scalar, order, order> Q1 = M1 * Z1 * M;
	Q1.colwise().reverseInPlace();
	Q1.rowwise().reverseInPlace();
	result = Q * P;
	std::cout << result << "\n\n";
	result = Q1 * P;
	std::cout << result << "\n\n";
	//std::cout << "Q\n" << Q << "\n\n";
	//std::cout << "Q1\n" << Q1 << "\n\n";
	//std::cout << "Z\n" << Z << "\n\n";
	//std::cout << "Z1\n" << Z1 << "\n\n";
	//std::cout << "Q * P\n" << result << "\n\n";


	//}


}

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
	TestMatrix();

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

	//using Point = Eigen::Vector3d;
	//Point control_point;
	//std::vector<Point> control_points;
	//Eigen::Matrix3Xd C;
	//size_t col = 0;
	//while (std::cin >> control_point.x() >> control_point.y() >> control_point.z()) {
	//	C.conservativeResize(Eigen::NoChange, col + 1);
	//	control_points.push_back(control_point);
	//	C.block(0, col, 3, 1) = control_point;
	//	++col;
	//}
	//size_t columns = C.cols();
	//std::cout << std::fixed << std::setprecision(6);
	//std::cout << C.transpose() << '\n';
	////for (auto i : control_points) {
	////	std::cout << i.transpose() << '\n';
	////}
	//Eigen::Vector3d coordinate;
	//Eigen::Vector3d tangent;
	//Eigen::Vector3d normal;
	//std::cout << std::fixed << std::setprecision(14);
	//for (auto i = 1; i <= control_points.size() / 4; ++i) {
	//	coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
	//	std::cout << coordinate.transpose() << '\n';
	//	tangent = bezier::CalculateTangent<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
	//	tangent.normalize();
	//	tangent += coordinate;
	//	std::cout << tangent.transpose() << '\n';
	//	normal = bezier::CalculateNormal<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
	//	normal.normalize();
	//	normal += coordinate;
	//	std::cout << normal.transpose() << '\n';
	//}

	//return 0;
}