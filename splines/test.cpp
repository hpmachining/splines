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

template <typename RealScalar, size_t rows, size_t cols>
void WriteFile(const std::string& file, const Eigen::Matrix<RealScalar, rows, cols>& out) {
	std::ofstream out_file(file);
	file << out;
}

void TestDegreeElevation() {
	using Point = Eigen::Vector2d;
  size_t dimension = 2; // 2d or 3d points
	const size_t degree = 3;
	const size_t order = degree + 1;

	// Read control points from file
	Point point;
	std::vector<Point> points;
	while (std::cin >> point.x() >> point.y()) {// >> point.z()) {
		points.push_back(point);
	}
  std::cout << std::setprecision(3) << std::fixed;
  std::cout << "Original control points\n";
  for (auto i : points) {
    std::cout << i[0] << " " << i[1] << "\n";
  }
  size_t segment_count = points.size() / order;
  std::vector<double> all_coefficients;
  for (size_t i = 0; i < segment_count; ++i) {
    std::vector<double> coefficients = bezier::GetCoefficients<double>(points, i, degree, dimension);
    all_coefficients.insert(std::end(all_coefficients), std::begin(coefficients), std::end(coefficients));
  }
  std::cout << "\nCoefficients\n";
  for (auto i : all_coefficients) {
    std::cout << i << '\n';
  }
  std::cout << "\nElevated control points\n";
	std::vector<Point> elevated_points = bezier::ElevateDegree<double, bezier::k2d, degree>(points);
  size_t elevated_degree = degree + 1;
  size_t elevated_order = order + 1;
  for (auto i : elevated_points) {
		std::cout << i[0] << " " << i[1] << "\n";
	}

  all_coefficients.clear();
  
  for (size_t i = 0; i < segment_count; ++i) {
    std::vector<double> elevated_coefficients = bezier::GetCoefficients<double>(elevated_points, i, elevated_degree, dimension);
    all_coefficients.insert(std::end(all_coefficients), std::begin(elevated_coefficients), std::end(elevated_coefficients));
  }
  std::cout << "\nElevated Coefficients\n";
  for (auto i : all_coefficients) {
    std::cout << i << '\n';
  }

  std::cout << "\nControl points from coefficients\n";
  for (size_t i = 0; i < segment_count; ++i) {
    std::vector<Point> cp_coeff = bezier::GetControlPoints<double, Point>(all_coefficients, i, elevated_degree, dimension);
    for (auto i : cp_coeff) {
      std::cout << i[0] << " " << i[1] << "\n";
    }
  }


  //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> test_invert = bezier::GetPowerCoefficients<double>(6);
  //test_invert.colwise().reverseInPlace();
  //std::cout << "\nBinomial Coefficients degree 4\n" << test_invert.inverse();
	//for (auto i = 0.0; i <= 1.0; i += .25) {
	//	Point coord = bezier::GetPosition<double>(elevated_points, i, 0, 3, dimension);
	//	Point first = bezier::GetFirstDerivative<double>(elevated_points, i, 0, 3, dimension);
	//	Point second = bezier::GetSecondDerivative<double>(elevated_points, i, 0, 3, dimension);
	//	std::cout << "t = " << i << "\n";
	//	std::cout << "Coordinate = x: " << coord.x() << " y: " << coord.y() << "\n";
	//	std::cout << "first derivative = x: " << first.x() << " y: " << first.y() << "\n";
	//	std::cout << "second derivative = x: " << second.x() << " y: " << second.y() << "\n\n";
	//}
	//for (size_t i = 0; i < points.size() / order; ++i) {
	//	const size_t segment_index = i * order;

	//	// Create and fill Eigen::Matrix with control points for specified segment
	//	Eigen::Matrix<double, order, dimension> P = Eigen::Map<Eigen::Matrix<double, order, dimension, Eigen::RowMajor>>
	//		(points[segment_index].data(), order, dimension);

	//	// First point stays the same so add to elevated points
	//	std::vector<Point> elevated_points;
	//	elevated_points.push_back(P.row(0).transpose());

	//	// Calculate the new control points
	//	Eigen::Matrix<double, degree, 1> M1(Eigen::Matrix<double, degree, 1>::LinSpaced(degree, (1.0 / order),
	//		static_cast<double>(degree) / order));
	//	Eigen::Matrix<double, degree, dimension> Q =
	//		(M1.asDiagonal() * P.topRows<degree>()) + (M1.reverse().asDiagonal() * P.bottomRows<degree>());

	//	// Add new control points to elevated points
	//	for (size_t j = 0; j < degree; ++j) {
	//		elevated_points.push_back(Q.row(j).transpose());
	//	}

	//	// Last point stays the same so add to elevated points
	//	elevated_points.push_back(P.row(degree).transpose());

	//	// Output elevated points
	//	size_t k = 0;
	//	for (auto j : elevated_points) {
	//		std::cout << j.transpose() << '\n';
	//		++k;
	//		if (k % (order + 1) == 0) {
	//			std::cout << '\n';
	//		}
	//	}
	//	// Get coefficients
	//	std::vector<double> coeff = bezier::GetCoefficients<double>(elevated_points, 0, 3, 2);
	//	for (auto c : coeff) {
	//		std::cout << c << "\n";
	//	}
	//}
}

void TestMatrix() {
	using Point = Eigen::Vector3d;
	const size_t degree = 3;
	const size_t order = degree + 1;
	const size_t dimension = 3;
	double t = .75;
	Eigen::Matrix4d lh;
	Eigen::Matrix4d rh;

	Eigen::Matrix<double, order, order> test = bezier::GetPowerCoefficients<double>(degree);
	std::cout << test << '\n';

	Point point;
	std::vector<Point> points;
	while (std::cin >> point.x() >> point.y() >> point.z()) {
		points.push_back(point);
	}
	std::vector<Point, Eigen::aligned_allocator<Point>> all_segments;
	for (size_t i = 0; i < points.size() / order; ++i) {
		std::vector<Point> split_segments = bezier::SplitSegment<double>(points, .25, i + 1, degree, dimension);
		all_segments.insert(std::end(all_segments), std::begin(split_segments), std::end(split_segments));
	}

	size_t j = 0;
	for (auto i : all_segments) {
		std::cout << i.transpose() << '\n';
		++j;
		if (j % order == 0) {
			std::cout << '\n';
		}
	}
}

void TestTriangle() {
	using Point = Eigen::Vector3d;
	const size_t degree = 3;
	const size_t order = degree + 1;
	const size_t dimension = 3;
	double t = .75;
	Eigen::Matrix4d lh;
	Eigen::Matrix4d rh;

	// Read control points from file
	Point point;
	std::vector<Point, Eigen::aligned_allocator<Point>> points;
	while (std::cin >> point.x() >> point.y() >> point.z()) {
		points.push_back(point);
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	Eigen::Matrix<double, order, dimension> P;
	//for (auto i = 0; i < points.size() / order; ++i) {
	size_t i = 0;
	const size_t segment_index = (i)* order;
	for (size_t j = segment_index, p = 0; j < segment_index + order; ++j, ++p) {
		for (size_t k = 0; k < dimension; ++k) {
			P(p, k) = points[j][k];
		}
	}

	// Create and fill the power basis (t) matrix
	Eigen::Matrix<double, order, order> Z;
	Z.setZero();
	Eigen::Matrix<double, 1, order> power;
	for (size_t i = 0; i < order; ++i) {
		power(0, i) = std::pow(t, i);
		Z.diagonal()[i] = power(0, i);
	}

	Eigen::Matrix<double, order, order> M = bezier::GetPowerCoefficients<double>(degree);
	Eigen::Matrix<double, order, order> M1 = M.inverse();

	//std::cout <<
	//	"T\n" << power << "\n\n" <<
	//	"M1\n" << M1 << "\n\n" <<
	//	"Z\n" << Z << "\n\n" <<
	//	"M\n" << M << "\n\n" <<
	//	"P\n" << P << "\n\n";

	Eigen::Matrix<double, order, dimension> result;
	Eigen::Matrix<double, order, order> Q = M1 * Z * M;
	Eigen::Matrix<double, order, order> Q1;
	Q1.setZero();
	for (size_t i = 0; i < order; ++i) {
		Q1.block(order - (i + 1), order - (i + 1), 1, i + 1) = Q.block(i, 0, 1, i + 1);
	}
	result = Q * P;
	std::cout << result << "\n\n";
	result = Q1 * P;
	std::cout << result << "\n\n";
	//std::cout << "Q\n" << Q << "\n\n";
	//std::cout << "Q1\n" << Q1 << "\n\n";
	//std::cout << "Z\n" << Z << "\n\n";
	//std::cout << "Z1\n" << Z1 << "\n\n";
	//std::cout << "Q * P\n" << result << "\n\n";

}

void TestSecondDerivative() {
	const size_t degree = 5;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coefficients;
	coefficients.resize((degree - 1), (degree - 1));
	coefficients = bezier::GetSecondDerivativeCoefficients<double>(degree);
	std::cout << "2nd Derivative Coefficients\n" << coefficients << "\n\n";
}

void TestFirstDerivative() {
	using Eigen::Dynamic;
	const size_t degree = 3;
	double t = .5;

	// Create and fill the power basis matrix
	Eigen::Matrix<double, 1, Dynamic> parameter;
	parameter.resize(Eigen::NoChange, degree);
	for (size_t i = 0; i < degree; ++i) {
		parameter(0, i) = std::pow(t, i);
	}

	Eigen::Matrix<double, Dynamic, Dynamic> basis;
	basis.resize(degree, degree);
	basis << bezier::GetPowerCoefficients<double>(degree - 1);
	basis *= (degree * 2);

	std::cout << "parameter\n" << parameter << "\n\n";
	std::cout << "basis\n" << basis << "\n\n";

}

void GnuPlotPoints() {
  const double G = 120.0;
  const double ER = 180.0;
  double T;
  std::cout << std::fixed << std::setprecision(6);
  for (double x = 150.0; x <= 200.; x += 1.0) {
    for (double y = -70.0; y <= -50.0; y += .5) {
      T = (G / 2.)*(y + (G / 2.)) / (ER*ER - 2.*ER*x + G*y + (G*G / 4.) + x*x + y*y);
      std::cout << '\n' << T;
    }
    std::cout << '\n';
  }

}

void TestLength() {
  using Point = Eigen::Vector2d;
  size_t dimension = 2; // 2d or 3d points
  const size_t degree = 3;
  const size_t order = degree + 1;

  // Read control points from file
  Point point;
  std::vector<Point> points;
  while (std::cin >> point.x() >> point.y()) {// >> point.z()) {
    points.push_back(point);
  }
  //Eigen::Matrix<double, 24, 1> A;
  //A <<
  //  -0.0640568928626056,
  //  0.0640568928626056,
  //  -0.1911188674736163,
  //  0.1911188674736163,
  //  -0.3150426796961634,
  //  0.3150426796961634,
  //  -0.4337935076260451,
  //  0.4337935076260451,
  //  -0.5454214713888396,
  //  0.5454214713888396,
  //  -0.6480936519369755,
  //  0.6480936519369755,
  //  -0.7401241915785544,
  //  0.7401241915785544,
  //  -0.8200019859739029,
  //  0.8200019859739029,
  //  -0.8864155270044011,
  //  0.8864155270044011,
  //  -0.9382745520027328,
  //  0.9382745520027328,
  //  -0.9747285559713095,
  //  0.9747285559713095,
  //  -0.9951872199970213,
  //  0.9951872199970213;

  //Eigen::Matrix<double, 24, 1> W;
  //W <<
  //  0.1279381953467522,
  //  0.1279381953467522,
  //  0.1258374563468283,
  //  0.1258374563468283,
  //  0.1216704729278034,
  //  0.1216704729278034,
  //  0.1155056680537256,
  //  0.1155056680537256,
  //  0.1074442701159656,
  //  0.1074442701159656,
  //  0.0976186521041139,
  //  0.0976186521041139,
  //  0.0861901615319533,
  //  0.0861901615319533,
  //  0.0733464814110803,
  //  0.0733464814110803,
  //  0.0592985849154368,
  //  0.0592985849154368,
  //  0.0442774388174198,
  //  0.0442774388174198,
  //  0.0285313886289337,
  //  0.0285313886289337,
  //  0.0123412297999872,
  //  0.0123412297999872;

  Eigen::Matrix<double, 32, 1> A;
  A <<
    -0.0483076656877383,
    0.0483076656877383,
    -0.1444719615827965,
    0.1444719615827965,
    -0.2392873622521371,
    0.2392873622521371,
    -0.3318686022821277,
    0.3318686022821277,
    -0.4213512761306353,
    0.4213512761306353,
    -0.5068999089322294,
    0.5068999089322294,
    -0.5877157572407623,
    0.5877157572407623,
    -0.6630442669302152,
    0.6630442669302152,
    -0.7321821187402897,
    0.7321821187402897,
    -0.7944837959679424,
    0.7944837959679424,
    -0.8493676137325700,
    0.8493676137325700,
    -0.8963211557660521,
    0.8963211557660521,
    -0.9349060759377397,
    0.9349060759377397,
    -0.9647622555875064,
    0.9647622555875064,
    -0.9856115115452684,
    0.9856115115452684,
    -0.9972638618494816,
    0.9972638618494816;
    //std::cout << std::setprecision(16) << std::fixed << A << '\n';
  Eigen::Matrix<double, 32, 1> W;
  W <<
    0.0965400885147278,
    0.0965400885147278,
    0.0956387200792749,
    0.0956387200792749,
    0.0938443990808046,
    0.0938443990808046,
    0.0911738786957639,
    0.0911738786957639,
    0.0876520930044038,
    0.0876520930044038,
    0.0833119242269467,
    0.0833119242269467,
    0.0781938957870703,
    0.0781938957870703,
    0.0723457941088485,
    0.0723457941088485,
    0.0658222227763618,
    0.0658222227763618,
    0.0586840934785355,
    0.0586840934785355,
    0.0509980592623762,
    0.0509980592623762,
    0.0428358980222267,
    0.0428358980222267,
    0.0342738629130214,
    0.0342738629130214,
    0.0253920653092621,
    0.0253920653092621,
    0.0162743947309057,
    0.0162743947309057,
    0.0070186100094701,
    0.0070186100094701;

  //std::cout << W;

  auto rows = W.rows();
  double z = .5;
  double sum = 0.0;
  for (auto i = 0; i < rows; ++i) {
    double a = A(i, 0);
    double w = W(i, 0);
    double t = z * a + z;
    Point derived = bezier::GetFirstDerivative(points, t, 0, degree, dimension);
    //std::cout << "\nt = " << t << '\n';
    //std::cout << "derived = " << derived.transpose() << '\n';
    derived = derived.array().square();
    double derived_sum = derived.sum();
    derived_sum = std::sqrt(derived_sum);
    derived_sum *= w;
    sum += derived_sum;
  }
  double length = z * sum;
  std::cout << std::setprecision(14) << std::fixed << "\nlength = " << length << '\n';
  return;
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
	using bezier::k2d;	
	using bezier::k3d;
	//TestMatrix();
	//TestFirstDerivative();
	//TestDegreeElevation();
  //GnuPlotPoints();
  TestLength();

	//Eigen::Vector3d controlPoint;
	//std::vector<Eigen::Vector3d> controlPoints;
	//while (std::cin >> controlPoint.x() >> controlPoint.y() >> controlPoint.z()) {
	//	controlPoints.push_back(controlPoint);
	//}
	//std::vector<double> coefficients = bezier::GetCoefficients<double, k3d, kCubic>(controlPoints);
	//for (auto i : coefficients) {
	//	std::cout << i << '\n';
	//}

	//Eigen::Vector2d controlPoint;
	//std::vector<Eigen::Vector2d> controlPoints;
	//while (std::cin >> controlPoint.x() >> controlPoint.y()) {
	//	controlPoints.push_back(controlPoint);
	//}
	//std::vector<double> coefficients = bezier::GetCoefficients<double, k2d, kQuadratic>(controlPoints, 1);
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
	//for (auto i = 0; i < controlPoints.size() / 3; ++i) {
	//	coordinate = bezier::GetPosition<double, k2d, kQuadratic>(controlPoints, i, .25);
	//	std::cout.precision(5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::GetPosition<double, k2d, kQuadratic>(controlPoints, i, .5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::GetPosition<double, k2d, kQuadratic>(controlPoints, i, .75);
	//	std::cout << coordinate.transpose() << '\n';
	//}

	//Eigen::Vector3d controlPoint;
	//std::vector<Eigen::Vector3d> controlPoints;
	//while (std::cin >> controlPoint.x() >> controlPoint.y() >> controlPoint.z()) {
	//	controlPoints.push_back(controlPoint);
	//}
	//Eigen::Vector3d coordinate;
	//for (auto i = 0; i < controlPoints.size() / 4; ++i) {
	//	coordinate = bezier::GetPosition<double, k3d, kCubic>(controlPoints, i, .25);
	//	std::cout.precision(5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::GetPosition<double, k3d, kCubic>(controlPoints, i, .5);
	//	std::cout << coordinate.transpose() << '\n';
	//	coordinate = bezier::GetPosition<double, k3d, kCubic>(controlPoints, i, .75);
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
	//for (auto i = 0; i < control_points.size() / 4; ++i) {
	//	coordinate = bezier::GetPosition<double, bezier::k3d, 3>(control_points, i, .5);
	//	std::cout << coordinate.transpose() << '\n';
	//	tangent = bezier::GetFirstDerivative<double, bezier::k3d, 3>(control_points, i, .5);
	//	tangent.normalize();
	//	tangent += coordinate;
	//	std::cout << tangent.transpose() << '\n';
	//	normal = bezier::GetNormal<double, bezier::k3d, 3>(control_points, i, .5);
	//	normal.normalize();
	//	normal += coordinate;
	//	std::cout << normal.transpose() << '\n';
	//}

	//return 0;
}