
#include "stdafx.h"
#include "CppUnitTest.h"
#include "splines.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <fstream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/StdVector>


using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{		
	TEST_CLASS(UnitTest1)

	{
	public:
		double max_double = std::numeric_limits<double>::max();

		TEST_METHOD(Cubic2dCoefficient)
		{
			std::ifstream in_file("cubic2d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("cubic2dcoeff.points");
			Eigen::Vector2d control_point;
			std::vector<Eigen::Vector2d> control_points;
			while (in_file >> control_point.x() >> control_point.y()) {
				control_points.push_back(control_point);
			}
			std::vector<double> coefficients;
			for (size_t i = 1; i <= control_points.size() / 4; ++i) {
				coefficients = bezier::CalculateCoefficients<double, bezier::k2d, bezier::kCubic>(control_points, i);
				Assert::IsTrue(coefficients.size() == 8);
				out_file << std::fixed << std::setprecision(14);
				for (auto j : coefficients) {
					out_file << j << '\n';
				}
			}
		}

		TEST_METHOD(Cubic3dCoefficient)
		{
			std::ifstream in_file("cubic3d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("cubic3dcoeff.points");
			Eigen::Vector3d control_point;
			std::vector<Eigen::Vector3d> control_points;
			while (in_file >> control_point.x() >> control_point.y() >> control_point.z()) {
				control_points.push_back(control_point);
			}
			std::vector<double> coefficients;
			for (size_t i = 1; i <= control_points.size() / 4; ++i) {
				coefficients = bezier::CalculateCoefficients<double, bezier::k3d, bezier::kCubic>(control_points, i);
				Assert::IsTrue(coefficients.size() == 12);
				out_file << std::fixed << std::setprecision(14);
				for (auto j : coefficients) {
					out_file << j << '\n';
				}
			}
		}

		TEST_METHOD(Quad2dCoefficient)
		{
			std::ifstream in_file("quad2d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quad2dcoeff.points");
			Eigen::Vector2d control_point;
			std::vector<Eigen::Vector2d> control_points;
			while (in_file >> control_point.x() >> control_point.y()) {
				control_points.push_back(control_point);
			}
			std::vector<double> coefficients;
			for (size_t i = 1; i <= control_points.size() / 3; ++i) {
				coefficients = bezier::CalculateCoefficients<double, bezier::k2d, bezier::kQuadratic>(control_points, i);
				Assert::IsTrue(coefficients.size() == 6);
				out_file << std::fixed << std::setprecision(14);
				for (auto j : coefficients) {
					out_file << j << '\n';
				}
			}
		}

		TEST_METHOD(Quad3dCoefficient)
		{
			std::ifstream in_file("quad3d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quad3dcoeff.points");
			Eigen::Vector3d control_point;
			std::vector<Eigen::Vector3d> control_points;
			while (in_file >> control_point.x() >> control_point.y() >> control_point.z()) {
				control_points.push_back(control_point);
			}
			std::vector<double> coefficients;
			for (size_t i = 1; i <= control_points.size() / 3; ++i) {
				coefficients = bezier::CalculateCoefficients<double, bezier::k3d, bezier::kQuadratic>(control_points, i);
				Assert::IsTrue(coefficients.size() == 9);
				out_file << std::fixed << std::setprecision(14);
				for (auto j : coefficients) {
					out_file << j << '\n';
				}
			}
		}

		TEST_METHOD(Cubic2dCoordinate)
		{
			std::ifstream in_file("cubic2d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("cubic2dcoord.points");
			Eigen::Vector2d control_point;
			std::vector<Eigen::Vector2d> control_points;
			while (in_file >> control_point.x() >> control_point.y()) {
				control_points.push_back(control_point);
			}

			Eigen::Vector2d coordinate;
			Eigen::Vector2d empty;
			empty << max_double, max_double;
			coordinate = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 4; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kCubic>(control_points, i, .25);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kCubic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kCubic>(control_points, i, .75);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
			}

		}

		TEST_METHOD(Cubic3dCoordinate)
		{
			std::ifstream in_file("cubic3d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("cubic3dcoord.points");
			Eigen::Vector3d control_point;
			std::vector<Eigen::Vector3d> control_points;
			while (in_file >> control_point.x() >> control_point.y() >> control_point.z()) {
				control_points.push_back(control_point);
			}

			Eigen::Vector3d coordinate;
			Eigen::Vector3d empty;
			empty << max_double, max_double, max_double;
			coordinate = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 4; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(control_points, i, .25);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(control_points, i, .75);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
			}

		}

		TEST_METHOD(Quad2dCoordinate)
		{
			std::ifstream in_file("quad2d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quad2dcoord.points");
			Eigen::Vector2d control_point;
			std::vector<Eigen::Vector2d> control_points;
			while (in_file >> control_point.x() >> control_point.y()) {
				control_points.push_back(control_point);
			}
			Eigen::Vector2d coordinate;
			Eigen::Vector2d empty;
			empty << max_double, max_double;
			coordinate = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 3; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kQuadratic>(control_points, i, .25);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kQuadratic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kQuadratic>(control_points, i, .75);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
			}

		}

		TEST_METHOD(Quad3dCoordinate)
		{
			std::ifstream in_file("quad3d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quad3dcoord.points");
			Eigen::Vector3d control_point;
			std::vector<Eigen::Vector3d> control_points;
			while (in_file >> control_point.x() >> control_point.y() >> control_point.z()) {
				control_points.push_back(control_point);
			}
			Eigen::Vector3d coordinate;
			Eigen::Vector3d empty;
			empty << max_double, max_double, max_double;
			coordinate = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 3; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kQuadratic>(control_points, i, .25);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kQuadratic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kQuadratic>(control_points, i, .75);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				coordinate = empty;
			}

		}

		TEST_METHOD(Cubic2dTangent)
		{
			std::ifstream in_file("cubic2d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("cubic2dtan.points");
			Eigen::Vector2d control_point;
			std::vector<Eigen::Vector2d> control_points;
			while (in_file >> control_point.x() >> control_point.y()) {
				control_points.push_back(control_point);
			}

			Eigen::Vector2d coordinate;
			Eigen::Vector2d tangent;
			Eigen::Vector2d empty;
			empty << max_double, max_double;
			coordinate = empty;
			tangent = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 4; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kCubic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				tangent = bezier::CalculateTangent<double, bezier::k2d, bezier::kCubic>(control_points, i, .5);
				Assert::IsFalse(tangent == empty);
				tangent.normalize();
				tangent += coordinate;
				out_file << tangent.transpose() << '\n';
				coordinate = empty;
				tangent = empty;
			}
		}

		TEST_METHOD(Cubic3dTangent)
		{
			using Point = Eigen::Vector3d;
			//using Point = std::array<double, 3>;
			std::ifstream in_file("cubic.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("cubic3dtan.points");
			Point control_point;
			std::vector<Point> control_points;
			//while (in_file >> control_point.x() >> control_point.y() >> control_point.z()) {
			//	control_points.push_back(control_point);
			//}
			while (in_file >> control_point[0] >> control_point[1] >> control_point[2]) {
				control_points.push_back(control_point);
			}

			Point coordinate;
			Point tangent;
			Point normal;
			Point empty;
			//empty << max_double, max_double, max_double;
			for (size_t i = 0; i < 3; ++i) {
				empty[i] = max_double;
			}
			coordinate = empty;
			normal = empty;
			tangent = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 4; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				Eigen::Vector3d eigen_coordinate;
				eigen_coordinate << coordinate[0], coordinate[1], coordinate[2];
				out_file << eigen_coordinate.transpose() << '\n';
				//out_file << coordinate[0] << ' ' << coordinate[1] << ' ' << coordinate[2] << '\n';
				tangent = bezier::CalculateTangent<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
				Assert::IsFalse(tangent == empty);
				Eigen::Vector3d normal_tangent;
				normal_tangent << tangent[0], tangent[1], tangent[2];
				normal_tangent.normalize();
				normal_tangent += eigen_coordinate;
				out_file << normal_tangent.transpose() << '\n';
				//out_file << tangent.transpose() << '\n';
				normal = bezier::CalculateNormal<double, bezier::k3d, bezier::kCubic>(control_points, i, .5);
				Assert::IsFalse(normal == empty);
				Eigen::Vector3d eigen_normal;
				eigen_normal << normal[0], normal[1], normal[2];
				eigen_normal.normalize();
				eigen_normal += eigen_coordinate;
				out_file << eigen_normal.transpose() << '\n';
				coordinate = empty;
				tangent = empty;
			}
		}

		TEST_METHOD(Quad2dTangent)
		{
			std::ifstream in_file("quad2d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quad2dtan.points");
			Eigen::Vector2d control_point;
			std::vector<Eigen::Vector2d> control_points;
			while (in_file >> control_point.x() >> control_point.y()) {
				control_points.push_back(control_point);
			}

			Eigen::Vector2d coordinate;
			Eigen::Vector2d tangent;
			Eigen::Vector2d normal;
			Eigen::Vector2d empty;
			empty << max_double, max_double;
			coordinate = empty;
			tangent = empty;
			normal = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 3; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kQuadratic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				tangent = bezier::CalculateTangent<double, bezier::k2d, bezier::kQuadratic>(control_points, i, .5);
				Assert::IsFalse(tangent == empty);
				tangent.normalize();
				tangent += coordinate;
				out_file << tangent.transpose() << '\n';
				//normal = bezier::CalculateNormal<double, bezier::k2d, bezier::kQuadratic>(control_points, i, .5);
				//Assert::IsFalse(normal == empty);
				//normal.normalize();
				//normal += coordinate;
				//out_file << normal.transpose() << '\n';
				coordinate = empty;
				tangent = empty;
			}
		}

		TEST_METHOD(Quad3dTangent)
		{
			using Point = Eigen::Vector3d;
			std::ifstream in_file("quad3d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quad3dtan.points");
			Point control_point;
			std::vector<Point> control_points;
			while (in_file >> control_point.x() >> control_point.y() >> control_point.z()) {
				control_points.push_back(control_point);
			}

			Point coordinate;
			Point tangent;
			Point normal;
			Point empty;
			empty << max_double, max_double, max_double;
			coordinate = empty;
			tangent = empty;
			normal = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 3; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kQuadratic>(control_points, i, .75);
				Assert::IsFalse(coordinate == empty);
				out_file << coordinate.transpose() << '\n';
				tangent = bezier::CalculateTangent<double, bezier::k3d, bezier::kQuadratic>(control_points, i, .75);
				Assert::IsFalse(tangent == empty);
				tangent.normalize();
				tangent += coordinate;
				out_file << tangent.transpose() << '\n';
				normal = bezier::CalculateNormal<double, bezier::k3d, bezier::kQuadratic>(control_points, i, .75);
				Assert::IsFalse(normal == empty);
				normal.normalize();
				normal += coordinate;
				out_file << normal.transpose() << '\n';
				coordinate = empty;
				tangent = empty;
			}
		}

		TEST_METHOD(Quartic3dTangent)
		{
			using Point = Eigen::Vector3d;
			//using Point = std::array<double, 3>;
			std::ifstream in_file("quartic3d.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quartic3dtan.points");
			Point control_point;
			std::vector<Point> control_points;
			//while (in_file >> control_point.x() >> control_point.y() >> control_point.z()) {
			//	control_points.push_back(control_point);
			//}
			while (in_file >> control_point[0] >> control_point[1] >> control_point[2]) {
				control_points.push_back(control_point);
			}

			Point coordinate;
			Point tangent;
			Point normal;
			Point empty;
			//empty << max_double, max_double, max_double;
			for (size_t i = 0; i < 3; ++i) {
				empty[i] = max_double;
			}
			coordinate = empty;
			normal = empty;
			tangent = empty;
			out_file << std::fixed << std::setprecision(14);
			for (size_t i = 1; i <= control_points.size() / 4; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kQuartic>(control_points, i, .5);
				Assert::IsFalse(coordinate == empty);
				Eigen::Vector3d eigen_coordinate;
				eigen_coordinate << coordinate[0], coordinate[1], coordinate[2];
				out_file << eigen_coordinate.transpose() << '\n';
				//out_file << coordinate[0] << ' ' << coordinate[1] << ' ' << coordinate[2] << '\n';
				tangent = bezier::CalculateTangent<double, bezier::k3d, bezier::kQuartic>(control_points, i, .5);
				Assert::IsFalse(tangent == empty);
				Eigen::Vector3d normal_tangent;
				normal_tangent << tangent[0], tangent[1], tangent[2];
				normal_tangent.normalize();
				normal_tangent += eigen_coordinate;
				out_file << normal_tangent.transpose() << '\n';
				normal = bezier::CalculateNormal<double, bezier::k3d, bezier::kQuartic>(control_points, i, .5);
				Assert::IsFalse(normal == empty);
				Eigen::Vector3d eigen_normal;
				eigen_normal << normal[0], normal[1], normal[2];
				eigen_normal.normalize();
				eigen_normal += eigen_coordinate;
				out_file << eigen_normal.transpose() << '\n';
				coordinate = empty;
				tangent = empty;
			}
		}

		TEST_METHOD(PascalRow) {
			std::ofstream out_file("pascalRow.dat");
			out_file << std::fixed << std::setprecision(0);
			for (int row = 0; row < 100; ++row) {
				Eigen::Matrix<double, Eigen::Dynamic, 1> pascal_row = bezier::GetPascalRow(row);
				out_file << pascal_row.transpose() << '\n';
			}
		}

	};
}