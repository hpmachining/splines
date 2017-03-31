#include "stdafx.h"
#include "CppUnitTest.h"
#include "splines.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/StdVector>


using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		double max_double = std::numeric_limits<double>::max();

		TEST_METHOD(Cubic2dCoordinate)
		{
			// TODO: Your test code here
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
			for (auto i = 1; i <= control_points.size() / 4; ++i) {
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
			// TODO: Your test code here
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
			for (auto i = 1; i <= control_points.size() / 4; ++i) {
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
			// TODO: Your test code here
			std::ifstream in_file("parabola.dat");
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
			for (auto i = 1; i <= control_points.size() / 3; ++i) {
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
			std::vector<double> coefficients = bezier::CalculateCoefficients<double, bezier::k3d, bezier::kCubic>(control_points);
			Assert::IsTrue(coefficients.size() == (control_points.size() / 4) * 12);
			out_file << std::fixed << std::setprecision(14);
			for (auto i : coefficients) {
				out_file << i << '\n';
			}
		}

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
			std::vector<double> coefficients = bezier::CalculateCoefficients<double, bezier::k2d, bezier::kCubic>(control_points);
			Assert::IsTrue(coefficients.size() == (control_points.size() / 4) * 8);
			out_file << std::fixed << std::setprecision(14);
			for (auto i : coefficients) {
				out_file << i << '\n';
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
			for (auto i = 1; i <= control_points.size() / 4; ++i) {
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

	};
}