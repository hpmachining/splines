#include "stdafx.h"
#include "CppUnitTest.h"
#include "splines.h"
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
		
		TEST_METHOD(Cubic3dCoordinate)
		{
			// TODO: Your test code here
			std::ifstream in_file("ctrl.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("cubic3dcoord.points");
			Eigen::Vector3d controlPoint;
			std::vector<Eigen::Vector3d> controlPoints;
			while (in_file >> controlPoint.x() >> controlPoint.y() >> controlPoint.z()) {
				controlPoints.push_back(controlPoint);
			}
			Eigen::Vector3d coordinate;
			for (auto i = 1; i <= controlPoints.size() / 4; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(controlPoints, i, .25);
				out_file << coordinate.transpose() << '\n';
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(controlPoints, i, .5);
				out_file << coordinate.transpose() << '\n';
				coordinate = bezier::CalculateCoordinate<double, bezier::k3d, bezier::kCubic>(controlPoints, i, .75);
				out_file << coordinate.transpose() << '\n';
			}

		}

		TEST_METHOD(Quad2dCoordinate)
		{
			// TODO: Your test code here
			std::ifstream in_file("parabola.dat");
			Assert::IsTrue(in_file.is_open());
			std::ofstream out_file("quad2dcoord.points");
			Eigen::Vector2d controlPoint;
			std::vector<Eigen::Vector2d> controlPoints;
			while (in_file >> controlPoint.x() >> controlPoint.y()) {
				controlPoints.push_back(controlPoint);
			}
			Eigen::Vector2d coordinate;
			for (auto i = 1; i <= controlPoints.size() / 3; ++i) {
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kQuadratic>(controlPoints, i, .25);
				out_file << coordinate.transpose() << '\n';
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kQuadratic>(controlPoints, i, .5);
				out_file << coordinate.transpose() << '\n';
				coordinate = bezier::CalculateCoordinate<double, bezier::k2d, bezier::kQuadratic>(controlPoints, i, .75);
				out_file << coordinate.transpose() << '\n';
			}

		}

	};
}