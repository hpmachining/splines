/** @mainpage
* @author Paul J. Hentschel
*
* A collection of templated functions for working with Cubic Bézier curves.
*/

#pragma once

#include <vector>
#include <cmath>
#include <exception>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace hpm {

enum BezierDegree {
	kQuad = 2,
	kCubic,
};

/**	@brief	Calculate a vector of coefficients derived from a vector of control points
*			of a cubic Bézier spline.
*
*			Size of vector of control points must be in multiples of 4.
*
*	@param	points Vector of control points. This function will process the vector elements
*			in groups of 4.
*	@param	is3D Boolean value to determine if 2D or 3D spline coefficients are returned.
*	@return	A vector of coefficients calculated from the control points.
*			Coefficients are stored in the order
*			{\f$A_x, B_x, C_x, D_x, A_y, B_y, C_y, D_y\f$} and if the points are 3d
*			{\f$A_z, B_z, C_z, D_z\f$}
*/

template <typename Scalar, typename Point>
std::vector<Scalar> 
CalculateCoefficientsFromPoints(const BezierDegree degree, const std::vector<Point>& points, bool is3d) {
	std::vector<Scalar> coefficients;
	//if (degree != kQuadratic || degree != kCubic || ) {

	//}
	size_t increment = degree + 1;
	size_t num_control_points = points.size() / increment;

	if (num_control_points) {

		Scalar Ax = 0.0, Bx = 0.0, Cx = 0.0, Dx = 0.0;
		Scalar Ay = 0.0, By = 0.0, Cy = 0.0, Dy = 0.0;
		Scalar Az = 0.0, Bz = 0.0, Cz = 0.0, Dz = 0.0;

		Point cp[4];
		for (size_t i = 0; i < num_control_points; i += 4) {
			for (size_t j = 0; j < 4; ++j) {
				cp[j] = points[i + j];
			}

			Ax = cp[3].x() - 3.0 * cp[2].x() + 3.0 * cp[1].x() - cp[0].x();
			Ay = cp[3].y() - 3.0 * cp[2].y() + 3.0 * cp[1].y() - cp[0].y();
			Bx = 3.0 * cp[2].x() - 6.0 * cp[1].x() + 3.0 * cp[0].x();
			By = 3.0 * cp[2].y() - 6.0 * cp[1].y() + 3.0 * cp[0].y();
			Cx = 3.0 * (cp[1].x() - cp[0].x());
			Cy = 3.0 * (cp[1].y() - cp[0].y());
			Dx = cp[0].x();
			Dy = cp[0].y();

			coefficients.push_back(Ax);
			coefficients.push_back(Bx);
			coefficients.push_back(Cx);
			coefficients.push_back(Dx);

			coefficients.push_back(Ay);
			coefficients.push_back(By);
			coefficients.push_back(Cy);
			coefficients.push_back(Dy);

			if (is3d) {
				Az = cp[3].z() - 3.0 * cp[2].z() + 3.0 * cp[1].z() - cp[0].z();
				Bz = 3.0 * cp[2].z() - 6.0 * cp[1].z() + 3.0 * cp[0].z();
				Cz = 3.0 * (cp[1].z() - cp[0].z());
				Dz = cp[0].z();

				coefficients.push_back(Az);
				coefficients.push_back(Bz);
				coefficients.push_back(Cz);
				coefficients.push_back(Dz);
			}
		}
	}
	return coefficients;
}


template <typename Scalar, typename Point>
Point
CalculateCoordinate(const BezierDegree degree,
	const long segment_id,
	const Scalar t,
	const std::vector<Point>& points,
	const bool is_3d) {
	
	// Determine the size of a control point set for one Bézier curve segment,
	// then check if the input data is valid.
	const size_t control_size = degree + 1;
	const size_t num_segments = points.size() / control_size;
	Point coordinate;
	if (num_segments == 0 || segment_id <= 0 || segment_id > num_segments) {
		return coordinate;
	}

	size_t row_size = points.size();
	size_t column_size = 2;
	if (is_3d) {
		column_size = 3;
	}

	// Copy std::vector of control points to an Eigen::Matrix
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> control_points;
	control_points.resize(row_size, column_size);
	if (is_3d) {
		for (auto i = 0; i < row_size; ++i) {
			control_points(i, 0) = points[i].x();
			control_points(i, 1) = points[i].y();
			control_points(i, 2) = points[i].z();
		}
	}
	else {
		for (auto i = 0; i < row_size; ++i) {
			control_points(i, 0) = points[i].x();
			control_points(i, 1) = points[i].y();
		}
	}

	// Create and fill matrix with control points for specified segment
	size_t segment_index = (segment_id - 1) * control_size;
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> P;
	P.resize(control_size, column_size);
	P.topLeftCorner(control_size, column_size).setZero();
	P =	control_points.block(segment_index, 0, control_size, column_size);
	
	// Create and fill the time matrix
	Eigen::Matrix<Scalar, 1, Eigen::Dynamic> time;
	time.resize(1, control_size);
	for (auto i = 0; i < control_size; ++i) {
		time(0, i) = std::pow(t, i);
	}
	
	// Create and fill the basis matrix
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> basis;
	basis.resize(control_size, control_size);
	switch (degree) {
	case kCubic: 	
		basis <<	 
			 1,  0,  0, 0,
			-3,  3,  0, 0,
			 3, -6,  3, 0,
			-1,  3, -3, 1;

		break;
	case kQuad: 
		basis <<
			 1,  0, 0,
			-2,  2, 0,
			 1, -2, 1;
		break;
	}

	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> result;
	result.resize(control_size, 1);
	result = time * basis * P;
	coordinate.x() = result(0, 0);
	coordinate.y() = result(0, 1);
	if (is_3d) {
		coordinate.z() = result(0, 2);
	}
	std::cout << "\n\nTime Matrix\n" << time << "\n\nBasis Matrix\n" << basis 
		<< "\n\nControl Points\n" << P << "\n\nCalculated Coordinate\n" << coordinate;
	return coordinate;
}

}	// end namespace hpm
