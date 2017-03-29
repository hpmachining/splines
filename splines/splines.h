/** @mainpage
* @author Paul J. Hentschel
*
* A collection of templated functions for working with Cubic B�zier curves.
*/

#pragma once

#include <vector>
#include <cmath>
#include <exception>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace bezier {

/** Used to specify the degree of the curve */
enum BezierDegree {
	/** Degree 2 */
	kQuadratic = 2,
	/** Degree 3 */
	kCubic,
};

enum Dimension {
	k2d = 2,
	k3d,
};

template <typename Point>
std::vector<double> CalculateCoefficientsFromPoints(const BezierDegree degree,
	const std::vector<Point>& points, bool is_3d);

template <typename Scalar, Dimension dimension, BezierDegree degree, typename Point>
Point CalculateCoordinate(const long segment_id, const Scalar t, const std::vector<Point>& points);

/**	@brief	Calculate a vector of coefficients derived from a vector of control points
*			of a cubic B�zier spline.
*
*			Size of vector of control points must be in multiples of 4.
*
*	@param	degree Degree of the B�zier curve. Currently supports Quadratic(2) and Cubic(3) curves.
*	@param	points Vector of control points. This function will process the vector elements
*			in groups of 4.
*	@param	is_3d Boolean value to determine if 2D or 3D spline coefficients are returned.
*	@return	A vector of coefficients calculated from the control points.
*			Coefficients are stored in the order
*			{\f$A_x, B_x, C_x, D_x, A_y, B_y, C_y, D_y\f$} and if the points are 3d
*			{\f$A_z, B_z, C_z, D_z\f$}
*/

template <typename Point>
std::vector<double> 
CalculateCoefficientsFromPoints(const BezierDegree degree, const std::vector<Point>& points, bool is_3d) {
	std::vector<double> coefficients;
	//if (degree != kQuadratic || degree != kCubic || ) {

	//}
	size_t increment = degree + 1;
	size_t num_control_points = points.size() / increment;

	if (num_control_points) {

		double Ax = 0.0, Bx = 0.0, Cx = 0.0, Dx = 0.0;
		double Ay = 0.0, By = 0.0, Cy = 0.0, Dy = 0.0;
		double Az = 0.0, Bz = 0.0, Cz = 0.0, Dz = 0.0;

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

			if (is_3d) {
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


/**	@brief	Calculate a position on a Quadratic or Cubic B�zier curve segment
			using matrices

	@tparam Scalar Type of data being passed in. Valid types are float and double.
	@tparam	degree Degree of the B�zier curve. Valid options are bezier::kQuadratic(2) and Cubic(3) curves.
	@tparam Point Container type for the points. Must have [] accessor. [0] = x coordinate,
			[1] = y coordinate and if 3d, [2] = z coordinate.
	@param	segment_id Segment number to use for calculation.
	@param	t Time parameter in the interval \f$t \in [0,1]\f$ for a point on the curve segemnt.
			Values outside of this range may be used to calculate a coordinate on a natural
			extension of the curve.
	@param	points Vector of control points. Number of control points 
			for each segment is \f$ degree + 1 \f$
	@param	is_3d Boolean value indicating if the coordinates provided for the
			control points are 3d \f$(x, y, z)\f$. Set to false for 2d \f$(c,y)\f$
	@return	A vector of coordinates of type Point for the calculated position.
*/

template <typename Scalar, Dimension dimension, BezierDegree degree, typename Point>
Point CalculateCoordinate(const long segment_id, const Scalar t, const std::vector<Point>& points) {
	using Eigen::Dynamic;
	const bool is_3d = (dimension == k3d);

	// Determine the size of a control point set for one B�zier curve segment,
	// then check if the input data is valid.
	const size_t control_size = degree + 1;
	const size_t num_segments = points.size() / control_size;	// Omits last set of control points if not a complete set.
	Point coordinate;
	if (num_segments == 0 || segment_id <= 0 || segment_id > num_segments) {
		return coordinate;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * control_size;
	Eigen::Matrix<Scalar, control_size, 3> control_points;
	if (is_3d) {
		for (size_t i = 0, p = segment_index; p < segment_index + control_size; ++i, ++p) {
			control_points(i, 0) = points[p][0];
			control_points(i, 1) = points[p][1];
			control_points(i, 2) = points[p][2];
		}
	}
	else {
		for (size_t i = 0, p = segment_index; p < segment_index + control_size; ++i, ++p) {
			control_points(i, 0) = points[p][0];
			control_points(i, 1) = points[p][1];
		}
	}
	
	// Create and fill the time matrix
	Eigen::Matrix<Scalar, 1, control_size> time;
	for (auto i = 0; i < control_size; ++i) {
		time(0, i) = std::pow(t, i);
	}
	
	// Create and fill the basis matrix
	Eigen::Matrix<Scalar, control_size, control_size> basis;
	switch (degree) {
	case kCubic:
		basis <<
			1,  0,  0, 0,
			-3,  3,  0, 0,
			3, -6,  3, 0,
			-1,  3, -3, 1;
		break;
	case kQuadratic:
		basis <<
			1, 0, 0,
			-2, 2, 0,
			1, -2, 1;
		break;
	}

	// Calculate the coordinate
	Eigen::Matrix<Scalar, control_size, 1> result;
	result = time * basis * control_points;
	coordinate[0] = result(0, 0);
	coordinate[1] = result(1, 0);
	if (is_3d) {
		coordinate[2] = result(2, 0);
	}
	return coordinate;
}

}	// end namespace bezier
