/** @mainpage
* @author Paul J. Hentschel
*
* A collection of templated functions for working with Cubic Bézier curves.
*/

#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace bezier {

/** Used to specify the degree of the curve */
enum Degree {
	/** Degree = 2 */
	kQuadratic = 2,
	/** Degree = 3 */
	kCubic,
};

/** Used to specify the coordinate dimension */
enum Dimension {
	/** 2d coordinates \f$(x, y)\f$ */
	k2d = 2,
	/** 3d coordinates \f$(x, y, z)\f$ */
	k3d,
};

template <typename Scalar, Dimension, Degree, typename Point>
std::vector<Scalar> CalculateCoefficients(const std::vector<Point>& points);

template <typename Scalar, Dimension, Degree, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const long segment_id, const Scalar t);

/**	@brief	Calculate a vector of coefficients derived from a vector of control points
*			of a cubic Bézier spline.
*
*			Size of vector of control points must be in multiples of 4.
*
*	@param	degree Degree of the Bézier curve. Currently supports Quadratic(2) and Cubic(3) curves.
*	@param	points Vector of control points. This function will process the vector elements
*			in groups of 4.
*	@param	is_3d Boolean value to determine if 2D or 3D spline coefficients are returned.
*	@return	A vector of coefficients calculated from the control points.
*			Coefficients are stored in the order
*			{\f$A_x, B_x, C_x, D_x, A_y, B_y, C_y, D_y\f$} and if the points are 3d
*			{\f$A_z, B_z, C_z, D_z\f$}
*/

template <typename Scalar, Dimension dimension, Degree degree, typename Point>
std::vector<Scalar> CalculateCoefficients(const std::vector<Point>& points) {
	//const bool is_3d = (dimension == k3d);
	const size_t control_size = degree + 1;
	//Eigen::Matrix<Scalar, control_size, control_size> basis;
	//basis <<
	//	-1, 3, -3, 1,
	//	3, -6, 3, 0,
	//	-3, 3, 0, 0,
	//	1, 0, 0, 0;
	Eigen::Matrix<Scalar, control_size, control_size> basis;
	switch (degree) {
	case kCubic:
		basis <<
			-1, 3, -3, 1,
			3, -6, 3, 0,
			-3, 3, 0, 0,
			1, 0, 0, 0;
		break;
	case kQuadratic:
		basis <<
			1, -2, 1,
			-2, 2, 0,
			1, 0, 0;
		break;
	}

	Eigen::Matrix<Scalar, control_size, dimension> P;
	Eigen::Matrix<Scalar, control_size, dimension> result;
	Eigen::Matrix<Scalar, dimension, control_size> output;
	std::vector<Scalar> coefficients;
	const size_t num_segments = points.size() / control_size;	// Omits last set of control points if not a complete set.
	if (num_segments) {
		for (auto i = 0; i < num_segments * control_size; i += control_size) {
			for (auto j = 0; j < control_size; ++j) {
				P.block<1, dimension>(j, 0) = points[i + j];
			}
			result = basis * P;
			output = result.transpose();
			for (auto p = 0; p < output.rows(); ++p) {
				for (auto q = 0; q < output.cols(); ++q) {
					coefficients.push_back(output(p, q));
				}
			}
		}
	}
	//std::cout << basis << "\n\n";
	//basis.colwise().reverseInPlace();
	//std::cout << basis << "\n\n";
	return coefficients;
}


/**	@brief	Calculate a position on a Quadratic or Cubic Bézier curve segment
			using matrices

	@tparam Scalar Type of data being passed in. Valid types are float and double.
	@tparam Dimension Dimension of control point coordinates.
			Valid options are bezier::k2d \f$(x, y)\f$ or bezier::k3d \f$(x, y, z)\f$
	@tparam	Degree Degree of the Bézier curve. 
			Valid options are bezier::kQuadratic (2) and bezier::kCubic (3).
	@tparam Point Container type for the points. Must have [] accessor. [0] = x coordinate,
			[1] = y coordinate and if 3d, [2] = z coordinate.
	@param	segment_id Segment number to use for calculation.
	@param	t Time parameter in the interval \f$t \in [0,1]\f$ for a point on the curve segemnt.
			Values outside of this range may be used to calculate a coordinate on a natural
			extension of the curve.
	@param	points std::vector of control points. Number of control points 
			for each segment must be \f$ (degree + 1) \f$
	@return	A vector of coordinates of type \<**Point**\> for the calculated position.
*/

template <typename Scalar, Dimension dimension, Degree degree, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const long segment_id, const Scalar t) {
	using Eigen::Dynamic;
	const bool is_3d = (dimension == k3d);

	// Determine the size of a control point set for one Bézier curve segment,
	// then check if the input data is valid.
	const size_t control_size = degree + 1;
	const size_t num_segments = points.size() / control_size;	// Omits last set of control points if not a complete set.
	Point coordinate;
	if (num_segments == 0 || segment_id <= 0 || segment_id > num_segments) {
		return coordinate;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * control_size;
	Eigen::Matrix<Scalar, control_size, dimension> control_points;
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
	Eigen::Matrix<Scalar, dimension, 1> result;
	result = time * basis * control_points;
	coordinate[0] = result(0, 0);
	coordinate[1] = result(1, 0);
	if (is_3d) {
		coordinate[2] = result(2, 0);
	}
	return coordinate;
}

}	// end namespace bezier
