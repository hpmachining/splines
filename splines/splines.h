/**
@mainpage
A collection of templated functions for working with Bézier curves.



*/

#pragma once

#include "SplineHelper.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace bezier {

/** Specifies the coordinate dimension */
enum Dimension {
	/** 2d coordinates \f$(x, y)\f$ */
	k2d = 2,
	/** 3d coordinates \f$(x, y, z)\f$ */
	k3d
};


template <typename Scalar, Dimension dimension, size_t degree, typename Point>
std::vector<Scalar> CalculateCoefficients(const std::vector<Point>& points, const size_t segment_id);

template <typename Scalar, Dimension dimension, size_t degree, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const size_t segment_id, const Scalar t);

template <typename Scalar, Dimension dimension, size_t degree, typename Point>
Point CalculateTangent(const std::vector<Point>& points, const size_t segment_id, const Scalar t);

template <typename Scalar, Dimension dimension, size_t degree, typename Point>
Point CalculateNormal(const std::vector<Point>& points, const size_t segment_id, const Scalar t);

template <typename Scalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> SplitSegment(const std::vector<Point>& points, const size_t segment_id, const Scalar t);

/**
@brief	Calculate the coefficients derived from a segment of a composite Bézier curve.

		This function takes a set of 2d or 3d control points and a segment number of a 
		composite Bézier curve and returns the coefficients. Control point sets should be in 
		multiples of \f$degree + 1\f$.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	degree Degree of the Bézier curve.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		starts at 1.
@return	Coefficients calculated from the control points. Coefficients are stored in the order
		\f${C_{1_{x,y,(z)}},C_{2_{x,y,(z)}}, ..., C_{n_{x,y,(z)}}}\f$ where \f$n = degree + 1\f$.
		3d curves will have \f$n \cdot 3\f$ coefficients and 2d curves will have \f$n \cdot 2\f$.
*/
template <typename Scalar, Dimension dimension, size_t degree, typename Point>
std::vector<Scalar> CalculateCoefficients(const std::vector<Point>& points, const size_t segment_id) {
	const size_t order = degree + 1;
	std::vector<Scalar> coefficients;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return coefficients;
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<Scalar, order, order> basis;
	basis << GetPowerCoefficients<Scalar, degree>();
	basis.colwise().reverseInPlace();

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * order;
	Eigen::Matrix<Scalar, order, dimension> control_points;
	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}

	// Calculate the coefficients
	Eigen::Matrix<Scalar, order, dimension> result;
	result = basis * control_points;
	for (auto p = 0; p < result.rows(); ++p) {
		for (auto q = 0; q < result.cols(); ++q) {
			coefficients.push_back(result(p, q));
		}
	}

	return coefficients;
}

/**	
@brief	Calculate a position on a segment of a composite Bézier curve.

		This function calculates the coordinate at parameter \f$t\f$ on a Bézier curve segment.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	Degree Degree of the Bézier curve.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
		outside of this range may be used to calculate a coordinate on a natural extension of the 
		curve segment.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@return	Coordinate of type *Point* for the calculated position.
*/
template <typename Scalar, Dimension dimension, size_t degree, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const size_t segment_id, const Scalar t) {
	const size_t order = degree + 1;
	Point coordinate;
	//if (num_segments == 0 || segment_id < 1 || segment_id > num_segments) {
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return coordinate;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * order;
	Eigen::Matrix<Scalar, order, dimension> control_points;
	//for (size_t i = segment_index, j = 0; i < segment_index + order; ++i, ++j) {
	//	control_points.block<1, dimension>(j, 0) = points[i];
	//}

	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}

	// Create and fill the power basis (t) matrix
	Eigen::Matrix<Scalar, 1, order> parameter;
	for (size_t i = 0; i < order; ++i) {
		parameter(0, i) = std::pow(t, i);
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<Scalar, order, order> basis;
	basis << GetPowerCoefficients<Scalar, degree>();

	// Calculate the coordinate
	Eigen::Matrix<Scalar, dimension, 1> result;
	result = parameter * basis * control_points;
	for (size_t i = 0; i < dimension; ++i) {
		coordinate[i] = result(i, 0);
	}

	return coordinate;
}

/**
@brief	Calculate a tangent vector on a segment of a composite Bézier curve.

		This function calculates the tangent vector at parameter \f$t\f$ on a Bézier curve segment.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	Degree Degree of the Bézier curve. 
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values 
		outside of this range may be used to calculate a coordinate on a natural extension of the 
		curve segment.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename Scalar, Dimension dimension, size_t degree, typename Point>
Point CalculateTangent(const std::vector<Point>& points, const size_t segment_id, const Scalar t) {
	const size_t order = degree + 1;
	Point tangent;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return tangent;
	}

	// Get coeffecients of specified segment
	std::vector<Scalar> coefficients;
	coefficients = CalculateCoefficients<Scalar, dimension, degree>(points, segment_id);

	// Create and fill Eigen matrix with coefficients
	Eigen::Matrix<Scalar, degree, dimension> C;
	for (size_t i = 0, p = 0; p < degree; i += dimension, ++p) {
		for (size_t q = 0; q < dimension; ++q) {
			C(p, q) = coefficients[i + q];
		}
	}

	// Create and fill the power basis matrix
	Eigen::Matrix<Scalar, 1, degree> parameter;
	for (size_t i = 0; i < degree; ++i) {
		parameter(0, i) = std::pow(t, i);
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<Scalar, degree, degree> basis;
	basis << GetTangentCoefficients<Scalar, degree>();

	Eigen::Matrix<Scalar, dimension, 1> result;
	result = parameter * basis * C;
	for (size_t i = 0; i < dimension; ++i) {
		tangent[i] = result(i, 0);
	}

	return tangent;
}

/**
@brief	Calculate a normal vector on a segment of a composite Bézier curve.

This function calculates the normal vector at parameter \f$t\f$ on a Bézier curve segment.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	Degree Degree of the Bézier curve.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
outside of this range may be used to calculate a coordinate on a natural extension of the
curve segment.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename Scalar, Dimension dimension, size_t degree, typename Point>
Point CalculateNormal(const std::vector<Point>& points, const size_t segment_id, const Scalar t) {
	const size_t order = degree + 1;
	Point normal;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return normal;
	}

	// Get coordinate and tangent of curve at specified t parameter
	//Scalar next_t = t + .0001;
	Scalar next_t = t - .0001;
	Point coord_1 = CalculateCoordinate<Scalar, dimension, degree>(points, segment_id, t);
	Point coord_2 = CalculateCoordinate<Scalar, dimension, degree>(points, segment_id, next_t);
	Point tan_1 = CalculateTangent<Scalar, dimension, degree>(points, segment_id, t);
	Point tan_2 = CalculateTangent<Scalar, dimension, degree>(points, segment_id, next_t);

	// Copy to Eigen matrices
	Eigen::Matrix<Scalar, dimension, 1> tangent;
	Eigen::Matrix<Scalar, dimension, 1> next_tangent;
	Eigen::Matrix<Scalar, dimension, 1> coordinate;
	Eigen::Matrix<Scalar, dimension, 1> next_coordinate;

	for (size_t i = 0; i < dimension; ++i) {
		tangent(i, 0) = tan_1[i];
		next_tangent(i, 0) = tan_2[i];
		coordinate(i, 0) = coord_1[i];
		next_coordinate(i, 0) = coord_2[i];
	}
	Eigen::Matrix<Scalar, dimension, 1> offset = coordinate - next_coordinate;
	next_tangent += offset;
	tangent.normalize();
	next_tangent.normalize();
	//Eigen::Matrix<Scalar, dimension, 1> z_axis = next_tangent.cross(tangent);
	Eigen::Matrix<Scalar, dimension, 1> z_axis = tangent.cross(next_tangent);
	
	// Create rotation matrix
	Eigen::Transform<Scalar, dimension, Eigen::Affine> transform;
	transform.setIdentity();
	Eigen::AngleAxis<Scalar> rotation(M_PI_2, z_axis);
	transform.rotate(rotation);
	Eigen::Matrix<Scalar, dimension, 1> result;
	result = transform * tangent;
	result.normalize();
	for (size_t i = 0; i < dimension; ++i) {
		normal[i] = result(i, 0);
	}

	return normal;
}

/**
@brief	Divide a Bézier curve segment into 2 smaller segments.

This function divides a Bézier curve segment at parameter \f$t\f$ and returns 2 sets of control
points representing the subdivided segments.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	Degree Degree of the Bézier curve.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ which indicates the position to split the curve segment.
@param	points Control points of a composite Bézier. Number of control points for each segment 
		should be \f$degree + 1\f$
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename Scalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> SplitSegment(const std::vector<Point>& points, const size_t segment_id, const Scalar t) {
	const size_t order = degree + 1;
	std::vector<Point> split_segments;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return split_segments;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * order;
	Eigen::Matrix<Scalar, order, dimension> P;
	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			P(p, j) = points[i][j];
		}
	}

	// Create and fill the power basis matrices
	Eigen::Matrix<Scalar, order, order> Z;
	Z.setZero();
	Eigen::Matrix<Scalar, 1, order> power;
	for (size_t i = 0; i < order; ++i) {
		power(0, i) = std::pow(t, i);
		Z.diagonal()[i] = power(0, i);
	}

	Eigen::Matrix<Scalar, order, order> M = bezier::GetPowerCoefficients<Scalar, degree>();
	Eigen::Matrix<Scalar, order, order> M1 = M.inverse();
	Eigen::Matrix<Scalar, order, order> Q = M1 * Z * M;
	Eigen::Matrix<Scalar, order, order> Q1;
	Q1.setZero();
	for (size_t i = 0; i < order; ++i) {
		Q1.block(order - (i + 1), order - (i + 1), 1, i + 1) = Q.block(i, 0, 1, i + 1);
	}

	Eigen::Matrix<Scalar, order, dimension> result = Q * P;
	Point point;
	for (size_t i = 0; i < order; ++i) {
		for (size_t j = 0; j < dimension; ++j) {
			point[j] = result(i, j);
		}
		split_segments.push_back(point);
	}

	result = Q1 * P;
	for (size_t i = 0; i < order; ++i) {
		for (size_t j = 0; j < dimension; ++j) {
			point[j] = result(i, j);
		}
		split_segments.push_back(point);
	}

	return split_segments;
}

}	// end namespace bezier
