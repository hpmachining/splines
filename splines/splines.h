/**
@mainpage
A collection of templated functions for working with Bézier curves.



*/

#pragma once

#include "SplineHelper.h"
#include <cmath>
#include <limits>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace bezier {


template <typename RealScalar, typename Point>
std::vector<RealScalar> CalculateCoefficients(const std::vector<Point>& points,
	const size_t segment_id = 1, const size_t degree = 3, const size_t dimension = k3d);

template <typename RealScalar, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const RealScalar t, 
	const size_t segment_id = 1, const size_t degree = 3, const size_t dimension= k3d);

template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
Point CalculateTangent(const std::vector<Point>& points, const size_t segment_id, const RealScalar t);

template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
Point CalculateNormal(const std::vector<Point>& points, const size_t segment_id, const RealScalar t);

template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> SplitSegment(const std::vector<Point>& points, const size_t segment_id, const RealScalar t);

template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> ElevateDegree(const std::vector<Point>& points);

/**
@brief	Calculate the coefficients derived from a segment of a composite Bézier curve.

		This function takes a set of 2d or 3d control points and a segment number of a 
		composite Bézier curve and returns the coefficients. Control point sets should be in 
		multiples of \f$degree + 1\f$.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		starts at 1.
@return	Coefficients calculated from the control points. Coefficients are stored in the order
		\f${C_{1_{x,y,(z)}},C_{2_{x,y,(z)}}, ..., C_{n_{x,y,(z)}}}\f$ where \f$n = degree + 1\f$.
		3d curves will have \f$n \cdot 3\f$ coefficients and 2d curves will have \f$n \cdot 2\f$.
*/
template <typename RealScalar, typename Point>
std::vector<RealScalar> CalculateCoefficients(const std::vector<Point>& points, 
	const size_t segment_id, const size_t degree, const size_t dimension) {
	const size_t order = degree + 1;
	std::vector<RealScalar> coefficients;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return coefficients;
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis;
	basis.resize(order, order);
	basis << GetPowerCoefficients<RealScalar>(degree);
	basis.colwise().reverseInPlace();

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * order;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> control_points;
	control_points.resize(order, dimension);
	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}

	// Calculate the coefficients
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> result;
	result.resize(order, dimension);
	result = basis * control_points;
	size_t row_count = static_cast<size_t>(result.rows());
	size_t col_count = static_cast<size_t>(result.cols());
	for (size_t p = 0; p < row_count; ++p) {
		for (size_t q = 0; q < col_count; ++q) {
			coefficients.push_back(result(p, q));
		}
	}

	return coefficients;
}

/**	
@brief	Calculate a position on a segment of a composite Bézier curve.

		This function calculates the coordinate at parameter \f$t\f$ on a Bézier curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	Degree Degree of the Bézier curve.
@param	Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
		outside of this range may be used to calculate a coordinate on a natural extension of the
		curve segment.
@return	Coordinate of type *Point* for the calculated position.
*/
template <typename RealScalar, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const RealScalar t, 
	const size_t segment_id, const size_t degree, const size_t dimension){
	const size_t order = degree + 1;
	Point coordinate;
	coordinate[0] = std::numeric_limits<RealScalar>::quiet_NaN();
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return coordinate;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * order;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> control_points;
	control_points.resize(order, dimension);
	//for (size_t i = segment_index, j = 0; i < segment_index + order; ++i, ++j) {
	//	control_points.block<1, dimension>(j, 0) = points[i];
	//}

	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}

	// Create and fill the power basis (t) matrix
	Eigen::Matrix<RealScalar, 1, Dynamic> parameter;
	parameter.resize(Eigen::NoChange, order);
	for (size_t i = 0; i < order; ++i) {
		parameter(0, i) = std::pow(t, i);
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis;
	basis.resize(order, order);
	basis << GetPowerCoefficients<RealScalar>(degree);

	// Calculate the coordinate
	Eigen::Matrix<RealScalar, Dynamic, 1> result;
	result.resize(dimension, Eigen::NoChange);
	result = parameter * basis * control_points;
	for (size_t i = 0; i < dimension; ++i) {
		coordinate[i] = result(i, 0);
	}

	return coordinate;
}

/**
@brief	Calculate a tangent vector on a segment of a composite Bézier curve.

		This function calculates the tangent vector at parameter \f$t\f$ on a Bézier curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
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
template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
Point CalculateTangent(const std::vector<Point>& points, const size_t segment_id, const RealScalar t) {
	const size_t order = degree + 1;
	Point tangent;
	tangent[0] = std::numeric_limits<RealScalar>::quiet_NaN();
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return tangent;
	}

	// Get coeffecients of specified segment
	std::vector<RealScalar> coefficients;
	coefficients = CalculateCoefficients<RealScalar>(points, degree, dimension, segment_id);

	// Create and fill Eigen matrix with coefficients
	Eigen::Matrix<RealScalar, degree, dimension> C;
	for (size_t i = 0, p = 0; p < degree; i += dimension, ++p) {
		for (size_t q = 0; q < dimension; ++q) {
			C(p, q) = coefficients[i + q];
		}
	}

	// Create and fill the power basis matrix
	Eigen::Matrix<RealScalar, 1, degree> parameter;
	for (size_t i = 0; i < degree; ++i) {
		parameter(0, i) = std::pow(t, i);
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<RealScalar, degree, degree> basis;
	basis << GetTangentCoefficients<RealScalar>(degree);

	Eigen::Matrix<RealScalar, dimension, 1> result;
	result = parameter * basis * C;
	for (size_t i = 0; i < dimension; ++i) {
		tangent[i] = result(i, 0);
	}

	return tangent;
}

/**
@brief	Calculate a normal vector on a segment of a composite Bézier curve.

This function calculates the normal vector at parameter \f$t\f$ on a Bézier curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
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
template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
Point CalculateNormal(const std::vector<Point>& points, const size_t segment_id, const RealScalar t) {
	const size_t order = degree + 1;
	Point normal;
	normal[0] = std::numeric_limits<RealScalar>::quiet_NaN();
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return normal;
	}

	// Get coordinate and tangent of curve at specified t parameter
	//RealScalar next_t = t + .0001;
	RealScalar next_t = t - .0001;
	Point coord_1 = CalculateCoordinate<RealScalar>(points, t, segment_id, degree, dimension);
	Point coord_2 = CalculateCoordinate<RealScalar>(points, next_t, segment_id, degree, dimension);
	Point tan_1 = CalculateTangent<RealScalar, dimension, degree>(points, segment_id, t);
	Point tan_2 = CalculateTangent<RealScalar, dimension, degree>(points, segment_id, next_t);

	// Copy to Eigen matrices
	Eigen::Matrix<RealScalar, dimension, 1> tangent;
	Eigen::Matrix<RealScalar, dimension, 1> next_tangent;
	Eigen::Matrix<RealScalar, dimension, 1> coordinate;
	Eigen::Matrix<RealScalar, dimension, 1> next_coordinate;

	for (size_t i = 0; i < dimension; ++i) {
		tangent(i, 0) = tan_1[i];
		next_tangent(i, 0) = tan_2[i];
		coordinate(i, 0) = coord_1[i];
		next_coordinate(i, 0) = coord_2[i];
	}
	Eigen::Matrix<RealScalar, dimension, 1> offset = coordinate - next_coordinate;
	next_tangent += offset;
	tangent.normalize();
	next_tangent.normalize();
	//Eigen::Matrix<RealScalar, dimension, 1> z_axis = next_tangent.cross(tangent);
	Eigen::Matrix<RealScalar, dimension, 1> z_axis = tangent.cross(next_tangent);
	
	// Create rotation matrix
	Eigen::Transform<RealScalar, dimension, Eigen::Affine> transform;
	transform.setIdentity();
	Eigen::AngleAxis<RealScalar> rotation(M_PI_2, z_axis);
	transform.rotate(rotation);
	Eigen::Matrix<RealScalar, dimension, 1> result;
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

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	Degree Degree of the Bézier curve.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ which indicates the position to split the curve segment.
@param	points Control points of a composite Bézier. Number of control points for each segment 
		should be \f$degree + 1\f$
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> SplitSegment(const std::vector<Point>& points, const size_t segment_id, const RealScalar t) {
	const size_t order = degree + 1;
	std::vector<Point> split_segments;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return split_segments;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * order;
	Eigen::Matrix<RealScalar, order, dimension> P;
	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			P(p, j) = points[i][j];
		}
	}

	// Create and fill the power basis matrices
	Eigen::Matrix<RealScalar, order, order> Z;
	Z.setZero();
	Eigen::Matrix<RealScalar, 1, order> power;
	for (size_t i = 0; i < order; ++i) {
		power(0, i) = std::pow(t, i);
		Z.diagonal()[i] = power(0, i);
	}

	Eigen::Matrix<RealScalar, order, order> M = bezier::GetPowerCoefficients<RealScalar>(degree);
	Eigen::Matrix<RealScalar, order, order> M1 = M.inverse();
	Eigen::Matrix<RealScalar, order, order> Q = M1 * Z * M;
	Eigen::Matrix<RealScalar, order, order> Q1;
	Q1.setZero();
	for (size_t i = 0; i < order; ++i) {
		Q1.block(order - (i + 1), order - (i + 1), 1, i + 1) = Q.block(i, 0, 1, i + 1);
	}

	Eigen::Matrix<RealScalar, order, dimension> result = Q * P;
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

template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> ElevateDegree(const std::vector<Point>& points) {
	const size_t order = degree + 1;

	// Create and fill Eigen::Matrix with control points for specified segment
	std::vector<Point> elevated_points;
	Eigen::Matrix<double, order, dimension> P;
	for (size_t i = 0; i < points.size() / order; ++i) {
		const size_t segment_index = i * order;
		//P = Eigen::Map<Eigen::Matrix<double, order, dimension, Eigen::RowMajor>>
		//	(points[segment_index].data(), order, dimension);
		for (size_t j = segment_index, p = 0; j < segment_index + order; ++j, ++p) {
			for (size_t k = 0; k < dimension; ++k) {
				P(p, k) = points[j][k];
			}
		}

		// First point stays the same so add to elevated points
		//Point point = P.block(0, 0, 1, dimension).transpose();
		Point point;
		for (size_t j = 0; j < dimension; ++j) {
			point[j] = P(0, j);
		}
		elevated_points.push_back(point);

		// Calculate the new control points
		Eigen::Matrix<double, degree, 1> M1(Eigen::Matrix<double, degree, 1>::LinSpaced(degree, (1.0 / order),
			static_cast<double>(degree) / order));
		Eigen::Matrix<double, degree, dimension> Q =
			(M1.asDiagonal() * P.topRows<degree>()) + (M1.reverse().asDiagonal() * P.bottomRows<degree>());

		// Add new control points to elevated points
		for (size_t j = 0; j < degree; ++j) {
			//point = Q.block(j, 0, 1, dimension).transpose();
			for (size_t k = 0; k < dimension; ++k) {
				point[k] = Q(j, k);
			}
			elevated_points.push_back(point);
		}

		// Last point stays the same so add to elevated points
		for (size_t j = 0; j < dimension; ++j) {
			point[j] = P(degree, j);
		}
		//point = P.block(degree, 0, 1, dimension).transpose();
		elevated_points.push_back(point);
	}
	return elevated_points;
}

}	// end namespace bezier
