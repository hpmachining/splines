/**
@mainpage
A collection of templated functions for working with Quadratic and Cubic Bézier curves.



*/

#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace bezier {

/** Specifies the degree of the curve */
enum Degree {
	/** Degree = 2 */
	kQuadratic = 2,
	/** Degree = 3 */
	kCubic,
};

/** Specifies the coordinate dimension */
enum Dimension {
	/** 2d coordinates \f$(x, y)\f$ */
	k2d = 2,
	/** 3d coordinates \f$(x, y, z)\f$ */
	k3d,
};

template <typename Scalar>
Eigen::Matrix<Scalar, 4, 4> GetCubicBasis();

template <typename Scalar>
Eigen::Matrix<Scalar, 3, 3> GetQuadraticBasis();

template <typename Scalar, Dimension, Degree, typename Point>
std::vector<Scalar> CalculateCoefficients(const std::vector<Point>& points, const size_t segment_id);

template <typename Scalar, Dimension, Degree, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const size_t segment_id, const Scalar t);

template <typename Scalar, Dimension, Degree, typename Point>
Point CalculateTangent(const std::vector<Point>& points, const size_t segment_id, const Scalar t);

template <typename Scalar, Dimension, Degree, typename Point>
Point CalculateNormal(const std::vector<Point>& points, const size_t segment_id, const Scalar t);


/** 
@return	A 4 x 4 matrix of the coefficients of the cubic power basis.
*/
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 4> GetCubicBasis() {
	Eigen::Matrix<Scalar, 4, 4> basis;
	basis <<
		1, 0, 0, 0,
		-3, 3, 0, 0,
		3, -6, 3, 0,
		-1, 3, -3, 1;
	return basis;
}

/**
@return	A 3 x 3 matrix of the coefficients of the quadratic power basis.
*/
template <typename Scalar>
Eigen::Matrix<Scalar, 3, 3> GetQuadraticBasis() {
	Eigen::Matrix<Scalar, 3, 3> basis;
	basis <<
		1, 0, 0,
		-2, 2, 0,
		1, -2, 1;
	return basis;
}

/**	
@brief	Calculate the coefficients derived from a segment of a composite Bézier curve.

		This function takes a set of 2d or 3d control points and a segment number of a 
		quadratic or cubic composite Bézier curve and returns the coefficients.
		Control point sets must be in multiples of \f$degree + 1\f$.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam dimension Dimension of control point coordinates.
		Valid options are bezier::k2d or bezier::k3d
@tparam	degree Degree of the Bézier curve.
		Valid options are bezier::kQuadratic and bezier::kCubic.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment must be \f$degree + 1\f$
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		starts at 1.
@return	Coefficients calculated from the control points. Coefficients are stored in the order
		\f${A_{x,y,(z)},B_{x,y,(z)}, C_{x,y,(z)}, D_{x,y,(z)}}\f$. 3d curves will have 12
		coefficients and 2d curves will have 8.
*/
template <typename Scalar, Dimension dimension, Degree degree, typename Point>
std::vector<Scalar> CalculateCoefficients(const std::vector<Point>& points, const size_t segment_id) {
	// Determine the size of a control point set for one Bézier curve segment,
	// then check if the input data is valid.
	const size_t control_size = degree + 1;
	// Omits the last set of control points if not a complete set.
	const size_t num_segments = points.size() / control_size;	 
	std::vector<Scalar> coefficients;
	if (num_segments == 0 || segment_id < 1 || segment_id > num_segments) {
		return coefficients;
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<Scalar, control_size, control_size> basis;
	switch (degree) {
	case kCubic:
		basis << GetCubicBasis<Scalar>();
		basis.colwise().reverseInPlace();
		break;
	case kQuadratic:
		basis << GetQuadraticBasis<Scalar>();
		basis.colwise().reverseInPlace(); 
		break;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * control_size;
	Eigen::Matrix<Scalar, control_size, dimension> control_points;
	//for (size_t i = segment_index, j = 0; i < segment_index + control_size; ++i, ++j) {
	//	control_points.block<1, dimension>(j, 0) = points[i];
	//}

	for (size_t i = segment_index, p = 0; i < segment_index + control_size; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}

	// Calculate the coefficients
	Eigen::Matrix<Scalar, control_size, dimension> result;
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

		This function calculates the coordinate at parameter \f$t\f$ on a quadratic or cubic
		Bézier curve segment. It is an implementation of the equations
		\f$B(t) = (1-t)^3P_0 + 3(1-t)^2tP_1 + 3(1-t)t^2P_2 + t^3P_3,\; 0 \le t \le 1\f$ for cubic 
		curves and \f$B(t) = (1-t)^2P_0 + 2(1-t)tP_1 + t^2P_2,\; 0 \le t \le 1\f$ for quadratic 
		curves.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	Degree Degree of the Bézier curve. Valid options are bezier::kQuadratic and bezier::kCubic.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
		outside of this range may be used to calculate a coordinate on a natural extension of the 
		curve segment.
@param	points Control points. Number of control points for each segment must be \f$degree + 1\f$
@return	Coordinate of type *Point* for the calculated position.
*/
template <typename Scalar, Dimension dimension, Degree degree, typename Point>
Point CalculateCoordinate(const std::vector<Point>& points, const size_t segment_id, const Scalar t) {
	// Determine the size of a control point set for one Bézier curve segment, then check if the
	// input data is valid.
	const size_t control_size = degree + 1;
	// Omits last set of control points if not a complete set.
	const size_t num_segments = points.size() / control_size;
	Point coordinate;
	if (num_segments == 0 || segment_id < 1 || segment_id > num_segments) {
		return coordinate;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = (segment_id - 1) * control_size;
	Eigen::Matrix<Scalar, control_size, dimension> control_points;
	//for (size_t i = segment_index, j = 0; i < segment_index + control_size; ++i, ++j) {
	//	control_points.block<1, dimension>(j, 0) = points[i];
	//}

	for (size_t i = segment_index, p = 0; i < segment_index + control_size; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}

	// Create and fill the power basis (t) matrix
	Eigen::Matrix<Scalar, 1, control_size> parameter;
	for (size_t i = 0; i < control_size; ++i) {
		parameter(0, i) = std::pow(t, i);
	}

	// Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<Scalar, control_size, control_size> basis;
	switch (degree) {
	case kCubic:
		basis << GetCubicBasis<Scalar>();
		break;
	case kQuadratic:
		basis << GetQuadraticBasis<Scalar>();
		break;
	}

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

		This function calculates the tangent vector at parameter \f$t\f$ on a quadratic or cubic
		Bézier curve segment.

@tparam Scalar Type of data being passed in. Valid types are float and double.
@tparam Dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	Degree Degree of the Bézier curve. Valid options are bezier::kQuadratic and bezier::kCubic.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	segment_id Segment number to use for calculation.
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values 
		outside of this range may be used to calculate a coordinate on a natural extension of the 
		curve segment.
@param	points Control points. Number of control points for each segment must be \f$degree + 1\f$
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename Scalar, Dimension dimension, Degree degree, typename Point>
Point CalculateTangent(const std::vector<Point>& points, const size_t segment_id, const Scalar t) {
	// Determine the size of a control point set for one Bézier curve segment, then check if the
	// input data is valid.
	const size_t control_size = degree + 1;
	// Omits last set of control points if not a complete set.
	const size_t num_segments = points.size() / control_size;
	Point tangent;
	if (num_segments == 0 || segment_id < 1 || segment_id > num_segments) {
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
	switch (degree) {
	case kCubic:
		basis <<
			0, 0, 1,
			0, 2, 0,
			3, 0, 0;
		break;
	case kQuadratic:
		basis <<
			0, 1,
			2, 0;
		break;
	}

	Eigen::Matrix<Scalar, dimension, 1> result;
	result = parameter * basis * C;
	for (size_t i = 0; i < dimension; ++i) {
		tangent[i] = result(i, 0);
	}

	return tangent;
}

template <typename Scalar, Dimension dimension, Degree degree, typename Point>
Point CalculateNormal(const std::vector<Point>& points, const size_t segment_id, const Scalar t) {
	// Determine the size of a control point set for one Bézier curve segment, then check if the
	// input data is valid.
	const size_t control_size = degree + 1;
	// Omits last set of control points if not a complete set.
	const size_t num_segments = points.size() / control_size;
	Point normal;
	if (num_segments == 0 || segment_id < 1 || segment_id > num_segments) {
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
	//normal = transform * tangent;
	//normal.normalize();
	//switch (degree) {
	//case kCubic:
	//	//transform.rotate(Eigen::AngleAxis<Scalar>(.5 * M_PI, Eigen::Matrix<Scalar, dimension, 1>::UnitZ()));
	//	//normal = transform * tangent;
	//	break;
	//case kQuadratic:
	//	break;
	//}

	return normal;
}

}	// end namespace bezier
