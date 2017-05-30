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
std::vector<RealScalar> GetCoefficients(const std::vector<Point>& points,
	const size_t segment_id = 0, const size_t degree = 3, const size_t dimension = k3d);

template <typename RealScalar, typename Point>
std::vector<Point> GetControlPoints(const std::vector<RealScalar>& coefficients,
  const size_t segment_id = 0, const size_t degree = 3, const size_t dimension = k3d);

template <typename RealScalar, typename Point>
Point GetPosition(const std::vector<Point>& points, const RealScalar t, 
	const size_t segment_id = 0, const size_t degree = 3, const size_t dimension = k3d);

template <typename RealScalar, typename Point>
Point GetFirstDerivative(const std::vector<Point>& points, const RealScalar t, 
	const size_t segment_id = 0, const size_t degree = 3, const size_t dimension = k3d);

template <typename RealScalar, typename Point>
Point GetSecondDerivative(const std::vector<Point>& points, const RealScalar t,
	const size_t segment_id = 0, const size_t degree = 3, const size_t dimension = k3d);

template <typename RealScalar, typename Point>
Point GetNormal(const std::vector<Point>& points, const RealScalar t,
	const size_t segment_id = 0, const size_t degree = 3, const size_t dimension = k3d);

template <typename RealScalar, typename Point>
Point Get2dNormal(const std::vector<Point>& points, const RealScalar t,
  const size_t segment_id = 0, const size_t degree = 3);

template <typename RealScalar, typename Point>
Point Get3dNormal(const std::vector<Point>& points, const RealScalar t,
  const size_t segment_id = 0, const size_t degree = 3);

template <typename RealScalar, typename Point>
std::vector<Point> SplitSegment(const std::vector<Point>& points, const RealScalar t, 
	const size_t segment_id = 0, size_t degree = 3, size_t dimension = k3d);

template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> ElevateDegree(const std::vector<Point>& points);

/**
@brief	Calculate the coefficients derived from a segment of a composite Bézier curve.

		This function takes a set of 2d or 3d control points and a segment number of a 
		composite Bézier curve and returns the coefficients. Control point sets should be in 
		multiples of \f$degree + 1\f$.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. 
        [0] = x, [1] = y, and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be 
        \f$degree + 1\f$
@param	segment_id Indicates which segment of the composite Bézier curve to process. 
        Numbering starts at 0.
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of control point coordinates. 
        Valid options are bezier::k2d or bezier::k3d
@return	Coefficients calculated from the control points. Order of stored coefficients is
		    \f$\left[C_{1,1}, C_{2,1}, ..., C_{n,1};\; C_{1,2}, C_{2,2}, ..., 
        C_{n,2};\; C_{1,m}, ..., C_{n,m}\right]\f$ where \f$n = degree + 1\f$ 
        and \f$m = dimension\f$. Curves will have \f$n \cdot m\f$ coefficients.
*/
template <typename RealScalar, typename Point>
std::vector<RealScalar> GetCoefficients(const std::vector<Point>& points, 
	const size_t segment_id, const size_t degree, const size_t dimension) {
	const size_t order = degree + 1;

  std::vector<RealScalar> coefficients;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return coefficients;
	}
	
  // Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = segment_id * order;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> control_points(order, dimension);
	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}
	
  // Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis(order, order);
	basis << GetPowerCoefficients<RealScalar>(degree);
	basis.colwise().reverseInPlace();
	
  // Calculate the coefficients
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> result(order, dimension);
	result = basis * control_points;
	
  // Add coefficients to std::vector
	size_t row_count = static_cast<size_t>(result.rows());
	size_t col_count = static_cast<size_t>(result.cols());
	for (size_t p = 0; p < col_count; ++p) {
		for (size_t q = 0; q < row_count; ++q) {
			coefficients.push_back(result(q, p));
		}
	}
	
  return coefficients;
}

/**
@brief	Calculate the control points derived from a segment of a composite Bézier curve.

This function takes a set of 2d or 3d coefficients and a segment number of a
composite Bézier curve and returns the control points for that segment. 
Coefficient sets should be in multiples of \f$\left(degree + 1\right) * dimension\f$ 
and in the order \f$\left[C_{1,1}, C_{2,1}, ..., C_{n,1};\; C_{1,2}, C_{2,2}, ...,
C_{n,2};\; C_{1,m}, ..., C_{n,m}\right]\f$ where \f$n = degree + 1\f$ and \f$m = dimension\f$.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the  control points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	coefficients Coefficients of the curve. Each segment should have 
        \f$\left(degree + 1\right) * dimension\f$ coefficients.
@param	segment_id Indicates which segment of the composite Bézier curve to process. 
        Numbering starts at 0.
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of coefficients. Valid options are bezier::k2d or bezier::k3d
@return	Control points calculated from the coefficients.
*/
template <typename RealScalar, typename Point>
std::vector<Point> GetControlPoints(const std::vector<RealScalar>& coefficients,
  const size_t segment_id, const size_t degree, const size_t dimension) {
  std::vector<Point> control_points;
  const size_t order = degree + 1;
  const size_t segment_size = order * dimension;
  const size_t segment_count = coefficients.size() / segment_size;
  
  if (segment_id >= segment_count) {
    return control_points;
  }
  const size_t segment_index = segment_id * segment_size;

  Eigen::Map<const Eigen::Matrix<RealScalar, Dynamic, Dynamic>>
    coeff_map(coefficients.data() + segment_index, order, dimension);

  // Create and fill the coefficients of the power basis matrix
  Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis = GetPowerCoefficients<RealScalar>(degree);
  basis.colwise().reverseInPlace();
  Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis_inverse = basis.inverse();
  
  // Calculate the control points
  Eigen::Matrix<RealScalar, Dynamic, Dynamic> result(order, dimension);
  result = basis_inverse * coeff_map;

  // Add control points to std::vector
  size_t row_count = static_cast<size_t>(result.rows());
  size_t col_count = static_cast<size_t>(result.cols());
  Point control_point;
  for (size_t i = 0; i < row_count; ++i) {
    for (size_t j = 0; j < col_count; ++j) {
      control_point[j] = result(i, j);
    }
    control_points.push_back(control_point);
  }

  return control_points;
}

/**
@brief	Calculate a position on a segment of a composite Bézier curve.

		This function calculates the coordinate at parameter \f$t\f$ on a Bézier curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y, and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
		outside of this range may be used to calculate a coordinate on a natural extension of the
		curve segment.
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		starts at 0.
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@return	Coordinate of type *Point* for the calculated position.
*/
template <typename RealScalar, typename Point>
Point GetPosition(const std::vector<Point>& points, const RealScalar t, 
	const size_t segment_id, const size_t degree, const size_t dimension){
	const size_t order = degree + 1;
	
  Point coordinate;
	coordinate[0] = std::numeric_limits<RealScalar>::quiet_NaN();
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return coordinate;
	}
	
  // Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = segment_id * order;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> control_points(order, dimension);
	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			control_points(p, j) = points[i][j];
		}
	}
	
  // Create and fill the power basis (t) matrix
	Eigen::Matrix<RealScalar, 1, Dynamic> parameter(1, order);
	for (size_t i = 0; i < order; ++i) {
		parameter(0, i) = std::pow(t, i);
	}
	
  // Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis(order, order);
	basis << GetPowerCoefficients<RealScalar>(degree);
	
  // Calculate the coordinate
	Eigen::Matrix<RealScalar, Dynamic, 1> result;
	result = parameter * basis * control_points;
	for (size_t i = 0; i < dimension; ++i) {
		coordinate[i] = result(i, 0);
	}
	
  return coordinate;
}

/**
@brief	Calculate a tangent vector on a segment of a composite Bézier curve.

		This function calculates the first derivative (tangent) at parameter \f$t\f$ on a Bézier 
		curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
		outside of this range may be used to calculate a coordinate on a natural extension of the
		curve segment.
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		starts at 0.
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename RealScalar, typename Point>
Point GetFirstDerivative(const std::vector<Point>& points, const RealScalar t,
	const size_t segment_id, const size_t degree, const size_t dimension) {
	const size_t order = degree + 1;
	
  Point tangent;
	tangent[0] = std::numeric_limits<RealScalar>::quiet_NaN();
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return tangent;
	}
	
  // Get coeffecients of specified segment
	std::vector<RealScalar> coefficients;
	coefficients = GetCoefficients<RealScalar>(points, segment_id, degree, dimension);
	
  // Create and fill Eigen matrix with coefficients
  //Eigen::Map<const Eigen::Matrix<RealScalar, Dynamic, Dynamic>>
  //  C(coefficients.data(), order, dimension);
  Eigen::Matrix<RealScalar, Dynamic, Dynamic> C;
  C.resize(order, dimension);
  C = Eigen::Map<const Eigen::Matrix<RealScalar, Dynamic, Dynamic>>(coefficients.data(), order, dimension);
  C = C.block(0, 0, degree, dimension).eval();
  // Create and fill the power basis matrix
	Eigen::Matrix<RealScalar, 1, Dynamic> parameter(1, degree);
	for (size_t i = 0; i < degree; ++i) {
		parameter(0, i) = std::pow(t, i);
	}

  // Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis(degree, degree);
	basis << GetTangentCoefficients<RealScalar>(degree);

  // Calculate first derivative (tangent)
	Eigen::Matrix<RealScalar, Dynamic, 1> result;
	result = parameter * basis * C;
	for (size_t i = 0; i < dimension; ++i) {
		tangent[i] = result(i, 0);
	}
	
  return tangent;
}

/**
@brief	Calculate the curvature on a segment of a composite Bézier curve.

		This function calculates the second derivative (curvature) at parameter \f$t\f$ on a Bézier 
		curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
		outside of this range may be used to calculate a coordinate on a natural extension of the
		curve segment.
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		starts at 0.
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@return	Coordinate of type *Point* for the calculated second derivative.
*/
template <typename RealScalar, typename Point>
Point GetSecondDerivative(const std::vector<Point>& points, const RealScalar t,
	const size_t segment_id, const size_t degree, const size_t dimension) {
	const size_t order = degree + 1;
	const size_t matrix_size = degree - 1;
	
  Point curvature;
	curvature[0] = std::numeric_limits<RealScalar>::quiet_NaN();
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return curvature;
	}
	
  // Get coeffecients of specified segment
	std::vector<RealScalar> coefficients;
	coefficients = GetCoefficients<RealScalar>(points, segment_id, degree, dimension);
	
  // Create and fill Eigen matrix with coefficients
  Eigen::Map<const Eigen::Matrix<RealScalar, Dynamic, Dynamic>>
    C(coefficients.data(), degree, dimension);
	
  // Create and fill the power basis matrix
	Eigen::Matrix<RealScalar, 1, Dynamic> parameter(1, matrix_size);
	for (size_t i = 0; i < matrix_size; ++i) {
		parameter(0, i) = std::pow(t, i);
	}
	
  // Create and fill the coefficients of the power basis matrix
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis(matrix_size, matrix_size);
	basis << GetSecondDerivativeCoefficients<RealScalar>(degree);
	
  // Calculate the second derivative (curvature)
	Eigen::Matrix<RealScalar, Dynamic, 1> result;
	result = parameter * basis * C.topRows(matrix_size);
	for (size_t i = 0; i < dimension; ++i) {
		curvature[i] = result(i, 0);
	}
	
  return curvature;
}

/**
@brief	Calculate a normal vector on a segment of a composite Bézier curve.

		This function calculates the normal vector at parameter \f$t\f$ on a Bézier curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
		outside of this range may be used to calculate a coordinate on a natural extension of the
		curve segment.
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		starts at 0.
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename RealScalar, typename Point>
Point GetNormal(const std::vector<Point>& points, const RealScalar t,
	const size_t segment_id, const size_t degree, const size_t dimension) {
	const size_t order = degree + 1;
	
  Point normal;
	normal[0] = std::numeric_limits<RealScalar>::quiet_NaN();
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return normal;
	}

  switch (dimension) {
  case k2d:
    normal = Get2dNormal(points, t, segment_id, degree);
    break;
  case k3d:
    normal = Get3dNormal(points, t, segment_id, degree);
  }

	return normal;
}

/**
@brief	Get the 3d normal vector on a segment of a composite Bézier curve.

This function returns the 3d normal vector at parameter \f$t\f$ on a Bézier curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and [2] = z.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
outside of this range may be used to calculate a coordinate on a natural extension of the
curve segment.
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
starts at 0.
@param	degree Degree of the Bézier curve.
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename RealScalar, typename Point>
Point Get3dNormal(const std::vector<Point>& points, const RealScalar t,
  const size_t segment_id, const size_t degree) {
  const size_t dimension = bezier::k3d;
  const size_t order = degree + 1;

  Point normal;
  normal[0] = std::numeric_limits<RealScalar>::quiet_NaN();
  if (!IsSegmentDataValid(points, order, segment_id)) {
    return normal;
  }

  // Get coordinate and tangent of curve at specified t parameter
  Point coord_1 = GetPosition<RealScalar>(points, t, segment_id, degree, dimension);
  Point tan_1 = GetFirstDerivative<RealScalar>(points, t, segment_id, degree, dimension);
  
  // Get coordinate and tangent close to t parameter
  RealScalar next_t = t - .0001;
  Point coord_2 = GetPosition<RealScalar>(points, next_t, segment_id, degree, dimension);
  Point tan_2 = GetFirstDerivative<RealScalar>(points, next_t, segment_id, degree, dimension);
  
  // Copy to Eigen matrices
  Eigen::Matrix<RealScalar, 3, 1> tangent = Eigen::Matrix<RealScalar, 3, 1>::Zero();
  Eigen::Matrix<RealScalar, 3, 1> next_tangent = Eigen::Matrix<RealScalar, 3, 1>::Zero();
  Eigen::Matrix<RealScalar, 3, 1> coordinate = Eigen::Matrix<RealScalar, 3, 1>::Zero();
  Eigen::Matrix<RealScalar, 3, 1> next_coordinate = Eigen::Matrix<RealScalar, 3, 1>::Zero();
  for (size_t i = 0; i < dimension; ++i) {
    tangent(i, 0) = tan_1[i];
    next_tangent(i, 0) = tan_2[i];
    coordinate(i, 0) = coord_1[i];
    next_coordinate(i, 0) = coord_2[i];
  }

  // Calculate rotation axis
  Eigen::Matrix<RealScalar, 3, 1> offset;
  offset = coordinate - next_coordinate;
  next_tangent += offset;
  tangent.normalize();
  next_tangent.normalize();
  Eigen::Matrix<RealScalar, 3, 1> z_axis;
  z_axis = tangent.cross(next_tangent);
  z_axis.normalize();

  // Create rotation matrix and rotate tangent vector
  Eigen::Transform<RealScalar, 3, Eigen::Affine> transform;
  transform.setIdentity();
  Eigen::AngleAxis<RealScalar> rotation(M_PI_2, z_axis);
  transform.rotate(rotation);
  Eigen::Matrix<RealScalar, 3, 1> result;
  result = transform * tangent;
  result.normalize();
  for (size_t i = 0; i < dimension; ++i) {
    normal[i] = result(i, 0);
  }

  return normal;
}

/**
@brief	Get the 2d normal vector on a segment of a composite Bézier curve.

This function returns the 2d normal vector at parameter \f$t\f$ on a Bézier curve segment.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y.
@param	points Control points. Number of control points for each segment should be \f$degree + 1\f$
@param	t Parameter in the range \f$0 \le t \le 1\f$ for a point on the curve segment. Values
outside of this range may be used to calculate a coordinate on a natural extension of the
curve segment.
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
starts at 0.
@param	degree Degree of the Bézier curve.
@return	Coordinate of type *Point* for the calculated tangent vector.
*/
template <typename RealScalar, typename Point>
Point Get2dNormal(const std::vector<Point>& points, const RealScalar t,
  const size_t segment_id, const size_t degree) {
  const size_t dimension = bezier::k2d;
  const size_t order = degree + 1;

  Point normal;
  normal[0] = std::numeric_limits<RealScalar>::quiet_NaN();
  if (!IsSegmentDataValid(points, order, segment_id)) {
    return normal;
  }

  // Get coordinate and tangent of curve at specified t parameter
  Point coord_1 = GetPosition<RealScalar>(points, t, segment_id, degree, dimension);
  Point tan_1 = GetFirstDerivative<RealScalar>(points, t, segment_id, degree, dimension);

  // Copy to Eigen matrices
  Eigen::Matrix<RealScalar, 2, 1> tangent;
  Eigen::Matrix<RealScalar, 2, 1> coordinate;
  for (size_t i = 0; i < dimension; ++i) {
    tangent(i, 0) = tan_1[i];
    coordinate(i, 0) = coord_1[i];
  }
  tangent.normalize();

  // Rotate tangent vector
  Eigen::Matrix<RealScalar, 2, 1> result = Eigen::Rotation2D<RealScalar>::Rotation2D(M_PI_2) * tangent;
  result.normalize();
  for (size_t i = 0; i < dimension; ++i) {
    normal[i] = result(i, 0);
  }

  return normal;
}

/**
@brief	Divide a Bézier curve segment into 2 smaller segments.

This function divides a Bézier curve segment at parameter \f$t\f$ and returns 2 sets of control
points representing the divided segments.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	points Control points of a composite Bézier. Number of control points for each segment
		    should be \f$degree + 1\f$
@param	t Parameter in the range \f$0 \le t \le 1\f$ which indicates the position to split the curve segment.
@param	segment_id Indicates which segment of the composite Bézier curve to process. Numbering
		    starts at 0.
@param	degree Degree of the Bézier curve.
@param	dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@return	Two sets of control points representing the divided segments.
*/
template <typename RealScalar, typename Point>
std::vector<Point> SplitSegment(const std::vector<Point>& points, const RealScalar t,
	const size_t segment_id, size_t degree, size_t dimension) {
	const size_t order = degree + 1;
	
  std::vector<Point> split_segments;
	if (!IsSegmentDataValid(points, order, segment_id)) {
		return split_segments;
	}

	// Create and fill Eigen::Matrix with control points for specified segment
	const size_t segment_index = segment_id * order;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> P;
	P.resize(order, dimension);
	for (size_t i = segment_index, p = 0; i < segment_index + order; ++i, ++p) {
		for (size_t j = 0; j < dimension; ++j) {
			P(p, j) = points[i][j];
		}
	}
	
  // Create and fill the power basis matrices
	Eigen::Matrix<RealScalar, 1, Dynamic> Z;
	Z.resize(Eigen::NoChange, order);
	for (size_t i = 0; i < order; ++i) {
		Z(0, i) = std::pow(t, i);
	}

	Eigen::Matrix<RealScalar, Dynamic, Dynamic> M;
	M.resize(order, order);
	M = bezier::GetPowerCoefficients<RealScalar>(degree);
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> Q;
	Q.resize(order, order);
	Q = M.inverse() * Z.asDiagonal() * M;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> Q1;
	Q1.resize(order, order);
	Q1.setZero();
	for (size_t i = 0; i < order; ++i) {
		Q1.block(order - (i + 1), order - (i + 1), 1, i + 1) = Q.block(i, 0, 1, i + 1);
	}

	Eigen::Matrix<RealScalar, Dynamic, Dynamic> result;
	result.resize(order, dimension);
	result = Q * P;
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

/**
@brief	Elevate degree of a composite Bézier curve by one.

@tparam RealScalar Type of data being passed in. Valid types are float and double.
@tparam	dimension Dimension of control point coordinates. Valid options are bezier::k2d or bezier::k3d
@tparam	degree Degree of the inputed composite Bézier curve.
@tparam Point Container type for the points. Must have [] accessor. [0] = x, [1] = y and if 3d, [2] = z.
@param	points Control points of a composite Bézier. Number of control points for each segment
        should be \f$degree + 1\f$
@return	Control points of the elevated composite curve.
*/
template <typename RealScalar, Dimension dimension, size_t degree, typename Point>
std::vector<Point> ElevateDegree(const std::vector<Point>& points) {
	const size_t order = degree + 1;
	const size_t segment_count = points.size() / order;

	std::vector<Point> elevated_points;
	Eigen::Matrix<double, order, dimension> P;
	for (size_t i = 0; i < segment_count; ++i) {
		const size_t segment_index = i * order;

    // Fill Eigen::Matrix with control points for current segment
		for (size_t j = segment_index, p = 0; j < segment_index + order; ++j, ++p) {
			for (size_t k = 0; k < dimension; ++k) {
				P(p, k) = points[j][k];
			}
		}
		// First point stays the same so add to elevated points
		Point point;
		for (size_t j = 0; j < dimension; ++j) {
			point[j] = P(0, j);
		}
		elevated_points.push_back(point);

    // Calculate the new elevated control points
		Eigen::Matrix<double, degree, 1> M1(Eigen::Matrix<double, degree, 1>::LinSpaced(degree, (1.0 / order),
			static_cast<double>(degree) / order));
		Eigen::Matrix<double, degree, dimension> Q =
			(M1.asDiagonal() * P.topRows<degree>()) + (M1.reverse().asDiagonal() * P.bottomRows<degree>());

    // Add new control points to elevated points
		for (size_t j = 0; j < degree; ++j) {
			for (size_t k = 0; k < dimension; ++k) {
				point[k] = Q(j, k);
			}
			elevated_points.push_back(point);
		}

    // Last point stays the same so add to elevated points
		for (size_t j = 0; j < dimension; ++j) {
			point[j] = P(degree, j);
		}

    elevated_points.push_back(point);
	}
	
  return elevated_points;
}

}	// end namespace bezier
