#pragma once

#include <vector>
#include <iostream>
//#include <limits>
//#include <complex>
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

using Eigen::Dynamic;

template <typename Point>
bool IsSegmentDataValid(const std::vector<Point>& points, size_t order, const size_t segment_id);

template <typename RealScalar>
Eigen::Matrix<RealScalar, Dynamic, Dynamic> GetTangentCoefficients(const size_t degree);

template <typename RealScalar>
Eigen::Matrix<RealScalar, Dynamic, Dynamic> GetSecondDerivativeCoefficients(const size_t degree);
  
template <typename RealScalar>
Eigen::Matrix<RealScalar, Dynamic, Dynamic> GetPowerCoefficients(const size_t degree);

template <typename RealScalar>
Eigen::DiagonalMatrix<RealScalar, Dynamic> GetBinomialCoefficients(const size_t degree);

template<typename RealScalar>
std::vector<RealScalar> ConvertCoefficientLayoutToKC(const std::vector<RealScalar>& coefficients,
  const size_t degree, const size_t dimension);

template<typename RealScalar>
int SolveQuadratic(const std::vector<RealScalar>& coefficients, std::vector<RealScalar>& solutions);

template<typename RealScalar>
RealScalar SolveLinear(const std::vector<RealScalar>& coefficients);

template <typename Point>
bool IsSegmentDataValid(const std::vector<Point>& points, const size_t order, const size_t segment_id) {
	const size_t num_segments = points.size() / order;
	return  (num_segments > 0 && segment_id < num_segments);
}

template <typename RealScalar>
Eigen::Matrix<RealScalar, Dynamic, Dynamic> GetTangentCoefficients(const size_t degree) {
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis;
	basis.resize(degree, degree);
	Eigen::DiagonalMatrix<RealScalar, Dynamic> coefficients;
	coefficients.resize(degree);
	for (size_t i = 0; i < degree; ++i) {
		coefficients.diagonal()[i] = static_cast<RealScalar>(i + 1);
	}
	basis = coefficients;
	basis.rowwise().reverseInPlace();
	return basis;
}

template <typename RealScalar>
Eigen::Matrix<RealScalar, Dynamic, Dynamic> GetSecondDerivativeCoefficients(const size_t degree) {
	const size_t matrix_size = degree - 1;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> basis;
	basis.resize(matrix_size, matrix_size);
	Eigen::DiagonalMatrix<RealScalar, Dynamic> coefficients;
	coefficients.resize(matrix_size);
	for (size_t i = 0; i < degree - 1; ++i) {
		coefficients.diagonal()[i] = static_cast<RealScalar>((degree - i) * (degree - (i + 1)));
	}
	basis = coefficients;
	basis.rowwise().reverseInPlace();
	return basis;
}

template <typename RealScalar>
Eigen::Matrix<RealScalar, Dynamic, Dynamic> GetPowerCoefficients(const size_t degree) {
	const size_t order = degree + 1;
	Eigen::Matrix<RealScalar, Dynamic, Dynamic> coefficients;
	coefficients.resize(order, order);
	coefficients = GetBinomialCoefficients<RealScalar>(degree);
	for (size_t i = 1; i < order; ++i) {
		for (size_t j = 0; j < order - i; ++j) {
			RealScalar element = -(coefficients(i + j, j + 1) / i * (j + 1));
			coefficients(i + j, j) = element;
		}
	}
  return coefficients;
}

template<typename RealScalar>
Eigen::DiagonalMatrix<RealScalar, Dynamic> GetBinomialCoefficients(const size_t degree) {
	const size_t order = degree + 1;
	Eigen::DiagonalMatrix<RealScalar, Dynamic> coefficients;
	coefficients.resize(order);
	coefficients.diagonal()[0] = 1;
	for (size_t i = 0; i < degree; ++i) {
		coefficients.diagonal()[i + 1] = coefficients.diagonal()[i] * (degree - i) / (i + 1);
	}
  return coefficients;
}

template<typename RealScalar>
std::vector<RealScalar> ConvertCoefficientLayoutToKC(const std::vector<RealScalar>& coefficients,
	const size_t degree, const size_t dimension) {
	std::vector<RealScalar> sorted;
	const size_t order = degree + 1;
	const size_t segment_size = (degree + 1) * dimension;
	if (coefficients.size() % segment_size != 0) {
		return sorted;
	}
	for (size_t i = 0; i < coefficients.size() / segment_size; ++i) {
		Eigen::Matrix<RealScalar, Dynamic, Dynamic> segment_coefficients;
		segment_coefficients.resize(order, dimension);
		for (size_t j = 0; j < segment_size / dimension; ++j) {
			for (size_t k = 0; k < dimension; ++k) {
				segment_coefficients(j, k) = coefficients[(i * segment_size) + (j * dimension) + k];
			}
		}
		for (size_t p = 0; p < dimension; ++p) {
			for (size_t q = 0; q < order; ++q) {
				sorted.push_back(segment_coefficients(q, p));
			}
		}
	}
  return sorted;
}

template<typename RealScalar>
int SolveQuadratic(const std::vector<RealScalar>& coefficients, std::vector<RealScalar>& solutions) {
  if (coefficients.size() != 3) {
    return 1;
  }
  RealScalar a = coefficients[0];
  RealScalar b = coefficients[1];
  RealScalar c = coefficients[2];
  RealScalar square_root = std::sqrt(b * b - 4 * a * c);
  RealScalar t = (-b + square_root) / (2 * a);
  solutions.push_back(t);
  t = (-b - square_root) / (2 * a);
  solutions.push_back(t);
  return 0;
}

template<typename RealScalar>
RealScalar SolveLinear(const std::vector<RealScalar>& coefficients) {
  RealScalar a = coefficients[0];
  RealScalar b = coefficients[1];
  RealScalar t = -b / a;
  return t;
}

} // end namespace bezier