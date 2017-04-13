#pragma once

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace bezier {

//template <typename Point>
//bool IsSegmentDataValid(const std::vector<Point>& points, size_t order, const size_t segment_id);
//
//template <typename Scalar, size_t degree>
//Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> GetPowerCoefficients();
//
//template <typename Scalar, size_t degree>
//Eigen::Matrix<Scalar, degree, degree> GetTangentCoefficients();
//
//template <typename Scalar, size_t degree>
//Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> GetBinomialCoefficients();


template <typename Point>
bool IsSegmentDataValid(const std::vector<Point>& points, size_t order, const size_t segment_id) {
	// Omits the last set of control points if not a complete set.
	const size_t num_segments = points.size() / order;
	return  (num_segments > 0 && segment_id > 0 && segment_id <= num_segments);
}

template <typename Scalar, size_t degree>
Eigen::Matrix<Scalar, degree, degree> GetTangentCoefficients() {
	Eigen::Matrix<Scalar, degree, degree> basis;
	Eigen::DiagonalMatrix<Scalar, degree> coefficients;
	for (auto i = 0; i < degree; ++i) {
		coefficients.diagonal()[i] = i + 1;
	}
	basis = coefficients;
	basis.rowwise().reverseInPlace();
	return basis;
}

template <typename Scalar, size_t degree>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> GetPowerCoefficients() {
	const size_t order = degree + 1;
	Eigen::Matrix<Scalar, order, order> coefficients = GetBinomialCoefficients<Scalar, degree>();
	for (int i = 1; i < order; ++i) {
		for (int j = 0; j < order - i; ++j) {
			Scalar element = -(coefficients(i + j, j + 1) / i * (j + 1));
			coefficients(i + j, j) = element;
		}
	}
	return coefficients;
}

template<typename Scalar, size_t degree>
Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> GetBinomialCoefficients() {
	const size_t order = degree + 1;
	Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> coefficients;
	coefficients.resize(order);
	coefficients.diagonal()[0] = 1;
	for (int i = 0; i < degree; ++i) {
		coefficients.diagonal()[i + 1] = coefficients.diagonal()[i] * (degree - i) / (i + 1);
	}
	return coefficients;
}
} // end namespace bezier