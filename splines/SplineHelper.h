#pragma once

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

template<typename Scalar, Dimension dimension, size_t degree>
std::vector<Scalar> ConvertCoefficientLayoutToKC(const std::vector<Scalar>& coefficients) {
	std::vector<Scalar> sorted;
	const size_t order = degree + 1;
	const size_t segment_size = (degree + 1) * dimension;
	if (coefficients.size() % segment_size != 0) {
		return sorted;
	}
	//sorted.resize(coefficients.size());
	for (size_t i = 0; i < coefficients.size() / segment_size; ++i) {
		Eigen::Matrix<Scalar, order, dimension> segment_coefficients;
		for (size_t j = 0; j < segment_size / dimension; ++j) {
			for (size_t k = 0; k < dimension; ++k) {
				segment_coefficients(j, k) = coefficients[(i * segment_size) + (j * dimension) + k];
			}
		}
		for (size_t p = 0; p < segment_coefficients.cols(); ++p) {
			for (size_t q = 0; q < segment_coefficients.rows(); ++q) {
				sorted.push_back(segment_coefficients(q, p));
			}
		}
	}
	return sorted;
}
} // end namespace bezier