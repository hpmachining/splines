// Implementation for kc_test
//

#include "stdafx.h"
#include <cmath>
#include "funcs.h"
#include "ck_sdk.h"
#include "KCSdkUtilities.h"
#include "SMask.h"
#include "TestUtil.h"
#include "splines.h"

const double MIN_TOL = .00005;

int TestSplineLibrary() {
	int status = CKNoError;
	CKPart part = CKGetActivePart();
	if (!part.IsValid()) {
		return CK_NO_PART;
	}
	CKSCoordArray control_points;   // Holds the 4 control points for entire segment
	CKSEntity spline = SplineSelect(part);
	if (!spline.IsValid()) {
		MessageBox(nullptr, _T("Error with spline selection"), _T("Spline Data"), MB_OK_STOP);
		status = CKError;
	}
	else {
		std::vector <double> coeffs;
		bool is3D = false;
		bool isClosed = false;
		HPMatrix splineMatrix;
		CKSCoordArray nodePoints;
		CKSCoord startVector, endVector;
		part.GetSpline(spline, NULL, coeffs, is3D, isClosed, NULL, &splineMatrix);
		part.GetSpline(spline, NULL, nodePoints, &startVector, &endVector, is3D, isClosed, NULL, &splineMatrix);
		size_t blockSize = 0;
		if (is3D) {
			blockSize = 12;
		}
		else {
			blockSize = 8;
		}
		double param[3][4] = { 0.0 };   // Array to hold current segment coefficients
		for (size_t i = 0; i < coeffs.size(); i += blockSize) {
			if (is3D) {
				for (size_t j = 0; j < 3; ++j) {
					for (size_t k = 0; k < 4; ++k) {
						param[j][k] = coeffs[i + (j * 4) + k];
					}
				}
			}
			else {  // 2D Spline
				for (size_t j = 0; j < 2; ++j) {
					for (size_t k = 0; k < 4; ++k) {
						param[j][k] = coeffs[i + (j * 4) + k];
					}
				}
			}
			GetSplineControlPoints(param, control_points);
		}

		//// Test 2D spline points
		//CKSCoord2D cp_2d;
		//std::vector<CKSCoord2D> points_2d;
		//for (auto i : control_points) {
		//	cp_2d.m_dX = i.m_dX;
		//	cp_2d.m_dY = i.m_dY;
		//	points_2d.push_back(cp_2d);
		//}
		//const size_t degree = 3;

		//// Test elevate degree 2d. Works
		//std::vector<CKSCoord2D> elevated_points = bezier::ElevateDegree<double, bezier::Dimension::k2d, degree>(points_2d);
		//for (auto i: elevated_points) {
		//	part.AddPoint(i.m_dX, i.m_dY, 0.0);
		//}
		//part.NoteState();

		//// Test segment split 2d. Works.
		//std::vector<CKSCoord2D> split_points = bezier::SplitSegment<double>(points_2d, .5, 1, degree, 2);
		//for (auto i : split_points) {
		//	part.AddPoint(i.m_dX, i.m_dY, 0.0);
		//}
		//part.NoteState();

		//// Create spline from calculated coefficients 2d
		//for (size_t i = 0; i < 2; ++i) {
		//	std::vector<double> new_coeff = bezier::CalculateCoefficients<double>(split_points, i + 1, degree, bezier::k2d);
		//	std::vector<double> kc_coeff = bezier::ConvertCoefficientLayoutToKC<double>(new_coeff, degree, bezier::k2d);
		//	part.AddSpline(false, false, kc_coeff);
		//}
		//part.NoteState();

		//// Test tangent and normal functions 2d. Works
		//CKSCoord2D coordinate;
		//CKSCoord2D tangent;
		//CKSCoord2D normal;
		//for (size_t i = 1; i <= points_2d.size() / 4; ++i) {
		//	coordinate = bezier::CalculateCoordinate<double>(points_2d, .25, i, degree, bezier::k2d);
		//	part.AddPoint(coordinate.m_dX, coordinate.m_dY, 0.0);
		//	tangent = bezier::CalculateTangent<double>(points_2d, .25, i, degree, bezier::k2d);
		//	normal = bezier::CalculateNormal<double>(points_2d, .25, i, degree, bezier::k2d);
		//	CKSMatrix temp;
		//	CKSCoord v1(coordinate.m_dX, coordinate.m_dY, 0.0);
		//	CKSCoord v2(tangent.m_dX, tangent.m_dY, 0.0);
		//	CKSCoord v3(normal.m_dX, normal.m_dY, 0.0);
		//	CKSMath::MatrixVector(v1, v1 + v2, temp);
		//	CKEntityAttrib attrib;
		//	attrib.m_ucColorNumber = 7;
		//	part.AddVector(1.0, &temp, &attrib);
		//	attrib.m_ucColorNumber = 10;
		//	CKSMath::MatrixVector(v1, v1 + v3, temp);
		//	part.AddVector(1.0, &temp, &attrib);
		//}
		//part.NoteState();

		//// Test elevate degree 3d. Works
		//CKSCoordArray elevated_points = bezier::ElevateDegree<double, bezier::Dimension::k3d, 3>(control_points);
		//for (size_t i = 0; i < elevated_points.size(); ++i) {
		//	part.AddPoint(elevated_points[i]);
		//}
		//part.NoteState();

		// Test segment split 3d. Works.
		const size_t degree = 3;
		CKSCoordArray split_points = bezier::SplitSegment<double>(control_points, .5, 1);
		for (auto i : split_points) {
			part.AddPoint(i);
		}
		part.NoteState();

		//// Create spline from calculated coefficients 3d. Works
		//for (size_t i = 0; i < 2; ++i) {
		//	std::vector<double> new_coeff = bezier::CalculateCoefficients<double>(split_points, i + 1);
		//	std::vector<double> kc_coeff = bezier::ConvertCoefficientLayoutToKC<double>(new_coeff, degree, bezier::k3d);
		//	part.AddSpline(true, false, kc_coeff);
		//}
		//part.NoteState();

		//// Test tangent and normal functions 3d. Works
		//CKSCoord coordinate;
		//CKSCoord tangent;
		//CKSCoord normal;
		//for (size_t i = 1; i <= control_points.size() / 4; ++i) {
		//	coordinate = bezier::CalculateCoordinate<double>(control_points, .25, i, 3, bezier::k3d);
		//	part.AddPoint(coordinate);
		//	tangent = bezier::CalculateTangent<double>(control_points, .25, i, 3, bezier::k3d);
		//	normal = bezier::CalculateNormal<double>(control_points, .25, i, 3, bezier::k3d);
		//	tangent.Normalize();
		//	normal.Normalize();
		//	CKSMatrix temp;
		//	CKSMath::MatrixVector(coordinate, coordinate + tangent, temp);
		//	CKEntityAttrib attrib;
		//	attrib.m_ucColorNumber = 7;
		//	part.AddVector(1.0, &temp, &attrib);
		//	attrib.m_ucColorNumber = 10;
		//	CKSMath::MatrixVector(coordinate, coordinate + normal, temp);
		//	part.AddVector(1.0, &temp, &attrib);
		//}
		//part.NoteState();

		WriteData("nodes.dat", nodePoints);
		WriteCoefficients("coeff.dat", coeffs);
		WriteControlPoints("ctrl.dat", control_points);
		CKSCoordArray fitPoints;
		SplineToPoints(part, spline, fitPoints);
		WriteData("spline.dat", fitPoints);
		fitPoints.clear();
	}
	return status;
}
