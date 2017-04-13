// Implementation for kc_test
//

#include "stdafx.h"
#include "funcs.h"
#include "ck_sdk.h"
#include "KCSdkUtilities.h"
#include "SMask.h"
#include "TestUtil.h"
#include "splines.h"

const double MIN_TOL = .00005;

int SplineData() {
	int status = CKNoError;
	CKPart part = CKGetActivePart();
	if (!part.IsValid()) {
		return CK_NO_PART;
	}
	CKSCoordArray ctrlPoints;   // Holds the 4 control points for entire segment
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
			GetSplineControlPoints(param, ctrlPoints);
		}
		CKSCoordArray elevated_points = bezier::ElevateDegree<double, bezier::Dimension::k3d, 3>(ctrlPoints);
		for (size_t i = 0; i < elevated_points.size(); ++i) {
			part.AddPoint(elevated_points[i]);
		}
		CKSCoordArray split_points = bezier::SplitSegment<double, bezier::k3d, 3>(ctrlPoints, 1, .75);
		std::vector<double> point_coefficients = bezier::CalculateCoefficients<double, bezier::k3d, 3>(ctrlPoints, 1);
		std::vector<double> new_coeff(12);
		for (size_t i = 0; i < point_coefficients.size() / 4; i += 4) {
			new_coeff[i] = point_coefficients[i];
			new_coeff[i + 1] = point_coefficients[i + 3];
			new_coeff[i + 2] = point_coefficients[i + 6];
			new_coeff[i + 3] = point_coefficients[i + 9];
		}
		part.AddSpline(true, false, new_coeff);
		point_coefficients = bezier::CalculateCoefficients<double, bezier::k3d, 3>(ctrlPoints, 2);
		for (size_t i = 0; i < point_coefficients.size() / 4; i += 4) {
			new_coeff[i] = point_coefficients[i];
			new_coeff[i + 1] = point_coefficients[i + 3];
			new_coeff[i + 2] = point_coefficients[i + 6];
			new_coeff[i + 3] = point_coefficients[i + 9];
		}
		part.AddSpline(true, false, new_coeff);
		part.NoteState();
		WriteData("nodes.dat", nodePoints);
		WriteCoefficients("coeff.dat", coeffs);
		WriteControlPoints("ctrl.dat", ctrlPoints);
		CKSCoordArray fitPoints;
		SplineToPoints(part, spline, fitPoints);
		WriteData("spline.dat", fitPoints);
		fitPoints.clear();
	}
	return status;
}
