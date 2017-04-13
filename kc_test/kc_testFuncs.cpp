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

		//// Works
		//CKSCoordArray elevated_points = bezier::ElevateDegree<double, bezier::Dimension::k3d, 3>(control_points);
		//for (size_t i = 0; i < elevated_points.size(); ++i) {
		//	part.AddPoint(elevated_points[i]);
		//}

		//// Works with calculate coefficients reversed
		//CKSCoordArray split_points = bezier::SplitSegment<double, bezier::k3d, 3>(control_points, 1, .5);
		//for (auto i : split_points) {
		//	part.AddPoint(i);
		//}
		//for (size_t i = 0; i < 2; ++i) {
		//	std::vector<double> new_coeff = bezier::CalculateCoefficients<double, bezier::k3d, 3>(split_points, i + 1);
		//	part.AddSpline(true, false, new_coeff);
		//}

		// Test tangent and normal functions. Works
		CKSCoord coordinate;
		CKSCoord tangent;
		CKSCoord normal;
		for (size_t i = 1; i <= control_points.size() / 4; ++i) {
			coordinate = bezier::CalculateCoordinate<double, bezier::k3d, 3>(control_points, i, .25);
			part.AddPoint(coordinate);
			tangent = bezier::CalculateTangent<double, bezier::k3d, 3>(control_points, i, .25);
			normal = bezier::CalculateNormal<double, bezier::k3d, 3>(control_points, i, .25);
			tangent.Normalize();
			normal.Normalize();
			CKSMatrix temp;
			CKSMath::MatrixVector(coordinate, coordinate + tangent, temp);
			CKEntityAttrib attrib;
			attrib.m_ucColorNumber = 7;
			part.AddVector(1.0, &temp, &attrib);
			attrib.m_ucColorNumber = 10;
			CKSMath::MatrixVector(coordinate, coordinate + normal, temp);
			part.AddVector(1.0, &temp, &attrib);
			//part.AddLine(coordinate, coordinate + tangent);
			//part.AddLine(coordinate, coordinate + normal);
		}
		part.NoteState();
		//WriteData("nodes.dat", nodePoints);
		//WriteCoefficients("coeff.dat", coeffs);
		//WriteControlPoints("ctrl.dat", control_points);
		//CKSCoordArray fitPoints;
		//SplineToPoints(part, spline, fitPoints);
		//WriteData("spline.dat", fitPoints);
		//fitPoints.clear();
	}
	return status;
}
