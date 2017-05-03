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
#include <fstream>
#include <chrono>
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>

const double MIN_TOL = .00005;
using namespace Eigen;
using timer = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
typedef Matrix<double, -1, -1> Points;
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

		//// Test old normal function using 2nd derivative
		//CKSCoord old_normal;
		//CKSCoord old_tangent;
		//CKSCoord old_position;

		//CKSMath::Evaluate(spline, NULL, false, false, false, 1, .0625, &splineMatrix, &old_position,
		//	&old_normal, &old_tangent);
		////GetSplineNormal(coeffs, 1, .5, old_normal);
		//if (old_normal.Magnitude() > .001) {
		//	old_normal.Normalize();
		//}
		//if (old_tangent.Magnitude() > .001) {
		//	old_tangent.Normalize();
		//}

		////GetSplineCoord(coeffs, 1, .5, old_position);
		//old_normal = old_position + old_normal;
		//old_tangent = old_position + old_tangent;
		//part.AddPoint(old_position);
		//part.AddPoint(old_normal);
		//part.AddPoint(old_tangent);
		//part.NoteState();

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
		//std::vector<CKSCoord2D> split_points;
		//for (size_t i = 0; i < points_2d.size() / (degree + 1); ++i) {
		//	split_points = bezier::SplitSegment<double>(points_2d, .5, i, degree, 2);
		//	for (auto j : split_points) {
		//		part.AddPoint(j.m_dX, j.m_dY, 0.0);
		//	}
		//}
		//part.NoteState();

		//// Create spline from calculated coefficients 2d
		//for (size_t i = 0; i < 2; ++i) {
		//	std::vector<double> new_coeff = bezier::CalculateCoefficients<double>(split_points, i, degree, bezier::k2d);
		//	std::vector<double> kc_coeff = bezier::ConvertCoefficientLayoutToKC<double>(new_coeff, degree, bezier::k2d);
		//	part.AddSpline(false, false, kc_coeff);
		//}
		//part.NoteState();

		//// Test tangent and normal functions 2d. Works
		//CKSCoord2D coordinate;
		//CKSCoord2D tangent;
		//CKSCoord2D normal;
		//for (size_t i = 0; i < points_2d.size() / 4; ++i) {
		//	coordinate = bezier::CalculateCoordinate<double>(points_2d, .25, i, degree, bezier::k2d);
		//	part.AddPoint(coordinate.m_dX, coordinate.m_dY, 0.0);
		//	tangent = bezier::CalculateFirstDerivative<double>(points_2d, .25, i, degree, bezier::k2d);
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
		//std::vector<CKSCoord> split_points;
		//for (size_t i = 0; i < control_points.size() / (degree + 1); ++i) {
		//	split_points = bezier::SplitSegment<double>(control_points, .5, i, degree, 3);
		//	for (auto j : split_points) {
		//		part.AddPoint(j);
		//	}
		//}
		//part.NoteState();

		// Create spline from calculated coefficients 3d. Works
		//for (size_t i = 0; i < 2; ++i) {
		//	std::vector<double> new_coeff = bezier::CalculateCoefficients<double>(split_points, i);
		//	std::vector<double> kc_coeff = bezier::ConvertCoefficientLayoutToKC<double>(new_coeff, degree, bezier::k3d);
		//	part.AddSpline(true, false, kc_coeff);
		//}
		//part.NoteState();

		// Test tangent and normal functions 3d. Works
		CKSCoord coordinate;
		CKSCoord tangent;
		CKSCoord normal;
		CKSCoord curvature;
		for (size_t i = 0; i < control_points.size() / 4; ++i) {
			for (auto j = 0; j < 5; ++j) {
				double t = j / 4.0;
				coordinate = bezier::CalculateCoordinate<double>(control_points, t, i, 3, bezier::k3d);
				part.AddPoint(coordinate);
				tangent = bezier::CalculateFirstDerivative<double>(control_points, t, i, 3, bezier::k3d);
				normal = bezier::CalculateNormal<double>(control_points, t, i, 3, bezier::k3d);
				curvature = bezier::CalculateSecondDerivative<double>(control_points, t, i, 3, bezier::k3d);
				tangent.Normalize();
				normal.Normalize();
				curvature.Normalize();
				CKSMatrix temp;
				CKSMath::MatrixVector(coordinate, coordinate + tangent, temp);
				CKEntityAttrib attrib;
				attrib.m_ucColorNumber = 7;
				part.AddVector(1.0, &temp, &attrib);
				attrib.m_ucColorNumber = 10;
				CKSMath::MatrixVector(coordinate, coordinate + normal, temp);
				part.AddVector(1.0, &temp, &attrib);
				attrib.m_ucColorNumber = 2;
				CKSMath::MatrixVector(coordinate, coordinate + curvature, temp);
				part.AddVector(1.0, &temp, &attrib);
			}
		}
		part.NoteState();

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

int SplineHelix()
{
  CWnd* pWnd = AfxGetMainWnd();
  int status = CKNoError;
  CKPart part = CKGetActivePart();
  if (!part.IsValid())
  {
    return CK_NO_PART;
  }

  CKSMatrix worldMat;
  CKSMask mask;
  mask.AddEntity(CKMaskLine);
  mask.AddEntity(CKMaskArc);
  mask.AddEntity(CKMaskSpline);
  mask.AddEntity(CKMaskNURBSpline);
  mask.AddEntity(CKMaskPolyline);
  mask.AddEntity(CKMaskEllipse);
  mask.AddEntity(CKMaskParabola);
  mask.AddEntity(CKMaskHyperbola);

  // Test create helical spline along a 3D curve
  Events keyCheck;
  double diameter = 1.0;
  double pitch = 0.25;
  CRegistry reg;
  if (reg.KeyExists(_T("Software\\HPM\\HPMTools\\SplineHelix")))
  {
    reg.SetKey(_T("Software\\HPM\\HPMTools\\SplineHelix"), FALSE);
    diameter = reg.ReadFloat(_T("Diameter"), 1.0);
    pitch = reg.ReadFloat(_T("Pitch"), .25);
  }
  else
  {
    reg.CreateKey(_T("Software\\HPM\\HPMTools\\SplineHelix"));
    reg.SetKey(_T("Software\\HPM\\HPMTools\\SplineHelix"), FALSE);
    reg.WriteFloat(_T("Diameter"), 1.0);
    reg.WriteFloat(_T("Pitch"), .25);
  }
  int step = 0;
  while (true)
  {
    switch (step)
    {
    case 0:
    {
      keyCheck = ck_get_input(_T("Enter diameter: "), _T(""), diameter, true, 0, CKS::GreaterThan, 0.0);
      switch (keyCheck)
      {
      case CKBackup:
      case CKEscape:
        return keyCheck;
      case CKNoError:
      {
        reg.WriteFloat(_T("Diameter"), diameter);
        step++;
        break;
      }
      default:
        return keyCheck;
      }
    }
    case 1:
    {
      keyCheck = ck_get_input(_T("Enter pitch: "), _T(""), pitch);
      switch (keyCheck)
      {
      case CKBackup:
      {
        step--;
        continue;
      }
      case CKEscape:
        return keyCheck;
      case CKNoError:
      {
        reg.WriteFloat(_T("Pitch"), pitch);
        step++;
        break;
      }
      default:
        return keyCheck;
      }
      if (CKSMath::CompareToZero(pitch, .01) <= 0)
      {
        pWnd->MessageBox(_T("The pitch entered is too small"), MB_TITLE, MB_OK_INFO);
        step--;
        continue;
      }
    }
    case 2:
    {
      CKSEntityArray driveCurves;
      status = part.GenSel(_T("Select the sweep path chain of curves"), driveCurves);
      switch (status)
      {
      case CKNoError:
        break;
      case CKBackup:
      {
        step--;
        continue;
      }
      default:
        if ((status < CKMenu1) || (status >= CKEscape))
          return status;
        //case CKEscape:
        //case CK_NO_PART:
        //  return status;
      }
      CKSCoordArray helixPnts;
      CKSCoord startVec, endVec;
      status = GetHelicalSplinePoints(part, driveCurves, helixPnts, startVec, endVec, diameter, pitch, .0001);
      if (helixPnts.size())
      {
        std::ofstream out_file("helix.points");
        for (auto p : helixPnts) {
          out_file << p.m_dX << '\t' << p.m_dY << '\t' << p.m_dZ << '\n';
        }
        time_point<timer> start_time = timer::now();
        size_t point_count = helixPnts.size();
        Points points(3, point_count);
        for (size_t i = 0; i < point_count; ++i) {
          points(0, i) = helixPnts[i][0];
          points(1, i) = helixPnts[i][1];
          points(2, i) = helixPnts[i][2];
        }
        Spline3d testspline = SplineFitting<Spline3d>::Interpolate(points, 3);
        size_t knots_size = testspline.knots().size();
        std::vector<double> knots(knots_size);
        for (size_t j = 0; j < knots_size; ++j) {
          knots[j] = testspline.knots().data()[j];
        }
        Matrix<double, -1, -1> control_points = testspline.ctrls();
        size_t ctrl_size = control_points.cols();
        std::vector<CKSCoord> ctrl_points(ctrl_size);
        for (size_t j = 0; j < ctrl_size; ++j) {
          ctrl_points[j].m_dX = control_points(0, j);
          ctrl_points[j].m_dY = control_points(1, j);
          ctrl_points[j].m_dZ = control_points(2, j);
        }

        CKSEntity helicalSpline;
        CKEntityAttrib attrib;
        attrib.m_ucColorNumber = 7;
        std::vector<double> weights(ctrl_size, 1.0);
        helicalSpline = part.AddNURBSpline(3, true, false, knots, ctrl_points, weights, &attrib);
        time_point<timer> end_time = timer::now();
        milliseconds elapsed = duration_cast<milliseconds>(end_time - start_time);
        CString eigen_time;
        eigen_time.Format(_T("Eigen Spline creation time: %d ms\n"), elapsed.count());
        part.NoteState();
        //helicalSpline = part.AddSpline(true, false, true, true, startVec, endVec, helixPnts, NULL, &worldMat);
        //helicalSpline = part.AddSpline(true, false, false, false, startVec, endVec, helixPnts, NULL, &worldMat);
        start_time = timer::now();
        attrib.m_ucColorNumber = 9;
        helicalSpline = part.AddNURBSpline(3, true, false, helixPnts, &attrib, &worldMat);
        end_time = timer::now();
        elapsed = duration_cast<milliseconds>(end_time - start_time);
        CString kc_time;
        kc_time.Format(_T("KC spline time: %d ms\n"), elapsed.count());
        if (!helicalSpline.IsValid()) {
          pWnd->MessageBox(_T("Error creating helical spline"), MB_TITLE, MB_OK_STOP);
          return CKError;
        }
        part.NoteState();
        CString all_time = eigen_time + kc_time;
        pWnd->MessageBox(all_time, MB_TITLE, MB_OK_STOP);
      }
    }
    }
  }
  return status;
}
