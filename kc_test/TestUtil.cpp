#include "stdafx.h"

//#include "TestUtil.h"
#include "ck_sdk.h"
#include "SMask.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>

CKSEntity SplineSelect(CKPart& part) {
    CKSEntity entity;
    CKSDrawInst nullInst;
    CKSMask mask;
    mask.AddEntity(CKMaskSpline);
    mask.AddEntity(CKMaskNURBSpline);
    part.GetEnt(_T("Select Spline"), entity, nullInst);
    return entity;
}

CKSEntityArray CurvesSelect(CKPart& part) {
  CKSEntityArray curves;
  CKSDrawInst nullInst;
  CKSMask mask;
  mask.AddEntity(CKMaskLine);
  mask.AddEntity(CKMaskArc);
  mask.AddEntity(CKMaskSpline);
  mask.AddEntity(CKMaskNURBSpline);
  mask.AddEntity(CKMaskPolyline);
  mask.AddEntity(CKMaskEllipse);
  mask.AddEntity(CKMaskParabola);
  mask.AddEntity(CKMaskHyperbola);
  part.GenSel(_T("Select the sweep path chain of curves"), curves);
  mask.Clear();
  return curves;
}

void WriteData(const std::string& dataFile, const std::vector<CKSCoord>& points) {
  std::ofstream data(dataFile);
  if (!data.is_open()) return;
  data << std::setprecision(15);
  for (auto coord : points) {
    data << coord.m_dX << " " << coord.m_dY << " " << coord.m_dZ << '\n';
  }
}

void WriteControlPoints(const std::string& dataFile, const std::vector<CKSCoord>& points) {
  std::ofstream data(dataFile);
  if (!data.is_open()) return;
  data << std::setprecision(15);
  int counter = 1;
  for (auto coord : points) {
    data << coord.m_dX << " " << coord.m_dY << " " << coord.m_dZ << '\n';
    if (counter == 4) {
      data << '\n';
      counter = 0;
    }
    ++counter;
  }
}

void WriteCoefficients(const std::string& dataFile, const std::vector<double>& points) {
  std::ofstream data(dataFile);
  if (!data.is_open()) return;
  data << std::setprecision(15);
  for (auto coeff : points) {
    data << coeff << '\n';
  }
}
void WriteCoefficients(const std::string& dataFile, const std::vector<std::complex<double>>& points) {
  std::ofstream data(dataFile);
  if (!data.is_open()) return;
  data << std::setprecision(15);
  for (auto coeff : points) {
    data << coeff << '\n';
  }
}
