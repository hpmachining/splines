#include "stdafx.h"

//#include "TestUtil.h"
#include "ck_sdk.h"
#include "SMask.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

CKSEntity SplineSelect(CKPart& part) {
    CKSEntity entity;
    CKSDrawInst nullInst;
    CKSMask mask;
    mask.AddEntity(CKMaskSpline);
    part.GetEnt(_T("Select Spline"), entity, nullInst);
    return entity;
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
