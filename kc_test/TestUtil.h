#pragma once

#include "ck_sdk.h"

#define MB_TITLE (_T("Hello!"))
#define MB_OK_STOP (MB_OK | MB_ICONSTOP)
#define MB_OK_WARN (MB_OK | MB_ICONWARNING)
#define MB_OK_INFO (MB_OK | MB_ICONASTERISK)
#define MB_OK_QUES (MB_OK | MB_ICONQUESTION)


CKSEntity SplineSelect(CKPart& part);
CKSEntityArray CurvesSelect(CKPart& part);
void WriteData(const std::string& dataFile, const std::vector<CKSCoord>& points);
void WriteControlPoints(const std::string& dataFile, const std::vector<CKSCoord>& points);
void WriteCoefficients(const std::string& dataFile, const std::vector<double>& points);