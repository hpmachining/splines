#include "splines.h"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

// @func Calculate a vector of coefficients derived from a vector of control points
// of a cubic Bezier spline.
// Size of vector of control points must be in multiples of 4. Coefficients are
// stored in the order Ax, Bx, Cx, Dx, Ay, By, Cy, Dy, and if 3d, Az, Bz, Cz, Dz.
// @parm Vector of control points. Size of vector must be a multiple of 4.
// @parm Boolean value to determine if 2D or 3D spline coefficients are returned.
// @flag TRUE   | Default value. Will return 3D spline coefficients. (X, Y, & Z coefficients for each 4 control points)
// @flag FALSE  | Will return 2D spline coefficients. (X & Y coefficients for each 4 control points)
// @rdesc Returns one of the following:
// @flag CK_NO_DATA | Point array not valid size for control points.
// @flag CKNoError  | Function completed succesfully.
std::vector<double> GetCoeffFromPoints(std::vector<Eigen::Vector3d> &points, bool is3D)
{
	std::vector <double> coeff;
	if (points.size() < 4 || points.size() % 4)
		return coeff;

	double Ax = 0.0, Bx = 0.0, Cx = 0.0, Dx = 0.0;
	double Ay = 0.0, By = 0.0, Cy = 0.0, Dy = 0.0;
	double Az = 0.0, Bz = 0.0, Cz = 0.0, Dz = 0.0;

	Eigen::Vector3d cp[4];
	for (size_t i = 0; i < points.size(); i += 4)
	{
		for (size_t j = 0; j < 4; ++j)
		{
			cp[j] = points[i + j];
		}

		Ax = cp[3].x() - 3.0 * cp[2].x() + 3.0 * cp[1].x() - cp[0].x();
		Ay = cp[3].y() - 3.0 * cp[2].y() + 3.0 * cp[1].y() - cp[0].y();
		Az = cp[3].z() - 3.0 * cp[2].z() + 3.0 * cp[1].z() - cp[0].z();

		Bx = 3.0 * cp[2].x() - 6.0 * cp[1].x() + 3.0 * cp[0].x();
		By = 3.0 * cp[2].y() - 6.0 * cp[1].y() + 3.0 * cp[0].y();
		Bz = 3.0 * cp[2].z() - 6.0 * cp[1].z() + 3.0 * cp[0].z();

		Cx = 3.0 * (cp[1].x() - cp[0].x());
		Cy = 3.0 * (cp[1].y() - cp[0].y());
		Cz = 3.0 * (cp[1].z() - cp[0].z());

		Dx = cp[0].x();
		Dy = cp[0].y();
		Dz = cp[0].z();

		coeff.push_back(Ax);
		coeff.push_back(Bx);
		coeff.push_back(Cx);
		coeff.push_back(Dx);

		coeff.push_back(Ay);
		coeff.push_back(By);
		coeff.push_back(Cy);
		coeff.push_back(Dy);

		if (is3D)
		{
			coeff.push_back(Az);
			coeff.push_back(Bz);
			coeff.push_back(Cz);
			coeff.push_back(Dz);
		}
	}
	return coeff;
}

