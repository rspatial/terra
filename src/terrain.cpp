/* Robert Hijmans, October 2011 */

using namespace std;
#include <vector>
#include <algorithm>
#include <math.h>
#include "distance.h"

double get_TRI (double v[], int n, double x) {
	double s = 0;
	for (int i=0; i<n; i++) {
		s = s + fabs(v[i] - x);
	}
	s = s / n;
	return(s);
}

double get_TPI (double v[], int n, double x) {
	double s = 0;
	if (n > 0) {
		for (int i=0; i<n; i++) {
			s = s + v[i];
		}
		s = s / n;
		s = x - s;
	}
	return(s);
}

double get_roughness (double v[], int n) {
	double min = v[0];
	double max = v[0];
	for (int i=1; i < n; i++) {
		if (v[i] < min) {
			min = v[i];
		} else if (v[i] > max) {
			max = v[i];
		}
	}
	return(max - min);
}


std::vector<std::vector<double> > terrain_indices(std::vector<std::vector<double> > &m, std::vector<std::vector<bool> > f, std::string option) {

	int nrow = m.size();
	int ncol = m[0].size();
	std::vector<std::vector<double>> v (nrow, std::vector<double>(ncol, NAN));

	int opt;
	if (option == "TPI") {
		// TPI (Topographic Position Index)
		// difference between the value of a cell and the mean value of its 8 surrounding cells.
		opt = 0;
	} else if (option == "TRI") {
		// TRI (Terrain Ruggedness Index)
		//  mean of the absolute differences between the value of a cell and the value of its 8 surrounding cells
		opt = 1;
	} else { //(option == 'roughness')
		// Roughness
		// difference between the maximum and the minimum value of a cell and its 8 surrounding cells.
		opt = 2;
	}

	int r, c, n;
	double va[9];
	for (int row = 0; row < (nrow-1); row++) {
		for (int col = 0; col < (ncol-1); col++) {
			n = 0;
			for(int i = -2; i < 3; i++) {
				r = i + row;
				if (r < 0 || r > (nrow-1)) {
					continue;
				}
				for(int j = -2; j < 3; j++) {
					c = j + col;
					if (c < 0 || c > (ncol-1)) {
						continue;
					}
					// center not included for TPI and TRI
					if (r==0 && c==0 && opt != 2) {
						continue;
					}
					if (! isnan(m[r][c]) ) {
						va[n] = m[r][c];
						n++;
					}
				}
			}
			if (opt==0) {
				v[row][col] = get_TPI(va, n, m[row][col]);
			} else if (opt==1) {
				v[row][col] = get_TRI(va, n, m[row][col]);
			} else { // rough
				v[row][col] = get_roughness(va, n);
			}
		}
	}
	return(v);
}




std::vector<std::vector<double> > slope4lonlat(std::vector<std::vector<double> > d, double dx, double dy, double ymax, int unit) {

	double a = 6378137;
	double f = 1/298.257223563;
    double PI = 3.14159265358979323846;
	double ddy = distance_lonlat(0, -dy, 0, dy, a, f) / 2 ;
	double yw = 1 / (2 * ddy);
	double xw;
	double ddx;

	int nrow = d.size();
	int ncol = d[0].size();
	std::vector<std::vector<double>> v (nrow, std::vector<double>(ncol, NAN));

	double y, zx, zy;

	for (int row = 1; row < (nrow-2); row++) {
		y = ymax - row * dy;
		ddx = distance_lonlat(-dx, y, dx, y, a, f) / 2 ;
		xw = 1 / (-2 * ddx);

		for (int col = 1; col < (ncol-2); col++) {
			zx = -xw * d[row][col-1] + d[row][col+1] * xw;
			zy = -yw * d[row-1][col] + d[row+1][col] * yw;
			v[row][col] = sqrt( pow(zy, 2) + pow(zx, 2) ) ;
		}
	}

	if (unit == 0) {  // degrees
		double adj = 180 / PI;
		for (int row = 1; row < (nrow-2); row++) {
			for (int col = 1; col < (ncol-2); col++) {
				v[row][col] = atan(v[row][col]) * adj;
			}
		}
	} else if (unit == 1) { // radians
		for (int row = 1; row < (nrow-2); row++) {
			for (int col = 1; col < (ncol-2); col++) {
				v[row][col] = atan(v[row][col]);
			}
		}
	}  // tangent

	return(v);
}


std::vector<std::vector<double> > slope4plane(std::vector<std::vector<double> > d, double dx, double dy, int unit) {

	double xw = 1 / (-2 * dx);
	double yw = 1 / (2 * dy);
    double PI = 3.14159265358979323846;

	int nrow = d.size();
	int ncol = d[0].size();
	std::vector<std::vector<double>> v (nrow, std::vector<double>(ncol, NAN));

	double zx, zy;
	for (int row = 1; row < (nrow-2); row++) {
		for (int col = 1; col < (ncol-2); col++) {
			zx = -xw * d[row][col-1] + d[row][col+1] * xw;
			zy = -yw * d[row-1][col] + d[row+1][col] * yw;
			v[row][col] = sqrt( pow(zy, 2) + pow(zx, 2) ) ;
		}
	}

	if (unit == 0) {
		double adj = 180 / PI;
		for (int row = 1; row < (nrow-2); row++) {
			for (int col = 1; col < (ncol-2); col++) {
				v[row][col] = atan(v[row][col]) * adj;
			}
		}
	} else if (unit == 1) {
		for (int row = 1; row < (nrow-2); row++) {
			for (int col = 1; col < (ncol-2); col++) {
				v[row][col] = atan(v[row][col]);
			}
		}
	}

	return(v);
}

