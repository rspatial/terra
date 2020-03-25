// Copyright (c) 2018-2020  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

/* Robert Hijmans, October 2011 */


#include "spatRaster.h"
#include "math_utils.h"
#include "distance.h"
#include <cmath>
#include <string>
//#include <vector>
#include <algorithm>


#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif



double dmod(double x, double n) {
	return(x - n * std::floor(x/n));
}


double distPlane(double x1, double y1, double x2, double y2) {
	return( sqrt(pow((x2-x1),2) + pow((y2-y1), 2)) );
}

double distHav(double lon1, double lat1, double lon2, double lat2) {
	double r = 6378137;
	double dLat, dLon, a;
	lon1 = toRad(lon1);
	lon2 = toRad(lon2);
	lat1 = toRad(lat1);
	lat2 = toRad(lat2);

	dLat = lat2-lat1;
	dLon = lon2-lon1;
	a = sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2.) * sin(dLon/2.);
	return 2. * atan2(sqrt(a), sqrt(1.-a)) * r;
}



double TRI (std::vector<double> v) {
	double s = 0;
	size_t n = 0;
	if (std::isnan(v[4])) return(NAN);
	for (size_t i=0; i<9; i++) {
		if (! (i==4 || std::isnan(v[i]))) {
			s += fabs(v[i] - v[4]);
			n++;
		}
	}
	if (n > 0) {
		s = s / n;
	} else {
		s = NAN;
	}
	return(s);
}

double TPI (std::vector<double> v) {
	double s = 0;
	size_t n = 0;
	for (size_t i=0; i<9; i++) {
		if (! (i==4 || std::isnan(v[i]))) {
			s += v[i];
			n++;
		}
	}
	if (n > 0) {
		s = s / n;
		s = v[4] - s;
	} else {
		s = NAN;
	}
	return(s);
}

double roughness (std::vector<double> v) {
	double vmin, vmax;
	minmax(v.begin(), v.end(), vmin, vmax);
	return(vmax - vmin);
}



std::vector<double> do_slope(std::vector<double> d, unsigned ngb, unsigned nrow, unsigned ncol, double dx, double dy, bool geo, std::vector<double> gy) {
	size_t n = nrow * ncol;
	
	std::vector<double> ddx;
	if (geo) {
		ddx.resize(nrow);	
		for (size_t i=0; i<nrow; i++) {
			ddx[i] = distHav(-dx, gy[i], dx, gy[i]) / 2 ;
		}
	} 

	double zy, zx; 
	std::vector<double> val(n);

	
	if (ngb == 4) {
		if (geo) {
			int q;
			double xwi[2] = {-1,1};
			double xw[2] = {0,0};
			double yw[2] = {-1,1};

			for (size_t i=0; i<2; i++) {
				yw[i] = yw[i] / (2 * dy);
			}			
			for (size_t i = ncol; i < (ncol * (nrow-1)-1); i++) {
				if (i % ncol == 0) {
					q = i / ncol;
					for (size_t k=0; k<2; k++) {
						xw[k] = xwi[k] / (-2 * ddx[q]);
					}
				}
				zx = d[i-1] * xw[0] + d[i+1] * xw[1];
				zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
				val[i] = sqrt( pow(zy, 2) + pow(zx, 2) ) ;
			}
		} else {
			
			double xw[2] = {-1,1};
			double yw[2] = {-1,1};
			for (size_t i=0; i<2; i++) {
				xw[i] = xw[i] / (-2 * dx);
				yw[i] = yw[i] / (2 * dy);
			}
			for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
				zx = d[i-1] * xw[0] + d[i+1] * xw[1];
				zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
				val[i] = sqrt( pow(zy, 2) + pow(zx, 2)  );
			}
		}
	} else {
		
		if (geo) {
			int q;
			double xwi[6] = {-1,-2,-1,1,2,1};
			double xw[6] = {0,0,0,0,0,0};
			double yw[6] = {-1,1,-2,2,-1,1};
			
			for (size_t i=0; i<6; i++) {
				yw[i] = yw[i] / (8 * dy);
			}
						
			for (size_t i = ncol; i < (ncol * (nrow-1)-1); i++) {
				if (i % ncol == 0) {
					q = i / ncol;
					for (size_t k=0; k<6; k++) {
						xw[k] = xwi[k] / (8 * ddx[q]);
					}
				}
				zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
						+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
				zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2] 
						+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
				val[i] = sqrt( pow(zy, 2) + pow(zx, 2)  );
								
			}
			
		} else {
		
			double xw[6] = {-1,-2,-1,1,2,1};
			double yw[6] = {-1,1,-2,2,-1,1};
			for (size_t i=0; i<6; i++) {
				xw[i] = xw[i] / (-8 * dx);
				yw[i] = yw[i] / (8 * dy);
			}
			for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
				zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
						+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
				zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2] 
						+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
				val[i] = sqrt( pow(zy, 2) + pow(zx, 2) );
			}
		}
	} 	
				
	return val;
}


void to_degrees(std::vector<double>& x) { 
	double adj = 180 / M_PI;
	for (double& d : x) d = atan(d)*adj;
}

void to_radians(std::vector<double>& x) { 
	for (double& d : x) d = atan(d);
}

SpatRaster SpatRaster::slope(unsigned neighbors, bool degrees, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (nlyr() > 1) {
		out.setError("provide a single layer object");
		return out;
	}
	std::vector<double> d = getValues();
	std::vector<double> y;
	bool lonlat = could_be_lonlat();
	double yr = yres();
	if (lonlat) {
		std::vector<unsigned> rows(nrow());
		std::iota(rows.begin(), rows.end(), 0);
		y = yFromRow(rows);
		yr = distHav(0, 0, 0, yr);
	}
	std::vector<double> val = do_slope(d, neighbors, nrow(), ncol(), xres(), yr, lonlat, y);
	if (degrees) { 
		to_degrees(val);
	} else {
		to_radians(val);
	}
	out.setValues(val);
	return out; 
}


/*
std::vector<std::vector<double> > terrain_indices(std::vector<std::vector<double> > &m, std::vector<std::vector<bool> > f, std::string option) {

	int nrows = m.size();
	int ncols = m[0].size();
	std::vector<std::vector<double>> v (nrows, std::vector<double>(ncols, NAN));


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
	std::vector<double> va(9);
	for (int row = 0; row < (nrows-1); row++) {
		for (int col = 0; col < (ncols-1); col++) {
			n = 0;
			for(int i = -2; i < 3; i++) {
				r = i + row;
				if (r < 0 || r > (nrows-1)) {
					continue;
				}
				for(int j = -2; j < 3; j++) {
					c = j + col;
					if (c < 0 || c > (ncols-1)) {
						continue;
					}
					// center not included for TPI and TRI
					if (r==0 && c==0 && opt != 2) {
						continue;
					}
					if (! std::isnan(m[r][c]) ) {
						va[n] = m[r][c];
						n++;
					}
				}
			}
			if (opt==0) {
				v[row][col] = TPI(va);
			} else if (opt==1) {
				v[row][col] = TRI(va);
			} else { // rough
				v[row][col] = roughness(va);
			}
		}
	}
	return(v);
}




std::vector<std::vector<double> > slope4lonlat(std::vector<std::vector<double> > d, double dx, double dy, double ymax, int unit) {

	double a = 6378137;
	double f = 1/298.257223563;
	double ddy = distance_lonlat(0, -dy, 0, dy, a, f) / 2 ;
	double yw = 1 / (2 * ddy);
	double xw;
	double ddx;

	int nrows = d.size();
	int ncols = d[0].size();
	std::vector<std::vector<double>> v (nrows, std::vector<double>(ncols, NAN));

	double y, zx, zy;

	for (int row = 1; row < (nrows-2); row++) {
		y = ymax - row * dy;
		ddx = distance_lonlat(-dx, y, dx, y, a, f) / 2 ;
		xw = 1 / (-2 * ddx);

		for (int col = 1; col < (ncols-2); col++) {
			zx = -xw * d[row][col-1] + d[row][col+1] * xw;
			zy = -yw * d[row-1][col] + d[row+1][col] * yw;
			v[row][col] = sqrt( pow(zy, 2) + pow(zx, 2) ) ;
		}
	}

	if (unit == 0) {  // degrees
		double adj = 180 / M_PI;
		for (int row = 1; row < (nrows-2); row++) {
			for (int col = 1; col < (ncols-2); col++) {
				v[row][col] = atan(v[row][col]) * adj;
			}
		}
	} else if (unit == 1) { // radians
		for (int row = 1; row < (nrows-2); row++) {
			for (int col = 1; col < (ncols-2); col++) {
				v[row][col] = atan(v[row][col]);
			}
		}
	}  // tangent

	return(v);
}


std::vector<std::vector<double> > slope4plane(std::vector<std::vector<double> > d, double dx, double dy, int unit) {

	double xw = 1 / (-2 * dx);
	double yw = 1 / (2 * dy);

	int nrows = d.size();
	int ncols = d[0].size();
	std::vector<std::vector<double>> v (nrows, std::vector<double>(ncols, NAN));

	double zx, zy;
	for (int row = 1; row < (nrows-2); row++) {
		for (int col = 1; col < (ncols-2); col++) {
			zx = -xw * d[row][col-1] + d[row][col+1] * xw;
			zy = -yw * d[row-1][col] + d[row+1][col] * yw;
			v[row][col] = sqrt( pow(zy, 2) + pow(zx, 2) ) ;
		}
	}

	if (unit == 0) {
		double adj = 180 / M_PI;
		for (int row = 1; row < (nrows-2); row++) {
			for (int col = 1; col < (ncols-2); col++) {
				v[row][col] = atan(v[row][col]) * adj;
			}
		}
	} else if (unit == 1) {
		for (int row = 1; row < (nrows-2); row++) {
			for (int col = 1; col < (ncols-2); col++) {
				v[row][col] = atan(v[row][col]);
			}
		}
	}

	return(v);
}




SpatRaster SpatRaster::terrain(std::string option, std::string unit, SpatOptions &opt) {

	SpatRaster out=geometry();
    std::function<std::vector<double>(std::vector<double>&, double&)> terrainFun;

	if (option == "TPI") {
		terrainFun = TPI;
	} else if (option == "TRI") {
		terrainFun = TRI;
	} else { //(option == 'roughness')
		terrainFun = roughness;
	}

	std::vector<double> va(9, NAN);
	unsigned nc = ncol();
	unsigned nr = nrow();
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		unsigned startrow = std::max(0, bs.row[i]-1);
		unsigned nrows = std::min(nrow()-startrow, bs.nrows[i]+1);
		std::vector<double> vv;
		if (startrow = 0) {
			vv.insert(vv.begin(), std::vector<double>(nc, NAN))
		}
        std::vector<double> v = readValues(startrow, nrows, 0, nc);
		vv.insert(vv.begin(), v);
		if ((startrow+nrows) > nr) {
			vv.insert(vv.begin(), std::vector<double>(nc, NAN))		
		}
		
		std::vector<double> a = readBlock(out.bs, i);
		for (size_t j=0; j<a.size(); j++) {
			
			
			
		}
	}
	
	for (int row = 0; row < (nrows-1); row++) {
		for (int col = 0; col < (ncols-1); col++) {
			n = 0;
			for(int i = -2; i < 3; i++) {
				r = i + row;
				if (r < 0 || r > (nrows-1)) {
					continue;
				}
				for(int j = -2; j < 3; j++) {
					c = j + col;
					if (c < 0 || c > (ncols-1)) {
						continue;
					}
					// center not included for TPI and TRI
					if (r==0 && c==0 && opt != 2) {
						continue;
					}
					if (! std::isnan(m[r][c]) ) {
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
	return(

	return out;
}

*/

