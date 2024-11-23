// Copyright (c) 2018-2023  Robert J. Hijmans
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


#include "spatRaster.h"
#include "distance.h"

/* Compute surface area, using method from:
Jeff S. Jenness, 2004. Calculating Landscape Surface Area from Digital Elevation Models. Wildlife Society Bulletin 32(3):829-839. http://www.jstor.org/stable/3784807

With edge adjustments. 

From C code in R package "sp" by Barry Rowlingson 2010 <b.rowlingson@lancaster.ac.uk>

Adapted for terra by Robert Hijmans (C++, support for lonlat) 
*/

inline double height(const std::vector<double> &heights, const long &ncols, const size_t &row, const long &col) {
	size_t ncr = ncols*row;
	if (col < 0) {
		return heights[ncr]; // + 0
	} else if (col == ncols) {
		return heights[ncr + (ncols-1)];
	} else {
		return heights[ncr + col];
	}
}

inline double triarea(const double &a, const double &b, const double &c) {
	// triangle area given side lengths
	double s = (a + b + c) / 2.0;
	return sqrt(s * (s-a) * (s-b) * (s-c));
}


void sarea(std::vector<double> &heights, const size_t &nrow, const long &ncol, const std::vector<double> &w, const double &h, bool lonlat, std::vector<double> &sa) {

// given an nx by ny matrix of heights with single-cell edge border, compute the surface area.

// offsets to neighbours
	std::vector<int> dyv = {-1, -1, -1, 0, 1, 1, 1, 0, -1};
	std::vector<int> dxv = {-1, 0, 1, 1, 1, 0, -1, -1, -1};

// triangle diagonal length
	double s2 = sqrt((w[0]*w[0]) + (h*h));
// radial side lengths
	std::vector<double> side = {s2, h, s2, w[0], s2, h, s2, w[0], s2};
// outer edges lengths
	std::vector<double> l3v = {w[0], w[0], h, h, w[0], w[0], h, h};

	sa = std::vector<double>(heights.size() - 2*ncol, NAN);

	size_t cell = 0;
	for (size_t i=1; i<(nrow-1); i++){
		if (lonlat) {
			size_t k = i - 1;
			s2 = sqrt((w[k]*w[k]) + (h*h));
			side = {s2, h, s2, w[k], s2, h, s2, w[k], s2};
			l3v = {w[k], w[k], h, h, w[k], w[k], h, h};
		}

		for (long j=0; j<ncol; j++) {
			double z1 = height(heights, ncol, i, j);
			if (!std::isnan(z1)) {
				sa[cell] = 0;
				for (size_t tri=0; tri<8; tri++){
					double z2 = height(heights, ncol, i+dxv[tri], j+dyv[tri]);
					// replace missing adjacent values with the current cell value
					if (std::isnan(z2)) z2=z1;
					double z3 = height(heights, ncol, i+dxv[tri+1], j+dyv[tri+1]);
					if (std::isnan(z3)) z3=z1;
					double l1 = 0.5 * sqrt(side[tri] * side[tri] + (z1-z2) * (z1-z2));
					double l2 = 0.5 * sqrt(side[tri+1] * side[tri+1] + (z1-z3) * (z1-z3));
					double l3 = 0.5 * sqrt(l3v[tri] * l3v[tri] + (z2-z3) * (z2-z3));
					sa[cell] += triarea(l1, l2, l3);
				}
			}
			cell++;
		}	
	}
}


SpatRaster SpatRaster::surfaceArea(SpatOptions &opt) {

	SpatRaster out = geometry(1, false);
	
	if (!hasValues()) {
		out.setError("cannot compute surfaceArea for a raster with no values");
		return out;
	}

	size_t nl = nlyr();
	if (nl != 1) {
		out.setError("can only compute surfaceArea for a single raster layer");
		return out;		
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	BlockSize cbs = out.bs;
	for (size_t i = 0; i < (cbs.n-1); i++) {
		cbs.nrows[i] += 1;
	}	
	for (size_t i = 1; i < cbs.n; i++) {
		cbs.row[i] -= 1;
		cbs.nrows[i] += 1;
	}
	
	size_t nc = ncol();	
	bool lonlat = is_lonlat();
	double xr = xres();	
	std::vector<double> resx = { xr };	
	double resy = yres();
	if (lonlat) {
		resy = distance_lonlat(0, 0, 0, resy);		
	} 
	std::vector<int_64> rows;

	for (size_t i = 0; i < cbs.n; i++) {
		std::vector<double> v;
		readBlock(v, cbs, i);
		if (i==0) {
			v.insert(v.begin(), v.begin(), v.begin()+nc);
		}
		if (i == (cbs.n - 1)) {
			v.insert(v.end(), v.end()-nc, v.end());
		}
		if (lonlat) {
			rows.resize(out.bs.nrows[i]);
			std::iota(rows.begin(), rows.end(), out.bs.row[i]);
			std::vector<double> y = yFromRow(rows);
			resx = distance_lon(xr, y);
		}
		std::vector<double> sa;
		sarea(v, out.bs.nrows[i]+2, nc, resx, resy, lonlat, sa);
		if (!out.writeBlock(sa, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}

