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

/* Compute surface area, using method from:
Jeff S. Jenness, 2004. Calculating Landscape Surface Area from Digital Elevation Models. Wildlife Society Bulletin 32(3):829-839. http://www.jstor.org/stable/3784807

With edge adjustments. 

Adapted by Robert Hijmans from code in R package "sp" by Barry Rowlingson 2010 <b.rowlingson@lancaster.ac.uk>
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


void sarea(std::vector<double> &heights, const size_t &nrow, const long &ncol, const double &w, const double &h, 
	std::vector<double> &sa) {

// given an nx by ny matrix of heights with single-cell edge border, compute the surface area.

// point values
	double z1, z2, z3;
// side lengths
	double l1, l2, l3; 
// diagonal length
	double s2 = sqrt((w*w)+(h*h));

// offsets to neighbours
	std::vector<int> dyv = {-1, -1, -1, 0, 1, 1, 1, 0, -1};
	std::vector<int> dxv = {-1, 0, 1, 1, 1, 0, -1, -1, -1};

// triangle side lengths
// first the radial sides
	double side[] = {s2, h, s2, w, s2, h, s2, w, s2};
// outer edges
	double l3v[] = {w, w, h, h, w, w, h, h};

	size_t outsize = heights.size() - 2*ncol;
	sa = std::vector<double>(outsize, NAN);
	size_t cell = 0;

	for (size_t i=1; i<(nrow-1); i++){
		for (long j=0; j<ncol; j++) {
			z1 = height(heights, ncol, i, j);
			if (!std::isnan(z1)) {
				double cellArea = 0;
				for (size_t tri=0; tri<8; tri++){
					z2 = height(heights, ncol, i+dxv[tri], j+dyv[tri]);
					// replace missing adjacent values with the current cell value
					if (std::isnan(z2)) z2=z1;
					z3 = height(heights, ncol, i+dxv[tri+1], j+dyv[tri+1]);
					if (std::isnan(z3)) z3=z1;
					l1 = 0.5 * sqrt(side[tri] * side[tri] + (z1-z2) * (z1-z2));
					l2 = 0.5 * sqrt(side[tri+1] * side[tri+1] + (z1-z3) * (z1-z3));
					l3 = 0.5 * sqrt(l3v[tri] * l3v[tri] + (z2-z3) * (z2-z3));
					cellArea += triarea(l1, l2, l3);
				}
				sa[cell] = cellArea;
			}
			cell++;
		}	
	}
}



SpatRaster SpatRaster::surfaceArea(SpatOptions &opt) {

	SpatRaster out = geometry(1, false);
	if (is_lonlat()) {
		out.setError("not yet implemented for lonlat data");
		return out;
	}
	
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
	std::vector<double> wh = resolution();
	for (size_t i = 0; i < cbs.n; i++) {
		std::vector<double> v;
		readBlock(v, cbs, i);
		if (i==0) {
			v.insert(v.begin(), v.begin(), v.begin()+nc);
		}
		if (i == (cbs.n - 1)) {
			v.insert(v.end(), v.end()-nc, v.end());
		}
		std::vector<double> sa;
		sarea(v, out.bs.nrows[i]+2, nc, wh[0], wh[1], sa);
		if (!out.writeBlock(sa, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}

