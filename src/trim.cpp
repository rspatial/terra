// Copyright (c) 2018-2019  Robert J. Hijmans
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

SpatRaster SpatRaster::trim(unsigned padding, SpatOptions &opt) {

	long nrl = nrow() * nlyr();
	long ncl = ncol() * nlyr();

	std::vector<double> v;
	unsigned r;
	for (r=0; r<nrow(); r++) {
		v = readValues(r, 1, 0, ncol());
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
			break;
		}
	}

	if (r == nrow()) { //stop('only NA values found')
	}
	unsigned firstrow = std::min(std::max(r - padding, unsigned(0)), nrow());

	for (r=nrow()-1; r>firstrow; r--) {
		v = readValues(r, 1, 0, ncol());
		if (std::count_if(v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
			break;
		}
	}

	unsigned lastrow = std::max(std::min(r+padding, nrow()), unsigned(0));

	unsigned tmp;
	if (lastrow < firstrow) {
		tmp = firstrow;
		firstrow = lastrow;
		lastrow = tmp;
	}
	unsigned c;
	for (c=0; c<ncol(); c++) {
		v = readValues(0, nrow(), c, 1);
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
			break;
		}
	}
	unsigned firstcol = std::min(std::max(c-padding, unsigned(0)), ncol());


	for (c=ncol()-1; c>firstcol; c--) {
		v = readValues(0, nrow(), c, 1);
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
			break;
		}
	}
	unsigned lastcol = std::max(std::min(c+padding, ncol()), unsigned(0));

	if (lastcol < firstcol) {
		tmp = firstcol;
		firstcol = lastcol;
		lastcol = tmp;
	}

	std::vector<double> res = resolution();
	double xr = res[0];
	double yr = res[1];
	SpatExtent e = SpatExtent(xFromCol(firstcol)-0.5*xr, xFromCol(lastcol)+0.5*xr, yFromRow(lastrow)-0.5*yr, yFromRow(firstrow)+0.5*yr);

	return( crop(e, "near", opt) ) ;
}

