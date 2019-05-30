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


SpatRaster SpatRaster::crop(SpatExtent e, std::string snap, SpatOptions &opt) {

	SpatRaster out = geometry();

	e.intersect(out.getExtent());

/*	if ( !e.valid() ) {
		return NULL;
		stop("extents do not overlap")
	} */

	out.setExtent(e, true, snap);

	if (!source[0].hasValues ) {
		return(out);
	}

	double xr = xres();
	double yr = yres();

	unsigned col1 = colFromX(out.extent.xmin + 0.5 * xr);
	unsigned col2 = colFromX(out.extent.xmax - 0.5 * xr);
	unsigned row1 = rowFromY(out.extent.ymax - 0.5 * yr);
	unsigned row2 = rowFromY(out.extent.ymin + 0.5 * yr);
	if ((row1==0) && (row2==nrow()-1) && (col1==0) && (col2==ncol()-1)) {
		return(out);
	}

	unsigned ncols = out.ncol();

 	if (!out.writeStart(opt)) { return out; }
	readStart();
	std::vector<double> v;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(row1+out.bs.row[i], out.bs.nrows[i], col1, ncols);
		if (!out.writeValues(v, out.bs.row[i])) return out;
	}
	out.writeStop();
	readStop();

	return(out);
}
