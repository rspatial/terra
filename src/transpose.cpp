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

#include <vector>
#include "spatRaster.h"

SpatRaster SpatRaster::transpose(SpatOptions &opt) {

	SpatRaster out = geometry();
	SpatExtent eold = getExtent(); 
	SpatExtent enew = getExtent(); 
	enew.xmin = eold.ymin;
	enew.xmax = eold.ymax;
	enew.ymin = eold.xmin;
	enew.ymax = eold.xmax;
	out.setExtent(enew, false, "");
	out.source[0].ncol = nrow();
	out.source[0].nrow = ncol();
	if (!hasValues()) return out;
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i=0; i < out.bs.n; i++) {
		std::vector<double> v = readValues(0, nrow(), out.bs.row[i], out.bs.nrows[i]);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}

