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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURP0OSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include "spatRaster.h"


SpatRaster SpatRaster::collapse(SpatRaster x, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!out.compare_geom(x, true, true)) {
		out.setError("dimensions and/or extent do not match");
		return(out);
	}	
	if (!hasValues()) return(out);
	if (!x.hasValues()) {
		out.setError("index raster has no values");
		return out;
	}
	
	if (x.nlyr() > 1) {
		SpatOptions ops;
		std::vector<unsigned> lyr = {0};
		x = x.subset(lyr, ops);
	}

	int nl = nlyr();
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	x.readStart();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		std::vector<double> idx = x.readBlock(out.bs, i);
		size_t is = idx.size();
		std::vector<double> vv(is, NAN);
		size_t ncell = out.bs.nrows[i] * ncol();
		for (size_t j=0; j<is; j++) {
			int index = idx[j] - 1;
			if ((index >= 0) && (index < nl)) {
				vv[j] = v[j + index * ncell];
			}				
		}
		if (!out.writeValues(vv, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	readStop();
	x.readStop();
	out.writeStop();
	return(out);
}

