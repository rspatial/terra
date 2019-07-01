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

SpatRaster SpatRasterCollection::merge(SpatOptions &opt) {

	SpatRaster out;
	unsigned n = size();

	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	if (n == 1) {
		out = x[0].deepCopy();
		return(out);
	}

	out = x[0].geometry();
	SpatExtent e = x[0].getExtent();
	for (size_t i=1; i<n; i++) {
		// for now, must have same nlyr; but should be easy to recycle.
		if (!x[0].compare_geom(x[i], true, true, false, false, false, true)) {
			return(out);
		}
		e.unite(x[i].getExtent());
	}
	out.setExtent(e, true);
 //   out.setResolution(xres(), yres());
 	if (!out.writeStart(opt)) { return out; }
	if (out.source[0].driver == "gdal") {
		out.fillValuesGDAL(NAN);
	}

	for (size_t i=0; i<n; i++) {
		SpatRaster r = x[i];
		if (!r.hasValues()) continue;
		BlockSize bs = r.getBlockSize(4);
		r.readStart();
		for (size_t j=0; j<bs.n; j++) {
            std::vector<double> v = r.readValues(bs.row[j], bs.nrows[j], 0, r.ncol());
            unsigned row1 = out.rowFromY(r.yFromRow(bs.row[j]));
            unsigned row2 = out.rowFromY(r.yFromRow(bs.row[j]+bs.nrows[j]-1));
            unsigned col1 = out.colFromX(r.xFromCol(0));
            unsigned col2 = out.colFromX(r.xFromCol(r.ncol()-1));
            if (!out.writeValues(v, row1, row2-row1+1, col1, col2-col1+1)) return out;
		}
		r.readStop();
	}

	out.writeStop();
	return(out);
}

