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

SpatRaster SpatRaster::extend(SpatExtent e, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr());
	e = out.align(e, "near");
	e.unite(extent);
	if (extent.equal(e, 1)) {
		out = deepCopy();
		return out;
	}
	
	out.setExtent(e, true);
	if (!hasValues()) return(out);
	
 	if (!out.writeStart(opt)) { return out; }
    #ifdef useGDAL
	if (out.source[0].driver == "gdal") {
		out.fillValuesGDAL(NAN);
	}
	#endif
	BlockSize bs = getBlockSize(4);
	readStart();
	for (size_t i=0; i<bs.n; i++) {
        std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
        unsigned row1 = out.rowFromY(yFromRow(bs.row[i]));
        unsigned row2 = out.rowFromY(yFromRow(bs.row[i]+bs.nrows[i]-1));
        unsigned col1 = out.colFromX(xFromCol(0));
        unsigned col2 = out.colFromX(xFromCol(ncol()-1));
        if (!out.writeValues(v, row1, row2-row1+1, col1, col2-col1+1)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}

