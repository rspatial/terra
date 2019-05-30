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
#include "math_utils.h"


bool SpatRaster::compare_geom(SpatRaster x, bool lyrs, bool crs, bool warncrs, bool ext, bool rowcol, bool res) {
	
	
	if (ext) {
		bool e1 = is_equal(x.extent.xmax, extent.xmax, 1);
		bool e2 = is_equal(x.extent.xmin, extent.xmin, 1);
		bool e3 = is_equal(x.extent.ymax, extent.ymax, 1);
		bool e4 = is_equal(x.extent.ymin, extent.ymin, 1);
		if (!(e1 && e2 && e3 && e4)) {
			setError("extents do not match");
			return false;
		}
	}
	if (rowcol) {
		if (! ((nrow() == x.nrow()) && (ncol() == x.ncol())) ) {
			setError("number of rows and/or columns do not match");
			return false;
		}
	}
	if (res) {
		if (! ((x.xres() == xres()) & (x.yres() == yres()))) {
			setError("resolution does not match");
			return false;
		}
	}

	if (lyrs) {
		if (!(nlyr() == x.nlyr())) {
			setError("number of layers does not match");
			return false;
		}
	}

	if (crs) {
		if (!(getCRS() == x.getCRS())) {
			if (warncrs) {
				addWarning("CRS do not match");
			} else {
				setError("CRS do not match");
				return false;
			}
		}
	}
	return true;
}

