// Copyright (c) 2018  Robert J. Hijmans
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
#include "recycle.h"


SpatRaster SpatRaster::rotate(SpatOptions &opt) {

	SpatOptions optint(opt);
	double xhalf = (extent.xmax - extent.xmin) / 2;
	SpatExtent ext1 = SpatExtent(extent.xmin, extent.xmin + xhalf, extent.ymin, extent.ymax);
	SpatExtent ext2 = SpatExtent(extent.xmin + xhalf, extent.xmax, extent.ymin, extent.ymax);
	SpatRaster r1 = crop(ext1, "", optint);
	SpatRaster r2 = crop(ext2, "", optint);
	// change extents
	SpatRaster out = r1.merge(r2, opt);
	
	return(out);
}

