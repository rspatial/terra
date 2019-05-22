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


SpatRaster SpatRaster::merge(SpatRaster x, SpatOptions &opt) {

	SpatRaster out;
	
	if (!compare_geom(x, true, true, false, false, false, true)) {	
		return(out);
	}
	
	SpatExtent e = extent;
	e.unite(x.extent);
	out.setExtent(e, true);
	
	return(out);
}
