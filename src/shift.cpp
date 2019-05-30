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

SpatRaster SpatRaster::shift(double x, double y, SpatOptions &opt) {

	SpatRaster out = deepCopy();
	out.extent.xmin = out.extent.xmin + x;
	out.extent.xmax = out.extent.xmax + x;
	out.extent.ymin = out.extent.ymin + y;
	out.extent.ymax = out.extent.ymax + y;
	return out;
	
}




