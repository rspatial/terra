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

/*
#include "spatRaster.h"

bool SpatRaster::isLonLat() {
  return (crs.find("longlat") != std::string::npos);
}

bool SpatRaster::couldBeLonLat() {
	if (crs != "") return isLonLat();
	if ((extent.xmin > -181) & (extent.xmax < 181) &
		(extent.ymin > -91 ) & (extent.ymax < 91 )) {
		return true;
	} 
	return false;
}

bool equal(double a, double b, double epsilon) {
    return fabs(a - b) < epsilon;
}

bool SpatRaster::isGlobalLonLat() {
	double tolerance = 0.1 * xres();
	if (equal(extent.xmin, -180, tolerance) & 
		equal(extent.xmax, 180, tolerance)) {
		if (couldBeLonLat()) {
			return true;
		}
	}
	return false;
}


*/