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
#include "string_utils.h"
#include "math_utils.h"

void SpatRaster::setExtent(SpatExtent ext, bool keepRes, std::string snap) {

	if (snap != "") {
		ext = align(ext, snap);
	}
	
	if (keepRes) {
		std::vector<double> res = resolution();
		double xrs = res[0];
		double yrs = res[1];
		unsigned nc = std::max(1.0, round( (ext.xmax - ext.xmin) / xrs ));
		unsigned nr = std::max(1.0, round( (ext.ymax - ext.ymin) / yrs ));
		source[0].ncol = nc;
		source[0].nrow = nr;
		ext.xmax = ext.xmin + nc * xrs;
		ext.ymax = ext.ymin + nr * yrs;
	}
	
	extent = ext;
}


SpatExtent SpatRaster::align(SpatExtent e, std::string snap) {

	snap = is_in_set_default(snap, std::vector<std::string> {"near", "in", "out"}, "near", true);
	std::vector<double> res = resolution();
	std::vector<double> orig = origin(); 
	
	// snap points to pixel boundaries
	double xmn, xmx, ymn, ymx;
	if (snap == "near") {
		xmn = round((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = round((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = round((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = round((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
	} else if (snap == "out") {
		xmn = floor((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = ceil((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = floor((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = ceil((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
	} else { //if (snap == "in") {
		xmn = ceil((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = floor((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = ceil((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = floor((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
	}
	
	if (xmn == xmx) {
		if (xmn < e.xmin) {
			xmx = xmx + res[0];
		} else {
			xmn = xmn - res[0];	
		}
	}
	if (ymn == ymx) {
		if (ymn < e.ymin) {
			ymx = ymx + res[1];
		} else {
			ymn = ymn - res[1];	
		}
	}
	return SpatExtent(xmn, xmx, ymn, ymx);
}


std::vector<double> SpatRaster::origin() {
	std::vector<double> r = resolution();
	double x = extent.xmin - r[0] * (round(extent.xmin / r[0]));
	double y = extent.ymax - r[1] * (round(extent.ymax / r[1]));
	if (is_equal((r[0] + x), abs(x))) {
		x = fabs(x);
	}
	if (is_equal((r[1] + y), abs(y))) {
		y = fabs(y);
	}
	std::vector<double> out {x, y};
	return out;
}

