// Copyright (c) 2018-2020  Robert J. Hijmans
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


bool extent_operator(std::string oper) {
	std::vector<std::string> f {"==", "!=", ">", "<", ">=", "<="};
	return (std::find(f.begin(), f.end(), oper) != f.end());
}

bool SpatExtent::compare(SpatExtent e, std::string oper, double tolerance) {

	if (!extent_operator(oper)) {
		return false;  // not very useful
	}

	//double xr = (xmax - xmin) / tolerance;
	//double yr = (ymax - ymin) / tolerance;

	bool e1 = fabs(xmax - e.xmax) <= tolerance;
	bool e2 = fabs(xmin - e.xmin) <= tolerance;
	bool e3 = fabs(ymax - e.ymax) <= tolerance;
	bool e4 = fabs(ymin - e.ymin) <= tolerance;
	bool equal = (e1 && e2 && e3 && e4);
	if (oper == "==") {
		return equal;
	} else if (oper == "!=") {
		return (!equal);
	}
	if (oper == "<" || oper == "<=") {
		bool c1 = xmax < e.xmax;
		bool c2 = xmin > e.xmin;
		bool c3 = ymax < e.ymax;
		bool c4 = ymin > e.ymin;
		bool smaller = (c1 && c2 && c3 && c4);
		if (oper == "<") {
			return smaller;
		} else {
			return (equal || smaller);
		}
	} 
	if (oper == ">" || oper == ">=") {
		bool c1 = xmax > e.xmax;
		bool c2 = xmin < e.xmin;
		bool c3 = ymax > e.ymax;
		bool c4 = ymin < e.ymin;
		bool larger = (c1 && c2 && c3 && c4);
		if (oper == ">") {
			return larger;
		} else {
			return (equal || larger);
		}
	}
	return false;
}

SpatExtent SpatExtent::round(int n) {
	double xn = roundn(xmin, n);
	double xx = roundn(xmax, n);
	double yn = roundn(ymin, n);
	double yx = roundn(ymax, n);
	SpatExtent e(xn, xx, yn, yx);
	return e;
}


SpatExtent SpatExtent::floor() {
	double xn = std::floor(xmin);
	double xx = std::ceil(xmax);
	double yn = std::floor(ymin);
	double yx = std::ceil(ymax);
	SpatExtent e(xn, xx, yn, yx);
	return e;
}

SpatExtent SpatExtent::ceil() {
	double xn = std::ceil(xmin);
	double xx = std::floor(xmax);
	double yn = std::ceil(ymin);
	double yx = std::floor(ymax);
	SpatExtent e(xn, xx, yn, yx);
	return e;
}

SpatExtent SpatRaster::getExtent() { 
	if (source.size() > 0) {
		return source[0].extent;
	} else {
		SpatExtent e;
		return e;
	}
}



void SpatRaster::setExtent(SpatExtent e) { 
	for (size_t i=0; i<nsrc(); i++) {
		source[i].extent = e;
	}
}


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
		source[0].extent = ext;
	}

	for (size_t i=0; i<nsrc(); i++) {
		source[i].extent = ext;
		//source[i].nrow = source[0].nrow;
		//source[i].ncol = source[0].ncol;
	}
}


SpatExtent SpatExtent::align(double d, std::string snap) {
    std::vector<double> e = asVector();
	if (d == 0) {
		SpatExtent out = *this;
		return(out);
	}
	d = d < 0 ? -d : d;


	for (size_t i=0; i<4; i++) {
		double x = d * trunc(e[i] / d);
		if ((i == 0) | (i == 2)) {
			if (x > e[i]) {
				x -= d;
			} 
		} else {
			if (x < e[i]) {
				x += d;
			}
		}
		e[i] = x;
	}
	SpatExtent out(e[0], e[1], e[2], e[3]);
	return(out)	;
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
		xmn = std::floor((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = std::ceil((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = std::floor((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = std::ceil((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
	} else { //if (snap == "in") {
		xmn = std::ceil((e.xmin-orig[0]) / res[0]) * res[0] + orig[0];
		xmx = std::floor((e.xmax-orig[0]) / res[0]) * res[0] + orig[0];
		ymn = std::ceil((e.ymin-orig[1]) / res[1]) * res[1] + orig[1];
		ymx = std::floor((e.ymax-orig[1]) / res[1]) * res[1] + orig[1];
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
	SpatExtent extent = getExtent();
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



bool SpatRaster::compare_geom(SpatRaster x, bool lyrs, bool crs, bool warncrs, bool ext, bool rowcol, bool res) {


	if (ext) {
		SpatExtent extent = getExtent();
		if (extent.compare(x.getExtent(), "!=", 100)) {
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
		if (! ((is_equal_relative(x.xres(), xres(), 0.0001)) & (is_equal_relative(x.yres(), yres(), 0.0001)))) {
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
		if (!source[0].srs.is_equal(x.source[0].srs)) {
			if (warncrs) {
				addWarning("SRS do not match");
			} else {
				setError("SRS do not match");
				return false;
			}
		}
	}
	return true;
}


