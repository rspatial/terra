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
#include <math.h>
#include "GeographicLib_geodesic.h"
#include "spatRaster.h"
#include "distance.h"

double area_polygon_lonlat(std::vector<double> lon, std::vector<double> lat, double a, double f) {
	struct geod_geodesic g;
	struct geod_polygon p;
	geod_init(&g, a, f);
	geod_polygon_init(&p, 0);
	size_t n = lat.size();
	for (size_t i=0; i < n; i++) {
		geod_polygon_addpoint(&g, &p, lat[i], lon[i]);
	}
	double area, P;
	geod_polygon_compute(&g, &p, 0, 1, &area, &P);
	return(area < 0 ? -area : area);
}



double area_polygon_plane(std::vector<double> x, std::vector<double> y) {
// based on http://paulbourke.net/geometry/polygonmesh/source1.c
	size_t n = x.size();
	double area = x[n-1] * y[0];
	area -= y[n-1] * x[0];
	for (size_t i=0; i < (n-1); i++) {
		area += x[i] * y[i+1];
		area -= x[i+1] * y[i];
	}
	area /= 2;
	return(area < 0 ? -area : area);
}



double SpatGeom::area_lonlat(double a, double f) {
	double area = 0;
	if (gtype != polygons) return area;
	for (size_t i=0; i<parts.size(); i++) {
		area += area_polygon_lonlat(parts[i].x, parts[i].y, a, f);
		for (size_t j=0; j<parts[i].holes.size(); j++) {
			area -= area_polygon_lonlat(parts[i].holes[j].x, parts[i].holes[j].y, a, f);
		}
	}
	return area;
}


double SpatGeom::area_plane() {
	double area = 0;
	if (gtype != polygons) return area;
	for (size_t i=0; i<parts.size(); i++) {
		area += area_polygon_plane(parts[i].x, parts[i].y);
		for (size_t j=0; j<parts[i].holes.size(); j++) {
			area -= area_polygon_plane(parts[i].holes[j].x, parts[i].holes[j].y);
		}
	}
	return area;
}


std::vector<double> SpatVector::area() {

	size_t s = size();
	std::vector<double> ar;
	ar.reserve(s);
	if (could_be_lonlat()) {
		double a = 6378137;
		double f = 1 / 298.257223563;
		for (size_t i=0; i<s; i++) {
			ar.push_back(lyr.geoms[i].area_lonlat(a, f));
		}
	} else {
		for (size_t i=0; i<s; i++) {
			ar.push_back(lyr.geoms[i].area_plane());
		}
	}
	return ar;
}




double length_line_lonlat(std::vector<double> lon, std::vector<double> lat, double a, double f) {
	struct geod_geodesic g;
	geod_init(&g, a, f);
	size_t n = lat.size();
	double length = 0;
	for (size_t i=1; i < n; i++) {
		length += distance_lonlat(lon[i-1], lat[i-1], lon[i], lat[i], a, f);
	}
	return(length);
}



double length_line_plane(std::vector<double> x, std::vector<double> y) {
	size_t n = x.size();
	double length = 0;
	for (size_t i=1; i<n; i++) {
		length += sqrt(pow(x[i-1] - x[i], 2) + pow(y[i-1] - y[i], 2));
	}
	return(length);
}



double SpatGeom::length_lonlat(double a, double f) {
	double length = 0;
	if (gtype == points) return length;
	for (size_t i=0; i<parts.size(); i++) {
		length += length_line_lonlat(parts[i].x, parts[i].y, a, f);
		for (size_t j=0; j<parts[i].holes.size(); j++) {
			length += length_line_lonlat(parts[i].holes[j].x, parts[i].holes[j].y, a, f);
		}
	}
	return length;
}


double SpatGeom::length_plane() {
	double length = 0;
	if (gtype == points) return length;
	for (size_t i=0; i<parts.size(); i++) {
		length += length_line_plane(parts[i].x, parts[i].y);
		for (size_t j=0; j<parts[i].holes.size(); j++) {
			length += length_line_plane(parts[i].holes[j].x, parts[i].holes[j].y);
		}
	}
	return length;
}


std::vector<double> SpatVector::length() {

	size_t s = size();
	std::vector<double> r;
	r.reserve(s);

	if (could_be_lonlat()) {
		double a = 6378137;
		double f = 1 / 298.257223563;
		for (size_t i=0; i<s; i++) {
			r.push_back(lyr.geoms[i].length_lonlat(a, f));
		}
	} else {
		for (size_t i=0; i<s; i++) {
			r.push_back(lyr.geoms[i].length_plane());
		}
	}
	return r;
}

SpatRaster SpatRaster::area(SpatOptions &opt) {

	SpatRaster out = geometry(1);
  	if (!out.writeStart(opt)) { return out; }
	if (could_be_lonlat()) {
		SpatExtent e = {extent.xmin, extent.xmin+xres(), extent.ymin, extent.ymax};
		SpatOptions optint(opt);
		SpatRaster onecol = out.crop(e, "near", optint);
		SpatVector p = onecol.as_polygons(false, false);
		std::vector<double> a = p.area();
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			for (size_t j=0; j<out.bs.nrows[i]; j++) {
				size_t r = out.bs.row[i] + j;
				v.insert(v.end(), ncol(), a[r]);
			}
			if (!out.writeValues(v, out.bs.row[i])) return out;
		}
		out.writeStop();
	} else {
		double a = xres() * yres();
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v(out.bs.nrows[i]*ncol(), a);
			if (!out.writeValues(v, out.bs.row[i])) return out;
		}
		out.writeStop();
	}
	return(out);
}

