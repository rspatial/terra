// Copyright (c) 2018-2021  Robert J. Hijmans
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
#include "ggeodesic.h"
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

	double m = srs.to_meter();
	m = std::isnan(m) ? 1 : m;
	// could_be_lonlat() no more
	if (m == 0) {
		double a = 6378137;
		double f = 1 / 298.257223563;
		for (size_t i=0; i<s; i++) {
			ar.push_back(geoms[i].area_lonlat(a, f));
		}
	} else {
		for (size_t i=0; i<s; i++) {
			ar.push_back(geoms[i].area_plane() * m * m);
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
		length += distance_lonlat(lon[i-1], lat[i-1], lon[i], lat[i]);
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

	double m = srs.to_meter();
	m = std::isnan(m) ? 1 : m;

//	if (could_be_lonlat()) {
	if (m == 0) {
		double a = 6378137;
		double f = 1 / 298.257223563;
		for (size_t i=0; i<s; i++) {
			r.push_back(geoms[i].length_lonlat(a, f));
		}
	} else {
		for (size_t i=0; i<s; i++) {
			r.push_back(geoms[i].length_plane() * m);
		}
	}
	return r;
}

SpatRaster SpatRaster::rst_area(bool adjust, SpatOptions &opt) {

	double m = source[0].srs.to_meter();
	m = std::isnan(m) ? 1 : m;

	SpatRaster out = geometry(1);
  	if (!out.writeStart(opt)) { return out; }

	if (m == 0) { //lonlat
		SpatExtent extent = getExtent();
		SpatExtent e = {extent.xmin, extent.xmin+xres(), extent.ymin, extent.ymax};
		SpatOptions optint(opt);
		SpatRaster onecol = out.crop(e, "near", optint);
		SpatOptions popt(opt);
		SpatVector p = onecol.as_polygons(false, false, false, false, popt);
		if (p.hasError()) {
			out.setError(p.getError());
			return out;
		}
		std::vector<double> a = p.area();
		size_t nc = ncol();
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			for (size_t j=0; j<out.bs.nrows[i]; j++) {
				size_t r = out.bs.row[i] + j;
				v.insert(v.end(), nc, a[r]);
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		}

	} else if (adjust) {
		SpatExtent extent = getExtent();
		double dy = yres() / 2;
		SpatOptions popt(opt);
		for (size_t i = 0; i < out.bs.n; i++) {
			double ymax = yFromRow(out.bs.row[i]) + dy;
			double ymin = yFromRow(out.bs.row[i] + out.bs.nrows[i]-1) - dy;
			SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
			SpatRaster onechunk = out.crop(e, "near", popt);
			SpatVector p = onechunk.as_polygons(false, false, false, false, popt);
			//std::vector<double> cells(onechunk.ncell());
			//std::iota (cells.begin(), cells.end(), 0);
			//onechunk.setValues(cells);
			//SpatVector p = onechunk.as_polygons(false, true, false, false, popt);
			p = p.project("EPSG:4326");
			std::vector<double> v = p.area();
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		}
	} else {
		double a = xres() * yres() * m * m;
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v(out.bs.nrows[i]*ncol(), a);
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		}
	}
	out.writeStop();
	return(out);
}


std::vector<double> SpatRaster::sum_area(bool adjust, SpatOptions &opt) {

	std::vector<double> out(nlyr(), 0);

	if (adjust) { //avoid very large polygon objects
		opt.set_memfrac(std::max(0.1, opt.get_memfrac()/2));
	}
	BlockSize bs = getBlockSize(opt);
	if (!readStart()) {
		std::vector<double> err(nlyr(), -1);
		return(err);
	}

	double m = source[0].srs.to_meter();
	m = std::isnan(m) ? 1 : m;

	if (m == 0) {
		SpatRaster x = geometry(1);
		SpatExtent extent = x.getExtent();
		SpatExtent e = {extent.xmin, extent.xmin+xres(), extent.ymin, extent.ymax};
		SpatOptions opt;
		SpatRaster onecol = x.crop(e, "near", opt);
		SpatVector p = onecol.as_polygons(false, false, false, false, opt);
		std::vector<double> ar = p.area();
		size_t nc = ncol();
		if (!hasValues()) {
			out.resize(1);
			for (size_t i=0; i<ar.size(); i++) {
				out[0] += ar[i] * nc;
			}
		} else {
			for (size_t i=0; i<bs.n; i++) {
				std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
				size_t blockoff = bs.nrows[i] * nc;
				for (size_t lyr=0; lyr<nlyr(); lyr++) {
					size_t lyroff = lyr * blockoff;
					for (size_t j=0; j<bs.nrows[i]; j++) {
						size_t row = bs.row[i] + j;
						size_t offset = lyroff + row * nc;
						size_t n = offset + nc;
						for (size_t k=offset; k<n; k++) {
							if (!std::isnan(v[k])) out[lyr] += ar[row];
						}
					}
				}
			}
		}
	} else if (adjust) {
		SpatRaster x = geometry(1);
		SpatExtent extent = x.getExtent();
		double dy = x.yres() / 2;
		SpatOptions popt(opt);
		if (!hasValues()) {
			out.resize(1);
			for (size_t i=0; i<bs.n; i++) {
				double ymax = x.yFromRow(bs.row[i]) + dy;
				double ymin = x.yFromRow(bs.row[i] + bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = x.crop(e, "near", popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, popt);
				p = p.project("EPSG:4326");
				std::vector<double> v = p.area();
				out[0] += accumulate(v.begin(), v.end(), 0);
			}
		} else {
			for (size_t i=0; i<bs.n; i++) {
				double ymax = x.yFromRow(bs.row[i]) + dy;
				double ymin = x.yFromRow(bs.row[i] + bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = x.crop(e, "near", popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, popt);
				p = p.project("EPSG:4326");
				std::vector<double> par = p.area();
			
				std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
				unsigned off = bs.nrows[i] * ncol() ;
				for (size_t lyr=0; lyr<nlyr(); lyr++) {
					unsigned offset = lyr * off;
					unsigned n = offset + off;
					for (size_t j=offset; j<n; j++) {
						if (!std::isnan(v[j])) {
							out[lyr] += par[j - offset];
						}
					}
				}
			}
		
		}
	} else {
		double ar = xres() * yres() * m * m;
		if (!hasValues()) {
			out.resize(1);
			out[0] = ncell() * ar;
		} else {
			for (size_t i=0; i<bs.n; i++) {
				std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
				unsigned off = bs.nrows[i] * ncol() ;
				for (size_t lyr=0; lyr<nlyr(); lyr++) {
					unsigned offset = lyr * off;
					unsigned n = offset + off;
					for (size_t j=offset; j<n; j++) {
						if (!std::isnan(v[j])) out[lyr]++;
					}
				}
			}
			for (size_t lyr=0; lyr<nlyr(); lyr++) {
				out[lyr] = out[lyr] * ar;
			}
		}
	}
	readStop();
	return(out);
}



//layer<value-area
std::vector<std::vector<double>> SpatRaster::area_by_value(SpatOptions &opt) {

	double m = source[0].srs.to_meter();
	m = std::isnan(m) ? 1 : m;

	if (m != 0) {
		double ar = xres() * yres() * m * m;
		std::vector<std::vector<double>> f = freq(true, false, 0, opt);
		for (size_t i=0; i<f.size(); i++) {
			size_t fs = f[i].size();
			for (size_t j=fs/2; j<fs; j++) {
				f[i][j] *= ar;
			}
		}
		return f;
	} else {
		// to do
		// combine freq and area to get area by latitudes
	
		std::vector<std::vector<double>> out(nlyr());
		return out;
	}
}

