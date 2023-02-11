// Copyright (c) 2018-2023  Robert J. Hijmans
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

#include "spatVector.h"
#include "string_utils.h"
#include "vecmath.h"
#include "recycle.h"

#include "gdal_alg.h"
#include "ogrsf_frmts.h"

/*
std::vector<bool> SpatVector::is_valid() {
	std::vector<bool> out;
	out.reserve(nrow());
	GDALDataset* src;
	if (!write_ogr(src, "", "layer", "Memory", false, false, true, std::vector<std::string>())) {
		if (src != NULL) GDALClose( src );
		setError("cannot do it");
		return false;
	}
	OGRLayer *inLayer = src->GetLayer(0);
	inLayer->ResetReading();
	OGRFeature *inFeature;
	while( (inFeature = inLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = inFeature->GetGeometryRef();
		out.push_back(poGeometry->IsValid());
		OGRFeature::DestroyFeature( inFeature );
	}
	return out;
}


SpatVector SpatVector::make_valid() {
	SpatVector out;
	GDALDataset* src = write_ogr("", "layer", "Memory", false, false, true, std::vector<std::string>());
	OGRLayer *inLayer = src->GetLayer(0);
	inLayer->ResetReading();
	OGRFeature *inFeature;
	while( (inFeature = inLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = inFeature->GetGeometryRef();
		//OGRGeometry *poGeom = poGeometry->MakeValid();
		if (inFeature->SetGeometry( poGeometry ) != OGRERR_NONE) {
			out.setError("cannot set geometry");
			return out;
		}
		if (inLayer->SetFeature( inFeature ) != OGRERR_NONE) {
			out.setError("cannot set feature");
			return out;
		}
		OGRFeature::DestroyFeature( inFeature );
	}
	std::vector<double> fext;
	SpatVector fvct;
	out.read_ogr(src, "", "", fext, fvct, false, "");
	GDALClose(src);
	return out;
}
*/


SpatVector SpatVector::disaggregate(bool segments) {

	SpatVector out;
	out.srs = srs;
	out.df = df.skeleton();
	if (nrow() == 0) {
		return out;
	}
	size_t n=0;
	for (size_t i=0; i<nrow(); i++) {
		n += geoms[i].parts.size();
	}
	out.reserve(n);

	for (size_t i=0; i<nrow(); i++) {
		SpatGeom g = getGeom(i);
		SpatDataFrame row = df.subset_rows(i);
		for (size_t j=0; j<g.parts.size(); j++) {
			SpatGeom gg = SpatGeom(g.parts[j], g.gtype);
			out.addGeom(gg);
			if (!out.df.rbind(row)) {
				out.setError("cannot add row");
				return out;
			}
		}
	}
	if (segments && (type() != "points")) {
		SpatVector x;
		x.srs = srs;
		x.df = df.skeleton();

		for (size_t i=0; i<out.nrow(); i++) {
			SpatGeom g = out.getGeom(i);
			SpatDataFrame row = out.df.subset_rows(i);
			size_t n = g.parts[0].x.size() - 1;
			for (size_t j=0; j<n; j++) {
				std::vector<double> sx = {g.parts[0].x[j], g.parts[0].x[j+1]};
				std::vector<double> sy = {g.parts[0].y[j], g.parts[0].y[j+1]};
				SpatPart p(sx, sy);
				SpatGeom gg = SpatGeom(p, lines);
				x.addGeom(gg);
				if (!x.df.rbind(row)) {
					x.setError("cannot add row");
					return x;
				}
			}
		}
		return x;
	}

	return out;

}


SpatVector SpatVector::aggregate(std::string field, bool dissolve) {

	SpatVector out;
	int i = where_in_vector(field, get_names(), false);
	if (i < 0) {
		out.setError("cannot find field: " + field);
		return out;
	}
	SpatDataFrame uv;
	std::vector<int> idx = df.getIndex(i, uv);
	out.reserve(uv.nrow());
	for (size_t i=0; i<uv.nrow(); i++) {
		SpatGeom g;
		g.gtype = geoms[0].gtype;
		for (size_t j=0; j<idx.size(); j++) {
			if (i == (size_t)idx[j]) {
				g.unite( getGeom(j) );
			}
		}
		out.addGeom(g);
	}
	if (dissolve) {
		out = out.unaryunion();
	}
	out.srs = srs;
	out.df  = uv;
	return out;
}



SpatVector SpatVector::aggregate(bool dissolve) {
	SpatVector out;
	SpatGeom g;
	g.gtype = geoms[0].gtype;
	for (size_t i=0; i<size(); i++) {
		g.unite( getGeom(i) );
	}
	out.addGeom(g);
	if (dissolve) {
		out = out.unaryunion();
	}
	out.srs = srs;
	return out;
}


#include "geodesic.h"

void extend_line(const double &x1, const double &y1, const double &x2, const double &y2, double &x, double &y, const bool &geo, const double &distance, bool plus) {
	if (geo) {
		double a = 6378137.0;
		double f = 1/298.257223563;
		double s12, azi1, azi2;
		struct geod_geodesic g;
		geod_init(&g, a, f);
		geod_inverse(&g, y1, x1, y2, x2, &s12, &azi1, &azi2);		
		geod_direct(&g, y2, x2, azi2, distance, &y, &x, &azi1);
	} else {
		double bearing;
		double dx = x2 - x1;
		if (dx == 0) {
			if (y2 > y1) {
				bearing = 0;
			} else {
				bearing = M_PI;				
			}
		} else {
			bearing = atan((y2-y1)/dx);
		}
		if (plus) {
			x = x2 + distance * sin(bearing);
			y = y2 + distance * cos(bearing);
		} else {
			x = x2 - distance * sin(bearing);
			y = y2 - distance * cos(bearing);			
		}
	}
}

SpatVector SpatVector::elongate(double length) {

	SpatVector out = *this;
	size_t n = size();
	if (n == 0) {
		return out;
	}
	if (geoms[0].gtype != lines) {
		out.setError("you can only elongate lines");
		return out;
	}
	if (length < 0) {
		out.setError("length must be > 0");
		return out;
	}
	if (length == 0) {
		return out;
	}

	double x, y;
	bool geo = is_lonlat();
	
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j < out.geoms[i].size(); j++) {
			SpatPart p = out.geoms[i].parts[j];
			size_t n = p.x.size();
			if (n < 2) continue;
			extend_line(p.x[1], p.y[1], p.x[0], p.y[0], x, y, geo, length, false);
			p.x.insert(p.x.begin(), x);
			p.y.insert(p.y.begin(), y);
			
			extend_line(p.x[n-1], p.y[n-1], p.x[n], p.y[n], x, y, geo, length, true);
			p.x.push_back(x);			
			p.y.push_back(y);
			
			out.geoms[i].parts[j] = p;
		}
		out.geoms[i].computeExtent();
	}
	out.computeExtent();
	return out;
}



SpatVectorCollection SpatVector::split(std::string field) {

	SpatVectorCollection out;

	int i = where_in_vector(field, get_names(), false);
	if (i < 0) {
		out.setError("cannot find field: " + field);
		return out;
	}
	SpatDataFrame uv;
	std::vector<int> idx = df.getIndex(i, uv);

	for (size_t i=0; i<uv.nrow(); i++) {
		SpatVector v;
		std::vector<unsigned> r;
		for (size_t j=0; j<idx.size(); j++) {
			if (i == (size_t)idx[j]) {
				v.addGeom( getGeom(j) );
				r.push_back(j);
			}
		}
		v.srs = srs;
		v.df = df.subset_rows(r);
		out.push_back(v);
	}
	return out;
}



SpatVector SpatVector::remove_holes() {

	SpatVector out = *this;

	size_t n = size();
	if (n == 0) {
		return out;
	}
	if (geoms[0].gtype != polygons) {
		return out;
	}

	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j < out.geoms[i].size(); j++) {
			SpatPart p = out.geoms[i].parts[j];
			if (p.hasHoles()) {
				p.holes.resize(0);
				out.geoms[i].parts[j] = p;
			}
		}
	}
	return out;
}


SpatVector SpatVector::get_holes() {

	SpatVector out;
	out.srs = srs;
	size_t n = size();
	if (n == 0) {
		return out;
	}
	if (geoms[0].gtype != polygons) {
		return out;
	}
	std::vector<unsigned> atts;

	for (size_t i=0; i<n; i++) {
		SpatGeom g;
		g.gtype = polygons;
		bool found = false;
		for (size_t j=0; j < geoms[i].size(); j++) {
			SpatPart p = geoms[i].parts[j];
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatPart h(p.holes[k].x, p.holes[k].y);
					g.addPart(h);
				}
				found = true;
			}
		}
		if (found) {
			out.addGeom(g);
			atts.push_back(i);
		}
	}
	out.df = df.subset_rows(atts);
	return out;
}


SpatVector SpatVector::set_holes(SpatVector x, size_t i) {

	SpatVector out;
	if (size() == 0) {
		out.setError("object has no geometries");
		return out;
	}
	if (i > size()) {
		out.setError("invalid index");
		return out;
	}
	if (x.type() != "polygons") {
		out.setError("holes must be polygons");
		return out;
	}
	if (out.geoms[i].size() > 1) {
		out.setError("selected object has multiple geometries");
	}

	x = x.unaryunion();
	SpatPart p  = out.geoms[i].parts[0];
	SpatGeom g =   x.geoms[0];
	for (size_t i=0; i<g.size(); i++) {
		p.addHole(g.parts[i].x, g.parts[i].y);
	}
	out = *this;
	out.geoms[i].parts[0] = p;
	return out;
}




/*
std::vector<OGRGeometry *> geoms_from_ds(GDALDataset* src, int field, int value) {
	std::vector<OGRGeometry *> g;
	OGRLayer *poLayer = src->GetLayer(0);
	poLayer->ResetReading();
	OGRFeature *poFeature;

	while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		g.push_back(poGeometry);
	}
	return g;
}
// create output dataset
	GDALDataset* dst;
// get unique values in field
// loop over unique values
	// for value in uvalues
	std::vector<OGRGeometry *> gvec = geoms_from_ds(src, field, value);
	OGRGeometry *geom;
	geom = (OGRGeometry *) gvec.data();
	OGRGeometry *gout;
	gout = geom->UnionCascaded();
// set geometry to output
   return dst;
*/




SpatVector SpatVector::shift(double x, double y) {

	SpatVector out = *this;

	for (size_t i=0; i < size(); i++) {
		for (size_t j=0; j < geoms[i].size(); j++) {
			for (size_t q=0; q < geoms[i].parts[j].x.size(); q++) {
				out.geoms[i].parts[j].x[q] += x;
				out.geoms[i].parts[j].y[q] += y;
			}
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					for (size_t q=0; q < geoms[i].parts[j].holes[k].x.size(); q++) {
						out.geoms[i].parts[j].holes[k].x[q] += x;
						out.geoms[i].parts[j].holes[k].y[q] += y;
					}
					out.geoms[i].parts[j].holes[k].extent.xmin += x;
					out.geoms[i].parts[j].holes[k].extent.xmax += x;
					out.geoms[i].parts[j].holes[k].extent.ymin += y;
					out.geoms[i].parts[j].holes[k].extent.ymax += y;
				}
			}
			out.geoms[i].parts[j].extent.xmin += x;
			out.geoms[i].parts[j].extent.xmax += x;
			out.geoms[i].parts[j].extent.ymin += y;
			out.geoms[i].parts[j].extent.ymax += y;
		}
		out.geoms[i].extent.xmin += x;
		out.geoms[i].extent.xmax += x;
		out.geoms[i].extent.ymin += y;
		out.geoms[i].extent.ymax += y;
	}
	out.extent.xmin += x;
	out.extent.xmax += x;
	out.extent.ymin += y;
	out.extent.ymax += y;
	return out;
}


void resc(double &value, const double &base, const double &f) {
	value = base + f * (value - base);
}


SpatVector SpatVector::rescale(double fx, double fy, double x0, double y0) {

	SpatVector out = *this;
	for (size_t i=0; i < size(); i++) {
		for (size_t j=0; j < geoms[i].size(); j++) {
			for (size_t q=0; q < geoms[i].parts[j].x.size(); q++) {
				resc(out.geoms[i].parts[j].x[q], x0, fx);
				resc(out.geoms[i].parts[j].y[q], y0, fy);
			}
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					for (size_t q=0; q < geoms[i].parts[j].holes[k].x.size(); q++) {
						resc(out.geoms[i].parts[j].holes[k].x[q], x0, fx);
						resc(out.geoms[i].parts[j].holes[k].y[q], y0, fy);
					}
					resc(out.geoms[i].parts[j].holes[k].extent.xmax, x0, fx);
					resc(out.geoms[i].parts[j].holes[k].extent.ymax, y0, fy);
				}
			}
			resc(out.geoms[i].parts[j].extent.xmin, x0, fx);
			resc(out.geoms[i].parts[j].extent.xmax, x0, fx);
			resc(out.geoms[i].parts[j].extent.ymin, y0, fy);
			resc(out.geoms[i].parts[j].extent.ymax, y0, fy);
		}
		resc(out.geoms[i].extent.xmin, x0, fx);
		resc(out.geoms[i].extent.xmax, x0, fx);
		resc(out.geoms[i].extent.ymin, y0, fy);
		resc(out.geoms[i].extent.ymax, y0, fy);
	}
	resc(out.extent.xmin, x0, fx);
	resc(out.extent.xmax, x0, fx);
	resc(out.extent.ymin, y0, fy);
	resc(out.extent.ymax, y0, fy);
	return out;
}

void dswap(double &a, double&b) {
	double tmp = a;
	a = b;
	b = tmp;
}

SpatVector SpatVector::transpose() {

	SpatVector out = *this;
	for (size_t i=0; i < size(); i++) {
		for (size_t j=0; j < geoms[i].size(); j++) {
			out.geoms[i].parts[j].x.swap(out.geoms[i].parts[j].y);
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					out.geoms[i].parts[j].holes[k].x.swap(out.geoms[i].parts[j].holes[k].y);

					dswap(out.geoms[i].parts[j].holes[k].extent.xmin,
						 out.geoms[i].parts[j].holes[k].extent.ymin);
					dswap(out.geoms[i].parts[j].holes[k].extent.xmax,
						 out.geoms[i].parts[j].holes[k].extent.ymax);
				}
			}
			dswap(out.geoms[i].parts[j].extent.xmin,
				 out.geoms[i].parts[j].extent.ymin);
			dswap(out.geoms[i].parts[j].extent.xmax,
				 out.geoms[i].parts[j].extent.ymax);
		}
		dswap(out.geoms[i].extent.xmin, out.geoms[i].extent.ymin);
		dswap(out.geoms[i].extent.xmax, out.geoms[i].extent.ymax);
	}
	dswap(out.extent.xmin, out.extent.ymin);
	dswap(out.extent.xmax, out.extent.ymax);
	return out;
}


void flipd(double &value, const double &base) {
	value = base - (value - base);
}

void flipv(std::vector<double> &v, const double &base) {
	for (double &d : v) d = base - (d - base);
}

SpatVector SpatVector::flip(bool vertical) {
	double x0 = extent.xmin;
	double y0 = extent.ymin;
	SpatVector out = *this;
	bool horizontal = !vertical;
	for (size_t i=0; i < size(); i++) {
		for (size_t j=0; j < geoms[i].size(); j++) {
			if (horizontal) {
				flipv(out.geoms[i].parts[j].x, x0);
				flipd(out.geoms[i].parts[j].extent.xmin, x0);
				flipd(out.geoms[i].parts[j].extent.xmax, x0);
				dswap(out.geoms[i].parts[j].extent.xmin, out.geoms[i].parts[j].extent.xmax);
			} else {
				flipv(out.geoms[i].parts[j].y, y0);
				flipd(out.geoms[i].parts[j].extent.ymin, y0);
				flipd(out.geoms[i].parts[j].extent.ymax, y0);
				dswap(out.geoms[i].parts[j].extent.ymin, out.geoms[i].parts[j].extent.ymax);
			}
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					if (horizontal) {
						flipv(out.geoms[i].parts[j].holes[k].x, x0);
						flipd(out.geoms[i].parts[j].holes[k].extent.xmin, x0);
						flipd(out.geoms[i].parts[j].holes[k].extent.xmax, x0);
						dswap(out.geoms[i].parts[j].holes[k].extent.xmin,
							  out.geoms[i].parts[j].holes[k].extent.xmax);
					} else {
						flipv(out.geoms[i].parts[j].holes[k].y, y0);
						flipd(out.geoms[i].parts[j].holes[k].extent.ymin, y0);
						flipd(out.geoms[i].parts[j].holes[k].extent.ymax, y0);
						dswap(out.geoms[i].parts[j].holes[k].extent.ymin,
							  out.geoms[i].parts[j].holes[k].extent.ymax);
					}
				}
			}
		}
		if (horizontal) {
			flipd(out.geoms[i].extent.xmin, x0);
			flipd(out.geoms[i].extent.xmax, x0);
			dswap(out.geoms[i].extent.xmin, out.geoms[i].extent.xmax);
		} else {
			flipd(out.geoms[i].extent.ymin, y0);
			flipd(out.geoms[i].extent.ymax, y0);
			dswap(out.geoms[i].extent.ymin, out.geoms[i].extent.ymax);
		}
	}
	if (horizontal) {
		flipd(out.extent.xmin, x0);
		flipd(out.extent.xmax, x0);
		dswap(out.extent.xmin, out.extent.xmax);
	} else {
		flipd(out.extent.ymin, y0);
		flipd(out.extent.ymax, y0);
		dswap(out.extent.ymin, out.extent.ymax);
	}
	return out;
}



void rotit(std::vector<double> &x, std::vector<double> &y, const double &x0, const double &y0, const double &cos_angle, const double &sin_angle) {
	for (size_t i=0; i<x.size(); i++) {
		double sx = x[i] - x0;
		double sy = y[i] - y0;
		x[i] = sx * cos_angle + sy * -sin_angle + x0;
		y[i] = sx * sin_angle + sy * cos_angle + y0;
	}
}



SpatVector SpatVector::rotate(double angle, std::vector<double> x0, std::vector<double> y0) {
	angle = -M_PI * angle / 180;
	size_t n = size();
	if ((x0.size() == 0) || (y0.size() == 0)) {
		SpatVector out;
		out.setError("no center of rotation provided");
		return out;
	}
	bool multi = true;
	double ix0, iy0; 
	if ((x0.size() == 1) && (y0.size() == 1)) {
		multi = false;
		ix0 = x0[0];
		iy0 = y0[0]; 
	} else {
		recycle(x0, n);
		recycle(y0, n);
	}
	double cos_angle = cos(angle);
	double sin_angle = sin(angle);
	SpatVector out = *this;
	for (size_t i=0; i < n; i++) {
		if (multi) {
			ix0 = x0[i];
			iy0 = y0[i];
		} 
		for (size_t j=0; j < geoms[i].size(); j++) {
			rotit(out.geoms[i].parts[j].x, out.geoms[i].parts[j].y, ix0, iy0, cos_angle, sin_angle);
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					rotit(out.geoms[i].parts[j].holes[k].x,
						  out.geoms[i].parts[j].holes[k].y, ix0, iy0, cos_angle, sin_angle);

					out.geoms[i].parts[j].holes[k].extent.xmin =
						vmin(out.geoms[i].parts[j].holes[k].x, true);
					out.geoms[i].parts[j].holes[k].extent.xmax =
						vmax(out.geoms[i].parts[j].holes[k].x, true);
					out.geoms[i].parts[j].holes[k].extent.ymin =
						vmin(out.geoms[i].parts[j].holes[k].y, true);
					out.geoms[i].parts[j].holes[k].extent.ymax =
						vmax(out.geoms[i].parts[j].holes[k].y, true);
				}
			}
			out.geoms[i].parts[j].extent.xmin = vmin(out.geoms[i].parts[j].x, true);
			out.geoms[i].parts[j].extent.xmax = vmax(out.geoms[i].parts[j].x, true);
			out.geoms[i].parts[j].extent.ymin = vmin(out.geoms[i].parts[j].y, true);
			out.geoms[i].parts[j].extent.ymax = vmax(out.geoms[i].parts[j].y, true);
			if (j==0) {
				out.geoms[i].extent = out.geoms[i].parts[j].extent;
			} else {
				out.geoms[i].extent.unite(out.geoms[i].parts[j].extent);
			}
		}
		if (i==0) {
			out.extent = out.geoms[i].extent;
		} else {
			out.extent.unite(out.geoms[i].extent);
		}
	}
	return out;
}


inline double cartdist(const double& x1, const double& y1, const double &x2, const double &y2) {
	return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}



bool thinnodes(std::vector<double> &x, std::vector<double> &y, const double &threshold, const size_t &mnsize) {
	std::vector<double> xout, yout;
	size_t n = x.size();
	xout.reserve(n);
	yout.reserve(n);
	n--;
	for (size_t i=0; i<n; i++) {
		if (cartdist(x[i], y[i], x[i+1], y[i+1]) <= threshold) {
			xout.push_back((x[i] + x[i+1])/2);
			yout.push_back((y[i] + y[i+1])/2);
		} else {
			xout.push_back(x[i]);
			yout.push_back(y[i]);
		}
	}
	if (cartdist(x[n], y[n], xout[0], yout[0]) <= threshold) {
		xout.push_back((x[n] + xout[0])/2);
		yout.push_back((y[n] + yout[0])/2);
		xout[0] = xout[n];
		yout[0] = xout[n];
	} else {
		xout.push_back(xout[0]);
		yout.push_back(yout[0]);
	}
	if (xout.size() == (n+1)) {
		return false;
	}
	if (xout.size() >= mnsize) {
		x = std::move(xout);
		y = std::move(yout);
		return true;
	}
	return false;
}



SpatVector SpatVector::thin(double threshold) {

	SpatVector out;
	if (threshold < 0) {
		out.setError("threshold must be a positive number");
		return out;
	}
	size_t mnode = 4;
	if (geoms[0].gtype == lines) {
		mnode = 3;
	} else if (geoms[0].gtype != polygons) {
		out.setError("can only thin lines or polygons");
		return out;
	}

	out = *this;
	bool objext = false;
	for (size_t i=0; i < size(); i++) {
		bool geomext = false;
		for (size_t j=0; j < out.geoms[i].size(); j++) {
			if (thinnodes(out.geoms[i].parts[j].x, out.geoms[i].parts[j].y, threshold, mnode)) {
				geomext = true;
			}
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					thinnodes(geoms[i].parts[j].holes[k].x, geoms[i].parts[j].holes[k].y, threshold, mnode);
				}
			}
		}
		if (geomext) {
			objext = true;
			geoms[i].computeExtent();
		}
	}
	if (objext) {
		computeExtent();
	}

	return out;
}






/*
SpatVector SpatVector::removeSlivers(double dthres, double athres, size_t n) {

	SpatVector out;
	if (geoms[0].gtype != polygons) {
		out.setError("can only remove slivers from polygons");
		return out;
	}
	if ((dthres < 0) || (athres < 0)) {
		out.setError("thresholds must be a positive number");
		return out;
	}
	if (n < 2)) {
		out.setError("n must be at least 2");
		return out;
	}

	out = *this;
	bool objext = false;
	for (size_t i=0; i < size(); i++) {
		bool geomext = false;
		for (size_t j=0; j < out.geoms[i].size(); j++) {
			if (remove_slivers(out.geoms[i].parts[j].x, out.geoms[i].parts[j].y, threshold, mnode)) {
				geomext = true;
			}
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					remove_slivers(geoms[i].parts[j].holes[k].x, geoms[i].parts[j].holes[k].y, threshold, mnode);
				}
			}
		}
		if (geomext) {
			objext = true;
			geoms[i].computeExtent();
		}
	}
	if (objext) {
		computeExtent();
	}
	return out;
}
*/

