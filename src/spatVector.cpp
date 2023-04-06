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
#include <numeric>
#include "string_utils.h"
#include "math_utils.h"
#include "vecmath.h"
#include "recycle.h"

#ifdef useGDAL
	#include "crs.h"
#endif


SpatHole::SpatHole() {}

SpatHole::SpatHole(std::vector<double> X, std::vector<double> Y) {
	x = X; y = Y;
	extent.xmin = *std::min_element(X.begin(), X.end());
	extent.xmax = *std::max_element(X.begin(), X.end());
	extent.ymin = *std::min_element(Y.begin(), Y.end());
	extent.ymax = *std::max_element(Y.begin(), Y.end());
}

bool SpatPart::addHole(std::vector<double> X, std::vector<double> Y) {
	SpatHole h(X, Y);
	holes.push_back(h);
	// check if inside pol?
	return true;
}


bool SpatPart::addHole(SpatHole h) {
	holes.push_back(h);
	// check if inside pol?
	return true;
}


SpatPart::SpatPart() {}

SpatPart::SpatPart(double X, double Y) {
	x.push_back(X);
	y.push_back(Y);
	extent.xmin = X;
	extent.xmax = X;
	extent.ymin = Y;
	extent.ymax = Y;
}

SpatPart::SpatPart(std::vector<double> X, std::vector<double> Y) {
	x = X; y = Y;
	extent.xmin = *std::min_element(X.begin(), X.end());
	extent.xmax = *std::max_element(X.begin(), X.end());
	extent.ymin = *std::min_element(Y.begin(), Y.end());
	extent.ymax = *std::max_element(Y.begin(), Y.end());
}


size_t SpatPart::ncoords() {
	size_t ncrds = x.size();
	size_t nh = holes.size();
	for (size_t k=0; k < nh; k++) {
		ncrds += holes[k].x.size();
	}
	return ncrds;
}


SpatGeom::SpatGeom() {}

SpatGeom::SpatGeom(SpatPart p, SpatGeomType type) {
	parts = {p};
	gtype = type;
	extent = p.extent;
}

SpatGeom::SpatGeom(SpatGeomType g) {
	gtype = g;
}

bool SpatGeom::unite(SpatGeom g) {
	if (parts.empty()) {
		parts = g.parts;
		extent = g.extent;
	} else {
		parts.insert(parts.end(), g.parts.begin(), g.parts.end());
		extent.unite(g.extent);
	}
	return true;
}


bool SpatGeom::addPart(SpatPart p) {
	parts.push_back(p);
	if (parts.size() > 1) {
		extent.unite(p.extent);
	} else {
		extent = p.extent;
	}
	return true;
}


bool SpatGeom::addHole(SpatHole h) {
	long i = parts.size()-1;
	if (i > -1) {
		parts[i].addHole(h);
		return true;
	} else {
		return false;
	}
}


bool SpatGeom::setPart(SpatPart p, unsigned i) {
	parts[i] = p;
	if (parts.size() > 1) {
		extent.unite(p.extent);
	} else {
		extent = p.extent;
	}
	return true;
}

bool SpatGeom::reSetPart(SpatPart p) {
	parts.resize(1);
	parts[0] = p;
	extent = p.extent;
	return true;
}



SpatPart SpatGeom::getPart(unsigned i) {
	return parts[i];
}


size_t SpatGeom::ncoords() {
	size_t ncrds = 0;
	size_t np = parts.size();
	for (size_t j=0; j<np; j++) {
		ncrds += parts[j].x.size();
		if (parts[j].hasHoles()) {
			size_t nh = parts[j].holes.size();
			for (size_t k=0; k < nh; k++) {
				ncrds += parts[j].holes[k].x.size();
			}
		}
	}
	return ncrds;
}


std::vector<std::vector<double>> SpatGeom::coordinates() {
	std::vector<std::vector<double>> out(2);
	size_t np = size();
	size_t ncrds = ncoords();
	out[0].reserve(ncrds);
	out[1].reserve(ncrds);
	for (size_t j=0; j<np; j++) {
		size_t nx = parts[j].x.size();
		for (size_t q=0; q < nx; q++) {
			out[0].insert(out[0].end(), parts[j].x.begin(), parts[j].x.end());
			out[1].insert(out[1].end(), parts[j].y.begin(), parts[j].y.end());
		}
		if (parts[j].hasHoles()) {
			size_t nh =  parts[j].nHoles();
			for (size_t k=0; k < nh; k++) {
				out[0].insert(out[0].end(), parts[j].holes[k].x.begin(), parts[j].holes[k].x.end());
				out[1].insert(out[1].end(), parts[j].holes[k].y.begin(), parts[j].holes[k].y.end());
			}
		}
	}
	return out;
}

void SpatGeom::computeExtent() {
	if (parts.empty()) return;
	extent.xmin = *std::min_element(parts[0].x.begin(), parts[0].x.end());
	extent.xmax = *std::max_element(parts[0].x.begin(), parts[0].x.end());
	extent.ymin = *std::min_element(parts[0].y.begin(), parts[0].y.end());
	extent.ymax = *std::max_element(parts[0].y.begin(), parts[0].y.end());
	for (size_t i=1; i<parts.size(); i++) {
		extent.xmin = std::min(extent.xmin, *std::min_element(parts[0].x.begin(), parts[0].x.end()));
		extent.xmax = std::max(extent.xmin, *std::max_element(parts[0].x.begin(), parts[0].x.end()));
		extent.ymin = std::min(extent.xmin, *std::min_element(parts[0].y.begin(), parts[0].y.end()));
		extent.ymax = std::max(extent.xmin, *std::max_element(parts[0].y.begin(), parts[0].y.end()));
	}
}



SpatVector::SpatVector() {
	extent.xmin = 0;
	extent.xmax = 0;
	extent.ymin = 0;
	extent.ymax = 0;

	srs.proj4 = "";
	srs.wkt = "";

//	srs.proj4 = "+proj=longlat +datum=WGS84";
//	srs.wkt = "GEOGCRS[\"WGS 84\", DATUM[\"World Geodetic System 1984\", ELLIPSOID[\"WGS 84\",6378137,298.257223563, LENGTHUNIT[\"metre\",1]]], PRIMEM[\"Greenwich\",0, ANGLEUNIT[\"degree\",0.0174532925199433]], CS[ellipsoidal,2], AXIS[\"geodetic latitude (Lat)\",north, ORDER[1], ANGLEUNIT[\"degree\",0.0174532925199433]], AXIS[\"geodetic longitude (Lon)\",east, ORDER[2], ANGLEUNIT[\"degree\",0.0174532925199433]], USAGE[ SCOPE[\"Horizontal component of 3D system.\"], AREA[\"World.\"], BBOX[-90,-180,90,180]], ID[\"EPSG\",4326]]";
}

SpatVector::SpatVector(SpatGeom g) {
	addGeom(g);
}



/*
SpatVector::SpatVector(const SpatVector &x) {
	srs = x.srs;
	df = SpatDataFrame(x.df);
}
*/

SpatVector::SpatVector(SpatExtent e, std::string crs) {
	std::vector<double> x = { e.xmin, e.xmin, e.xmax, e.xmax, e.xmin };
	std::vector<double> y = { e.ymin, e.ymax, e.ymax, e.ymin, e.ymin };
	SpatPart p(x, y);
	SpatGeom g(p, polygons);
	setGeom(g);
	setSRS( {crs});
}

SpatVector::SpatVector(std::vector<double> x, std::vector<double> y, SpatGeomType g, std::string crs) {
	if (x.empty()) return;

	if (g == points) {
		SpatPart p(x[0], y[0]);
		SpatGeom geom(p, g);
		setGeom(geom);
		for (size_t i=1; i<x.size(); i++) {
			SpatPart p(x[i], y[i]);
			geom.setPart(p, 0);
			addGeom(geom);
		}
	} else {
		SpatPart p(x, y);
		SpatGeom geom(p, g);
		setGeom(geom);
	}
	setSRS( {crs} );
}

std::vector<double> SpatVector::getDv(unsigned i) {
	return df.getD(i);
}

std::vector<long> SpatVector::getIv(unsigned i){
	return df.getI(i);
}

std::vector<std::string> SpatVector::getSv(unsigned i){
	return df.getS(i);
}

std::vector<unsigned> SpatVector::getItype(){
	return df.itype;
}

std::vector<unsigned> SpatVector::getIplace(){
	return df.iplace;
}


std::vector<std::string> SpatVector::get_names(){
	return df.get_names();
}

void SpatVector::set_names(std::vector<std::string> s){
	df.set_names(s);
}

unsigned SpatVector::ncol() {
	return df.ncol();
}

unsigned SpatVector::nrow() {
	return geoms.size();
}

size_t SpatVector::size() {
	return geoms.size();
}

bool SpatVector::empty() {
	return geoms.empty();
}

bool SpatVector::is_lonlat() {
	if (srs.is_lonlat()) {
		SpatExtent e = getExtent();
		if ((e.xmin < -181) || (e.xmax > 361) || (e.ymin < -90.001) || (e.ymax > 90.001)) {
			addWarning("coordinates are out of range for lon/lat");
		}
		return true;
	}
	return false;
}


bool SpatVector::could_be_lonlat() {
	if (srs.is_lonlat()) return true;
	SpatExtent e = getExtent();
	return srs.could_be_lonlat(e);
}


SpatExtent SpatVector::getExtent(){
	return extent;
}


/*
void SpatVector::setPRJ(std::string PRJ){
	crs[0] = PRJ;
}

std::string SpatVector::getPRJ(){
	return crs[0];
}
*/

std::string SpatVector::type(){
	if (size() == 0) {
		return "none";
	} else if (geoms[0].gtype == points) {
		return "points";
	} else if (geoms[0].gtype == lines) {
		return "lines";
	} else if (geoms[0].gtype == polygons) {
		return "polygons";
	} else if (geoms[0].gtype == null) {
		return("null");
	} else {
		return("unknown");
	}
}



SpatGeom SpatVector::getGeom(unsigned i) {
	return geoms[i];
}

bool SpatVector::addGeom(SpatGeom p) {
	geoms.push_back(p);
	if (geoms.size() > 1) {
		extent.unite(p.extent);
	} else {
		extent = p.extent;
	}
	return true;
}


bool SpatVector::setGeom(SpatGeom p) {
	geoms.resize(1);
	geoms[0] = p;
	extent = p.extent;
	return true;
}

void SpatVector::reserve(size_t n) {
	geoms.reserve(n);
}


void SpatVector::computeExtent() {
	if (geoms.empty()) return;
	extent = geoms[0].extent;
	for (size_t i=1; i<geoms.size(); i++) {	
		extent.unite(geoms[i].extent);
	}
}

std::vector<unsigned> SpatVector::nullGeoms(){
	std::vector<unsigned> ids;
	for (size_t i=0; i<geoms.size(); i++) {
		if ((geoms[i].gtype == null) || (geoms[i].gtype == unknown)) {
			ids.push_back(i);
		}
	}
	return ids;
}



bool SpatVector::replaceGeom(SpatGeom p, unsigned i) {
	if (i < geoms.size()) {
		if ((geoms[i].extent.xmin == extent.xmin) || (geoms[i].extent.xmax == extent.xmax) ||
			(geoms[i].extent.ymin == extent.ymin) || (geoms[i].extent.ymax == extent.ymax)) {
			geoms[i] = p;
			computeExtent();
		} else {
			geoms[i] = p;
		}
	} else {
		return false;
	}
	return true;
}


unsigned SpatVector::nxy() {
	unsigned n = 0;
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		if (g.empty()) {
			n++; // empty
		}
		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			n += p.x.size();
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					n += h.x.size();
				}
			}
		}
	}
	return n;
}

/*
std::vector<std::vector<double>> SpatVector::coordinates() {
	std::vector<std::vector<double>> out(2);
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			for (size_t q=0; q < p.x.size(); q++) {
				out[0].push_back( p.x[q] );
				out[1].push_back( p.y[q] );
			}
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					for (size_t q=0; q < h.x.size(); q++) {
						out[0].push_back( h.x[q] );
						out[1].push_back( h.y[q] );
					}
				}
			}
		}
	}
	return out;
}
*/


size_t SpatVector::ncoords() {
	size_t ncrds = 0;
	size_t ng = geoms.size();
	for (size_t i=0; i<ng; i++) {
		size_t np = geoms[i].parts.size();
		for (size_t j=0; j<np; j++) {
			ncrds += geoms[i].parts[j].x.size();
			if (geoms[i].parts[j].hasHoles()) {
				for (size_t k=0; k < geoms[i].parts[j].nHoles(); k++) {
					ncrds += geoms[i].parts[j].holes[k].x.size();
				}
			}
		}
	}
	return ncrds;
}

size_t SpatVector::nparts(bool holes) {
	size_t totnp = 0;
	size_t ng = geoms.size();
	for (size_t i=0; i<ng; i++) {
		size_t np = geoms[i].parts.size();
		totnp += np;
		if (holes) {
			for (size_t j=0; j<np; j++) {
				totnp += geoms[i].parts[j].nHoles();
			}
		}
	}
	return totnp;
}


std::vector<std::vector<double>> SpatVector::coordinates() {
	std::vector<std::vector<double>> out(2);
	size_t ncrds = ncoords();
	out[0].reserve(ncrds);
	out[1].reserve(ncrds);
	size_t ng = size();
	for (size_t i=0; i<ng; i++) {
		size_t np = geoms[i].size();
		for (size_t j=0; j<np; j++) {
			out[0].insert(out[0].end(), geoms[i].parts[j].x.begin(), geoms[i].parts[j].x.end());
			out[1].insert(out[1].end(), geoms[i].parts[j].y.begin(), geoms[i].parts[j].y.end());
			if (geoms[i].parts[j].hasHoles()) {
				size_t nh = geoms[i].parts[j].nHoles();
				for (size_t k=0; k < nh; k++) {
					out[0].insert(out[0].end(), geoms[i].parts[j].holes[k].x.begin(), geoms[i].parts[j].holes[k].x.end());
					out[1].insert(out[1].end(), geoms[i].parts[j].holes[k].y.begin(), geoms[i].parts[j].holes[k].y.end());
				}
			}
		}
	}
	return out;
}


SpatDataFrame SpatVector::getGeometryDF() {

	SpatDataFrame out;
	out.add_column(1, "geom");
	out.add_column(1, "part");
	out.add_column(0, "x");
	out.add_column(0, "y");
	out.add_column(1, "hole");

	unsigned n = nxy();
	out.resize_rows(n);

	size_t idx = 0;
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		if (g.empty()) { // empty
			out.iv[0][idx] = i+1;
			out.iv[1][idx] = 1;
			out.dv[0][idx] = NAN;
			out.dv[1][idx] = NAN;
			out.iv[2][idx] = 0;
			idx++;
		}

		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			for (size_t q=0; q < p.x.size(); q++) {
				out.iv[0][idx] = i+1;
				out.iv[1][idx] = j+1;
				out.dv[0][idx] = p.x[q];
				out.dv[1][idx] = p.y[q];
				out.iv[2][idx] = 0;
				idx++;
			}
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					for (size_t q=0; q < h.x.size(); q++) {
						out.iv[0][idx] = i+1;
						out.iv[1][idx] = j+1;
						out.dv[0][idx] = h.x[q];
						out.dv[1][idx] = h.y[q];
						out.iv[2][idx] = k+1;
						idx++;
					}
				}
			}
		}
	}
	return out;
}



std::vector<std::vector<double>> SpatVector::getGeometry() {

	unsigned n = nxy();
	std::vector<std::vector<double>> out(5);
	for (size_t i=0; i>out.size(); i++) {
		out[i].reserve(n);
	}
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		if (g.empty()) { // empty
			out[0].push_back(i+1);
			out[1].push_back(1);
			out[2].push_back(NAN);
			out[3].push_back(NAN);
			out[4].push_back(0);
		}

		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			for (size_t q=0; q < p.x.size(); q++) {
				out[0].push_back(i+1);
				out[1].push_back(j+1);
				out[2].push_back(p.x[q]);
				out[3].push_back(p.y[q]);
				out[4].push_back(0);
			}
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					for (size_t q=0; q < h.x.size(); q++) {
						out[0].push_back(i+1);
						out[1].push_back(j+1);
						out[2].push_back(h.x[q]);
						out[3].push_back(h.y[q]);
						out[4].push_back(k+1);
					}
				}
			}
		}
	}
	return out;
}

std::string nice_string(const double &x) {
	std::string s = std::to_string(x);
	s.erase(s.find_last_not_of('0') + 1, std::string::npos);
	s.erase(s.find_last_not_of('.') + 1, std::string::npos);
	return s;
}

std::vector<std::string> SpatVector::getGeometryWKT() {

	std::vector<std::string> out(size());
	std::string wkt;
	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		size_t n = g.size();
		if (g.gtype == points) {
			if (n > 1) {
				wkt = "MULTIPOINT ";
			} else {
				wkt = "POINT ";
			}
		} else if (g.gtype == lines) {
			if (n > 1) {
				wkt = "MULTILINESTRING ";
			} else {
				wkt = "LINESTRING ";
			}
		} else if (g.gtype == polygons) {
			if (n > 1) {
				wkt = "MULTIPOLYGON ";
			} else {
				wkt = "POLYGON ";
			}
		}

		if (n == 0) {
			wkt += "EMPTY";
			out[i] = wkt;
			continue;
		}

		if ((g.gtype == polygons) | (n > 1)) {
			wkt += "(";
		}

		for (size_t j=0; j < n; j++) {
			SpatPart p = g.getPart(j);
			if (j>0) wkt += ",";

			if ((g.gtype == polygons) & (n > 1)) {
				wkt += "(";
			}

			wkt += "(" + nice_string(p.x[0]) + " " + nice_string(p.y[0]);
			for (size_t q=1; q < p.x.size(); q++) {
				wkt += ", " + nice_string(p.x[q]) + " " + nice_string(p.y[q]);
			}
			wkt += ")";
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					wkt += ",(" + nice_string(h.x[0]) + " " + nice_string(h.y[0]);
					for (size_t q=1; q < h.x.size(); q++) {
						wkt += ", " + nice_string(h.x[q]) + " " + nice_string(h.y[q]);
					}
					wkt += ")";
				}
			}
			if ((g.gtype == polygons) & (n > 1)) {
				wkt += ")";
			}
		}
		if ((g.gtype == polygons) | (n > 1)) {
			wkt += ")";
		}
		out[i] = wkt;
	}
	return out;
}


SpatGeomType SpatVector::getGType(std::string &type) {
	if (type == "points") { return points; }
	else if (type == "lines") { return lines; }
	else if (type == "polygons") { return polygons; }
	else { return unknown; }
}


void SpatVector::setGeometry(std::string type, std::vector<unsigned> gid, std::vector<unsigned> part, std::vector<double> x, std::vector<double> y, std::vector<unsigned> hole) {

// it is assumed that values are sorted by gid, part, hole
	unsigned lastgeom = gid[0];
	unsigned lastpart = part[0];
	unsigned lasthole = hole[0];
	bool isHole = lasthole > 0;

	std::vector<double> X, Y;
	SpatGeom g;
	g.gtype = getGType(type);

	for (size_t i=0; i<gid.size(); i++) {
		if ((lastgeom != gid[i]) || (lastpart != part[i]) || (lasthole != hole[i])) {
			if (X.empty()) {
				SpatPart p(NAN, NAN);
				g.addPart(p);
			} else {
				if (g.gtype == polygons) {
					if ((X[0] != X[X.size()-1]) || (Y[0] != Y[Y.size()-1])) {
						X.push_back(X[0]);
						Y.push_back(Y[0]);
					}
					if (isHole) {
						SpatHole h(X, Y);
						g.addHole(h);
					} else {
						SpatPart p(X, Y);
						g.addPart(p);
					}
				} else {
					SpatPart p(X, Y);
					g.addPart(p);
				}
			}
			lastpart = part[i];
            lasthole = hole[i];
            isHole = lasthole > 0;
			X.resize(0);
			Y.resize(0);
			if (lastgeom != gid[i]) {
				addGeom(g);
				g.parts.resize(0);
				lastgeom = gid[i];
			}
		}
        if (!(std::isnan(x[i]) || std::isnan(y[i]))) {
			X.push_back(x[i]);
			Y.push_back(y[i]);
		}
	}

	if (X.empty()) {
		SpatPart p(NAN, NAN);
		g.addPart(p);
	} else {
		if (g.gtype == polygons) {
			if ((X[0] != X[X.size()-1]) || (Y[0] != Y[Y.size()-1])) {
				X.push_back(X[0]);
				Y.push_back(Y[0]);
			}
			if (isHole) {
				SpatHole h(X, Y);
				g.addHole(h);
			} else {
				SpatPart p(X, Y);
				g.addPart(p);
			}
        } else {
			SpatPart p(X, Y);
			g.addPart(p);
		}
	}
	addGeom(g);
}

/*
void SpatVector::setPointsGeometry(std::vector<double> &x, std::vector<double> &y) {
	size_t n = x.size();
	if (n == 0) return;
	reserve(n);
	SpatGeom g;
	g.gtype = points;
	SpatPart p(x[0],y[0]);
	g.addPart(p);;
	for (size_t i=0; i<n; i++) {
		g.parts[0].x[0] = x[i];
		g.parts[0].y[0] = y[i];
		g.extent.xmin = x[i];
		g.extent.xmax = x[i];
		g.extent.ymin = y[i];
		g.extent.ymax = y[i];
		geoms.push_back(g);
	}
	extent.xmin = vmin(x, true);
	extent.xmax = vmax(x, true);
	extent.ymin = vmin(y, true);
	extent.ymax = vmax(y, true);
}
*/

void SpatVector::setPointsGeometry(std::vector<double> &x, std::vector<double> &y) {
	size_t n = x.size();
	//reserve(n)
	if (n == 0) return;
	SpatGeom g;
	g.gtype = points;
	SpatPart p(x[0],y[0]);
	g.addPart(p);
	geoms.resize(n, g);
	for (size_t i=1; i<n; i++) {
		geoms[i].parts[0].x[0] = x[i];
		geoms[i].parts[0].y[0] = y[i];
		geoms[i].extent.xmin = x[i];
		geoms[i].extent.xmax = x[i];
		geoms[i].extent.ymin = y[i];
		geoms[i].extent.ymax = y[i];
	}
	extent.xmin = vmin(x, true);
	extent.xmax = vmax(x, true);
	extent.ymin = vmin(y, true);
	extent.ymax = vmax(y, true);
}



void SpatVector::setPointsDF(SpatDataFrame &x, std::vector<unsigned> geo, std::string crs, bool keepgeom) {
	if (x.nrow() == 0) return;
	if ((x.itype[geo[0]] != 0) || (x.itype[geo[1]] != 0)) {
		setError("coordinates must be numeric");
		return;
	}
	if (geo[0] == geo[1]) {
		setError("x and y coordinates are the same variable");
		return;
	}
	setPointsGeometry(x.dv[x.iplace[geo[0]]], x.dv[x.iplace[geo[1]]]);
	setSRS( {crs} );
	if (!keepgeom) {
		if (geo[0] > geo[1]) {
			x.remove_column(geo[0]);
			x.remove_column(geo[1]);
		} else {
			x.remove_column(geo[1]);
			x.remove_column(geo[0]);
		}
	}
	df = x;
}


SpatVector SpatVector::subset_rows(std::vector<int> range) {

	SpatVector out;
	int n = nrow();
	std::vector<unsigned> r;
	for (size_t i=0; i<range.size(); i++) {
		if ((range[i] >= 0) && (range[i] < n)) {
			r.push_back(range[i]);
		}
	}

	out.reserve(r.size());
	for (size_t i=0; i < r.size(); i++) {
		out.addGeom( geoms[r[i]] );
	}
	out.srs = srs;
	out.df = df.subset_rows(r);
	return out;
}


SpatVector SpatVector::subset_rows(std::vector<unsigned> range) {

	SpatVector out;
	unsigned n = nrow();
	std::vector<unsigned> r;
	out.reserve(r.size());
	for (size_t i=0; i<range.size(); i++) {
		if (range[i] < n) {
			r.push_back(range[i]);
		}
	}

	for (size_t i=0; i < r.size(); i++) {
		out.addGeom( geoms[r[i]] );
	}
	out.srs = srs;
	out.df = df.subset_rows(r);
	return out;
}




SpatVector SpatVector::subset_rows(int i) {
	std::vector<int> range(1, i);
	SpatVector out = subset_rows(range);
	return out;
}


SpatVector SpatVector::remove_rows(std::vector<unsigned> range) {

	std::sort(range.begin(), range.end());
	range.erase(std::unique(range.begin(), range.end()), range.end());
	std::reverse(range.begin(), range.end());
	std::vector<unsigned> id(size());
	std::iota(id.begin(), id.end(), 0);
	unsigned n = size();
	for (size_t i=0; i<range.size(); i++) {
		if (range[i] < n) {
			id.erase(id.begin()+range[i]);
		}
	}
	SpatVector out = subset_rows(id);
	return out;
}




SpatVector SpatVector::subset_cols(std::vector<int> range) {
	int nc = ncol();
	std::vector<unsigned> valid;
	valid.reserve(range.size());
	for (size_t i=0; i<range.size(); i++) {
		if ((range[i] >= 0) && (range[i] < nc)) {
			valid.push_back(range[i]);
		}
	}
	SpatVector out = *this;
	out.df = df.subset_cols(valid);
	return out;
}


SpatVector SpatVector::subset_cols(int i) {
	if (i < 0) {
		SpatVector out;
		out.geoms = geoms;
		out.srs = srs;
		return out;
	}
	std::vector<int> range = {i};
	return subset_cols(range);
}


SpatVector SpatVector::append(SpatVector x, bool ingnorecrs) {
	if (size() == 0) return x;
	if (x.empty()) return *this;

	SpatVector out;
	if (type() != x.type()) {
		out.setError("geom types do not match");
		return out;
	}

	if (!(ingnorecrs)) {
		if (!srs.is_same(x.srs, true)) {
			out.setError("append: crs does not match");
			return out;
		}
	}
	out = *this;
	out.reserve(out.size() + x.size());

	for (size_t i=0; i<x.size(); i++) {
		out.addGeom(x.getGeom(i));
	}
	if ((df.nrow() == 0) && (x.df.nrow() == 0)) {
		return out;
	}
	if ((df.nrow() > 0) && (x.df.nrow() > 0)) {
		out.df.rbind(x.df);
		return out;
	}
	if (x.df.nrow() == 0) {
		out.df.add_rows(x.size());
	} else {
		std::vector<unsigned> i;
		out.df = x.df.subset_rows(i);
		out.df.add_rows(size());
		out.df.rbind(x.df);
	}
	return out;
}



SpatVector SpatVector::cbind(SpatDataFrame d) {
	if (nrow() != d.nrow()) {
		SpatVector out;
		out.setError("nrow does not match");
		return out;
	}
	SpatVector out = *this;
	if (!out.df.cbind(d)) {
		out.setError("cbind failed");
	}
	return out;
}



SpatVector SpatVector::as_points(bool multi, bool skiplast) {
	if (geoms[0].gtype == points) {
		SpatVector v = *this;
		v.addWarning("returning a copy");
		return v;
	}

	SpatVector v = *this;

	if (geoms[0].gtype == lines) {
		for (size_t i=0; i < v.geoms.size(); i++) {
			SpatGeom g;
			g.gtype = points;
			for (size_t j=0; j<geoms[i].parts.size(); j++) {
				SpatPart p = geoms[i].parts[j];
				if (!p.empty()) {
					size_t n = p.size();
					for (size_t k=0; k<n; k++) {
						g.addPart(SpatPart(p.x[k], p.y[k]));
					}
				}
			}
			v.geoms[i] = g;
		}
	} else {
		size_t skip=skiplast;
		for (size_t i=0; i < v.geoms.size(); i++) {
			SpatGeom g;
			g.gtype = points;
			for (size_t j=0; j<geoms[i].parts.size(); j++) {
				SpatPart p = geoms[i].parts[j];
				if (!p.empty()) {
					size_t n = p.size() - skip;
					for (size_t k=0; k<n; k++) {
						g.addPart(SpatPart(p.x[k], p.y[k]));
					}
					if (p.hasHoles()) {
						size_t nh = p.nHoles();
						for (size_t h=0; h<nh; h++) {
							size_t n = p.holes[h].size()-skip;
							for (size_t k=0; k<n; k++) {
								g.addPart(SpatPart(p.holes[h].x[k], p.holes[h].y[k]));
							}
						}
					}
				}
			}
			v.geoms[i] = g;
		}
	}
	if (multi) {
		v.df = df;
	} else {
		v = v.disaggregate(false);
	}
	return(v);
}


SpatVector SpatVector::as_lines() {
	SpatVector v;

	if (geoms[0].gtype == lines) {
		return *this;
	}

	if (geoms[0].gtype == points) {
		std::vector<double> x, y;
		x.reserve(size());
		y.reserve(size());
		for (size_t i=0; i<size(); i++) {
			x.push_back(geoms[i].parts[0].x[0]);
			y.push_back(geoms[i].parts[0].y[0]);
		}
		SpatVector v;
		SpatPart p(x, y);
		SpatGeom g(p, lines)	;
		v.setGeom(g);
		v.srs = srs;
		return v;
	}
/*
		for (int i=geoms.size()-1; i >=0; i--) {
			if (v.geoms[i].parts.size() > 1) {
				for (size_t j=1; j<v.geoms[i].parts.size(); j++) {
					v.geoms[i].parts[0].x.push_back(v.geoms[i].parts[j].x[0]);
					v.geoms[i].parts[0].y.push_back(v.geoms[i].parts[j].y[0]);
					v.geoms[i].gtype = lines;
				}
				v.geoms[i].parts.resize(1);
			} else {
				v.geoms.erase(v.geoms.begin()+i);
			}
		}
		return v;
	}
*/
	// polygons, multipoints
	v = *this;
	for (size_t i=0; i<size(); i++) {
		for (size_t j=0; j < v.geoms[i].size(); j++) {
			SpatPart p = v.geoms[i].parts[j];
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					SpatPart pp(h.x, h.y);
					v.geoms[i].addPart(pp);
				}
				p.holes.resize(0);
				v.geoms[i].parts[j] = p;
			}
		}
		v.geoms[i].gtype = lines;
	}
	v.df = df;
	return(v);
}


void vecround(std::vector<double> &x, int digits) {
	for (double& d : x) d = roundn(d, digits);
}

void remove_duplicates(std::vector<double> &x, std::vector<double> &y, int digits) {
	if (digits > -1) {
		vecround(x, digits);
		vecround(y, digits);
	}
	size_t start = x.size() - 1;
	for (size_t i=start; i>0; i--) {
		if ((x[i] == x[i-1]) && (y[i] == y[i-1])) {
			x.erase(x.begin()+i);
			y.erase(y.begin()+i);
		}
	}
}


void SpatGeom::remove_duplicate_nodes(int digits) {
	size_t start = parts.size()-1;
	for (size_t i=start; i>0; i--) {
		remove_duplicates(parts[i].x, parts[i].y, digits);
		if (parts[i].x.size() < 4) {
			parts.erase(parts.begin()+i);
			continue;
		}
		if (parts[i].hasHoles()) {
			for (size_t j=0; j < parts[i].nHoles(); j++) {
				remove_duplicates(parts[i].holes[j].x, parts[i].holes[j].y, digits);
				if (parts[i].holes[j].x.size() < 4) {
					parts[i].holes.erase(parts[i].holes.begin()+j);
				}
			}
		}
	}
}


SpatVector SpatVector::remove_duplicate_nodes(int digits) {
	SpatVector v = *this;
	if (geoms[0].gtype == points) {
		v.addWarning("returning a copy");
		return v;
	}
	for (size_t i=0; i<size(); i++) {
		v.geoms[i].remove_duplicate_nodes(digits);
	}
	return(v);
}


SpatVector SpatVectorCollection::append() {
	SpatVector out;
	size_t n = size();
	if (n < 1) {
		out.setError("no data in collection");
		return out;
	}
	out = v[0];
	std::string gtype = out.type();
	for (size_t i=1; i<n; i++) {
		if (v[i].empty()) continue;
		if (out.empty()) {
			out = v[i];
			gtype = out.type();
			continue;
		}
		if (v[i].type() != gtype) {
			out.setError("all SpatVectors must have the same geometry type");
			return out;
		}
		//too much copying
		//out = out.append(v[i], true);

		out.geoms.insert(out.geoms.end(), v[i].geoms.begin(), v[i].geoms.end());
		out.extent.unite(v[i].extent);

		if ((out.df.nrow() == 0) && (v[i].df.nrow() == 0)) {
			continue;
		}
		if ((out.df.nrow() > 0) && (v[i].df.nrow() > 0)) {
			out.df.rbind(v[i].df);
			continue;
		}
		if (v[i].df.nrow() == 0) {
			out.df.add_rows(v[i].size());
		} else {
			std::vector<unsigned> r0;
			out.df = v[i].df.subset_rows(r0);
			out.df.add_rows(out.size()-v[i].size());
			out.df.rbind(v[i].df);
		}
	}
	return out;
}

bool SpatVectorCollection::setNames(std::vector<std::string> nms, bool make_valid) {
	recycle(nms, size());
	names = nms;
	return true;
}


SpatVector SpatVector::round(int digits) {
	SpatVector out = *this;
	size_t ng = out.size();
	for (size_t i=0; i<ng; i++) {
		size_t np = out.geoms[i].size();
		for (size_t j=0; j<np; j++) {
			vecround(out.geoms[i].parts[j].x, digits);
			vecround(out.geoms[i].parts[j].y, digits);
			if (out.geoms[i].parts[j].hasHoles()) {
				size_t nh = out.geoms[i].parts[j].holes.size();
				for (size_t k=0; k<nh; k++) {
					vecround(out.geoms[i].parts[j].holes[k].x, digits);
					vecround(out.geoms[i].parts[j].holes[k].y, digits);
				}
			}
		}
		out.geoms[i].computeExtent();
	}
	out.computeExtent();
	return(out);
}


SpatVector SpatVector::normalize_longitude() {
	SpatVector out = *this;
	SpatExtent e = {180, 361, -91, 91};
	SpatVector x = out.crop(e, false);
	if (x.nrow() > 0) {
		x = x.shift(-360, 0);
		SpatVector v(e, "");
		out = out.erase(v);
		out = out.append(x, true);
	}
	e = {-360, -180, -91, 91};
	x = out.crop(e, false);
	if (x.nrow() > 0) {
		x = x.shift(360, 0);
		SpatVector v(e, "");
		out = out.erase(v);
		out = out.append(x, true);
	}
	return out;
}


SpatVector SpatVector::rotate_longitude(double longitude, bool left) {
	SpatVector out = *this;
	size_t ng = out.size();
	for (size_t i=0; i<ng; i++) {
		size_t np = out.geoms[i].size();
		for (size_t j=0; j<np; j++) {
			size_t nx = out.geoms[i].parts[j].x.size();
			for (size_t k=0; k<nx; k++) {
				if (left) {
					if (out.geoms[i].parts[j].x[k] > longitude) {
						out.geoms[i].parts[j].x[k] = out.geoms[i].parts[j].x[k] - 360;
					}
				} else {
					if (out.geoms[i].parts[j].x[k] < longitude) {
						out.geoms[i].parts[j].x[k] = out.geoms[i].parts[j].x[k] + 360;
					}
				}
			}
			if (out.geoms[i].parts[j].hasHoles()) {
				size_t nh = out.geoms[i].parts[j].holes.size();
				for (size_t k=0; k<nh; k++) {
					size_t nx = out.geoms[i].parts[j].holes[k].x.size();
					for (size_t h=0; h<nx; h++) {
						if (left) {
							if (out.geoms[i].parts[j].holes[k].x[h] > longitude) {
								out.geoms[i].parts[j].holes[k].x[h] = out.geoms[i].parts[j].holes[k].x[h] - 360;
							} 
						} else {
							if (out.geoms[i].parts[j].holes[k].x[h] < longitude) {
								out.geoms[i].parts[j].holes[k].x[h] = out.geoms[i].parts[j].holes[k].x[h] + 360;
							} 
						}
					}
				}
			}
		}
		out.geoms[i].computeExtent();
	}
	out.computeExtent();
	return(out);
}


std::vector<std::vector<std::vector<double>>> SpatVector::linesList() {
	size_t ni = nrow();
	std::vector<std::vector<std::vector<double>>> out(ni);
	for (size_t i=0; i < ni; i++) {
		SpatGeom g = getGeom(i);
		size_t nj = g.size();
		if (nj == 0) { // empty
			continue;
		}
		out[i].resize(2);
		size_t ncr = g.ncoords()+nj-1;
		out[i][0].reserve(ncr);
		out[i][1].reserve(ncr);
		for (size_t j=0; j<nj; j++) {
			if (j > 0) {
				out[i][0].push_back(NAN);
				out[i][1].push_back(NAN);
			}
			out[i][0].insert(out[i][0].end(), g.parts[j].x.begin(), g.parts[j].x.end());
			out[i][1].insert(out[i][1].end(), g.parts[j].y.begin(), g.parts[j].y.end());
		}
	}
	return out;
}


std::vector<std::vector<double>> SpatVector::linesNA() {
	size_t ni = nrow();
	size_t n = ncoords() + ni;
	std::vector<std::vector<double>> out(2);
	out[0].reserve(n);
	out[1].reserve(n);
	for (size_t i=0; i < ni; i++) {
		SpatGeom g = getGeom(i);
		size_t nj = g.size();
		for (size_t j=0; j<nj; j++) {
			out[0].insert(out[0].end(), g.parts[j].x.begin(), g.parts[j].x.end());
			out[1].insert(out[1].end(), g.parts[j].y.begin(), g.parts[j].y.end());
			out[0].push_back(NAN);
			out[1].push_back(NAN);
			size_t nk = g.parts[j].nHoles();
			for (size_t k=0; k<nk ; k++) {
				out[0].insert(out[0].end(), g.parts[j].holes[k].x.begin(), g.parts[j].holes[k].x.end());
				out[1].insert(out[1].end(), g.parts[j].holes[k].y.begin(), g.parts[j].holes[k].y.end());
				out[0].push_back(NAN);
				out[1].push_back(NAN);
			}
		}
	}
	out[0].erase(out[0].end() - 1);
	out[1].erase(out[1].end() - 1);
	return out;
}


std::vector<std::vector<std::vector<std::vector<double>>>> SpatVector::polygonsList() {
	size_t ni = nrow();
	std::vector<std::vector<std::vector<std::vector<double>>>> out(ni);
	for (size_t i=0; i < ni; i++) {
		SpatGeom g = getGeom(i);
		size_t nj = g.size();
		if (nj == 0) { // empty
			continue;
		}
		out[i].resize(nj);
		for (size_t j=0; j<nj; j++) {
			out[i][j].resize(2);
			size_t nk = g.parts[j].nHoles();
			if (nk > 0) {
				size_t ncr = g.parts[j].ncoords()+nk;
				out[i][j][0].reserve(ncr);
				out[i][j][1].reserve(ncr);
				out[i][j][0].insert(out[i][j][0].end(), g.parts[j].x.begin(), g.parts[j].x.end());
				out[i][j][1].insert(out[i][j][1].end(), g.parts[j].y.begin(), g.parts[j].y.end());
				for (size_t k=0; k<nk ; k++) {
					out[i][j][0].push_back(NAN);
					out[i][j][1].push_back(NAN);
					out[i][j][0].insert(out[i][j][0].end(), g.parts[j].holes[k].x.begin(), g.parts[j].holes[k].x.end());
					out[i][j][1].insert(out[i][j][1].end(), g.parts[j].holes[k].y.begin(), g.parts[j].holes[k].y.end());
				}
			} else {
				out[i][j][0] = g.parts[j].x;
				out[i][j][1] = g.parts[j].y;
			}
		}
	}
	return out;
}


