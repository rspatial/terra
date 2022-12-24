// Copyright (c) 2018-2022  Robert J. Hijmans
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
#include "spatVector2.h"

SpatVector2::SpatVector2() {}


SpatVector SpatVector2::to_old() {
	SpatVector out;
	out.srs = srs;

	size_t ng = ngeoms();
	for (size_t i=0; i<ng; i++) { 
		SpatGeom geom;
		geom.gtype = gtype;
		for (size_t j=g[i]; j<g[i+1]; j++) { 
			std::vector<double> x = {xc.begin() + p[j], xc.begin() + p[j+1]};
			std::vector<double> y = {yc.begin() + p[j], yc.begin() + p[j+1]};
			if (gtype == polygons) {
				if (h[j] >= 0) {
					geom.parts[geom.parts.size()-1].addHole(x, y);
				} else {
					SpatPart prt(x, y);
					geom.addPart(prt);
				}
			} else {
				SpatPart prt(x, y);
				geom.addPart(prt);
			}
		}
		out.addGeom(geom);
	}
	return out;
}
		
SpatVector2 SpatVector2::from_old(SpatVector x) {
	SpatVector2 out;

	out.srs = x.srs;
	if (x.size() > 0) out.gtype = x.geoms[0].gtype;
	size_t nxy = x.nxy();
	out.xc.reserve(nxy);
	out.yc.reserve(nxy);
	size_t ng = x.size();
	out.g.reserve(ng);
	size_t np = x.nparts(true);
	out.p.reserve(np);
	if (x.type() == "polygons") {
		out.h.reserve(ng);
	}
	
	size_t pcnt = 0;
	size_t gcnt = 0;
	out.g.push_back(0); // so that we can use (j to j+1) for the first part
	out.p.push_back(0); // so that we can use (j to j+1) for the first part
	for (size_t i=0; i<ng; i++) {
		SpatGeom xg = x.getGeom(i);
		size_t np = xg.size();
		if (np == 0) { // empty
			pcnt++;
			out.p.push_back(pcnt);
			out.xc.push_back(NAN);
			out.yc.push_back(NAN);
			if (x.type() == "polygons") {
				out.h.push_back(-1);
			}
		}
		for (size_t j=0; j < np; j++) {
			SpatPart prt = xg.getPart(j);
			out.xc.insert(out.xc.end(), prt.x.begin(), prt.x.end());
			out.yc.insert(out.yc.end(), prt.y.begin(), prt.y.end());
			pcnt += prt.x.size();
			out.p.push_back(pcnt);
			out.h.push_back(-1);
			if (prt.hasHoles()) {
				for (size_t k=0; k < prt.nHoles(); k++) {
					SpatHole hle = prt.getHole(k);
					out.xc.insert(out.xc.end(), hle.x.begin(), hle.x.end());
					out.yc.insert(out.yc.end(), hle.y.begin(), hle.y.end());
					pcnt += hle.x.size();
					out.p.push_back(pcnt);
					out.h.push_back(j);
					gcnt++;
				}
			}
		}
		gcnt += np;
		out.g.push_back(gcnt);
	}
	return out;
}

size_t SpatVector2::ngeoms() {
	return (g.size() - 1);
}


	