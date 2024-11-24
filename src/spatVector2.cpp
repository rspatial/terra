// Copyright (c) 2018-2025  Robert J. Hijmans
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
		for (size_t j=G[i]; j<G[i+1]; j++) { 
			std::vector<double> x = {X.begin() + P[j], X.begin() + P[j+1]};
			std::vector<double> y = {Y.begin() + P[j], Y.begin() + P[j+1]};
			if (gtype == polygons) {
				if (H[j] >= 0) {
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
	if (!x.empty()) out.gtype = x.geoms[0].gtype;
	size_t nxy = x.nxy();
	out.X.reserve(nxy);
	out.Y.reserve(nxy);
	size_t ng = x.size();
	out.G.reserve(ng);
	size_t np = x.nparts(true);
	out.P.reserve(np);
	if (x.type() == "polygons") {
		out.H.reserve(ng);
	}
	
	size_t pcnt = 0;
	size_t gcnt = 0;
	out.G.push_back(0); // so that we can use (j to j+1) for the first part
	out.P.push_back(0); // so that we can use (j to j+1) for the first part
	for (size_t i=0; i<ng; i++) {
		SpatGeom xg = x.getGeom(i);
		size_t np = xg.size();
		if (np == 0) { // empty
			pcnt++;
			out.P.push_back(pcnt);
			out.X.push_back(NAN);
			out.Y.push_back(NAN);
			if (x.type() == "polygons") {
				out.H.push_back(-1);
			}
		}
		for (size_t j=0; j < np; j++) {
			SpatPart prt = xg.getPart(j);
			out.X.insert(out.X.end(), prt.x.begin(), prt.x.end());
			out.Y.insert(out.Y.end(), prt.y.begin(), prt.y.end());
			pcnt += prt.x.size();
			out.P.push_back(pcnt);
			out.H.push_back(-1);
			if (prt.hasHoles()) {
				for (size_t k=0; k < prt.nHoles(); k++) {
					SpatHole hle = prt.getHole(k);
					out.X.insert(out.X.end(), hle.x.begin(), hle.x.end());
					out.Y.insert(out.Y.end(), hle.y.begin(), hle.y.end());
					pcnt += hle.x.size();
					out.P.push_back(pcnt);
					out.H.push_back(j);
					gcnt++;
				}
			}
		}
		gcnt += np;
		out.G.push_back(gcnt);
	}
	return out;
}

size_t SpatVector2::ngeoms() {
	return (G.size() - 1);
}


*/	