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

#include "spatVector.h"

//see http://alienryderflex.com/polygon/
  
std::vector<bool> points_in_polygon(const std::vector<double> &polX, const std::vector<double> &polY, const std::vector<double> &pX, const std::vector<double> &pY) {

	unsigned nodes = polX.size();
	std::vector<double> constant(nodes);
	std::vector<double> multiple(nodes);
	std::vector<bool> result(nodes);
	
	size_t j = nodes-1 ;
	for(size_t i=0; i < nodes; i++) {
		if (polY[j] == polY[i]) {
			constant[i] = polX[i];
			multiple[i] = 0; 
		} else {
			constant[i] = polX[i]-(polY[i]*polX[j])/(polY[j]-polY[i])+(polY[i]*polX[i])/(polY[j]-polY[i]);
			multiple[i] = (polX[j]-polX[i])/(polY[j]-polY[i]); 
		}
		j=i;
	}

	j = nodes-1;
	for (size_t p=0; p<nodes; p++) {
		bool oddNodes = false;
		double x = pX[p];
		double y = pY[p];
		
		for (size_t i=0; i< nodes; i++) {
			if ((((polY[i]< y) && (polY[j]>=y)) || ((polY[j]< y) && (polY[i]>=y)))) {
				oddNodes ^= (y * multiple[i]+constant[i] < x); 
			}
			j=i; 
		}
		result[p] = oddNodes;
	}
	return result;
}



std::vector<int> pointsInPolygons(SpatVector pol, std::vector<double> pX, std::vector<double> pY) {

	unsigned n = pol.size();	
	std::vector<int> result(n, -1);
	
	for (size_t j = 0; j < n; j++) {
			
		SpatGeom geom = pol.getGeom(j);
		unsigned np = geom.size();
		std::vector<bool> inside;	
		for (size_t k = 0; k < np; k++) {
			SpatPart part = geom.getPart(k);
			if (part.hasHoles()) {
				inside = points_in_polygon(part.x, part.y, pX, pY);
				for (size_t h=0; h < part.nHoles(); h++) {
					std::vector<bool> inhole = points_in_polygon(part.x, part.y, pX, pY);
					for (size_t i=0; i<pX.size(); i++) {
						if (inhole[i]) inside[i] = false;
					}
					// remove inhole from inside
				}
			} else {
				inside = points_in_polygon(part.x, part.y, pX, pY);
			}
			for (size_t i=0; i<pX.size(); i++) {
				if (inside[i]) result[i] = j;
			}
		}
	}
	return result;
}
	