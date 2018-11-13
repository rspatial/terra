// Copyright (c) 2018  Robert J. Hijmans
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

// Copyright (c) 2018  Robert J. Hijmans
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

// Robert Hijmans, October 2011 
// This is an implementation of J. Ronald Eastman's pushbroom algorithm


#include "spatraster.h"	
#include <limits>
#include <cmath>

std::vector<double> broom_dist_planar(std::vector<double> &v, std::vector<double> &above, std::vector<double> res, std::vector<unsigned> dim, bool down) {
	
	double dx = res[0];
	double dy = res[1];
	double dxy = sqrt(dx * dx + dy *dy);	

	size_t n = v.size();
	unsigned nr = n / dim[0]; // must get entire rows
	unsigned nc = dim[1];

	double inf = std::numeric_limits<double>::infinity();
	std::vector<double> dist(n, 0);
	
	if (down) {	
		//left to right	
		if ( isnan(v[0])) { //first cell, no cell left of it
			dist[0] = above[0] + dy;
		} 
		for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
			if (isnan(v[i])) {
				dist[i] = std::min(std::min(above[i] + dy, above[i-1] + dxy), dist[i-1] + dx);
			} 
		}	
		for (size_t r=1; r<nr; r++) { //other rows
			size_t i=r*nc;
			if (isnan(v[i])) {
				dist[i] = dist[i-nc] + dy;
			} 
			for (size_t i=r*nc+1; i<((r+1)*nc); i++) {
				if (isnan(v[i])) {
					dist[i] = std::min(std::min(dist[i-1] + dx, dist[i-nc] + dy), dist[i-nc-1] + dxy);
				}
			}
		}

		//right to left
		if ( isnan(v[nc-1])) { //first cell
			dist[nc-1] = std::min(dist[nc-1], above[nc-1] + dy);
		} 				
		for (size_t i=(nc-2); i > -1; i--) { // other cells on first row
			if (isnan(v[i])) {
				dist[i] = std::min(std::min(std::min(dist[i], above[i] + dy), above[i+1] + dxy), dist[i+1] + dx);
			} 
		}
		for (size_t r=1; r<nr; r++) { // other rows
			size_t i=(r+1)*nc-1;
			if (isnan(v[i])) {
				dist[i] = std::min(dist[i], dist[i-nc] + dy);
			} 
			for (size_t i=(r+1)*nc-2; i>(r*nc-1); i--) {
				if (isnan(v[i])) {
					dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy), dist[i-nc+1] + dxy);		
				} 
			}
		}

	} else {  // bottom to top
		// left to right
		size_t r = nr-1; // first (=last) row
		size_t i = r*nc; // first cell
		if (isnan(v[i])) {
			dist[i] = std::min(dist[i], above[0] + dy);
		} 
		for (size_t i=(r*nc+1); i<n; i++) { // other cells on first row
			if (isnan(v[i])) {
				size_t j = i - r*nc;
				dist[i] = std::min(std::min(std::min(dist[i], above[j] + dy), above[j-1] + dxy),  dist[i-1] + dx);
			} 
		}
		for (size_t r=nr-2; r >= 0; r--) { // other rows
			size_t i=r*nc;
			if (isnan(v[i])) {
				dist[i] = std::min(dist[i], dist[i+nc] + dy);
			}  
			for (size_t i=(r*nc+1); i<((r+1)*nc); i++) {
				if (isnan(v[i])) {
					dist[i] = std::min(std::min(std::min(dist[i], dist[i-1] + dx), dist[i+nc] + dy), dist[i+nc-1] + dxy);
				} 
			}
		}
		
		// right to left
		if (isnan(v[n-1])) { // first cell
			dist[n-1] = std::min(dist[n-1], above[nc-1] + dy);
		} 
		r = nr-1; // other cells on first row
		for (size_t i=n-2; i > (r*nc-1); i--) {
			if (isnan(v[i])) {
				size_t j = i - r*nc;
				dist[i] = std::min(std::min(std::min(dist[i], above[j] + dx), above[j+1] + dxy), dist[i+1] + dx);
			} 
		}
		for (size_t r=nr-2; r >= 0; r--) { // other rows
			size_t i=(r+1)*nc-1;
			if (isnan(v[i])) {
				dist[i] = std::min(dist[i], dist[i+nc] + dy);
			} 
			for (size_t i=(r+1)*nc-2; i>(r*nc-1); i--) {
				if (isnan(v[i])) {
					dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i+nc] + dy), dist[i+nc+1] + dxy);
				} 
			}
		}
	}		
	return dist;
}

//std::vector<double> broom_dist_geo(std::vector<double> &v, std::vector<double> &above, std::vector<double> res, std::vector<unsigned> dim, bool down) {
//
//}

/*
SpatRaster SpatRaster::broomDistance(std::string filename, std::string format, std::string datatype, bool overwrite) {

	SpatRaster out = geometry();
	bool isgeo = out.islonlat
	std::vector<double> res = resolution();
	std::vector<double> dim = {nrow, ncol};

  	out.writeStart(filename, format, datatype, overwrite);
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		broom down
		out.writeValues(v, out.bs.row[i]);
	}
	readStop();
	out.writeStop();
	
	// broom up 
	
	// combine 
	return(out);
}

SpatRaster SpatRaster::broomCostDistance(std::string filename, std::string format, std::string datatype, bool overwrite) {

*/
