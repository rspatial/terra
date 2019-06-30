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

// Robert Hijmans, October 2011
// This is an implementation of J. Ronald Eastman's pushbroom algorithm


#include "spatRaster.h"
#include <limits>
#include <cmath>

std::vector<double> broom_dist_planar(std::vector<double> &v, std::vector<double> &above, std::vector<double> res, std::vector<unsigned> dim) {

	double dx = res[0];
	double dy = res[1];
	double dxy = sqrt(dx * dx + dy *dy);

	size_t n = v.size();
	unsigned nr = n / dim[0]; // must get entire rows
	unsigned nc = dim[1];

	std::vector<double> dist(n, 0);

	//top to bottom
    //left to right

	if ( std::isnan(v[0]) ) { //first cell, no cell left of it
		dist[0] = above[0] + dy;
	}
	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(above[i] + dy, above[i-1] + dxy), dist[i-1] + dx);
		}
	}
	for (size_t r=1; r<nr; r++) { //other rows
		size_t i=r*nc;
		if (std::isnan(v[i])) {
			dist[i] = dist[i-nc] + dy;
		}
		for (size_t i=r*nc+1; i<((r+1)*nc); i++) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(dist[i-1] + dx, dist[i-nc] + dy), dist[i-nc-1] + dxy);
			}
		}
	}
		//right to left
	if ( std::isnan(v[nc-1])) { //first cell
		dist[nc-1] = std::min(dist[nc-1], above[nc-1] + dy);
	}
	for (size_t i=(nc-1); i > 0; i--) { // other cells on first row
		if (std::isnan(v[i-1])) {
			dist[i] = std::min(std::min(std::min(dist[i-1], above[i-1] + dy), above[i] + dxy), dist[i] + dx);
		}
	}
	for (size_t r=1; r<nr; r++) { // other rows
		size_t i=(r+1)*nc-1;
		if (std::isnan(v[i])) {
			dist[i] = std::min(dist[i], dist[i-nc] + dy);
		}
		for (size_t i=(r+1)*nc-2; i>(r*nc-1); i--) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy), dist[i-nc+1] + dxy);
			}
		}
	}
	return dist;
}

/*


//std::vector<double> broom_dist_geo(std::vector<double> &v, std::vector<double> &above, std::vector<double> res, std::vector<unsigned> dim, bool down) {
//
//}
*/


SpatRaster SpatRaster::gridDistance(SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) {
		out.setError("cannot compute distance for a raster with no values");
		return out;
	}

	//bool isgeo = out.islonlat

	std::vector<double> res = resolution();
	std::vector<unsigned> dim = {nrow(), ncol()};

	SpatRaster first = out.geometry();

	std::string tempfile = "";
	std::vector<double> above(ncol(), std::numeric_limits<double>::infinity());
    std::vector<double> d, v, vv;
	readStart();
	std::string filename = opt.get_filename();
	opt.set_filename("");
 	if (!first.writeStart(opt)) { return first; }

	for (size_t i = 0; i < first.bs.n; i++) {
        v = readBlock(first.bs, i);
        d = broom_dist_planar(v, above, res, dim);
		if (!first.writeValues(d, first.bs.row[i], first.bs.nrows[i], 0, ncol())) return first;
	}
	first.writeStop();
	
  	first.readStart();
	opt.set_filename(filename);
	
  	if (!out.writeStart(opt)) { return out; }
	for (size_t i = out.bs.n; i>0; i--) {
        v = readBlock(out.bs, i-1);
		std::reverse(v.begin(), v.end());
        d = broom_dist_planar(v, above, res, dim);
		vv = first.readBlock(first.bs, i-1);
	    std::transform (d.rbegin(), d.rend(), vv.begin(), vv.begin(), [](double a, double b) {return std::min(a,b);});
		if (!out.writeValues(vv, out.bs.row[i-1], out.bs.nrows[i-1], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	first.readStop();
	return(out);
}


//SpatRaster SpatRaster::broomCostDistance(std::string filename, std::string format, std::string datatype, bool overwrite) {
