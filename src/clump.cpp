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

#include "spatRaster.h"
#include "file_utils.h"


std::vector<std::vector<double>> getRCL(std::vector<std::vector<size_t>> rcl, size_t n) {
	std::vector<std::vector<size_t>> rcl2(rcl[0].size());
	for (size_t i=0; i<rcl[0].size(); i++) {
		rcl2[i].push_back(rcl[0][i]);
		rcl2[i].push_back(rcl[1][i]);
	}
    std::sort(rcl2.begin(), rcl2.end());
    rcl2.erase(std::unique(rcl2.begin(), rcl2.end()), rcl2.end());
	std::vector<std::vector<double>> out(2);
	for (size_t i=0; i<rcl2.size(); i++) {
		out[0].push_back(rcl2[i][1]);
		out[1].push_back(rcl2[i][0]);
	}
	// from - to 
	// 3 - 1
	// 4 - 3
    // becomes
    // 3 - 1
    // 4 - 1
	for (size_t i=1; i<out[0].size(); i++) {
		for (size_t j=0; j<i; j++) {
			if (out[0][i] == out[1][j]) {
				out[1][j] = out[0][i];
			}
		}
	}

	std::vector<double> lost = out[0];
	lost.push_back(n);
	size_t sub = 0;
	for (size_t i=0; i<lost.size(); i++) {
		sub++;
		for (size_t j=lost[i]+1; j<lost[i+1]; j++) {
			out[0].push_back(j);
			out[1].push_back(j-sub);
		}
	}
	return out;
}


void replace(std::vector<double> &v, size_t n, const std::vector<double>& d, size_t cstart, std::vector<std::vector<size_t>>& rcl) {
	for (size_t i=0; i<n; i++) {
		for (size_t j=1; j<d.size(); j++) {
			if (v[i] == d[j]) {
				v[i] = d[0];
			}
		}
	}
	if (d[0] < cstart) {
		for (size_t j=1; j<d.size(); j++) {
			rcl[0].push_back(d[0]);
			rcl[1].push_back(d[j]);
		}
	}
}


void test(std::vector<double> &d) {
	d.erase(std::remove_if(d.begin(), d.end(),
		[](const double& v) { return std::isnan(v); }), d.end());
	std::sort(d.begin(), d.end());
	d.erase(std::unique(d.begin(), d.end()), d.end());
}

void broom_clumps(std::vector<double> &v, std::vector<double>& above, const size_t &dirs, size_t &ncps, const size_t &nr, const size_t &nc, std::vector<std::vector<size_t>> &rcl) {

	size_t nstart = ncps;

	bool d4 = dirs == 4;

	if ( !std::isnan(v[0]) ) { //first cell, no cell left of it
		if (std::isnan(above[0])) {
			v[0] = ncps;
			ncps++;
		} else {
			v[0] = above[0];
		}
	}

	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (!std::isnan(v[i])) {
			std::vector<double> d;
			if (d4) {
				d = {above[i], v[i-1]} ;
			} else {
				d = {above[i], above[i-1], v[i-1]} ;
			}
			test(d);
			if (d.size() > 0) {
				v[i] = d[0];
				if (d.size() > 1) {
					replace(v, i, d, nstart, rcl);
				}
			} else {
				v[i] = ncps;
				ncps++;
			}
		}
	}


	for (size_t r=1; r<nr; r++) { //other rows
		size_t i=r*nc;
		if (!std::isnan(v[i])) { // first cell
			if (std::isnan(v[i-nc])) {
				v[i] = ncps;
				ncps++;
			} else {
				v[i] = v[i-nc];
			}
		}
		for (size_t i=r*nc+1; i<((r+1)*nc); i++) { // other cells
			if (!std::isnan(v[i])) {
				std::vector<double> d;
				if (d4) {
					d = {v[i-nc], v[i-1]} ;
				} else {
					d = {v[i-nc], v[i-nc-1], v[i-1]} ;
				}
				test(d);
				if (d.size() > 0) {
					v[i] = d[0];
					if (d.size() > 1) {
						replace(v, i, d, nstart, rcl);
					}
				} else {
					v[i] = ncps;
					ncps++;
				}
			}
		}
	}
	size_t off = (nr-1) * nc;
	above = std::vector<double>(v.begin()+off, v.end());
}



SpatRaster SpatRaster::clumps(int directions, bool zeroAsNA, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (nlyr() > 1) {
		SpatOptions ops(opt);
		std::string filename = opt.get_filename();
		ops.set_filenames({""});
		for (size_t i=0; i<nlyr(); i++) {
			std::vector<unsigned> lyr = {(unsigned)i};
			SpatRaster x = subset(lyr, ops);
			x = x.clumps(directions, zeroAsNA, ops);
			out.addSource(x);
		}
		if (filename != "") {
			out = out.writeRaster(opt);
		}
		return out;
	}

	if (!(directions == 4 || directions == 8)) {
		out.setError("directions must be 4 or 8");
		return out;
	}
	if (!hasValues()) {
		out.setError("cannot compute clumps for a raster with no values");
		return out;
	}

	std::vector<size_t> dim = {nrow(), ncol()};

	std::string tempfile = "";
    std::vector<double> d, v, vv;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	std::string filename = opt.get_filename();
	if (filename != "") {
		bool overwrite = opt.get_overwrite();
		std::string errmsg;
		if (!can_write(filename, overwrite, errmsg)) {
			out.setError(errmsg + " (" + filename +")");
			return(out);
		}
	}

	opt.set_filenames({""});
 	if (!out.writeStart(opt)) { return out; }
	size_t nc = ncol();
	size_t ncps = 1;
	std::vector<double> above(nc, NAN);
	std::vector<std::vector<size_t>> rcl(2);
	for (size_t i = 0; i < out.bs.n; i++) {
        v = readBlock(out.bs, i);
		if (zeroAsNA) {
			std::replace(v.begin(), v.end(), 0.0, (double)NAN);
		}
        broom_clumps(v, above, directions, ncps, out.bs.nrows[i], nc, rcl);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, nc)) return out;
	}
	out.writeStop();
	readStop();

	opt.set_filenames({filename});
	if (rcl[0].size() > 0) {
		//for (size_t i=0; i<rcl[0].size(); i++) {
		//	Rcpp::Rcout << rcl[0][i] << " - " << rcl[1][i] << std::endl;
		//}
		//Rcpp::Rcout << std::endl;
		std::vector<std::vector<double>> rc = getRCL(rcl, ncps);
		//for (size_t i=0; i<rc[0].size(); i++) {
		//	Rcpp::Rcout << rc[0][i] << " - " << rc[1][i] << std::endl;
		//}
		out = out.reclassify(rc, 3, true, false, opt);
	} else if (filename != "") {
		out = out.writeRaster(opt);
	}
	return out;
}

