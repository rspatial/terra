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

#include "spatRaster.h"
#include "math_utils.h"

void patches_replace(std::vector<double> &v, size_t n, std::vector<double>& d, size_t cstart, std::vector<std::vector<size_t>>& rcl, size_t &ncps) {

	d.erase(std::remove_if(d.begin(), d.end(),
		[](const double& v) { return std::isnan(v); }), d.end());
	std::sort(d.begin(), d.end());
	d.erase(std::unique(d.begin(), d.end()), d.end());

	size_t nd = d.size();

	if (nd == 0) {
		v[n] = ncps;
		ncps++;
		return;
	} else if (nd == 1) {
		v[n] = d[0];
		return;
	}
	v[n] = d[0];
	for (size_t i=0; i<n; i++) {
		for (size_t j=1; j<nd; j++) {
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
	} else if (d[nd-1] == (ncps-1)) {  
		ncps--;
	}
}


void broom_patches(const std::vector<double> &vals, std::vector<double> &patches, std::vector<double>& above_p, std::vector<double>& above_v, const size_t &dirs, size_t &ncps, const size_t &nr, const size_t &nc, std::vector<std::vector<size_t>> &rcl, bool is_global) {

	size_t nstart = ncps;
	bool d4 = dirs == 4;
	size_t stopnc = nc-1;
	std::vector<double> d;
	
	auto check_neighbor = [&](double val, double neighbor_val, double neighbor_patch, std::vector<double>& candidates) {
		if (!std::isnan(neighbor_val) && !std::isnan(neighbor_patch)) {
			if (is_equal(val, neighbor_val)) {
				candidates.push_back(neighbor_patch);
			}
		}
	};

	//first row
	if ( !std::isnan(vals[0]) ) {
		d.clear();
		if (d4) {
			check_neighbor(vals[0], above_v[0], above_p[0], d);
		} else if (is_global) { 
			check_neighbor(vals[0], above_v[0], above_p[0], d);
			check_neighbor(vals[0], above_v[1], above_p[1], d);
			check_neighbor(vals[0], above_v[stopnc], above_p[stopnc], d);
		} else { 
			check_neighbor(vals[0], above_v[0], above_p[0], d);
			check_neighbor(vals[0], above_v[1], above_p[1], d);
		}
		patches_replace(patches, 0, d, nstart, rcl, ncps);
	}

	for (size_t i=1; i<stopnc; i++) {
		if (!std::isnan(vals[i])) {
			d.clear();
			check_neighbor(vals[i], above_v[i], above_p[i], d);
			check_neighbor(vals[i], vals[i-1], patches[i-1], d);
			if (!d4) {
				check_neighbor(vals[i], above_v[i-1], above_p[i-1], d);
				check_neighbor(vals[i], above_v[i+1], above_p[i+1], d);
			}
			patches_replace(patches, i, d, nstart, rcl, ncps);
		}
	}
	
	size_t i = stopnc;
	if (!std::isnan(vals[i])) {
		d.clear();
		check_neighbor(vals[i], above_v[i], above_p[i], d);
		check_neighbor(vals[i], vals[i-1], patches[i-1], d);
		if (is_global) {
			check_neighbor(vals[i], vals[0], patches[0], d);
			if (!d4) {
				check_neighbor(vals[i], above_v[i-1], above_p[i-1], d);
				check_neighbor(vals[i], above_v[0], above_p[0], d);
			}
		} else if (!d4) {
			check_neighbor(vals[i], above_v[i-1], above_p[i-1], d);
		}
		patches_replace(patches, i, d, nstart, rcl, ncps);
	}

	for (size_t r=1; r<nr; r++) {
		size_t start = r*nc;
		size_t i=start;
		if (!std::isnan(vals[i])) {
			d.clear();
			check_neighbor(vals[i], vals[i-nc], patches[i-nc], d);
			if (is_global) {
				if (!d4) {
					check_neighbor(vals[i], vals[i-nc+1], patches[i-nc+1], d);
					check_neighbor(vals[i], vals[i-1], patches[i-1], d);
				}
			} else if (!d4) {
				check_neighbor(vals[i], vals[i-nc+1], patches[i-nc+1], d);
			}
			patches_replace(patches, i, d, nstart, rcl, ncps);
		}

		size_t stop = start + stopnc;
		for (size_t i=(start+1); i<stop; i++) {
			if (!std::isnan(vals[i])) {
				d.clear();
				check_neighbor(vals[i], vals[i-1], patches[i-1], d);
				check_neighbor(vals[i], vals[i-nc], patches[i-nc], d);
				if (!d4) {
					check_neighbor(vals[i], vals[i-nc-1], patches[i-nc-1], d);
					check_neighbor(vals[i], vals[i-nc+1], patches[i-nc+1], d);
				}
				patches_replace(patches, i, d, nstart, rcl, ncps);
			}
		}

		i = stop;
		if (!std::isnan(vals[i])) {
			d.clear();
			check_neighbor(vals[i], vals[i-1], patches[i-1], d);
			check_neighbor(vals[i], vals[i-nc], patches[i-nc], d);
			if (is_global) {
				check_neighbor(vals[i], vals[start], patches[start], d);
				if (!d4) {
					check_neighbor(vals[i], vals[i-nc-1], patches[i-nc-1], d);
					check_neighbor(vals[i], vals[start-nc], patches[start-nc], d);
				}
			} else if (!d4) {
				check_neighbor(vals[i], vals[i-nc-1], patches[i-nc-1], d);
			}
			patches_replace(patches, i, d, nstart, rcl, ncps);
		}
	}
	size_t off = (nr-1) * nc;
	std::vector<double> last_row_p(patches.begin()+off, patches.end());
	std::vector<double> last_row_v(vals.begin()+off, vals.end());
	above_p = last_row_p;
	above_v = last_row_v;
}

std::vector<std::vector<double>> patches_getRCL(std::vector<std::vector<size_t>> rcl, size_t n) {
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

SpatRaster SpatRaster::patches(size_t directions, SpatOptions &opt) {

	SpatRaster out = geometry(1, false);

	if (nlyr() > 1) {
		SpatOptions ops(opt);
		std::vector<std::string> nms = getNames();
		if (ops.names.size() == nms.size()) {
			nms = opt.names;
		}
		for (size_t i=0; i<nlyr(); i++) {
			std::vector<size_t> lyr = {i};
			ops.names = {nms[i]};
			SpatRaster x = subset(lyr, ops);
			x = x.patches(directions, ops); 
			out.addSource(x, false, ops);
		}
		if (!opt.get_filename().empty()) {
			out = out.writeRaster(opt);
		}
		return out;
	}

	if (!hasValues()) {
		out.setError("cannot compute patches for a raster with no values");
		return out;
	}
	if (!((directions == 4) || (directions == 8))) {
		out.setError("directions should be 4 or 8");
		return out;		
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
	opt.set_filenames({""}); 
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

    size_t nc = ncol();
	size_t ncps = 1;
	std::vector<double> above_p(nc, NAN);
	std::vector<double> above_v(nc, NAN);
	std::vector<std::vector<size_t>> rcl(2);
	std::vector<double> v;
	
	bool is_global = is_global_lonlat();

	for (size_t i = 0; i < out.bs.n; i++) {
		readBlock(v, out.bs, i);
		std::vector<double> p(v.size(), NAN); 
		broom_patches(v, p, above_p, above_v, directions, ncps, out.bs.nrows[i], nc, rcl, is_global);
		if (!out.writeBlock(p, i)) return out;
	}
	
	readStop();
	out.writeStop();
	
	std::string filename = opt.get_filename();
	opt.set_filenames({filename});
	
	if (!rcl[0].empty()) {
		std::vector<std::vector<double>> rc = patches_getRCL(rcl, ncps);
		out = out.reclassify(rc, 3, true, false, 0.0, false, false, false, opt);
	} else if (!filename.empty()) {
		out = out.writeRaster(opt);
	}

	return(out);
}
