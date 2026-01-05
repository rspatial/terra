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

struct PatchUnionFind {
	std::vector<double> parents;
	size_t start_id;

	PatchUnionFind(size_t start) : start_id(start) {}

	void make_set(double id) {
		parents.push_back(id); 
	}

	double find(double id) {
		if (std::isnan(id)) return id;
		if (id < start_id) return id; 
		
		size_t idx = (size_t)(id - start_id);
		if (parents[idx] == id) return id;
		
		parents[idx] = find(parents[idx]); 
		return parents[idx];
	}

	void unite(double id1, double id2, std::vector<std::vector<size_t>>& rcl) {
		double root1 = find(id1);
		double root2 = find(id2);
		if (root1 == root2) return;

		bool g1 = root1 < start_id;
		bool g2 = root2 < start_id;

		if (g1 && g2) {
			rcl[0].push_back((size_t)root1);
			rcl[1].push_back((size_t)root2);
			return;
		}

		if (g1) {
			parents[(size_t)(root2 - start_id)] = root1;
		} else if (g2) {
			parents[(size_t)(root1 - start_id)] = root2;
		} else {
			if (root1 < root2)
				parents[(size_t)(root2 - start_id)] = root1;
			else
				parents[(size_t)(root1 - start_id)] = root2;
		}
	}
};

void broom_patches(const std::vector<double> &vals, std::vector<double> &patches, std::vector<double>& above_p, std::vector<double>& above_v, const size_t &dirs, size_t &ncps, const size_t &nr, const size_t &nc, std::vector<std::vector<size_t>> &rcl, bool is_global) {

	PatchUnionFind uf(ncps);
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

	auto assign_patch = [&](size_t idx, std::vector<double>& candidates) {
		if (candidates.empty()) {
			patches[idx] = ncps;
			uf.make_set(ncps);
			ncps++;
		} else {
			double chosen = candidates[0];
			patches[idx] = chosen;
			for (size_t k = 1; k < candidates.size(); ++k) {
				uf.unite(chosen, candidates[k], rcl);
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
		assign_patch(0, d);
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
			assign_patch(i, d);
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
		assign_patch(i, d);
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
			assign_patch(i, d);
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
				assign_patch(i, d);
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
			assign_patch(i, d);
		}
	}

	for (size_t k = 0; k < patches.size(); ++k) {
		if (!std::isnan(patches[k])) {
			patches[k] = uf.find(patches[k]);
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