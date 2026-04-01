// Copyright (c) 2018-2026  Robert J. Hijmans
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

// Dijkstra-based cost-distance.
// Replaces the iterative push-broom approach (costDistance) with a single-pass
// priority-queue algorithm that is exact and needs no convergence iterations.

#include "spatRaster.h"
#include "distance.h"
#include "file_utils.h"

#include <queue>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>


SpatRaster SpatRaster::costDistanceDijkstra(double target, double m, bool grid,
                                     SpatOptions &opt) {

	SpatRaster out = geometry(1);

	if (!hasValues()) {
		out.setError("cannot compute distance for a raster with no values");
		return out;
	}

	std::string filename = opt.get_filename();
	if (!filename.empty()) {
		if (file_exists(filename) && (!opt.get_overwrite())) {
			out.setError("output file exists. You can use 'overwrite=TRUE' to overwrite it");
			return out;
		}
	}

	SpatOptions ops(opt);
	if (nlyr() > 1) {
		std::vector<size_t> lyr = {0};
		out = subset(lyr, ops);
		out = out.costDistanceDijkstra(target, m, grid, opt);
		out.addWarning("distance computations are only done for the first input layer");
		return out;
	}

	if (source[0].srs.is_empty()) {
		out.addWarning("unknown CRS. Results can be wrong");
	}

	bool lonlat  = is_lonlat();
	bool global  = is_global_lonlat();
	int  polar   = ns_polar();
	bool npole   = (polar == 1) || (polar == 2);
	bool spole   = (polar == -1) || (polar == 2);

	double scale;
	if (!lonlat) {
		scale = source[0].srs.to_meter();
		scale = std::isnan(scale) ? 1 : scale;
		scale /= m;
	} else {
		scale = m;
	}

	size_t nr = nrow();
	size_t nc = ncol();
	size_t ncells = nr * nc;
	std::vector<double> res = resolution();

	// ── read friction surface ────────────────────────────────────────────
	std::vector<double> v;
	if (!readStart()) {
		out.setError(getError());
		return out;
	}
	readValues(v, 0, nr, 0, nc);
	readStop();

	// ── precompute per-row step distances ────────────────────────────────
	// cost mode : half-step (multiplied by sum of two frictions = avg × dist)
	// grid mode : full step (just geographic distance)
	double mult = grid ? 1.0 : 2.0;

	std::vector<double> rdx(nr), rdy(nr), rdxy(nr);

	if (!lonlat) {
		double dx  = res[0] * scale / mult;
		double dy  = res[1] * scale / mult;
		double dxy = std::sqrt(dx * dx + dy * dy);
		for (size_t r = 0; r < nr; r++) {
			rdx[r]  = dx;
			rdy[r]  = dy;
			rdxy[r] = dxy;
		}
	} else {
		for (size_t r = 0; r < nr; r++) {
			double lat = yFromRow((int64_t)r);
			rdx[r]  = distance_lonlat(0, lat, res[0], lat) / (mult * scale);
			rdy[r]  = distance_lonlat(0, 0, 0, res[1])     / (mult * scale);
			rdxy[r] = distance_lonlat(0, lat, res[0], lat - res[1])
			                                                / (mult * scale);
		}
	}

	// pole step: cost of reaching the virtual pole from the nearest row
	double pole_dy_n = 0, pole_dy_s = 0;
	if (npole) pole_dy_n = grid ? rdy[0]      : rdy[0];
	if (spole) pole_dy_s = grid ? rdy[nr - 1] : rdy[nr - 1];

	// ── initialise distances and priority queue ──────────────────────────
	const double INF = std::numeric_limits<double>::infinity();

	// extra slots for virtual pole nodes
	const size_t VPOLE_N = ncells;
	const size_t VPOLE_S = ncells + 1;
	size_t nnodes = ncells + (npole ? 1 : 0) + (spole ? 1 : 0);

	std::vector<double> dist(nnodes, INF);

	using Entry = std::pair<double, size_t>;
	std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

	// seed target cells
	bool found_target = false;
	for (size_t i = 0; i < ncells; i++) {
		if (std::isnan(v[i])) continue;
		if (v[i] == target) {
			dist[i] = 0;
			v[i] = 0;
			pq.push({0.0, i});
			found_target = true;
		} else if (!grid && v[i] < 0) {
			out.setError("negative friction values not allowed");
			return out;
		}
	}
	if (!found_target) {
		out.addWarning("no target cells found");
	}

	// ── Dijkstra ─────────────────────────────────────────────────────────
	const int drow[8] = {-1, -1, -1,  0, 0,  1,  1, 1};
	const int dcol[8] = {-1,  0,  1, -1, 1, -1,  0, 1};

	while (!pq.empty()) {
		double cost = pq.top().first;
		size_t idx  = pq.top().second;
		pq.pop();

		if (cost > dist[idx]) continue;  // stale entry

		// ── virtual pole nodes ───────────────────────────────────────
		if (idx == VPOLE_N && npole) {
			for (size_t c = 0; c < nc; c++) {
				if (std::isnan(v[c])) continue;
				double ec = grid ? pole_dy_n
				                 : v[c] * pole_dy_n;   // friction×half-step
				double nc2 = dist[idx] + ec;
				if (nc2 < dist[c]) {
					dist[c] = nc2;
					pq.push({nc2, c});
				}
			}
			continue;
		}
		if (idx == VPOLE_S && spole) {
			size_t base = (nr - 1) * nc;
			for (size_t c = 0; c < nc; c++) {
				size_t ci = base + c;
				if (std::isnan(v[ci])) continue;
				double ec = grid ? pole_dy_s
				                 : v[ci] * pole_dy_s;
				double nc2 = dist[idx] + ec;
				if (nc2 < dist[ci]) {
					dist[ci] = nc2;
					pq.push({nc2, ci});
				}
			}
			continue;
		}

		// ── regular cell ─────────────────────────────────────────────
		size_t r = idx / nc;
		size_t c = idx % nc;

		// 8-connected neighbours
		for (int d = 0; d < 8; d++) {
			int r2 = (int)r + drow[d];
			int c2 = (int)c + dcol[d];

			if (c2 < 0) {
				if (global) c2 += (int)nc; else continue;
			} else if (c2 >= (int)nc) {
				if (global) c2 -= (int)nc; else continue;
			}
			if (r2 < 0 || r2 >= (int)nr) continue;

			size_t nidx = (size_t)r2 * nc + (size_t)c2;
			if (std::isnan(v[nidx])) continue;

			bool is_diag  = (drow[d] != 0 && dcol[d] != 0);
			bool is_horiz = (drow[d] == 0);
			double step = is_diag ? rdxy[r] : (is_horiz ? rdx[r] : rdy[r]);

			double edge_cost = grid ? step : (v[idx] + v[nidx]) * step;

			double new_cost = dist[idx] + edge_cost;
			if (new_cost < dist[nidx]) {
				dist[nidx] = new_cost;
				pq.push({new_cost, nidx});
			}
		}

		// pole edges from first / last row
		if (npole && r == 0) {
			double ec = grid ? pole_dy_n : v[idx] * pole_dy_n;
			double nc2 = dist[idx] + ec;
			if (nc2 < dist[VPOLE_N]) {
				dist[VPOLE_N] = nc2;
				pq.push({nc2, VPOLE_N});
			}
		}
		if (spole && r == nr - 1) {
			double ec = grid ? pole_dy_s : v[idx] * pole_dy_s;
			double nc2 = dist[idx] + ec;
			if (nc2 < dist[VPOLE_S]) {
				dist[VPOLE_S] = nc2;
				pq.push({nc2, VPOLE_S});
			}
		}
	}

	// ── write result ─────────────────────────────────────────────────────
	dist.resize(ncells);   // drop virtual pole entries
	for (size_t i = 0; i < ncells; i++) {
		if (std::isinf(dist[i])) dist[i] = NAN;
	}

	if (!filename.empty()) {
		out.setValues(dist, ops);
		out = out.writeRaster(opt);
	} else {
		out.setValues(dist, ops);
	}
	return out;
}






inline double minCostDist(std::vector<double> &d) {
	d.erase(std::remove_if(d.begin(), d.end(),
		[](const double& v) { return std::isnan(v); }), d.end());
	std::sort(d.begin(), d.end());
	return d.empty() ? NAN : d[0];
}


inline void DxDxyCost(const double &lat, const size_t &row, double xres, double yres, const double &dir, double &dx,  double &dy, double &dxy, double distscale, const double mult=2) {
	double rlat = lat + (double)row * yres * dir;
	dx  = distance_lonlat(0, rlat, xres, rlat) / (mult * distscale);
	yres *= -dir;
	dy  = distance_lonlat(0, 0, 0, yres);
	dxy = distance_lonlat(0, rlat, xres, rlat+yres);
	dy = std::isnan(dy) ? NAN : dy / (mult * distscale);
	dxy = std::isnan(dxy) ? NAN : dxy / (mult * distscale);
}


void cost_dist(std::vector<double> &dist, std::vector<double> &dabove, std::vector<double> &v, std::vector<double> &vabove, std::vector<double> res, size_t nr, size_t nc, double lindist, bool geo, double lat, double latdir, bool global, bool npole, bool spole) {

	std::vector<double> cd;

	double dx, dy, dxy;
	if (geo) {
		DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist);
	} else {
		dx = res[0] * lindist / 2;
		dy = res[1] * lindist / 2;
		dxy = sqrt(dx*dx + dy*dy);
	}

	//top to bottom
    //left to right
	//first cell, no cell left of it
	if (!std::isnan(v[0])) {
		if (global) {
			cd = {dist[0], dabove[0] + (v[0]+vabove[0]) * dy,
				dist[nc-1] + (v[0] + v[nc-1]) * dx,
				dabove[nc-1] + dxy * (vabove[nc-1]+v[0])};
		} else {
			cd = {dist[0], dabove[0] + (v[0]+vabove[0]) * dy};
		}
		dist[0] = minCostDist(cd);
	}
	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (!std::isnan(v[i])) {
			cd = {dist[i], dabove[i]+(vabove[i]+v[i])*dy, dabove[i-1]+(vabove[i-1]+v[i])*dxy, dist[i-1]+(v[i-1]+v[i])*dx};
			dist[i] = minCostDist(cd);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}


	for (size_t r=1; r<nr; r++) { //other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist);
		size_t start=r*nc;
		if (!std::isnan(v[start])) {
			if (global) {
				cd = {dist[start-nc] + (v[start] + v[start-nc]) * dy, dist[start],
					dist[start+nc-1] + (v[start] + v[start+nc-1]) * dx,
					dist[start-1] + (v[start] + v[start-1]) * dxy};
			} else {
				cd = {dist[start-nc] + (v[start] + v[start-nc]) * dy, dist[start]};
			}
			dist[start] = minCostDist(cd);
		}
		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (!std::isnan(v[i])) {
				cd = {dist[i], dist[i-1]+(v[i]+v[i-1])*dx, dist[i-nc]+(v[i]+v[i-nc])*dy, dist[i-nc-1]+(v[i]+v[i-nc-1])*dxy};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//right to left
	// first row, no need for first (last) cell (unless is global)
	if (geo) DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist);
	if (global) {
		size_t i=(nc-1);
		cd = {dist[i],
			dist[0] + (v[0] + v[i]) * dx,
			dabove[0] + dxy * (vabove[0]+v[i])};
		dist[i] = minCostDist(cd);
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	for (int i=(nc-2); i > -1; i--) { // other cells on first row
		if (!std::isnan(v[i])) {
			cd = {dabove[i]+(vabove[i]+v[i])*dy, dabove[i+1]+(vabove[i+1]+v[i])*dxy, dist[i+1]+(v[i+1]+v[i])*dx, dist[i]};
			dist[i] = minCostDist(cd);
		}
	}

	for (size_t r=1; r<nr; r++) { // other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist);
		size_t start=(r+1)*nc-1;

		if (!std::isnan(v[start])) {
			if (global) {
				cd = { dist[start], dist[start-nc] + (v[start-nc]+v[start])* dy,
					dist[start-nc+1] + (v[start-nc+1] + v[start]) * dx,
					dist[start-(2*nc)+1] + (v[start-(2*nc)+1] + v[start]) * dxy
				};

			} else {
				cd = { dist[start], dist[start-nc] + (v[start-nc]+v[start])* dy };
			}
			dist[start] = minCostDist(cd);
		}

		size_t end=r*nc;
		start -= 1;
		for (size_t i=start; i>=end; i--) {
			if (!std::isnan(v[i])) {
				cd = { dist[i+1]+(v[i+1]+v[i])*dx, dist[i-nc+1]+(v[i]+v[i-nc+1])*dxy, dist[i-nc]+(v[i]+v[i-nc])*dy, dist[i]};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	size_t off = (nr-1) * nc;
	dabove = std::vector<double>(dist.begin()+off, dist.end());
	vabove = std::vector<double>(v.begin()+off, v.end());

}


void grid_dist(std::vector<double> &dist, std::vector<double> &dabove, std::vector<double> &v, std::vector<double> &vabove, std::vector<double> res, size_t nr, size_t nc, double lindist, bool geo, double lat, double latdir, bool global, bool npole, bool spole) {

	std::vector<double> cd;

	double dx, dy, dxy;
	if (geo) {
		DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
	} else {
		dx = res[0] * lindist;
		dy = res[1] * lindist;
		dxy = sqrt(dx*dx + dy*dy);
	}

	//top to bottom
    //left to right
	//first cell, no cell left of it
	if (!std::isnan(v[0])) {
		if (global) {
			cd = {dist[0], dabove[0] + dy,
				dist[nc-1] + dx,
				dabove[nc-1] + dxy};
		} else {
			cd = {dist[0], dabove[0] + dy};
		}
		dist[0] = minCostDist(cd);
	}
	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (!std::isnan(v[i])) {
			cd = {dist[i], dabove[i]+dy, dabove[i-1]+dxy, dist[i-1]+dx};
			dist[i] = minCostDist(cd);
		}
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}


	for (size_t r=1; r<nr; r++) { //other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
		size_t start=r*nc;
		if (!std::isnan(v[start])) {
			if (global) {
				cd = {dist[start-nc] + dy, dist[start],
					dist[start+nc-1] + dx,
					dist[start-1] + dxy};
			} else {
				cd = {dist[start-nc] + dy, dist[start]};
			}
			dist[start] = minCostDist(cd);
		}
		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (!std::isnan(v[i])) {
				cd = {dist[i], dist[i-1]+dx, dist[i-nc]+dy, dist[i-nc-1]+dxy};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	//right to left
	// first row, no need for first (last) cell (unless is global)
	if (geo) DxDxyCost(lat, 0, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
	if (global) {
		size_t i=(nc-1);
		cd = {dist[i],
			dist[0] + dx,
			dabove[0] + dxy};
		dist[i] = minCostDist(cd);
	}
	if (npole) {
		double minp = *std::min_element(dist.begin(), dist.begin()+nc);
		minp += dy;
		for (size_t i=0; i<nc; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	for (int i=(nc-2); i > -1; i--) { // other cells on first row
		if (!std::isnan(v[i])) {
			cd = {dabove[i]+dy, dabove[i+1]+dxy, dist[i+1]+dx, dist[i]};
			dist[i] = minCostDist(cd);
		}
	}

	for (size_t r=1; r<nr; r++) { // other rows
		if (geo) DxDxyCost(lat, r, res[0], res[1], latdir, dx, dy, dxy, lindist, 1);
		size_t start=(r+1)*nc-1;

		if (!std::isnan(v[start])) {
			if (global) {
				cd = { dist[start], dist[start-nc] + dy,
					dist[start-nc+1] + dx,
					dist[start-(2*nc)+1] + dxy
				};

			} else {
				cd = { dist[start], dist[start-nc] + dy };
			}
			dist[start] = minCostDist(cd);
		}

		size_t end=r*nc;
		start -= 1;
		for (size_t i=start; i>=end; i--) {
			if (!std::isnan(v[i])) {
				cd = { dist[i+1]+dx, dist[i-nc+1]+dxy, dist[i-nc]+dy, dist[i]};
				dist[i] = minCostDist(cd);
			}
		}
	}
	if (spole) {
		double minp = *std::min_element(dist.end()-nc, dist.end());
		minp += dy;
		size_t ds = dist.size();
		for (size_t i=ds-nc; i<ds; i++) {
			dist[i] = std::min(dist[i], minp);
		}
	}

	size_t off = (nr-1) * nc;
	dabove = std::vector<double>(dist.begin()+off, dist.end());
	vabove = std::vector<double>(v.begin()+off, v.end());
}

void block_is_same(bool& same, std::vector<double>& x,  std::vector<double>& y) {
	if (!same) return;
	for (size_t i=0; i<x.size(); i++) {
		if (!std::isnan(x[i]) && (x[i] != y[i])) {
			same = false;
			break;
		}
	}
}






SpatRaster SpatRaster::costDistanceRun(SpatRaster &old, bool &converged, double target, double m, bool lonlat, bool global, bool npole, bool spole, bool grid, SpatOptions &opt) {

	std::vector<double> res = resolution();

	SpatRaster first = geometry();
	SpatRaster second = first;
    std::vector<double> d, v, vv;
	if (!readStart()) {
		first.setError(getError());
		return(first);
	}
	opt.progressbar = false;
 	if (!first.writeStart(opt, filenames())) { return first; }

	size_t nc = ncol();
	std::vector<double> dabove(nc, NAN);
	std::vector<double> vabove(nc, 0);
	double lat = 0;
	if (old.hasValues()) {
		if (!old.readStart()) {
			first.setError(getError());
			return(first);
		}
		if (!first.writeStart(opt, filenames())) {
			readStop();
			old.readStop();
			return first;
		}

		for (size_t i = 0; i < first.bs.n; i++) {
			readBlock(v, first.bs, i);
			if (lonlat) {
				lat = yFromRow(first.bs.row[i]);
			}
			bool np = (i==0) && npole;
			bool sp = (i==first.bs.n-1) && spole;
			if (target != 0) {
				for (size_t j=0; j<v.size(); j++) {
					if (v[j] == target) {
						v[j] = 0;
					}
				}
			}
			old.readBlock(d, first.bs, i);
			if (grid) {
				grid_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			} else {
				cost_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			}
			if (!first.writeValues(d, first.bs.row[i], first.bs.nrows[i])) return first;
		}
	} else {
		converged = false;
		for (size_t i = 0; i < first.bs.n; i++) {
			if (lonlat) {
				lat = yFromRow(first.bs.row[i]);
			}
			bool np = (i==0) && npole;
			bool sp = (i==first.bs.n-1) && spole;
			readBlock(v, first.bs, i);
			d.clear();
			d.resize(v.size(), NAN);
			for (size_t j = 0; j < v.size(); j++) {
				if (v[j] == target) {
					v[j] = 0;
					d[j] = 0;
				} else if ((!grid) && (v[j] < 0)) {
					readStop();
					first.writeStop();
					first.setError("negative friction values not allowed");
					return first;
				}
			}
			if (grid) {
				grid_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			} else {
				cost_dist(d, dabove, v, vabove, res, first.bs.nrows[i], nc, m, lonlat, lat, -1, global, np, sp);
			}
			if (!first.writeValuesRect(d, first.bs.row[i], first.bs.nrows[i], 0, nc)) return first;
		}
	}
	first.writeStop();

	if (!first.readStart()) {
		return(first);
	}

	dabove = std::vector<double>(nc, NAN);
	vabove = std::vector<double>(nc, 0);
  	if (!second.writeStart(opt, filenames())) {
		readStop();
		first.readStop();
		return second;
	}
	for (int i = second.bs.n; i>0; i--) {
		if (lonlat) {
			lat = yFromRow(second.bs.row[i-1] + second.bs.nrows[i-1] - 1);
		}
		bool sp = (i==1) && spole; //! reverse order
		bool np = (i==(int)second.bs.n) && npole;
        readBlock(v, second.bs, i-1);
		if (target != 0) {
			for (size_t j=0; j<v.size(); j++) {
				if (v[j] == target) {
					v[j] = 0;
				}
			}
		}
		first.readBlock(d, second.bs, i-1);
		std::reverse(v.begin(), v.end());
		std::reverse(d.begin(), d.end());
		if (grid) {
			grid_dist(d, dabove, v, vabove, res, second.bs.nrows[i-1], nc, m, lonlat, lat, 1, global, np, sp);
		} else {
			cost_dist(d, dabove, v, vabove, res, second.bs.nrows[i-1], nc, m, lonlat, lat, 1, global, np, sp);
		}
		std::reverse(d.begin(), d.end());
		if (converged) {
			old.readBlock(v, second.bs, i-1);
			block_is_same(converged, d, v);
		}
		if (!second.writeValuesRect(d, second.bs.row[i-1], second.bs.nrows[i-1], 0, nc)) return second;
	}
	second.writeStop();
	first.readStop();
	if (old.hasValues()) {
		old.readStop();
	}
	readStop();
	return(second);
}

SpatRaster SpatRaster::costDistance(double target, double m, size_t maxiter, bool grid, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("cannot compute distance for a raster with no values");
		return out;
	}

	std::string filename = opt.get_filename();
	if (!filename.empty()) {
		if (file_exists(filename) && (!opt.get_overwrite())) {
			out.setError("output file exists. You can use 'overwrite=TRUE' to overwrite it");
			return(out);
		}
	}

	SpatOptions ops(opt);
	if (nlyr() > 1) {
		std::vector<size_t> lyr = {0};
		out = subset(lyr, ops);
		out = out.costDistance(target, m, maxiter, grid, opt);
		out.addWarning("distance computations are only done for the first input layer");
		return out;
	}

	if (source[0].srs.is_empty()) {
		out.addWarning("unknown CRS. Results can be wrong");
	}

	bool lonlat = is_lonlat();
	bool global = is_global_lonlat();
	int polar = ns_polar();
	bool npole = (polar == 1) || (polar == 2);
	bool spole = (polar == -1) || (polar == 2);

	double scale;
	if (!lonlat) {
		scale = source[0].srs.to_meter();
		scale = std::isnan(scale) ? 1 : scale;
		scale /= m;
	} else {
		scale = m;
	}

	std::vector<double> res = resolution();

	// if the raster fits in memory, use the exact single-pass Dijkstra
	SpatOptions memops(opt);
	memops.ncopies = 4;
	BlockSize membs = out.getBlockSize(memops);
	if (membs.nrows[0] >= nrow()) {
		return costDistanceDijkstra(target, m, grid, opt);
	}

	size_t i = 0;
	bool converged=false;
	while (i < maxiter) {
		out = costDistanceRun(out, converged, target, scale, lonlat, global, npole, spole, grid, ops);
		if (out.hasError()) return out;
		if (converged) break;
		converged = true;
		i++;
	}
	if (!filename.empty()) {
		out = out.writeRaster(opt);
	}
	if (i == maxiter) {
		out.addWarning("distance algorithm did not converge");
	}
	return(out);
}






