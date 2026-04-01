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


#include "spatRasterMultiple.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <numeric>


// helpers shared by blend and mosaic2

struct RastInfo {
	SpatExtent ext;
	size_t nrow, ncol;
	double xres, yres;
};

struct BlockRead {
	std::vector<double> vals;
	size_t nrow, ncol;
	double xmin, ymax;
	double xres, yres;
	bool active;
};

static bool setup_output(SpatRasterCollection &src, SpatRaster &out, SpatOptions &opt) {
	size_t n = src.ds.size();
	SpatExtent ue = src.ds[0].getExtent();
	for (size_t i = 1; i < n; i++) {
		if (!src.ds[0].compare_geom(src.ds[i], false, false, opt.get_tolerance(), false, false, false, true)) {
			out.setError(src.ds[0].getError());
			return false;
		}
		ue.unite(src.ds[i].getExtent());
	}
	out = src.ds[0].geometry(1, true);
	out.setExtent(ue, true, true, "");
	return true;
}

static std::vector<RastInfo> collect_info(SpatRasterCollection &src) {
	size_t n = src.ds.size();
	std::vector<RastInfo> ri(n);
	for (size_t i = 0; i < n; i++) {
		ri[i].ext  = src.ds[i].getExtent();
		ri[i].nrow = src.ds[i].nrow();
		ri[i].ncol = src.ds[i].ncol();
		ri[i].xres = src.ds[i].xres();
		ri[i].yres = src.ds[i].yres();
	}
	return ri;
}

static bool start_readers(SpatRasterCollection &src, SpatRaster &out) {
	for (size_t i = 0; i < src.ds.size(); i++) {
		if (!src.ds[i].readStart()) {
			out.setError(src.ds[i].getError());
			for (size_t j = 0; j < i; j++) src.ds[j].readStop();
			return false;
		}
	}
	return true;
}

static void stop_readers(SpatRasterCollection &src) {
	for (size_t i = 0; i < src.ds.size(); i++) {
		src.ds[i].readStop();
	}
}

static void read_block(SpatRasterCollection &src, std::vector<RastInfo> &ri, 
		std::vector<BlockRead> &br, double bymax, double bymin, size_t n) {
	
	for (size_t i = 0; i < n; i++) {
		br[i].active = false;
		SpatExtent& ei = ri[i].ext;

		if (bymin >= ei.ymax || bymax <= ei.ymin) continue;

		double y_top = std::min(bymax, ei.ymax);
		double y_bot = std::max(bymin, ei.ymin);

		size_t ir_first = (size_t)std::floor((ei.ymax - y_top) / ri[i].yres);
		size_t ir_last  = (size_t)std::ceil ((ei.ymax - y_bot) / ri[i].yres);
		if (ir_first >= ri[i].nrow) ir_first = ri[i].nrow - 1;
		if (ir_last  >  ri[i].nrow) ir_last  = ri[i].nrow;
		size_t ir_nrow = ir_last - ir_first;
		if (ir_nrow == 0) continue;

		src.ds[i].readValues(br[i].vals, ir_first, ir_nrow, 0, ri[i].ncol);

		br[i].nrow = ir_nrow;
		br[i].ncol = ri[i].ncol;
		br[i].xmin = ei.xmin;
		br[i].ymax = ei.ymax - ir_first * ri[i].yres;
		br[i].xres = ri[i].xres;
		br[i].yres = ri[i].yres;
		br[i].active = true;
	}
}

// Look up the value for a cell at (x,y) in block-read i.
// Returns false if the cell is outside the raster or the value is NA.
static bool lookup(std::vector<BlockRead> &br, std::vector<RastInfo> &ri,
                   size_t i, double x, double y, double &val) {
	if (!br[i].active) return false;
	SpatExtent& ei = ri[i].ext;
	if (x <= ei.xmin || x >= ei.xmax ||
	    y <= ei.ymin || y >= ei.ymax) return false;

	size_t ic = (size_t)((x - br[i].xmin) / br[i].xres);
	size_t ir = (size_t)((br[i].ymax - y)  / br[i].yres);
	if (ic >= br[i].ncol) ic = br[i].ncol - 1;
	if (ir >= br[i].nrow) ir = br[i].nrow - 1;

	val = br[i].vals[ir * br[i].ncol + ic];
	return !std::isnan(val);
}


SpatRaster SpatRasterCollection::blend(SpatOptions &opt) {

	size_t n = ds.size();
	SpatRaster out;

	if (n == 0) {
		out.setError("collection is empty");
		return out;
	}
	if (n == 1) {
		if (opt.get_filename() != "") {
			out = ds[0].writeRaster(opt);
		} else {
			out = ds[0].deepCopy();
		}
		return out;
	}

	if (!setup_output(*this, out, opt)) return out;

	double oxres = out.xres();
	double oyres = out.yres();
	SpatExtent oe = out.getExtent();
	size_t onc = out.ncol();

	std::vector<RastInfo> ri = collect_info(*this);

	if (!start_readers(*this, out)) return out;

	opt.ncopies = std::max(opt.ncopies, n + 1);
	if (!out.writeStart(opt, filenames())) {
		stop_readers(*this);
		return out;
	}

	for (size_t b = 0; b < out.bs.n; b++) {

		size_t brow   = out.bs.row[b];
		size_t bnrow  = out.bs.nrows[b];
		size_t bncell = bnrow * onc;

		double bymax = oe.ymax - brow * oyres;
		double bymin = oe.ymax - (brow + bnrow) * oyres;

		std::vector<BlockRead> br(n);
		read_block(*this, ri, br, bymax, bymin, n);

		std::vector<double> result(bncell,
		                           std::numeric_limits<double>::quiet_NaN());

		for (size_t r = 0; r < bnrow; r++) {
			double y = oe.ymax - (brow + r + 0.5) * oyres;

			for (size_t c = 0; c < onc; c++) {
				double x = oe.xmin + (c + 0.5) * oxres;

				double wsum = 0.0;
				double vsum = 0.0;

				for (size_t i = 0; i < n; i++) {
					double v;
					if (!lookup(br, ri, i, x, y, v)) continue;

					SpatExtent& ei = ri[i].ext;
					double d = std::min({x - ei.xmin, ei.xmax - x,
					                     y - ei.ymin, ei.ymax - y});
					if (d <= 0.0) d = 1e-10;

					wsum += d;
					vsum += v * d;
				}

				if (wsum > 0.0) {
					result[r * onc + c] = vsum / wsum;
				}
			}
		}

		if (!out.writeValues(result, brow, bnrow)) {
			stop_readers(*this);
			return out;
		}
	}

	out.writeStop();
	stop_readers(*this);
	return out;
}



SpatRaster SpatRasterCollection::mosaic(std::string fun, SpatOptions &opt) {

	SpatRaster out;

	std::vector<std::string> valid {"first", "last", "sum", "mean",
	                                "median", "min", "max"};
	if (std::find(valid.begin(), valid.end(), fun) == valid.end()) {
		out.setError(fun + "' is not a valid function name");
		return out;
	}

	size_t n = ds.size();

	if (n == 0) {
		out.setError("collection is empty");
		return out;
	}
	if (n == 1) {
		if (opt.get_filename() != "") {
			out = ds[0].writeRaster(opt);
		} else {
			out = ds[0].deepCopy();
		}
		return out;
	}

	if (fun == "first") return merge(true,  true, 1, "near", opt);
	if (fun == "last")  return merge(false, true, 1, "near", opt);

	if (!setup_output(*this, out, opt)) return out;

	double oxres = out.xres();
	double oyres = out.yres();
	SpatExtent oe = out.getExtent();
	size_t onc = out.ncol();

	std::vector<RastInfo> ri = collect_info(*this);

	if (!start_readers(*this, out)) return out;

	opt.ncopies = std::max(opt.ncopies, n + 1);
	if (!out.writeStart(opt, filenames())) {
		stop_readers(*this);
		return out;
	}

	for (size_t b = 0; b < out.bs.n; b++) {

		size_t brow   = out.bs.row[b];
		size_t bnrow  = out.bs.nrows[b];
		size_t bncell = bnrow * onc;

		double bymax = oe.ymax - brow * oyres;
		double bymin = oe.ymax - (brow + bnrow) * oyres;

		std::vector<BlockRead> br(n);
		read_block(*this, ri, br, bymax, bymin, n);

		std::vector<double> result(bncell, std::numeric_limits<double>::quiet_NaN());

		// small per-cell buffer to collect contributing values
		std::vector<double> cell_vals;
		cell_vals.reserve(n);

		for (size_t r = 0; r < bnrow; r++) {
			double y = oe.ymax - (brow + r + 0.5) * oyres;

			for (size_t c = 0; c < onc; c++) {
				double x = oe.xmin + (c + 0.5) * oxres;

				cell_vals.clear();
				for (size_t i = 0; i < n; i++) {
					double v;
					if (lookup(br, ri, i, x, y, v)) {
						cell_vals.push_back(v);
					}
				}

				size_t k = cell_vals.size();
				if (k == 0) continue;

				double val = NAN;

				if (fun == "mean") {
					double s = 0;
					for (size_t j = 0; j < k; j++) s += cell_vals[j];
					val = s / k;
				} else if (fun == "sum") {
					double s = 0;
					for (size_t j = 0; j < k; j++) s += cell_vals[j];
					val = s;
				} else if (fun == "min") {
					val = cell_vals[0];
					for (size_t j = 1; j < k; j++) {
						if (cell_vals[j] < val) val = cell_vals[j];
					}
				} else if (fun == "max") {
					val = cell_vals[0];
					for (size_t j = 1; j < k; j++) {
						if (cell_vals[j] > val) val = cell_vals[j];
					}
				} else if (fun == "median") {
					std::sort(cell_vals.begin(), cell_vals.end());
					if (k % 2 == 1) {
						val = cell_vals[k / 2];
					} else {
						val = (cell_vals[k / 2 - 1] +
						       cell_vals[k / 2]) / 2.0;
					}
				}
				result[r * onc + c] = val;
			}
		}

		if (!out.writeValues(result, brow, bnrow)) {
			stop_readers(*this);
			return out;
		}
	}

	out.writeStop();
	stop_readers(*this);
	return out;
}


