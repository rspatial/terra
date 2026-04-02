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

#include "recycle.h"
#include "string_utils.h"


// helpers shared by blend, mosaic, and merge

struct RastInfo {
	SpatExtent ext;
	size_t nrow, ncol, nlyr;
	double xres, yres;
};

struct BlockRead {
	std::vector<double> vals;
	size_t nrow, ncol;
	double xmin, ymax;
	double xres, yres;
	bool active;
};

static bool aligns(SpatRasterCollection &src, SpatRaster &out, bool &resample, std::string method, size_t &warn, SpatOptions &opt) {
	size_t n = src.ds.size();
	lrtrim(method);
	double tol = opt.get_tolerance();
	for (size_t i = 1; i < n; i++) {
		if (!src.ds[0].compare_geom(src.ds[i], false, false, tol, false, false, false, true)) {
			if (resample) {	
				if (warn == 0) {
					src = src.deepCopy();
				}
				if (method == "") {
					std::vector<bool> hascats = src.ds[i].hasCategories();
					method = hascats[0] ? "near" : "bilinear";
				}
				SpatOptions wopt(opt);
				SpatRaster tmp = out.geometry();
				tmp = tmp.crop(src.ds[i].getExtent(), "near", false, opt);
				src.ds[i] = src.ds[i].warper(tmp, "", method, false, false, true, wopt);
				if (src.ds[i].hasError()) {
					out.setError(src.ds[i].getError());
					return false;
				}
				warn++;
			} else {
				out.setError(src.ds[0].getError());
				return false;
			}
		}
	}
	return true;
}

static bool setup_output(SpatRasterCollection &src, SpatRaster &out, size_t &nl, SpatOptions &opt) {
	size_t n = src.ds.size();
	SpatExtent ue = src.ds[0].getExtent();
	nl = src.ds[0].nlyr();
	for (size_t i = 1; i < n; i++) {
		ue.unite(src.ds[i].getExtent());
		nl = std::max(nl, src.ds[i].nlyr());
	}
	out = src.ds[0].geometry(nl, true);
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
		ri[i].nlyr = src.ds[i].nlyr();
		ri[i].xres = src.ds[i].xres();
		ri[i].yres = src.ds[i].yres();
	}
	return ri;
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

// Look up the value for a cell at (x,y) in block-read i, layer lyr. Layers are recycled
// Returns false if the cell is outside the raster or the value is NA.
static bool lookup(std::vector<BlockRead> &br, std::vector<RastInfo> &ri,
                   size_t i, size_t lyr, double x, double y, double &val) {
	if (!br[i].active) return false;
	SpatExtent& ei = ri[i].ext;
	if (x <= ei.xmin || x >= ei.xmax ||
	    y <= ei.ymin || y >= ei.ymax) return false;

	size_t ic = (size_t)((x - br[i].xmin) / br[i].xres);
	size_t ir = (size_t)((br[i].ymax - y)  / br[i].yres);
	if (ic >= br[i].ncol) ic = br[i].ncol - 1;
	if (ir >= br[i].nrow) ir = br[i].nrow - 1;

	size_t rl = lyr % ri[i].nlyr;
	size_t cellsPerLyr = br[i].nrow * br[i].ncol;
	val = br[i].vals[rl * cellsPerLyr + ir * br[i].ncol + ic];
	return !std::isnan(val);
}


// Like lookup(), but returns true whenever the cell falls inside the
// raster's extent — even if the value is NaN.  Needed for narm=false.
static bool lookup_raw(std::vector<BlockRead> &br, std::vector<RastInfo> &ri,
                       size_t i, size_t lyr, double x, double y, double &val) {
	if (!br[i].active) return false;
	SpatExtent& ei = ri[i].ext;
	if (x <= ei.xmin || x >= ei.xmax ||
	    y <= ei.ymin || y >= ei.ymax) return false;

	size_t ic = (size_t)((x - br[i].xmin) / br[i].xres);
	size_t ir = (size_t)((br[i].ymax - y)  / br[i].yres);
	if (ic >= br[i].ncol) ic = br[i].ncol - 1;
	if (ir >= br[i].nrow) ir = br[i].nrow - 1;

	size_t rl = lyr % ri[i].nlyr;
	size_t cellsPerLyr = br[i].nrow * br[i].ncol;
	val = br[i].vals[rl * cellsPerLyr + ir * br[i].ncol + ic];
	return true;
}

// === main methods ===


SpatRaster SpatRasterCollection::blend(bool resample, std::string method, SpatOptions &opt) {

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

	size_t nl;
	if (!setup_output(*this, out, nl, opt)) return out;
	size_t warn=0;
	if (!aligns(*this, out, resample, method, warn, opt)) return out;


	double oxres = out.xres();
	double oyres = out.yres();
	SpatExtent oe = out.getExtent();
	size_t onc = out.ncol();

	std::vector<RastInfo> ri = collect_info(*this);

	if (!readStart()) {
		out.setError(getError());
		return out;
	}
	size_t total_lyrs = nl;
	for (size_t i = 0; i < n; i++) total_lyrs += ri[i].nlyr;
	opt.ncopies = std::max(opt.ncopies, total_lyrs);
	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t b = 0; b < out.bs.n; b++) {

		size_t brow   = out.bs.row[b];
		size_t bnrow  = out.bs.nrows[b];
		size_t lyrsize = bnrow * onc;

		double bymax = oe.ymax - brow * oyres;
		double bymin = oe.ymax - (brow + bnrow) * oyres;

		std::vector<BlockRead> br(n);
		read_block(*this, ri, br, bymax, bymin, n);

		std::vector<double> result(lyrsize * nl,
		                           std::numeric_limits<double>::quiet_NaN());

		for (size_t r = 0; r < bnrow; r++) {
			double y = oe.ymax - (brow + r + 0.5) * oyres;

			for (size_t c = 0; c < onc; c++) {
				double x = oe.xmin + (c + 0.5) * oxres;

				for (size_t l = 0; l < nl; l++) {
					double wsum = 0.0;
					double vsum = 0.0;

					for (size_t i = 0; i < n; i++) {
						double v;
						if (!lookup(br, ri, i, l, x, y, v)) continue;

						SpatExtent& ei = ri[i].ext;
						double d = std::min({x - ei.xmin, ei.xmax - x,
						                     y - ei.ymin, ei.ymax - y});
						if (d <= 0.0) d = 1e-10;

						wsum += d;
						vsum += v * d;
					}

					if (wsum > 0.0) {
						result[l * lyrsize + r * onc + c] = vsum / wsum;
					}
				}
			}
		}

		if (!out.writeValues(result, brow, bnrow)) {
			readStop();
			return out;
		}
	}

	out.writeStop();
	readStop();

	if (warn > 0) {
		out.addWarning(std::to_string(warn) + " raster(s) that did not share the base geometry of the first raster were resampled");
	}
	return out;
	
}



SpatRaster SpatRasterCollection::mosaic(std::string fun, bool resample, std::string method, SpatOptions &opt) {

	SpatRaster out;

	std::vector<std::string> valid {"first", "last", "blend", "sum", "mean", "median", "min", "max"};
	if (std::find(valid.begin(), valid.end(), fun) == valid.end()) {
		out.setError(fun + "' is not a valid function name");
		return out;
	}
	if (fun == "blend") return blend(resample, method, opt);
	if (fun == "first") return merge(true,  true, 1, resample, method, opt);
	if (fun == "last")  return merge(false, true, 1, resample, method, opt);

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

	size_t nl;
	if (!setup_output(*this, out, nl, opt)) return out;
	size_t warn=0;
	if (!aligns(*this, out, resample, method, warn, opt)) return out;

	double oxres = out.xres();
	double oyres = out.yres();
	SpatExtent oe = out.getExtent();
	size_t onc = out.ncol();

	std::vector<RastInfo> ri = collect_info(*this);

	if (!readStart()) {
		out.setError(getError());
		return out;
	}

	size_t total_lyrs = nl;
	for (size_t i = 0; i < n; i++) total_lyrs += ri[i].nlyr;
	opt.ncopies = std::max(opt.ncopies, total_lyrs);
	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t b = 0; b < out.bs.n; b++) {

		size_t brow   = out.bs.row[b];
		size_t bnrow  = out.bs.nrows[b];
		size_t lyrsize = bnrow * onc;

		double bymax = oe.ymax - brow * oyres;
		double bymin = oe.ymax - (brow + bnrow) * oyres;

		std::vector<BlockRead> br(n);
		read_block(*this, ri, br, bymax, bymin, n);

		std::vector<double> result(lyrsize * nl, std::numeric_limits<double>::quiet_NaN());

		std::vector<double> cell_vals;
		cell_vals.reserve(n);

		for (size_t r = 0; r < bnrow; r++) {
			double y = oe.ymax - (brow + r + 0.5) * oyres;

			for (size_t c = 0; c < onc; c++) {
				double x = oe.xmin + (c + 0.5) * oxres;

				for (size_t l = 0; l < nl; l++) {
					cell_vals.clear();
					for (size_t i = 0; i < n; i++) {
						double v;
						if (lookup(br, ri, i, l, x, y, v)) {
							cell_vals.push_back(v);
						}
					}

					size_t k = cell_vals.size();
					if (k == 0) continue;

					double val = NAN;

					if (fun == "mean") {
						val = 0;
						for (size_t j = 0; j < k; j++) val += cell_vals[j];
						val /= k;
					} else if (fun == "sum") {
						val = 0;
						for (size_t j = 0; j < k; j++) val += cell_vals[j];
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
					result[l * lyrsize + r * onc + c] = val;
				}
			}
		}

		if (!out.writeValues(result, brow, bnrow)) {
			readStop();
			return out;
		}
	}

	out.writeStop();
	readStop();
	if (warn > 0) {
		out.addWarning(std::to_string(warn) + " raster(s) that did not share the base geometry of the first raster were resampled");
	}
	
	return out;
}



bool write_part(SpatRaster& out, SpatRaster r, const double& hxr, size_t& nl, bool notfirstlyr, std::string method, size_t &warn, SpatOptions &opt) {

	BlockSize bs = r.getBlockSize(opt);
	SpatExtent re = r.getExtent();
	SpatRaster tmp = out.geometry();
	tmp = tmp.crop(r.getExtent(), "near", false, opt);

	if (!tmp.compare_geom(r, false, true, opt.get_tolerance(), false, true, true, false)) {
		std::vector<bool> hascats = r.hasCategories();
		if (method == "") method = hascats[0] ? "near" : "bilinear";
		//std::string filename = tempFile(opt.get_tempdir(), opt.tmpfile, ".tif");
		SpatOptions wopt(opt);
		r = r.warper(tmp, "", method, false, false, true, wopt);
		if (r.hasError()) {
			out.setError(r.getError());
			return false;
		}
		warn++;
		bs = r.getBlockSize(opt);
		re = r.getExtent();
	}

	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v, vout;
		size_t row1  = out.rowFromY(r.yFromRow(bs.row[i]));
		size_t row2  = out.rowFromY(r.yFromRow(bs.row[i]+bs.nrows[i]-1));
		size_t col1  = out.colFromX(re.xmin + hxr);
		size_t col2  = out.colFromX(re.xmax - hxr);
		size_t ncols = col2-col1+1;
		size_t nrows = row2-row1+1;

		if (!r.readStart()) {
			out.setError(r.getError());
			return false;
		}
		r.readBlock(v, bs, i);
		recycle(v, ncols * nrows * nl);

		if (notfirstlyr) {
			out.readValuesWhileWriting(vout, row1, nrows, col1, ncols);
			for (size_t j=0; j<v.size(); j++) {
				if (std::isnan(v[j])) {
					v[j] = vout[j];
				}
			}
		}
		if (!out.writeValuesRect(v, row1, nrows, col1, ncols)) return false;
	}
	r.readStop();
	return true;
}


SpatRaster SpatRasterCollection::merge(bool first, bool narm,  size_t algo, bool resample, std::string method, SpatOptions &opt) {

	SpatRaster out;
	size_t n = size();
	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	if (n == 1) {
		if (opt.get_filename() != "") {
			out = ds[0].writeRaster(opt);
		} else {
			out = ds[0].deepCopy();
		}
		return out;
	}
	if (algo == 1) {

		SpatExtent e = ds[0].getExtent();
		size_t nl = ds[0].nlyr();
		bool anyvals = false;
		for (size_t i=1; i<n; i++) {
			e.unite(ds[i].getExtent());
			if (ds[i].hasValues()) anyvals = true;
			nl = std::max(nl, ds[i].nlyr());
		}

		out = ds[0].geometry(nl, true);
		out.setExtent(e, true, true, "");

		if (!anyvals) return out;

		double hxr = out.xres()/2;

		std::vector<int> vt = getValueType(true);
			if (vt.size() == 1) {
			out.setValueType(vt[0]);
		}

		opt.ncopies = std::max(opt.ncopies, size() + nl);
		if (!out.writeStart(opt, filenames())) { return out; }

		std::vector<size_t> seq(n);
		if (first) {
			std::iota(seq.rbegin(), seq.rend(), 0);
		} else {
			std::iota(seq.begin(), seq.end(), 0);
		}

		SpatOptions topt(opt);

		size_t warn = 0;
		bool notfirst = false;
		for (size_t i=0; i<n; i++) {
			if (!ds[seq[i]].hasValues()) continue;
			if (narm) {
				notfirst = i > 0;
			}
			if (!write_part(out, ds[seq[i]], hxr, nl, notfirst, method, warn, topt)) {
				return out;
			}
		}
		out.writeStop();
		if (warn > 0) {
			out.addWarning(std::to_string(warn) + " raster(s) that did not share the base geometry of the first raster were resampled");
		}
		return(out);

	} else if (algo == 2) {

		size_t nl;
		if (!setup_output(*this, out, nl, opt)) return out;
		size_t warn=0;
		if (!aligns(*this, out, resample, method, warn, opt)) return out;

		double oxres = out.xres();
		double oyres = out.yres();
		SpatExtent oe = out.getExtent();
		size_t onc = out.ncol();

		std::vector<RastInfo> ri = collect_info(*this);

		if (!readStart()) {
			out.setError(getError());
			return out;
		}

		size_t total_lyrs = nl;
		for (size_t i = 0; i < n; i++) total_lyrs += ri[i].nlyr;
		opt.ncopies = std::max(opt.ncopies, total_lyrs);
		if (!out.writeStart(opt, filenames())) {
			readStop();
			return out;
		}

		std::vector<size_t> seq(n);
		if (first) {
			std::iota(seq.begin(), seq.end(), 0);
		} else {
			std::iota(seq.rbegin(), seq.rend(), 0);
		}

		for (size_t b = 0; b < out.bs.n; b++) {

			size_t brow   = out.bs.row[b];
			size_t bnrow  = out.bs.nrows[b];
			size_t lyrsize = bnrow * onc;

			double bymax = oe.ymax - brow * oyres;
			double bymin = oe.ymax - (brow + bnrow) * oyres;

			std::vector<BlockRead> br(n);
			read_block(*this, ri, br, bymax, bymin, n);

			std::vector<double> result(lyrsize * nl,
									   std::numeric_limits<double>::quiet_NaN());

			for (size_t r = 0; r < bnrow; r++) {
				double y = oe.ymax - (brow + r + 0.5) * oyres;

				for (size_t c = 0; c < onc; c++) {
					double x = oe.xmin + (c + 0.5) * oxres;

					for (size_t l = 0; l < nl; l++) {
						for (size_t s = 0; s < n; s++) {
							size_t i = seq[s];
							double v;
							if (narm) {
								if (lookup(br, ri, i, l, x, y, v)) {
									result[l * lyrsize + r * onc + c] = v;
									break;
								}
							} else {
								if (lookup_raw(br, ri, i, l, x, y, v)) {
									result[l * lyrsize + r * onc + c] = v;
									break;
								}
							}
						}
					}
				}
			}

			if (!out.writeValues(result, brow, bnrow)) {
				readStop();
				return out;
			}
		}

		out.writeStop();
		readStop();
		if (warn > 0) {
			out.addWarning(std::to_string(warn) + " raster(s) that did not share the base geometry of the first raster were resampled");
		}
		return out;


	} else if (algo==3) {

// narm is not used

		SpatExtent e = ds[0].getExtent();
		size_t nl = ds[0].nlyr();
		for (size_t i=1; i<n; i++) {
			if (nl != ds[i].nlyr()) {
				out.setError("you cannot use this algo with rasters with different numbers of layers");
				return out;
			}
			e.unite(ds[i].getExtent());
		}

		out = ds[0].geometry(nl, true);
		out.setExtent(e, true, true, "");

		SpatRaster tmp = out;
		for (size_t i=1; i<n; i++) {
//			SpatRaster tmp = out.crop(ds[i].getExtent(), "near", false, opt);
// check crs only
			tmp.compare_geom(ds[i], false, true, opt.get_tolerance(), false, false, false, false);
			if (tmp.hasError()) {
				return tmp;
			}
		}

		if (method == "") method = ds[0].hasCategories()[0] ? "nearest" : "bilinear";
		std::vector<std::string> options = {"-r", method};

		bool wvrt = false;
		std::string fout = opt.get_filename();
		if (!fout.empty()) {
			if (fout.size() > 4) {
				std::string ss = fout.substr(fout.size()-4, fout.size());
				lowercase(ss);
				wvrt = ss == ".vrt";
			}
		}

		std::vector<std::string> warnings;
		if (fout != "") {
			std::string fname;
			if (opt.names.empty()) {
				opt.names = ds[0].getNames();
			}

			if (wvrt) {
				fname = make_vrt(options, first, opt);
				warnings = opt.msg.warnings;
			} else {
				SpatOptions vopt(opt);
				fname = make_vrt(options, first, vopt);
				warnings = vopt.msg.warnings;
			}

			if (hasError()) {
				out.setError(getError());
				return out;
			}
			SpatRaster v(fname, {}, {}, {}, {}, false, false, {});
			if (warnings.size() > 0) {
				v.msg.warnings = warnings;
				v.msg.has_warning = true;
			}
			if (wvrt) {
				return v;
			} else {
				return v.writeRaster(opt);
			}
		}

		SpatOptions vopt(opt);
		std::string fname = make_vrt(options, first, vopt);
		SpatRaster v(fname, {}, {}, {}, {}, false, false, {});
		v.setNames(ds[0].getNames(), false);

		if (vopt.msg.warnings.size() > 0) {
			v.msg.warnings = warnings;
			v.msg.has_warning = true;
		}

		return v;

	} else {
		out.setError("invalid algo (should be 1, 2, or 3)");
		return out;
	}
}



SpatRaster SpatRasterCollection::morph(SpatRaster &x, SpatOptions &opt) {

	SpatRaster out;
	size_t n = size();
	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	std::string filename = opt.get_filename();
	opt.set_filenames({""});
	SpatExtent e = x.getExtent();

	out.source.resize(0);
	SpatRaster g = x.geometry();
	SpatOptions topt(opt);
	for (size_t i=0; i<n; i++) {
		if (g.compare_geom(ds[i], false, false, 0.01, false, true, true, false)) {
			out.source.insert(out.source.end(), ds[i].source.begin(), ds[i].source.end());
		} else {
			// should first consider whether going up or down in resolution
			// and perhaps use (dis) aggregate (first)
			std::vector<bool> hasCats = ds[i].hasCategories();
			// this should be done by layer
			bool call = true;
			for (size_t j=0; j<hasCats.size(); j++) {
				call = call && hasCats[j];
			}
			std::string method = call ? "near" : "bilinear";
			SpatRaster temp = ds[i].warper(g, "", method, false, false, false, topt);
			out.addSource(temp, false, topt);
		}
	}

	if (out.source.empty()) {
		out.setError("no data sources that overlap with x");
		return out;
	}

	out.setSRS(x.getSRS("wkt"));
	out.setExtent(e, false, true, "near");

	lrtrim(filename);
	if (!filename.empty()) {
		opt.set_filenames({filename});
		out.writeRaster(opt);
	}
	return(out);
}

/*
old merge algo2
	} else if (algo == 2) {

// narm is not used

		std::vector<size_t> use;
		use.reserve(n);
		if (ds[0].hasValues()) use.push_back(0);
		SpatExtent e = ds[0].getExtent();
		size_t nl = ds[0].nlyr();
		for (size_t i=1; i<n; i++) {
			e.unite(ds[i].getExtent());
			if (ds[i].hasValues()) use.push_back(i);
			nl = std::max(nl, ds[i].nlyr());
		}

		out = ds[0].geometry(nl, true);
		out.setExtent(e, true, true, "");
		if (use.empty()) return out;

		SpatOptions topt(opt);

		for (size_t i=0; i<use.size(); i++) {
			//  lyrs, crs, warncrs, ext, rowcol, res
			SpatRaster tmp = out.crop(ds[use[i]].getExtent(), "near", false, topt);
			if (!tmp.compare_geom(ds[use[i]], false, true, topt.get_tolerance(), true, true, true, false)) {
				out.setError("rasters geometries are not compatible");
				return(out);
			}
		}

		std::vector<int> vt = getValueType(true);
			if (vt.size() == 1) {
			out.setValueType(vt[0]);
		}

		opt.ncopies = std::max(opt.ncopies, size() + nl);
		if (!out.writeStart(opt, filenames())) { return out; }

		SpatOptions ropt(opt);

		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<std::vector<double>> v;
			readBlock(out, v, out.bs, i, use, ropt);
			if (hasError()) {
				out.writeStop();
				out.setError(getError());
				return out;
			}
			std::vector<size_t> sizes(v.size());

			size_t n = nl * out.bs.nrows[i] * out.ncol();
			size_t m = v.size();

			bool multi_sz = false;
			std::vector<size_t> sz(m);
			for (size_t j=0; j<m; j++) {
				sz[j] = v[j].size();
				if (v[j].size() < n) multi_sz = true;
				if (v[j].size() > n) {
				out.setError("something is not right. Exected: " + std::to_string(n) + " got: " + std::to_string(v[j].size()) + " values");
					return out;
				}
			}

			if (m == 1) {
				if (!out.writeBlock(v[0], i)) return out;
			} else if (first) {
				recycle(v[0], n);
				if (multi_sz) { // with recycling
					for (size_t j=0; j<n; j++) {
						for (size_t k=1; k<m; k++) {
							if (std::isnan(v[0][j])) {
								v[0][j] = v[k][j%sz[k]];
							} else {
								continue;
							}
						}
					}
				} else {
					for (size_t j=0; j<n; j++) {
						for (size_t k=1; k<m; k++) {
							if (std::isnan(v[0][j])) {
								v[0][j] = v[k][j];
							} else {
								continue;
							}
						}
					}
				}
				if (!out.writeBlock(v[0], i)) return out;
			} else { // last
				m -= 1;
				recycle(v[m], n);
				if (multi_sz) {
					for (size_t j=0; j<n; j++) {
						for (long k=(m-1); k >= 0; k--) {
							if (std::isnan(v[m][j])) {
								v[m][j] = v[k][j%sz[k]];
							} else {
								continue;
							}
						}
					}
				} else {
					for (size_t j=0; j<n; j++) {
						for (long k=(m-1); k >= 0; k--) {
							if (std::isnan(v[m][j])) {
								v[m][j] = v[k][j];
							} else {
								continue;
							}
						}
					}
				}
				if (!out.writeBlock(v[m], i)) return out;
			}
		}
		out.writeStop();
		return(out);
*/


/*
bool overlaps(const std::vector<unsigned>& r1, const std::vector<unsigned>& r2,
			  const std::vector<unsigned>& c1, const std::vector<unsigned>& c2) {
	size_t n = r1.size();
	for (size_t i=0; i<(n-1); i++) {
		for (size_t j=(i+1); j<n; j++) {
			if ((r1[i] <= r2[j]) && (r2[i] >= r1[j]) && (c1[i] <= c2[j]) && (c2[i] >= c1[j])) {
				return true;
			}
		}
	}
	return false;
}

// old 

SpatRaster SpatRasterCollection::mosaic2(std::string fun, SpatOptions &opt) {

	SpatRaster out;
	std::vector<std::string> f {"first", "last", "sum", "mean", "median", "min", "max", "modal"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("argument 'fun' is not a valid function name");
		return out;
	}

	if ((fun == "first") || (fun == "last")) {
		return merge(fun=="first", true, 1, "", opt);
	}

	size_t n = size();

	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	if (n == 1) {
		if (opt.get_filename() != "") {
			out = ds[0].writeRaster(opt);
		} else {
			out = ds[0].deepCopy();
		}
		return(out);
	}

	std::vector<bool> hvals(n);
	hvals[0] = ds[0].hasValues();
	SpatExtent e = ds[0].getExtent();
	size_t nl = ds[0].nlyr();
//std::vector<bool> resample(n, false);

	std::vector<size_t> use;
	use.reserve(n);
	if (hvals[0]) use.push_back(0);
	for (size_t i=1; i<n; i++) {
		SpatExtent ee = ds[i].getExtent();
									//  lyrs, crs, warncrs, ext, rowcol, res
		if (!ds[0].compare_geom(ds[i], false, false, opt.get_tolerance(), false, false, false, true)) {
			out.setError(ds[0].getError());
			return(out);
		}
		e.unite(ee);
		hvals[i] = ds[i].hasValues();
		if (hvals[i]) use.push_back(i);
		nl = std::max(nl, ds[i].nlyr());
	}
	out = ds[0].geometry(nl, false);
	out.setExtent(e, true, true, "");

	for (int i=(n-1); i>=0; i--) {
		if (!hvals[i]) {
			erase(i);
		}
	}

	n = size();
	if (size() == 0) {
		if (opt.get_filename() != "") {
			out = out.writeRaster(opt);
		} else {
			out = ds[0].deepCopy();
		}
		return(out);
	}

//	if (!overlaps(r1, r2, c1, c2)) {
//		return merge(true, true, opt);
//	}
	double ncl = 1000;
	if (n > 50) ncl = 500;
	if (n > 100) ncl = 250;
	double     ar = std::ceil(out.nrow() / ncl);
	size_t arow = std::ceil(out.nrow() / ar);
	double     ac = std::ceil(out.ncol() / ncl);
	size_t acol = std::ceil(out.ncol() / ac);

	SpatOptions sopt(opt);
	SpatExtent ae = out.getExtent();

	SpatRaster aout = out.aggregate({arow, acol}, "", true, sopt);
	SpatVector ve = aout.as_polygons(false, false, false, false, false, 0, sopt);
	SpatVector vcrp(out.getExtent(), "");

	ve = ve.intersect(vcrp, false);

	size_t nv = ve.nrow();
	bool warn = false;

 	if (!out.writeStart(opt, filenames())) { return out; }
	sopt.progressbar = false;

	SpatRasterStack s;

	for (size_t i=0; i<nv; i++) {
		SpatVector vi = ve.subset_rows(i);
		SpatExtent ce = vi.getExtent();
		SpatRasterCollection x = crop(ce, "near", true, use, sopt);
		if (x.empty()) {
			continue;
		}

//		Rcpp::Rcout << "ext e:" << ce.xmin << " " << ce.xmax << " " << ce.ymin << " " << ce.ymax << std::endl;
//		for (size_t j=0; j<x.size(); j++) {
//			Rcpp::Rcout << "ext " << j << ": " << x.ds[j].source[0].extent.xmin << " " << x.ds[j].source[0].extent.xmax
//				<<  " " << x.ds[j].source[0].extent.ymin <<  " " << x.ds[j].source[0].extent.ymax <<  " " << std::endl;
//		}
//		continue;

		s.ds = x.ds;
		SpatRaster r;
		if (nl == 1) {
		// work-around, need to investigate why this is needed for 1129
// see #1129
// if (i == 57 || i == 79 | i == 269) { // && (rcnt[i] == 6)) {
			r = s.collapse();
			r = r.summary(fun, true, sopt);
		} else {
			r = s.summary(fun, true, sopt);
		}
		if (r.hasError()) {
			return r;
		}
		if (!out.writeValuesRectRast(r, opt)) {
			return out;
		}
	}
	out.writeStop();

	if (warn) out.addWarning("rasters did not align and were resampled");
	return out;
}
*/
