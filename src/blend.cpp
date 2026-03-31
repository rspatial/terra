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

	// get output geometry
	SpatExtent ue = ds[0].getExtent();
	for (size_t i = 1; i < n; i++) {

									//  lyrs, crs, warncrs, ext, rowcol, res
		if (!ds[0].compare_geom(ds[i], false, false, opt.get_tolerance(), false, false, false, true)) {
			out.setError(ds[0].getError());
			return(out);
		}		
		ue.unite(ds[i].getExtent());
	}

	out = ds[0].geometry(1, true);
	out.setExtent(ue, true, true, "");

	double oxres = out.xres();
	double oyres = out.yres();
	SpatExtent oe = out.getExtent();
	size_t onc = out.ncol();

	// store per-raster info
	struct RastInfo {
		SpatExtent ext;
		size_t nrow, ncol;
		double xres, yres;
	};

	std::vector<RastInfo> ri(n);
	for (size_t i = 0; i < n; i++) {
		ri[i].ext  = ds[i].getExtent();
		ri[i].nrow = ds[i].nrow();
		ri[i].ncol = ds[i].ncol();
		ri[i].xres = ds[i].xres();
		ri[i].yres = ds[i].yres();
	}


	for (size_t i = 0; i < n; i++) {
		if (!ds[i].readStart()) {
			out.setError(ds[i].getError());
			for (size_t j = 0; j < i; j++) ds[j].readStop();
			return out;
		}
	}
	
	opt.ncopies = std::max(opt.ncopies, size() + nl);
	if (!out.writeStart(opt, filenames())) {
		for (size_t i = 0; i < n; i++) ds[i].readStop();
		return out;
	}

	for (size_t b = 0; b < out.bs.n; b++) {

		size_t brow  = out.bs.row[b];
		size_t bnrow = out.bs.nrows[b];
		size_t bncell = bnrow * onc;

		// y range of this block (top edge of first row, bottom of last row)
		double bymax = oe.ymax - brow * oyres;
		double bymin = oe.ymax - (brow + bnrow) * oyres;

		// Read from each input the rows that overlap this block.
		struct BlockRead {
			std::vector<double> vals;
			size_t nrow, ncol;
			double xmin, ymax; // of the read region
			double xres, yres;
			bool active;
		};

		std::vector<BlockRead> br(n);

		for (size_t i = 0; i < n; i++) {
			br[i].active = false;
			SpatExtent& ei = ri[i].ext;

			// y overlap with this block?
			if (bymin >= ei.ymax || bymax <= ei.ymin) continue;

			// rows that overlap this block
			double y_top = std::min(bymax, ei.ymax);
			double y_bot = std::max(bymin, ei.ymin);

			size_t ir_first = (size_t)std::floor((ei.ymax - y_top) / ri[i].yres);
			size_t ir_last  = (size_t)std::ceil ((ei.ymax - y_bot) / ri[i].yres);
			if (ir_first >= ri[i].nrow) ir_first = ri[i].nrow - 1;
			if (ir_last  >  ri[i].nrow) ir_last  = ri[i].nrow;
			size_t ir_nrow = ir_last - ir_first;
			if (ir_nrow == 0) continue;

			ds[i].readValues(br[i].vals, ir_first, ir_nrow, 0, ri[i].ncol);

			br[i].nrow = ir_nrow;
			br[i].ncol = ri[i].ncol;
			br[i].xmin = ei.xmin;
			br[i].ymax = ei.ymax - ir_first * ri[i].yres;
			br[i].xres = ri[i].xres;
			br[i].yres = ri[i].yres;
			br[i].active = true;
		}

		// Distance-weighted blend for every cell in the block
		std::vector<double> result(bncell, std::numeric_limits<double>::quiet_NaN());

		for (size_t r = 0; r < bnrow; r++) {
			double y = oe.ymax - (brow + r + 0.5) * oyres;

			for (size_t c = 0; c < onc; c++) {
				double x = oe.xmin + (c + 0.5) * oxres;

				double wsum = 0.0;
				double vsum = 0.0;

				for (size_t i = 0; i < n; i++) {
					if (!br[i].active) continue;

					SpatExtent& ei = ri[i].ext;

					// cell centre outside this raster's full extent?
					if (x <= ei.xmin || x >= ei.xmax ||
					    y <= ei.ymin || y >= ei.ymax) continue;

					// row/col in the block-read buffer
					size_t ic = (size_t)((x - br[i].xmin) / br[i].xres);
					size_t ir = (size_t)((br[i].ymax - y)  / br[i].yres);
					if (ic >= br[i].ncol) ic = br[i].ncol - 1;
					if (ir >= br[i].nrow) ir = br[i].nrow - 1;

					double v = br[i].vals[ir * br[i].ncol + ic];
					if (std::isnan(v)) continue;

					// weight = distance to nearest extent edge
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
			for (size_t i = 0; i < n; i++) ds[i].readStop();
			return out;
		}
	}

	out.writeStop();
	for (size_t i = 0; i < n; i++) {
		ds[i].readStop();
	}

	return out;
}
