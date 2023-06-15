// Copyright (c) 2018-2023  Robert J. Hijmans
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
#include "recycle.h"
#include "vecmath.h"
#include "vecmathse.h"
#include <cmath>
#include <functional>

#include "math_utils.h"
#include "file_utils.h"
#include "string_utils.h"
#include "sort.h"


/*
std::vector<double> flat(std::vector<std::vector<double>> v) {
    unsigned s1 = v.size();
    unsigned s2 = v[0].size();

	std::size_t s = s1 * s2;
    std::vector<double> result(s);
    for (size_t i=0; i<s1; i++) {
		for (size_t j=0; j<s2; j++) {
			result[i*s2+j] = v[i][j];
		}
	}
	return result;
}
*/

/*
SpatRaster SpatRaster::selectHighest(size_t n, bool low, SpatOptions &opt) {

	SpatVector out;

	if (nlyr() > 1) {
		SpatOptions ops(opt);
		out.addWarning("only processing the first layer");
		std::vector<unsigned> lyr = {0};
		*this = subset(lyr, ops);
	}
	if (!hasValues()) {
		return(out);
	}
	if (n >= ncell()) {
		return isnotnan(opt);
	}

	std::vector<double> sel;

	if (!readStart()) {
		return(out);
	}

	BlockSize bs = getBlockSize(opt);
	for (size_t i = 0; i < bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		for (size_t j=0; j<v.size(); j++) {
		}

	readStop();
	return(out);
}
*/



SpatExtent SpatRaster::ext_from_rc(int_64 r1, int_64 r2, int_64 c1, int_64 c2) {
	SpatExtent e = getExtent();
	double xrs = xres();
	double yrs = yres();
	int_64 nr = nrow();
	int_64 nc = ncol();
	c1 = std::min(std::max(c1, (int_64)0), nc);
	c2 = std::min(std::max(c2,  (int_64)0), nc);
	if (c1 > c2) {
		std::swap(c1, c2);
	}
	r1 = std::min(std::max(r1, (int_64)0), nr);
	r2 = std::min(std::max(r2, (int_64)0), nr);
	if (r1 > r2) {
		std::swap(r1, r2);
	}

	double xn = xFromCol(c1) - 0.5 * xrs;
	double xx = xFromCol(c2) + 0.5 * xrs;
	double yx = yFromRow(r1) + 0.5 * yrs;
	double yn = yFromRow(r2) - 0.5 * yrs;

	return SpatExtent(xn, xx, yn, yx);
}


SpatExtent SpatRaster::ext_from_cell(double cell) {
	std::vector<double> cells = {cell};
	std::vector<std::vector<int_64>> rc = rowColFromCell(cells);
	return ext_from_rc(rc[0][0], rc[0][0], rc[1][0], rc[1][0]);
}


std::vector<std::string> SpatRaster::make_tiles(SpatRaster x, bool expand, bool narm, std::string filename, SpatOptions &opt) {

	std::vector<std::string> ff;
	if (!hasValues()) {
		setError("input raster has no values");
		return ff;
	}
	x = x.geometry(1, false, false, false);
	SpatExtent e = getExtent();
	SpatOptions ops(opt);
	if (expand) {
		x = x.extend(e, "out", NAN, ops);
	}
	x = x.crop(e, "out", false, ops);

	std::vector<size_t> d(x.ncell());
	std::iota(d.begin(), d.end(), 1);

	std::string fext = getFileExt(filename);
	std::string f = noext(filename);
	ff.reserve(d.size());
	size_t nl = nlyr();
	bool overwrite = opt.get_overwrite();
	for (size_t i=0; i<d.size(); i++) {
		std::string fout = f + std::to_string(d[i]) + fext;
		if (file_exists(fout) && (!overwrite)) {
			ff.push_back(fout);			
			continue;
		}
		SpatExtent exi = x.ext_from_cell(i);
		opt.set_filenames({fout});
		SpatRaster out = crop(exi, "near", false, opt);
		if (out.hasError()) {
			setError(out.getError());
			return ff;
		}
		if ( out.hasValues() ) {
			if (narm) {
				std::vector<double> rmin = out.range_min();
				size_t cnt = 0;
				for (double &v : rmin) {
					if (std::isnan(v)) cnt++;
				}
				if (cnt == nl) {
					remove(fout.c_str());
					continue;
				}
			}
			ff.push_back(fout);
		}
	}
	return ff;
}



std::vector<std::string> SpatRaster::make_tiles_vect(SpatVector x, bool expand, bool narm, std::string filename, SpatOptions &opt) {

	std::vector<std::string> ff;
	if (!hasValues()) {
		setError("input raster has no values");
		return ff;
	}
	if (x.type() != "polygons") {
		setError("The SpatVector must have a polygons geometry");
		return ff;		
	}
	SpatExtent e = getExtent();
	SpatOptions ops(opt);
	std::vector<size_t> d(x.size());
	std::iota(d.begin(), d.end(), 1);

	std::string fext = getFileExt(filename);
	std::string f = noext(filename);
	ff.reserve(d.size());
	size_t nl = nlyr();
	bool overwrite = opt.get_overwrite();
	for (size_t i=0; i<d.size(); i++) {
		std::string fout = f + std::to_string(d[i]) + fext;
		if (file_exists(fout) && (!overwrite)) {
			ff.push_back(fout);			
			continue;
		}
		opt.set_filenames( {fout} );
		SpatRaster out;
		SpatExtent exi = x.geoms[i].extent;
		if (!e.intersects(exi)) continue;
		if (expand) {
			out = crop(exi, "near", false, ops);
			out = out.extend(exi, "out", NAN, opt);
		} else {
			out = crop(exi, "near", false, opt);
		}
		if (out.hasError()) {
			setError(out.getError());
			return ff;
		}
		if ( out.hasValues() ) {
			if (narm) {
				std::vector<double> rmin = out.range_min();
				size_t cnt = 0;
				for (double &v : rmin) {
					if (std::isnan(v)) cnt++;
				}
				if (cnt == nl) {
					remove(fout.c_str());
					continue;
				}
			}
			ff.push_back(fout);
		}
	}
	return ff;
}


bool SpatRaster::get_aggregate_dims(std::vector<unsigned> &fact, std::string &message ) {

	unsigned fs = fact.size();
	if ((fs > 3) | (fs == 0)) {
		message = "argument 'fact' should have length 1, 2, or 3";
		return false;
	}
	auto min_value = *std::min_element(fact.begin(),fact.end());
	if (min_value < 1) {
		message = "values in argument 'fact' should be > 0";
		return false;
	}
	auto max_value = *std::max_element(fact.begin(),fact.end());
	if (max_value == 1) {
		message = "all values in argument 'fact' are 1, nothing to do";
		return false;
	}

	fact.resize(6);
	if (fs == 1) {
		fact[1] = fact[0];
		fact[2] = 1;
	} else if (fs == 2) {
		fact[2] = 1;
	}
	// int dy = dim[0], dx = dim[1], dz = dim[2];
	fact[0] = fact[0] < 1 ? 1 : fact[0];
	fact[0] = fact[0] > nrow() ? nrow() : fact[0];
	fact[1] = fact[1] < 1 ? 1 : fact[1];
	fact[1] = fact[1] > ncol() ? ncol() : fact[1];
	fact[2] = std::max(unsigned(1), std::min(fact[2], nlyr()));
	// new dimensions: rows, cols, lays
	fact[3] = std::ceil(double(nrow()) / fact[0]);
	fact[4] = std::ceil(double(ncol()) / fact[1]);
	fact[5] = std::ceil(double(nlyr()) / fact[2]);
	return true;
}


std::vector<unsigned> SpatRaster::get_aggregate_dims2(std::vector<unsigned> fact) {
	// for use with R
	std::string message = "";
	get_aggregate_dims(fact, message);
	return(fact);
}


std::vector<std::vector<double>> SpatRaster::get_aggregates(std::vector<double> &in, size_t nr, std::vector<unsigned> dim) {

// dim 0, 1, 2, are the aggregations factors dy, dx, dz
// and 3, 4, 5 are the new nrow, ncol, nlyr

// adjust for chunk
	//dim[3] = std::ceil(double(nr) / dim[0]);
	//size_t bpC = dim[3];
	size_t bpC = std::ceil(double(nr) / dim[0]);

	size_t dy = dim[0], dx = dim[1], dz = dim[2];
	size_t bpR = dim[4];
	size_t bpL = bpR * bpC;

	// new number of layers
	size_t newNL = dim[5];

	// new number of rows, adjusted for additional (expansion) rows
	size_t adjnr = bpC * dy;

	// number of aggregates
	size_t nblocks = (bpR * bpC * newNL);
	// cells per aggregate
	size_t blockcells = dx * dy * dz;

	// output: each row is a block
	std::vector< std::vector<double> > a(nblocks, std::vector<double>(blockcells, std::numeric_limits<double>::quiet_NaN()));

    size_t nc = ncol();
    // size_t ncells = ncell();
    size_t ncells = nr * nc;
    size_t nl = nlyr();
    size_t lstart, rstart, cstart, lmax, rmax, cmax, f, lj, cell;

	for (size_t b = 0; b < nblocks; b++) {
		lstart = dz * (b / bpL);
		rstart = (dy * (b / bpR)) % adjnr;
		cstart = dx * (b % bpR);

		lmax = std::min(nl, (lstart + dz));
		rmax = std::min(nr, (rstart + dy));  // nrow -> nr
		cmax = std::min(nc, (cstart + dx));

		f = 0;
		for (size_t j = lstart; j < lmax; j++) {
			lj = j * ncells;
			for (size_t r = rstart; r < rmax; r++) {
				cell = lj + r * nc;
				for (size_t c = cstart; c < cmax; c++) {
					a[b][f] = in[cell + c];
					f++;
				}
			}
		}
	}
	return(a);
}


void compute_aggregates(const std::vector<double> &in, std::vector<double> &out, size_t nr, size_t nc, size_t nl, std::vector<unsigned> dim, std::function<double(std::vector<double>&, bool)> fun, bool narm) {

// dim 0, 1, 2, are the aggregations factors dy, dx, dz
// and 3, 4, 5 are the new nrow, ncol, nlyr

	size_t dy = dim[0], dx = dim[1], dz = dim[2];
//	size_t bpC = dim[3];
// adjust for chunk
//	size_t bpC = std::ceil(double(nr) / dim[0]);
	size_t bpC = std::ceil((double)nr / (double)dim[0]);

	size_t bpR = dim[4];
	size_t bpL = bpR * bpC;

	// new number of layers
	size_t newNL = dim[5];

	// new number of rows, adjusted for additional (expansion) rows
	size_t adjnr = bpC * dy;

	// number of aggregates
	size_t nblocks = (bpR * bpC * newNL);
	// cells per aggregate
	size_t blockcells = dx * dy * dz;

	// output: each row is a block
	out = std::vector<double>(nblocks, NAN);

//    size_t nl = nlyr();
//    size_t nc = ncol();
    size_t ncells = nr * nc;
    size_t lstart, rstart, cstart, lmax, rmax, cmax, f, lj, cell;

	for (size_t b = 0; b < nblocks; b++) {
		lstart = dz * (b / bpL);
		rstart = (dy * (b / bpR)) % adjnr;
		cstart = dx * (b % bpR);

		lmax = std::min(nl, (lstart + dz));
		rmax = std::min(nr, (rstart + dy));  // nrow -> nr
		cmax = std::min(nc, (cstart + dx));

		f = 0;
		std::vector<double> a(blockcells, NAN);
		for (size_t j = lstart; j < lmax; j++) {
			lj = j * ncells;
			for (size_t r = rstart; r < rmax; r++) {
				cell = lj + r * nc;
				for (size_t c = cstart; c < cmax; c++) {
					a[f] = in[cell + c];
					f++;
				}
			}
		}
		out[b] = fun(a, narm);
	}
}



SpatRaster SpatRaster::aggregate(std::vector<unsigned> fact, std::string fun, bool narm, SpatOptions &opt) {

	SpatRaster out;
	std::string message = "";
// fact 0, 1, 2, are the aggregation factors dy, dx, dz
// and  3, 4, 5 are the new nrow, ncol, nlyr
	if (!get_aggregate_dims(fact, message)) {
		if (message.substr(0,3) == "all") {
			std::string filename = opt.get_filename();
			if (filename.empty()) {
				out = *this;
				out.addWarning(message);
			} else {
				out = writeRaster(opt);
			}
		} else {
			out.setError(message);
		}
		return out;
	}

	SpatExtent extent = getExtent();
	double xmax = extent.xmin + fact[4] * fact[1] * xres();
	double ymin = extent.ymax - fact[3] * fact[0] * yres();
	SpatExtent e = SpatExtent(extent.xmin, xmax, ymin, extent.ymax);
	out = SpatRaster(fact[3], fact[4], fact[5], e, "");
	out.source[0].srs = source[0].srs;
	// there is much more. categories, time. should use geometry and then
	// set extent and row col
	if (fact[5] == nlyr()) {
		out.setNames(getNames());
	}

	if (!source[0].hasValues) {
		return out;
	}

	if (!haveFun(fun)) {
		out.setError("unknown function argument");
		return out;
	}

/*
	size_t ifun = std::distance(f.begin(), it);
	std::string gstring = "";
	if (ifun > 0) {
		std::vector<std::string> gf {"average", "min", "max", "med", "mode"};
		gstring = gf[ifun-1];
	}

#ifdef useGDAL
#if GDAL_VERSION_MAJOR >= 3
	if (gstring != "") {
		out = warper(out, "", gstring, opt);
		return out;
	}
#endif
#endif
*/
	std::function<double(std::vector<double>&, bool)> agFun = getFun(fun);

	//BlockSize bs = getBlockSize(4, opt.get_memfrac());
	opt.progress *= 300;
	BlockSize bs = getBlockSize(opt);
	//bs.n = floor(nrow() / fact[0]); # ambiguous on solaris
	bs.n = std::floor(static_cast <double> (nrow() / fact[0]));

	bs.nrows = std::vector<size_t>(bs.n, fact[0]);
	bs.row.resize(bs.n);
	for (size_t i =0; i<bs.n; i++) {
		bs.row[i] = i * fact[0];
	}
	size_t lastrow = bs.row[bs.n - 1] + bs.nrows[bs.n - 1]; // + 1;
	if (lastrow < nrow()) {

		bs.row.push_back(lastrow);
		bs.nrows.push_back(std::min(bs.nrows[bs.n-1], nrow()-lastrow));
		bs.n += 1;
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	opt.steps = bs.n;
	opt.minrows = fact[0];

	if (fun == "modal") {
		if (nlyr() == out.nlyr()) {
			out.source[0].hasColors = hasColors();
			out.source[0].cols = getColors();
			out.source[0].hasCategories = hasCategories();
			out.source[0].cats = getCategories();
		}
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	size_t nc = ncol();
	//size_t outnc = out.ncol();
	for (size_t i = 0; i < bs.n; i++) {
        std::vector<double> vin, v;
		readValues(vin, bs.row[i], bs.nrows[i], 0, nc);
		compute_aggregates(vin, v, bs.nrows[i], nc, nlyr(), fact, agFun, narm);
		if (!out.writeValues(v, i, 1)) return out;
		//if (!out.writeValuesRect(v, i, 1, 0, outnc)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}




SpatRaster SpatRaster::weighted_mean(SpatRaster w, bool narm, SpatOptions &opt) {
	SpatRaster out;
	if (nlyr() != w.nlyr()) {
		out.setError("nlyr of data and weights are different");
		return out;
	}

	SpatOptions topt(opt);
	out = arith(w, "*", false, topt);
	out = out.summary("sum", narm, topt);
	if (narm) {
		w = w.mask(*this, false, NAN, NAN, topt);
	}
	SpatRaster wsum = w.summary("sum", narm, topt);
	return out.arith(wsum, "/", false, opt);
}


SpatRaster SpatRaster::weighted_mean(std::vector<double> w, bool narm, SpatOptions &opt) {

	SpatRaster out;
	for (size_t i=0; i<w.size(); i++) {
		if (std::isnan(w[i]) || w[i] <= 0) {
			out.setError("all weights must be positive values");
			return out;
		}
	}

	unsigned nl = nlyr();
	//if (nl == 1) return *this; this is not consistent if weight is zero
	recycle(w, nl);

	if (narm) {
		if (!hasValues()) {
			out.setError("raster has no values");
			return out;
		}
		out = geometry(1);
		if (!readStart()) {
			out.setError(getError());
			return(out);
		}
		if (!out.writeStart(opt, filenames())) {
			readStop();
			return out;
		}
		unsigned nc = ncol();

		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			size_t off = out.bs.nrows[i] * nc;
			std::vector<double> wm(off, 0);
			std::vector<double> wv(off, 0);
			for (size_t j=0; j<nl; j++) {
				size_t start = j * off;
				size_t end = start + off;
				for (size_t k=start; k<end; k++) {
					if (!std::isnan(v[k])) {
						size_t kj = k - start;
						wm[kj] += v[k] * w[j];
						wv[kj] += w[j];
					}
				}
			}
			for (size_t k=0; k<wm.size(); k++) {
				if (wv[k] == 0) {
					wm[k] = NAN;
				} else {
					wm[k] /= wv[k];
				}
			}
			if (!out.writeBlock(wm, i)) return out;
		}
		out.writeStop();
		readStop();
		return(out);

	} else {
		SpatOptions topt(opt);
		out = arith(w, "*", false, false, topt);
		out = out.summary("sum", narm, topt);
		double wsum = vsum(w, narm);
		return out.arith(wsum, "/", false, false, opt);
	}
}


SpatRaster SpatRaster::separate(std::vector<double> classes, double keepvalue, double othervalue, bool round, int digits, SpatOptions &opt) {

	SpatRaster out;
	if (nlyr() > 1) {
		out.setError("input may only have one layer");
		return out;
	}
	if (classes.empty()) {
		SpatOptions topt(opt);
		std::vector<std::vector<double>> rc = unique(false, NAN, true, topt);
		classes = rc[0];
	}

	if (round) {
		for (size_t i=0; i<classes.size(); i++) {
			classes[i] = roundn(classes[i], digits);
		}
	}
	std::sort(classes.begin(), classes.end());
	classes.erase(std::unique(classes.begin(), classes.end()), classes.end());

	size_t n = classes.size();
	if (n == 0) {
		out.setError("no valid classes");
		return out;
	}
	out = geometry(n);
	std::vector<std::string> snms(n);
	for (size_t i=0; i<n; i++) {
		snms[i] = double_to_string(classes[i]);
	}
	out.setNames(snms);

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		if (round) {
			for(double& d : v) d = roundn(d, digits);
		}
		size_t nn = v.size();
		std::vector<double> vv(nn * n, NAN);
		for (size_t j=0; j<nn; j++) {
			if (!std::isnan(v[j])) {
				for (size_t k=0; k<classes.size(); k++) {
					if (v[j] == classes[k]) {
						if (keepvalue) {
							vv[j + k*nn] = classes[k];
						} else {
							vv[j + k*nn] = 1;	 // true
						}
					} else {
						vv[j + k*nn] = othervalue;
					}
				}

			}
		}
		if (!out.writeBlock(vv, i)) {
			readStop();
			return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::is_in(std::vector<double> m, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (m.empty()) {
		out.setError("no matches supplied");
		return(out);
	}
	if (!hasValues()) {
		out.setError("input has no values");
		return(out);
	}

	int hasNAN = 0;
	for (size_t i=0; i<m.size(); i++) {
		if (std::isnan(m[i])) {
			hasNAN = 1;
			m.erase(m.begin()+i);
			break;
		}
	}
	if (m.empty()) { // only NA
		return isnan(false, opt);
	}


	// if m is very long, perhaps first check if the value is in range?

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	out.setValueType(3);
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		std::vector<double> vv(v.size(), 0);
		for (size_t j=0; j<v.size(); j++) {
			if (std::isnan(v[j])) {
				vv[j] = hasNAN;
			} else {
				for (size_t k=0; k<m.size(); k++) {
					if (v[j] == m[k]) {
						vv[j] = 1;
						break;
					}
				}
			}
		}

		if (!out.writeBlock(vv, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}



std::vector<std::vector<double>> SpatRaster::is_in_cells(std::vector<double> m, SpatOptions &opt) {

	std::vector<std::vector<double>> out(nlyr());

	if (m.empty()) {
		return(out);
	}
	if (!hasValues()) {
		return(out);
	}
	bool hasNAN = false;
	for (size_t i=0; i<m.size(); i++) {
		if (std::isnan(m[i])) {
			hasNAN = true;
			m.erase(m.begin()+i);
			break;
		}
	}
//	if (m.size() == 0) { // only NA
		//nanOnly=true;
//	}

	if (!readStart()) {
		return(out);
	}

	BlockSize bs = getBlockSize(opt);
	size_t nc = ncol();
	for (size_t i = 0; i < bs.n; i++) {
		std::vector<double> v;
		readBlock(v, bs, i);
		size_t cellperlayer = bs.nrows[i] * nc;
		for (size_t j=0; j<v.size(); j++) {
			size_t lyr = j / cellperlayer;
			size_t cell = j % cellperlayer + bs.row[i] * nc;
			if (std::isnan(v[j])) {
				if (hasNAN)	out[lyr].push_back(cell);
			} else {
				for (size_t k=0; k<m.size(); k++) {
					if (v[j] == m[k]) {
						out[lyr].push_back(cell);
						break;
					}
				}
			}
		}
	}
	readStop();
	return(out);
}




SpatRaster SpatRaster::stretch(std::vector<double> minv, std::vector<double> maxv, std::vector<double> minq, std::vector<double> maxq, std::vector<double> smin, std::vector<double> smax, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return(out);

	size_t nl = nlyr();
	recycle(minv, nl);
	recycle(maxv, nl);
	recycle(minq, nl);
	recycle(maxq, nl);
	recycle(smin, nl);
	recycle(smax, nl);

	std::vector<std::vector<double>> q(nl);
	std::vector<bool> useS(nl, false);
	std::vector<double> mult(nl);

	for (size_t i=0; i<nl; i++) {
		if (minv[i] >= maxv[i]) {
			out.setError("maxv must be larger than minv");
			return out;
		}
		if ((!std::isnan(smin[i])) && (!std::isnan(smax[i]))) {
			if (smin[i] >= smax[i]) {
				out.setError("smax must be larger than smin");
				return out;
			}
			useS[i] = true;
			q[i] = {smin[i], smax[i]};
		} else {
			if (minq[i] >= maxq[i]) {
				out.setError("maxq must be larger than minq");
				return out;
			}
			if ((minq[i] < 0) || (maxq[i] > 1)) {
				out.setError("minq and maxq must be between 0 and 1");
				return out;
			}
		}
	}

	std::vector<bool> hR = hasRange();
	for (size_t i=0; i<nl; i++) {
		if (!useS[i]) {
			if ((minq[i]==0) && (maxq[i]==1) && hR[i]) {
				std::vector<double> rmn = range_min();
				std::vector<double> rmx = range_max();
				q[i] = {rmn[i], rmx[i]};
			} else {
				std::vector<double> probs = {minq[i], maxq[i]};
				SpatOptions xopt(opt);
				std::vector<double> v = getValues(i, xopt);
				q[i] = vquantile(v, probs, true);
			}
		}
		mult[i] = maxv[i] / (q[i][1]-q[i][0]);
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		size_t nc = out.bs.nrows[i] * ncol();
		for (size_t j=0; j<v.size(); j++) {
			size_t lyr = j / nc;
			v[j] = mult[lyr] * (v[j] - q[lyr][0]);
			if (v[j] < minv[lyr]) v[j] = minv[lyr];
			if (v[j] > maxv[lyr]) v[j] = maxv[lyr];
		}
		if (!out.writeBlock(v, i)) return out;
	}
	readStop();
	out.writeStop();

	return(out);
}



SpatRaster SpatRaster::apply(std::vector<unsigned> ind, std::string fun, bool narm, std::vector<std::string> nms, std::vector<int_64> time, std::string timestep, std::string timezone, SpatOptions &opt) {

	recycle(ind, nlyr());
	std::vector<unsigned> ui = vunique(ind);
	unsigned nl = ui.size();
	SpatRaster out = geometry(nl);
	recycle(nms, nl);
	out.setNames(nms);
	if (!time.empty()) {
		recycle(time, nl);
		if (!out.setTime(time, timestep, timezone)) {
			out.addWarning("could not set time");
		}
	}

	if (!haveFun(fun)) {
		out.setError("unknown function argument");
		return out;
	}

	if (!hasValues()) return(out);

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	out.bs = getBlockSize(opt);
//	#ifdef useRcpp
//	out.pbar = new Progress(out.bs.n+2, opt.show_progress(bs.n));
//	out.pbar->increment();
//	#endif

	std::vector<std::vector<double>> v(nl);
	std::vector<unsigned> ird(ind.size());
	std::vector<unsigned> jrd(ind.size());
	for (size_t i=0; i<nl; i++) {
		for (size_t j=0; j<ind.size(); j++) {
			if (ui[i] == ind[j]) {
				v[i].push_back(0);
				ird[j] = i;
				jrd[j] = v[i].size()-1;
			}
		}
	}

	std::function<double(std::vector<double>&, bool)> theFun = getFun(fun);

	for (size_t i=0; i<out.bs.n; i++) {
        std::vector<double> a;
		readBlock(a, out.bs, i);
		unsigned nc = out.bs.nrows[i] * ncol();
		std::vector<double> b(nc * nl);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<ird.size(); k++) {
				v[ird[k]][jrd[k]] = a[j+k*nc];
			}
			for (size_t k=0; k<ui.size(); k++) {
				size_t off = k * nc + j;
				b[off] = theFun(v[k], narm);
			}
		}
		if (!out.writeBlock(b, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::mask(SpatRaster &x, bool inverse, double maskvalue, double updatevalue, SpatOptions &opt) {

	unsigned nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl, true, true, true);

	if (!hasValues()) {
		out.setError("raster has no values");
		return out;
	}
	if (!x.hasValues()) {
		out.setError("mask raster has no values");
		return out;
	}

	if (!out.compare_geom(x, false, true, opt.get_tolerance(), true, true, true, false)) {
		return(out);
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}
	std::vector<int> vt = getValueType(true);
	if (vt.size() == 1) {
		out.setValueType(vt[0]);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		readValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		x.readValues(m, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(v, m);
		if (inverse) {
			if (std::isnan(maskvalue)) {
				for (size_t j=0; j < v.size(); j++) {
					if (!std::isnan(m[j])) {
						v[j] = updatevalue;
					}
				}
			} else {
				for (size_t j=0; j < v.size(); j++) {
					if (m[j] != maskvalue) {
						v[j] = updatevalue;
					}
				}
			}
		} else {
			if (std::isnan(maskvalue)) {
				for (size_t j=0; j < v.size(); j++) {
					if (std::isnan(m[j])) {
						v[j] = updatevalue;
					}
				}
			} else {
				for (size_t j=0; j < v.size(); j++) {
					if (m[j] == maskvalue) {
						v[j] = updatevalue;
					}
				}
			}
		}
		if (!out.writeBlock(v, i)) return out;

	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}



SpatRaster SpatRaster::mask(SpatRaster &x, bool inverse, std::vector<double> maskvalues, double updatevalue, SpatOptions &opt) {

	maskvalues = vunique(maskvalues);
	if (maskvalues.size() == 1) {
		return mask(x, inverse, maskvalues[0], updatevalue, opt);
	}

	unsigned nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl, true);

	if (!hasValues()) {
		out.setError("raster has no values");
		return out;
	}
	if (!x.hasValues()) {
		out.setError("mask raster has no values");
		return out;
	}

	if (!out.compare_geom(x, false, true, opt.get_tolerance(), true, true, true, false)) {
		return(out);
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}

	bool maskNA = false;
	for (int i = maskvalues.size()-1; i>=0; i--) {
		if (std::isnan(maskvalues[i])) {
			maskNA = true;
			maskvalues.erase(maskvalues.begin()+i);
		}
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		readValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		x.readValues(m, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(v, m);
		if (inverse) {
			for (size_t i=0; i < v.size(); i++) {
				if (maskNA && std::isnan(m[i])) {
					v[i] = updatevalue;
				} else {
					bool found = false;
					for (size_t j=0; j < maskvalues.size(); j++) {
						if (m[i] == maskvalues[j]) {
							found = true;
							break;
						}
					}
					if (!found) v[i] = updatevalue;
				}
			}
		} else {
			for (size_t i=0; i < v.size(); i++) {
				if (maskNA && std::isnan(m[i])) {
					v[i] = updatevalue;
				} else {
					for (size_t j=0; j < maskvalues.size(); j++) {
						if (m[i] == maskvalues[j]) {
							v[i] = updatevalue;
							break;
						}
					}
				}
			}
		}
		if (!out.writeBlock(v, i)) return out;

	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}

SpatRaster SpatRaster::mask(SpatOptions &opt) {
	SpatRaster out = geometry();

    if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	size_t nl = nlyr();
	size_t nc = ncol();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v;
		std::vector<bool> w;
		readBlock(v, out.bs, i);
		size_t off = out.bs.nrows[i] * nc;
		w.resize(off, false);
		for (size_t j=0; j<off; j++) {
			for (size_t k=0; k<nl; k++) {
				size_t cell = j + k * off;
				if (std::isnan(v[cell])) {
					w[j] = true;
					continue;
				}
			}
		}
		std::vector<size_t> koff;
		koff.reserve(nl);
		for (size_t k=0; k<nl; k++) {
			koff.push_back(( size_t)k * off );
		}
		for (size_t j=0; j<w.size(); j++) {
			if (w[j]) {
				for (size_t k=0; k<nl; k++) {
					v[j+koff[k]] = NAN;
				}
			}
		}
		if (!out.writeBlock(v, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);

}

SpatRaster SpatRaster::mask(SpatVector &x, bool inverse, double updatevalue, bool touches, SpatOptions &opt) {

	SpatRaster out;
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}
		
	if (inverse) {
		out = rasterize(x, "", {updatevalue}, NAN, touches, "", false, true, true, opt);
	} else {
		SpatOptions topt(opt);
		out = rasterize(x, "", {1.0}, 0, touches, "", false, false, false, topt);
		if (out.hasError()) {
			return out;
		}
		if (std::isnan(updatevalue)) {
			out = mask(out, false, 0, updatevalue, opt);
		} else {
			out = mask(out, false, 0, updatevalue, topt);
			out = out.mask(*this, false, NAN, NAN, opt);
		}
	}

	if (!source[0].srs.is_equal(x.srs)) {
		out.addWarning("CRS do not match");
	}

	return(out);
}




SpatRaster SpatRaster::transpose(SpatOptions &opt) {

	SpatRaster out = geometry(nlyr(), true);
	SpatExtent eold = getExtent();
	SpatExtent enew = getExtent();
	enew.xmin = eold.ymin;
	enew.xmax = eold.ymax;
	enew.ymin = eold.xmin;
	enew.ymax = eold.xmax;
	out.setExtent(enew, false, true, "");
	out.source[0].ncol = nrow();
	out.source[0].nrow = ncol();
	if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i=0; i < out.bs.n; i++) {
		unsigned nr = nrow();
		unsigned nc = out.bs.nrows[i];
		std::vector<double> v;
		readValues(v, 0, nr, out.bs.row[i], nc);
		std::vector<double> vv(v.size());
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			size_t off = lyr*ncell();
			for (size_t r = 0; r < nr; r++) {
				size_t rnc = off + r * nc;
				for (size_t c = 0; c < nc; c++) {
					vv[c*nr+r+off] = v[rnc+c];
				}
			}
		}
		if (!out.writeBlock(vv, i)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}





SpatRaster SpatRaster::trim1(double value, unsigned padding, SpatOptions &opt) {

	long nrl = nrow() * nlyr();
	long ncl = ncol() * nlyr();

	size_t r;
	size_t nr = nrow();
	bool rowfound = false;
	if (!readStart()) {
		SpatRaster out;
		out.setError(getError());
		return(out);
	}

	std::vector<double> v;
	size_t firstrow, lastrow, firstcol, lastcol;
	if (std::isnan(value)) {
		for (r=0; r<nr; r++) {
			readValues(v, r, 1, 0, ncol());
			if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
				rowfound = true;
				break;
			}
		}

		if (!rowfound) {
			SpatRaster out;
			out.setError("only cells with NA found");
			return out;
		}

		firstrow = std::max(r - padding, size_t(0));

		for (r=nrow()-1; r>firstrow; r--) {
			readValues(v, r, 1, 0, ncol());
			if (std::count_if(v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
				break;
			}
		}

		lastrow = std::max(std::min(r+padding, nrow()), size_t(0));

		if (lastrow < firstrow) {
			std::swap(firstrow, lastrow);
		}
		size_t c;
		for (c=0; c<ncol(); c++) {
			readValues(v, 0, nrow(), c, 1);
			if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
				break;
			}
		}
		firstcol = std::min(std::max(c-padding, size_t(0)), ncol());

		for (c=ncol()-1; c>firstcol; c--) {
			readValues(v, 0, nrow(), c, 1);
			if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
				break;
			}
		}
		lastcol = std::max(std::min(c+padding, ncol()), size_t(0));
	} else {
		for (r=0; r<nr; r++) {
			readValues(v, r, 1, 0, ncol());
			if (std::count( v.begin(), v.end(), value) < ncl) {
				rowfound = true;
				break;
			}
		}

		if (!rowfound) {
			SpatRaster out;
			out.setError("only cells with value: " + std::to_string(value) + " found");
			return out;
		}

		firstrow = std::max(r - padding, size_t(0));

		for (r=nrow()-1; r>firstrow; r--) {
			readValues(v, r, 1, 0, ncol());
			if (std::count( v.begin(), v.end(), value) < ncl) {
				break;
			}
		}

		lastrow = std::max(std::min(r+padding, nrow()), size_t(0));

		if (lastrow < firstrow) {
			std::swap(firstrow, lastrow);
		}
		size_t c;
		for (c=0; c<ncol(); c++) {
			readValues(v, 0, nrow(), c, 1);
			if (std::count( v.begin(), v.end(), value) < nrl) {
				break;
			}
		}
		firstcol = std::min(std::max(c-padding, size_t(0)), ncol());


		for (c=ncol()-1; c>firstcol; c--) {
			readValues(v, 0, nrow(), c, 1);
			if (std::count( v.begin(), v.end(), value) < nrl) {
				break;
			}
		}
		lastcol = std::max(std::min(c+padding, ncol()), size_t(0));

	}
	readStop();
	if (lastcol < firstcol) {
		std::swap(firstcol, lastcol);
	}

	std::vector<double> res = resolution();
	double xr = res[0];
	double yr = res[1];
	SpatExtent e = SpatExtent(xFromCol(firstcol)-0.5*xr, xFromCol(lastcol)+0.5*xr, yFromRow(lastrow)-0.5*yr, yFromRow(firstrow)+0.5*yr);

	return( crop(e, "near", false, opt) ) ;
}


void block_cols(const std::vector<double> &v, std::function<bool(double,double)> fun, const double &value, size_t &firstcol, size_t &lastcol, bool &firstcoldone, bool &lastcoldone, const size_t &firstrow, const size_t &lastrow, const size_t &nr, const size_t &nc, const size_t &nl, const size_t &padding) {

	size_t maxcol = nc - padding - 1;

	std::vector<size_t> loff(nl);
	for (size_t i=0; i<nl; i++) {
		loff[i] = i * nr * nc;
	}

	if (!firstcoldone) {
		for (size_t r=firstrow; r<lastrow; r++) {
			size_t roff = r * nc;
			for (size_t lyr=0; lyr<nl; lyr++) {
				for (size_t c=0; c<firstcol; c++) {
					if (fun(v[loff[lyr] + roff + c], value)) {
						firstcol = c;
						if (firstcol <= padding) firstcoldone = true;
						break;
					}
				}
				if (firstcoldone) break;
			}
			if (firstcoldone) break;
		}
	}
	if (!lastcoldone) {
		for (size_t r=firstrow; r<lastrow; r++) {
			size_t roff = r * nc;
			for (size_t lyr=0; lyr<nl; lyr++) {
				for (long c=(nc-1); c>=(long)lastcol; c--) {
					if (fun(v[loff[lyr] + roff + c], value)) {
						lastcol = c;
						if (lastcol >= maxcol) lastcoldone = true;
						break;
					}
				}
				if (lastcoldone) break;
			}
			if (lastcoldone) break;
		}
	}
}


inline bool trim_value(const double &x, const double &y) {
	return x != y;	
}

inline bool trim_nan(const double &x, const double &y) {
	return !std::isnan(x);	
}

SpatRaster SpatRaster::trim2(double value, unsigned padding, SpatOptions &opt) {

// check if opt.filename exists and overwrite=false?

	if (!readStart()) {
		SpatRaster out;
		out.setError(getError());
		return(out);
	}

	std::vector<double> v;
	BlockSize bs = getBlockSize(opt);
	size_t nl = nlyr();
	size_t nc = ncol();
	size_t nr = nrow();
	bool firstrowfound = false;
	bool lastrowfound = false;
	bool firstcolfound = false;
	bool lastcolfound = false;


	size_t bstart = 0;
	size_t bend = bs.n - 1;
	
	size_t firstrow = nr-1;
	size_t lastrow = 0;
	size_t firstcol = nc-1; 
	size_t lastcol = 0;

	if (padding >= nc) {
		firstcolfound = true;
		lastcolfound = true;
		firstcol = 0;
		lastcol = nc-1;
	}
	if (padding >= nr) {
		firstrowfound = true;
		lastrowfound = true;
		firstrow = 0;
		lastrow = nr-1;
	}

	std::function<bool(double,double)> fun;

	if (std::isnan(value)) {
		fun = trim_nan;
	} else {
		fun = trim_value;
	}
	
	bool rowfound = false;
	for (size_t i=0; i<bs.n; i++) {	
		if (firstrowfound) {
			rowfound = true;
			break;
		}
		bstart = i+1;
		readBlock(v, bs, i);
		std::vector<size_t> loff(nl);
		for (size_t j=0; j<nl; j++) {
			loff[j] = j * bs.nrows[i] * nc;
		}
		for (size_t r=0; r<bs.nrows[i]; r++) {
			size_t roff = r * nc;
			for (size_t lyr=0; lyr<nl; lyr++) {
				for (size_t j=0; j<nc; j++) {
					if (fun(v[loff[lyr] + roff + j], value)) {
						rowfound = true;
						firstrow = r; 
						break;
					}
				}
				if (rowfound) break;
			}
			if (rowfound) break;
		}
		if (rowfound) {
			block_cols(v, fun, value, firstcol, lastcol, firstcolfound, lastcolfound, firstrow, bs.nrows[i], bs.nrows[i], nc, nl, padding);
			break;
		}		
	}
	firstrow += bs.row[bstart-1];
			
	if (!rowfound) {
		SpatRaster out;
		out.setError("only cells with NA found");
		return out;
	}
	if (!rowfound) {
		SpatRaster out;
		out.setError("only cells with NA found");
		return out;
	}
	
	lastrow = firstrow;
	rowfound = false;
	if (!lastrowfound) {
		if (bstart == bs.n) { // no need to read v again
			bend = bstart;
			size_t i = bs.n - 1;
			std::vector<size_t> loff(nl);
			for (size_t j=0; j<nl; j++) {
				loff[j] = j * bs.nrows[i] * nc;
			}
			for (long r=(bs.nrows[i]-1); r>=0; r--) {
				size_t roff = r * nc;
				for (size_t lyr=0; lyr<nl; lyr++) {
					for (size_t j=0; j<nc; j++) {
						if (fun(v[loff[lyr] + roff + j], value)) {
							rowfound = true;
							lastrow = r;
							break;
						}
					}
					if (rowfound) break;
				}
				if (rowfound) break;
			}
			block_cols(v, fun, value, firstcol, lastcol, firstcolfound, lastcolfound, 0, lastrow, bs.nrows[i], nc, nl, padding);
			lastrow += bs.row[i];

		} else { // read blocks from bottom

			for (long i=(bs.n-1); i>=0; i--) {
				bend = i;
				readBlock(v, bs, i);
				std::vector<size_t> loff(nl);
				for (size_t j=0; j<nl; j++) {
					loff[j] = j * bs.nrows[i] * nc;
				}
				for (long r=(bs.nrows[i]-1); r>=0; r--) {
					size_t roff = r * nc;
					for (size_t lyr=0; lyr<nl; lyr++) {
						for (size_t j=0; j<nc; j++) {
							if (fun(v[loff[lyr] + roff + j], value)) {
								rowfound = true;
								lastrow = r;
								break;
							}
						}
						if (rowfound) break;
					}
					if (rowfound) break;
				}

				if (rowfound) {
					block_cols(v, fun, value, firstcol, lastcol, firstcolfound, lastcolfound, 0, lastrow, bs.nrows[i], nc, nl, padding);
					break;
				}
			}
			lastrow += bs.row[bend];
		}
	}
	for (size_t i=bstart; i<bend; i++) {
		if (firstcolfound && lastcolfound) break;
		readBlock(v, bs, i);
		block_cols(v, fun, value, firstcol, lastcol, firstcolfound, lastcolfound, 0, bs.nrows[i], bs.nrows[i], nc, nl, padding);
	}
	firstrow = std::max(firstrow - padding, size_t(0));
	lastrow = std::max(std::min(lastrow + padding, nr), size_t(0));
	if (lastrow < firstrow) {
		std::swap(firstrow, lastrow);
	}
	firstcol = std::min(std::max(firstcol-padding, size_t(0)), nc);
	lastcol = std::max(std::min(lastcol+padding, nc), size_t(0));
	if (lastcol < firstcol) {
		std::swap(firstcol, lastcol);
	}
	
	readStop();
	std::vector<double> res = resolution();
	double xr = res[0];
	double yr = res[1];
	SpatExtent e = SpatExtent(xFromCol(firstcol)-0.5*xr, xFromCol(lastcol)+0.5*xr, yFromRow(lastrow)-0.5*yr, yFromRow(firstrow)+0.5*yr);

	return( crop(e, "near", false, opt) ) ;
}



void clamp_vector(std::vector<double> &v, double low, double high, bool usevalue) {
	size_t n = v.size();
	if (usevalue) {
		for (size_t i=0; i<n; i++) {
			if ( v[i] < low ) {
				v[i] = low;
			} else if ( v[i] > high ) {
				v[i] = high;
			}
		}
	} else {
		for (size_t i=0; i<n; i++) {
			if ( (v[i] < low ) || (v[i] > high)) {
				v[i] = NAN;
			}
		}
	}
}



SpatRaster SpatRaster::clamp(std::vector<double> low, std::vector<double> high, bool usevalue, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr(), true);
	if (!hasValues()) {
		out.setError("cannot clamp a raster with no values");
		return out;
	}
	if (low.empty() || high.empty()) {
		out.setError("you must provide low and high clamp values");
		return out;
	}
	size_t nl = nlyr();
	if ((low.size() > nl) || (high.size() > nl)) {
		out.setError("there are more low and/or high values than layers");
		return out;
	}
	bool do_one = true;
	if ((low.size() > 1) || (high.size() > 1)) {
		do_one = false;
		recycle(low, nl);
		recycle(high, nl);
	}
	for (size_t i=0; i<low.size(); i++) {
		if (low[i] > high[i]) {
			out.setError("lower clamp value cannot be larger than the higher clamp value");
			return out;
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	if (do_one) {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			clamp_vector(v, low[0], high[0], usevalue);
			if (!out.writeBlock(v, i)) return out;
		}
	} else {
		size_t nc = ncol();
		for (size_t i = 0; i < out.bs.n; i++) {
			size_t off = out.bs.nrows[i] * nc;
			std::vector<double> v;
			readBlock(v, out.bs, i);
			if (usevalue) {
				for (size_t j=0; j<nl; j++) {
					size_t start = j * off;
					size_t end = start + off;
					for (size_t k=start; k<end; k++) {
						if (v[k] < low[j] ) {
							v[k] = low[j];
						} else if ( v[k] > high[j] ) {
							v[k] = high[j];
						}
					}
				}
			} else {
				for (size_t j=0; j<nl; j++) {
					size_t start = j * off;
					size_t end = start + off;
					for (size_t k=start; k<end; k++) {
						if ((v[k] < low[j] ) || (v[k] > high[j])) {
							v[k] = NAN;
						}
					}
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}

SpatRaster SpatRaster::clamp_raster(SpatRaster &x, SpatRaster &y, std::vector<double> low, std::vector<double> high, bool usevalue, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr(), true);
	if (!hasValues()) {
		out.setError("cannot clamp a raster with no values");
		return out;
	}
	size_t nl = nlyr();
	bool do_one = true;
	bool rA = false;
	bool rB = false;
	bool onex = true;
	bool oney = true;
	if (std::isnan(low[0])) {
		rA = true;
		if (!x.hasValues()) {
			out.setError("cannot clamp with raster that has no values");
			return out;
		}
		if (x.nlyr() > 1) {
			if (x.nlyr() != nl) {
				out.setError("clamp raster must have one layer or the same number of layers as x");
				return out;			
			} else {
				onex = false;
			}
		}
	} else {
		if (low.size() > nl) {
			out.setError("there are more low values than layers");
			return out;
		}
	}

	if (std::isnan(high[0])) {
		rB = true;
		if (!y.hasValues()) {
			out.setError("cannot clamp with raster that has no values");
			return out;
		}
		if (y.nlyr() > 1) {
			if (y.nlyr() != nl) {
				out.setError("clamp raster must have one layer or the same number of layers as x");
				return out;			
			} else {
				oney = false;
			}
		}
	} else {
		if (high.size() > nl) {
			out.setError("there are more high values than layers");
			return out;
		}
	}
	
	if ((low.size() > 1) || (high.size() > 1) || rA || rB) {
		do_one = false;
		recycle(low, nl);
		recycle(high, nl);
	}
	if (!(rA | rB)) {
		for (size_t i=0; i<low.size(); i++) {
			if (low[i] > high[i]) {
				out.setError("lower clamp value cannot be larger than the higher clamp value");
				return out;
			}
		}
	}
	
	if (rA) {
		if (!x.readStart()) {
			out.setError(x.getError());
			return(out);
		}		
	}
	if (rB) {
		if (!y.readStart()) {
			out.setError(y.getError());
			return(out);
		}		
	}
	
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
	opt.ncopies = (1 + oney + onex) * opt.ncopies; 
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	if (!(rA | rB)) {
		if (do_one) {
			for (size_t i = 0; i < out.bs.n; i++) {
				std::vector<double> v;
				readBlock(v, out.bs, i);
				clamp_vector(v, low[0], high[0], usevalue);
				if (!out.writeBlock(v, i)) return out;
			}
		} else {
			size_t nc = ncol();
			for (size_t i = 0; i < out.bs.n; i++) {
				size_t off = out.bs.nrows[i] * nc;
				std::vector<double> v;
				readBlock(v, out.bs, i);
				if (usevalue) {
					for (size_t j=0; j<nl; j++) {
						size_t start = j * off;
						size_t end = start + off;
						for (size_t k=start; k<end; k++) {
							if (v[k] < low[j] ) {
								v[k] = low[j];
							} else if ( v[k] > high[j] ) {
								v[k] = high[j];
							}
						}
					}
				} else {
					for (size_t j=0; j<nl; j++) {
						size_t start = j * off;
						size_t end = start + off;
						for (size_t k=start; k<end; k++) {
							if ((v[k] < low[j] ) || (v[k] > high[j])) {
								v[k] = NAN;
							}
						}
					}
				}
				if (!out.writeBlock(v, i)) return out;
			}
		}
	} else if (rA & rB) {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v, vx, vy;
			readBlock(v, out.bs, i);
			x.readBlock(vx, out.bs, i);
			y.readBlock(vy, out.bs, i);
			size_t ncl = vx.size();
			if (usevalue) {
				for (size_t j=0; j<v.size(); j++) {
					size_t kx = onex ? j % ncl : j;
					size_t ky = oney ? j % ncl : j;
					if (v[j] < vx[kx] ) {
						v[j] = vx[kx];
					} else if ( v[j] > vy[ky] ) {
						v[j] = vy[ky];
					}
				}
			} else {
				for (size_t j=0; j<v.size(); j++) {
					size_t kx = onex ? j % ncl : j;
					size_t ky = oney ? j % ncl : j;
					if (v[j] < vx[kx] ) {
						v[j] = NAN;
					} else if ( v[j] > vy[ky] ) {
						v[j] = NAN;
					}
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	} else if (rA) {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v, vx;
			readBlock(v, out.bs, i);
			x.readBlock(vx, out.bs, i);
			size_t ncl = vx.size();
			if (usevalue) {
				for (size_t j=0; j<v.size(); j++) {
					size_t k = onex ? j % ncl : j;
					size_t lyr = j / ncl;
					if (v[j] < vx[k] ) {
						v[j] = vx[k];
					} else if ( v[j] > high[lyr] ) {
						v[j] = high[lyr];
					}
				}
			} else {
				for (size_t j=0; j<v.size(); j++) {
					size_t k = onex ? j % ncl : j;
					size_t lyr = j / ncl;
					if (v[j] < vx[k] ) {
						v[j] = NAN;
					} else if ( v[j] > high[lyr] ) {
						v[j] = NAN;
					}
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	} else if (rB) {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v, vy;
			readBlock(v, out.bs, i);
			y.readBlock(vy, out.bs, i);
			size_t ncl = vy.size();
			if (usevalue) {
				for (size_t j=0; j<v.size(); j++) {
					size_t k = oney ? j % ncl : j;
					size_t lyr = j / ncl;
					if (v[j] < low[lyr]) {
						v[j] = low[lyr];
					} else if (v[j] > vy[k]) {
						v[j] = vy[k];
					}
				}
			} else {
				for (size_t j=0; j<v.size(); j++) {
					size_t k = oney ? j % ncl : j;
					size_t lyr = j / ncl;
					if (v[j] < low[lyr]) {
						v[j] = NAN;
					} else if (v[j] > vy[k]) {
						v[j] = NAN;
					}
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	}

	readStop();
	if (rA) x.readStop();
	if (rB) y.readStop();
	
	out.writeStop();	
	return(out);
}



std::vector<double> bip2bil(const std::vector<double> &v, size_t nl) {
	
	size_t n = v.size();
	size_t ncell = n / nl;
	
	std::vector<double> out(n);
	std::vector<size_t> offlyr(nl);
	for (size_t j=0; j<nl; j++) {
		offlyr[j] = j * ncell;
	}
	
	for (size_t i=0; i<ncell; i++) {
		size_t off = i * nl;
		for (size_t j=0; j<nl; j++) {
			out[offlyr[j] + i] = v[off + j];
		}
	}
	return out;
}



SpatRaster SpatRaster::clamp_ts(bool min, bool max, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr(), true);
	if (!hasValues()) {
		out.setError("cannot clamp a raster with no values");
		return out;
	}
	if (!(min || max)) {
		out.setError("min or max must be TRUE");
		return(out);		
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	size_t nl = nlyr();
	size_t nc = ncol();
	for (size_t i=0; i<out.bs.n; i++) {
		size_t ncells = out.bs.nrows[i] * nc;
		std::vector<double> v;
		readBlockIP(v, out.bs, i);
		
		for (size_t j=0; j<ncells; j++) {
			size_t start = j * nl;
			size_t end = start + nl;
			if (min) {
				double minv = min_se_rm(v, start, end);
				double wmin = whichmin_se_rm(v, start, end);
				for (size_t k=start; k<(start+wmin); k++) {
					v[k] = minv;
				}
			}
			if (max) {
				double maxv = max_se_rm(v, start, end);
				double wmax = whichmax_se_rm(v, start, end);
				for (size_t k=(start+wmax); k<end; k++) {
					v[k] = maxv;
				}
			}
		}
		v = bip2bil(v, nl);
		if (!out.writeBlock(v, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::selRange(SpatRaster x, int z, int recycleby, SpatOptions &opt) {

	int nl = nlyr();
	z = std::max(1, std::min(z, nl));
	size_t nrecs = 1;
	if (recycleby > 1 && recycleby < nl) {
		nrecs = nl / recycleby;
	} else {
		recycleby = 0;
	}
	SpatRaster out = geometry( z * nrecs );
	if (!out.compare_geom(x, false, false, opt.get_tolerance())) {
		return(out);
	}
	if (!hasValues()) return(out);

	if (x.nlyr() > 1) {
		out.setError("index raster must have only one layer");
		return out;
	}
	if (!x.hasValues()) {
		out.setError("index raster has no values");
		return out;
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v, idx;
		readBlock(v, out.bs, i);
		x.readBlock(idx, out.bs, i);
		size_t is = idx.size();
		std::vector<double> vv(is*z*nrecs, NAN);
		size_t ncell = out.bs.nrows[i] * ncol(); // same as is?

		for (size_t j=0; j<is; j++) {  //index.size (each cell)
			for (size_t k=0; k<nrecs; k++) {
				int start = idx[j] - 1 + k * recycleby;  // first layer for this cell
				int offbase = (k*z) * ncell;
				if ((start >= 0) && (start < nl)) {
					int zz = std::min(nl-start, z); // do not surpass the last layer
					for (int i=0; i<zz; i++){
						size_t offin = (start + i) * ncell + j;
						size_t offout = offbase + i * ncell + j;
						vv[offout] = v[offin];
					}
				}
			}
		}
		//for (size_t j=0; j<is; j++) {
		//	int index = idx[j] - 1;
		//	if ((index >= 0) && (index < nl)) {
		//		vv[j] = v[j + index * ncell];
		//	}
		//}
		if (!out.writeBlock(vv, i)) return out;
	}
	readStop();
	x.readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::roll(size_t n, std::string fun, std::string type, bool circular, bool narm, SpatOptions &opt) {
	
	SpatRaster out = geometry();
	if (!hasValues()) {
		out.setError("no values in input");
		return(out);
	}
	if (!haveFun(fun)) {
		out.setError("unknown function argument");
		return out;
	}
	if (n >= nlyr()) {
		out.setError("it makes no sense to use a rolling function with n >= nlyr(x)");
		return out;		
	}
	if (n <= 1) {
		out.setError("n should be > 1");
		return out;		
	}
	std::vector<std::string> types = {"around", "to", "from"};
	if (!is_in_vector(type, types)) {
		out.setError("unknown roll type, should be 'around', 'to', or 'from'");
		return out;					
	}

	// to do: use functions that iterate over vector instead of copying
	
	std::function<double(std::vector<double>&, bool)> theFun = getFun(fun);

	size_t nl = nlyr();
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (circular) {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlockIP(v, out.bs, i);
			size_t ncell = out.bs.nrows[i] * ncol();
			std::vector<double> vv(v.size(), NAN);
			if (type=="from") {
				for (size_t j=0; j<ncell; j++) {
					size_t offset = j*nl;
					for (size_t k=0; k<nl; k++) {
						std::vector<double> se;
						size_t start = offset + k;
						size_t end = k + n;
						if (end > nl) {
							size_t cend = end - nl;
							se = {v.begin()+offset, v.begin()+offset+cend};
							end = nl;
						}
						end += offset;
						se.insert(se.end(), v.begin()+start, v.begin()+end);
						vv[ncell * k + j] = theFun(se, narm);
					}
				}
			} else if (type=="around") {
				size_t halfn = n / 2;
				for (size_t j=0; j<ncell; j++) {
					size_t offset = j*nl;
					for (size_t k=0; k<nl; k++) {
						std::vector<double> se;
						size_t start, end;
						if (k < halfn) {
							start = 0;
							end = n + k - halfn;
							size_t cbegin = nl - (halfn - k);
							se = {v.begin()+offset+cbegin, v.begin()+offset+nl};
						} else {
							start = k - halfn;
							end = start + n;
						}
						if (end > nl) {
							end = nl;
							size_t cend = end - nl + 1;
							se = {v.begin()+offset, v.begin()+offset+cend};
						}
						start += offset;
						end += offset;
						se.insert(se.end(), v.begin()+start, v.begin()+end);
						vv[ncell * k + j] = theFun(se, narm);
					}
				}
			} else if (type=="to") {
				for (size_t j=0; j<ncell; j++) {
					size_t offset = j*nl;
					for (size_t k=0; k<nl; k++) {
						std::vector<double> se;
						size_t start;
						size_t end = offset + k + 1;
						if (k < (n-1)) {
							start = offset;
							size_t cbegin = nl - (n - k - 1);
							se = {v.begin()+offset+cbegin, v.begin()+offset+nl};
						} else {
							start = end - n;
						}
						se.insert(se.end(), v.begin()+start, v.begin()+end);
						vv[ncell * k + j] = theFun(se, narm);
					}
				}
			}	
			if (!out.writeBlock(vv, i)) return out;
		}
	} else { // not circular
		std::vector<double> se;
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlockIP(v, out.bs, i);
			size_t ncell = out.bs.nrows[i] * ncol();
			std::vector<double> vv(v.size(), NAN);
			if (type=="from") {
				for (size_t j=0; j<ncell; j++) {
					size_t offset = j*nl;
					for (size_t k=0; k<nl; k++) {
						size_t start = offset + k;
						size_t end = k + n;
						if (end > nl) {
							if (narm) {
								end = nl;
							} else {
								continue;
							}
						}
						end += offset;
						se = {v.begin()+start, v.begin()+end};
						vv[ncell * k + j] = theFun(se, narm);
					}
				}
			} else if (type=="around") {
				size_t halfn = n / 2;
				for (size_t j=0; j<ncell; j++) {
					size_t offset = j*nl;
					for (size_t k=0; k<nl; k++) {
						size_t start, end;
						if (k < halfn) {
							if (narm) {
								start = 0;
								end = n + k - halfn;
							} else {
								continue;	
							}
						} else {
							start = k - halfn;
							end = start + n;
						}
						if (end > nl) {
							if (narm) {
								end = nl;
							} else {
								continue;
							}
						}
						start += offset;
						end += offset;
						se = {v.begin()+start, v.begin()+end};
						vv[ncell * k + j] = theFun(se, narm);
					}
				}
			} else if (type=="to") {
				for (size_t j=0; j<ncell; j++) {
					size_t offset = j*nl;
					for (size_t k=0; k<nl; k++) {
						size_t start;
						size_t end = offset + k + 1;
						if (k < (n-1)) {
							if (narm) {
								start = offset;						
							} else {
								continue;
							}
						} else {
							start = end - n;
						}
						se = {v.begin()+start, v.begin()+end};
						vv[ncell * k + j] = theFun(se, narm);
					}
				}
			}	
			if (!out.writeBlock(vv, i)) return out;
		}
	}
	readStop();
	out.writeStop();	
	return out;	
}


SpatRaster SpatRaster::rapply(SpatRaster x, double first, double last, std::string fun, bool clamp, bool narm, bool circular, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!haveFun(fun)) {
		out.setError("unknown function argument");
		return out;
	}

	bool sval = !std::isnan(first);
	bool eval = !std::isnan(last);
	if (sval && eval) {
		out.setError("arguments `first` or `last` must be NA. See `app` for other cases");
		return out;
	}
	int start = sval ? first-1 : -99;
	int end = eval ? last-1 : -999;

	if (!out.compare_geom(x, false, false, opt.get_tolerance())) {
		return(out);
	}
	if (!x.hasValues()) {
		out.setError("index raster has no values");
		return out;
	}
	unsigned expnl = 2 - (sval + eval);
	if (x.nlyr() != expnl) {
		out.setError("index raster must have " + std::to_string(expnl) + "layer(s)");
		return out;
	}
	if (!hasValues()) {
		out.setError("no values in input");
		return(out);
	}

	std::function<double(std::vector<double>&, bool)> theFun = getFun(fun);

	int nl = nlyr();
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}

	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v, idx;
		readBlock(v, out.bs, i);
		x.readBlock(idx, out.bs, i);
		size_t ncell = out.bs.nrows[i] * ncol();
		std::vector<double> vv(ncell, NAN);
		for (size_t j=0; j<ncell; j++) {
			if (std::isnan(idx[j])) continue;
			if (sval) {
				end = idx[j] - 1;
			} else if (eval) {
				start = idx[j] - 1;
			} else {
				start = idx[j] - 1;
				double dend = idx[j+ncell]-1;
				if (std::isnan(dend)) continue;
				end   = dend;
			}

			if (clamp) {
				start = start < 0 ? 0 : start;
				end = end >= nl ? (nl-1) : end;
				if (circular) {
					end = end < 0 ? 0 : end;
					start = start >= nl ? (nl-1) : start;
				}
			}

			bool inrange = (start < nl) && (end < nl) && (start >= 0) && (end >= 0);
			bool circ = false;
			if (start > end) {
				if (circular) {
					circ = true;
				} else {
					inrange = false;
				}
			}

			if (inrange) {
				std::vector<double> se;
				if (circ) {
					se.reserve(end + nl - start + 1);
					for (int k = start; k<nl; k++){
						size_t off = k * ncell + j;
						se.push_back(v[off]);
					}
					for (int k = 0; k<=end; k++){
						size_t off = k * ncell + j;
						se.push_back(v[off]);
					}
					vv[j] = theFun(se, narm);
				} else {
					se.reserve(end-start+1);
					for (int k = start; k<=end; k++){
						size_t off = k * ncell + j;
						se.push_back(v[off]);
					}
					vv[j] = theFun(se, narm);
				}
			}
		}
		if (!out.writeBlock(vv, i)) return out;
	}
	readStop();
	x.readStop();
	out.writeStop();
	return(out);
}


std::vector<std::vector<double>> SpatRaster::rappvals(SpatRaster x, double first, double last, bool clamp, bool all, double fill, size_t startrow, size_t nrows, bool circular) {

	std::vector<std::vector<double>> r;

	bool sval = !std::isnan(first);
	bool eval = !std::isnan(last);
	if (sval && eval) {
		setError("first or last must be NA. See `app` for other cases");
		return r;
	}
	int start = sval ? first-1 : 0;
	int end = eval ? last-1 : 0;

	if (!compare_geom(x, false, false, 0.1)) {
		return(r);
	}
	if (!hasValues()) {
		return r;
	}
	if (!x.hasValues()) {
		setError("index raster has no values");
		return r;
	}
	unsigned expnl = 2 - (sval + eval);
	if (x.nlyr() != expnl) {
		setError("index raster must have " + std::to_string(expnl) + "layer(s)");
		return r;
	}

	int nl = nlyr();
	if (!readStart()) {
		return(r);
	}
	if (!x.readStart()) {
		setError(x.getError());
		return(r);
	}

	std::vector<double> v, idx;
	readValues(v, startrow, nrows, 0, ncol());
	x.readValues(idx, startrow, nrows, 0, ncol());
	size_t ncell = nrows * ncol();
	r.resize(ncell);

	for (size_t j=0; j<ncell; j++) {
		if (std::isnan(idx[j])) {
			if (all) {
				r[j].resize(nl, NAN);
			} else {
				r[j].push_back(NAN);
			}
			continue;
		}
		if (sval) {
			end = idx[j] - 1;
			//end = idx[j];
		} else if (eval) {
			start = idx[j] - 1;
		} else {
			start = idx[j] - 1;
			//double dend = idx[j+ncell];
			double dend = idx[j+ncell]-1;
			end = std::isnan(dend) ? -999 : (int) dend;
		}
		if (clamp) {
			start = start < 0 ? 0 : start;
			end = end >= nl ? (nl-1) : end;
			if (circular) {
				end = end < 0 ? 0 : end;
				start = start >= nl ? (nl-1) : start;
			}
		}

		bool inrange = (start < nl) && (end < nl) && (start >= 0) && (end >= 0);
		bool circ = false;
		if (start > end) {
			if (circular) {
				circ = true;
			} else {
				inrange = false;
			}
		}

		if (all) {
			if (inrange) {
				r[j].resize(nl, fill);
				if (circ) {
					for (int k=start; k<nl; k++){
						size_t off = k * ncell + j;
						r[j][k] = v[off];
					}
					for (int k=0; k<=end; k++){
						size_t off = k * ncell + j;
						r[j][k] = v[off];
					}
				} else {
					for (int k = start; k<=end; k++){
						size_t off = k * ncell + j;
						r[j][k] = v[off];
					}
				}
			} else {
				r[j].resize(nl, NAN);
			}
		} else if (inrange) {
			if (circ) {
				r[j].reserve(end + (nl-start) + 1);
				for (int k=start; k<nl; k++){
					size_t off = k * ncell + j;
					r[j].push_back(v[off]);
				}
				for (int k=0; k<=start; k++){
					size_t off = k * ncell + j;
					r[j].push_back(v[off]);
				}
			} else {
				r[j].reserve(end-start+1);
				for (int k=start; k<=end; k++){
					size_t off = k * ncell + j;
					r[j].push_back(v[off]);
				}
			}
		} else {
			r[j].push_back(NAN);
		}

	}
	readStop();
	x.readStop();
	return(r);
}


bool disaggregate_dims(std::vector<unsigned> &fact, std::string &message ) {

	unsigned fs = fact.size();
	if ((fs > 3) | (fs == 0)) {
		message = "argument 'fact' should have length 1, 2, or 3";
		return false;
	}
	auto min_value = *std::min_element(fact.begin(),fact.end());
	if (min_value < 1) {
		message = "values in argument 'fact' should be > 0";
		return false;
	}
	auto max_value = *std::max_element(fact.begin(),fact.end());
	if (max_value == 1) {
		message = "all values in argument 'fact' are 1, nothing to do";
		return false;
	}

	fact.resize(3);
	if (fs == 1) {
		fact[1] = fact[0];
	}
	fact[2] = 1;
	return true;
}



SpatRaster SpatRaster::disaggregate(std::vector<unsigned> fact, SpatOptions &opt) {

    SpatRaster out = geometry(nlyr(), true);
	std::string message = "";
	if (!disaggregate_dims(fact, message)) {
		if (message.substr(0,3) == "all") {
			out = *this;
			out.addWarning(message);
		} else {
			out.setError(message);
		}
		return out;
	}

    out.source[0].nrow = out.source[0].nrow * fact[0];
    out.source[0].ncol = out.source[0].ncol * fact[1];
    out.source[0].nlyr = out.source[0].nlyr * fact[2];

    if (!hasValues()) {
        return out;
    }

	opt.ncopies = 2*fact[0]*fact[1]*fact[2];
	BlockSize bs = getBlockSize(opt);
	opt.steps = bs.n;
	//opt.set_blocksizemp();
	unsigned nc = ncol();
	unsigned nl = nlyr();
	std::vector<double> newrow(nc*fact[1]);
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i = 0; i < bs.n; i++) {
		std::vector<double> v, vout;
		readValues(v, bs.row[i], bs.nrows[i], 0, nc);
		vout.reserve(v.size() * fact[0] * fact[1] * fact[2]);

		for (size_t lyr=0; lyr<nl; lyr++) {
			for (size_t row=0; row<bs.nrows[i]; row++) {
				unsigned rowoff = row*nc + lyr*nc*bs.nrows[i];
				// for each new column
				unsigned jfact = 0;
				for (size_t j=0; j<nc; j++) {
					unsigned coloff = rowoff + j;
					for (size_t k=0; k<fact[1]; k++) {
						newrow[jfact+k] = v[coloff];
					}
					jfact += fact[1];
				}
				// for each new row
				for (size_t j=0; j<fact[0]; j++) {
					vout.insert(vout.end(), newrow.begin(), newrow.end());
				}
			}
		}
		if (!out.writeValues(vout, bs.row[i]*fact[0], bs.nrows[i]*fact[0])) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::init(std::string value, bool plusone, SpatOptions &opt) {

	SpatRaster out = geometry(1);

	std::vector<std::string> f {"row", "col", "cell", "x", "y", "chess"};
	bool test = std::find(f.begin(), f.end(), value) == f.end();
	if (test) {
		out.setError("not a valid init option");
		return out;
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	size_t nc = ncol();
	std::vector<double> v;
	if (value == "row") {
		for (size_t i = 0; i < out.bs.n; i++) {
			v.resize(nc * out.bs.nrows[i]);
			for (size_t j = 0; j < out.bs.nrows[i]; j++) {
				size_t r = out.bs.row[i] + j + plusone;
				for (size_t k = 0; k < nc; k++) {
					v[j*nc+k] = r;
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
		//source[0].range_min.resize(1, 0 + plusone);
		//source[0].range_max.resize(1, nrow() - 1 + plusone);
		//source[0].hasRange.resize(1, true);
	} else if (value == "col") {
		std::vector<double> cnn(nc);
		double start = plusone ? 1 : 0;
		std::iota(cnn.begin(), cnn.end(), start);
		size_t oldnr = 0;
		for (size_t i = 0; i < out.bs.n; i++) {
			if (oldnr != out.bs.nrows[i]) {
				v = cnn;
				recycle(v, out.bs.nrows[i] * nc);
			}
			if (!out.writeBlock(v, i)) return out;
		}
		//source[0].range_min.resize(1, 0 + plusone);
		//source[0].range_max.resize(1, nc - 1 + plusone);
		//source[0].hasRange.resize(1, true);
	} else if (value == "cell") {
		for (size_t i = 0; i < out.bs.n; i++) {
			v.resize(nc * out.bs.nrows[i]);
			size_t firstcell = cellFromRowCol(out.bs.row[i], 0);
			firstcell = plusone ? firstcell + 1 : firstcell;
			std::iota(v.begin(), v.end(), firstcell);
			if (!out.writeBlock(v, i)) return out;
		}
		//source[0].range_min.resize(1, 0 + plusone);
		//source[0].range_max.resize(1, ncell() - 1 + plusone);
		//source[0].hasRange.resize(1, true);

	} else if (value == "x") {
		std::vector<int_64> col(nc);
		std::iota(col.begin(), col.end(), 0);
		std::vector<double> xcoords = xFromCol(col);
		size_t oldnr = 0;
		for (size_t i = 0; i < out.bs.n; i++) {
			if (oldnr != out.bs.nrows[i]) {
				v = xcoords;
				recycle(v, out.bs.nrows[i] * nc);
			}
			if (!out.writeBlock(v, i)) return out;
		}
		//source[0].range_min.resize(1, xcoords[0]);
		//source[0].range_max.resize(1, xcoords[nc-1]);
		//source[0].hasRange.resize(1, true);

	} else if (value == "y") {

		for (size_t i = 0; i < out.bs.n; i++) {
			v.resize(out.bs.nrows[i] * nc );
			for (size_t j = 0; j < out.bs.nrows[i]; j++) {
				double y = yFromRow(out.bs.row[i] + j);
				for (size_t k = 0; k < nc; k++) {
					v[j*nc+k] = y;
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
		//source[0].range_min.resize(1, yFromRow(0));
		//source[0].range_max.resize(1, yFromRow(nrow()-1));
		//source[0].hasRange.resize(1, true);

	} else if (value == "chess") {
		std::vector<double> a(nc);
		std::vector<double> b(nc);
		for (size_t i=0; i<nc; i++) {
			bool even = i%2 == 0;
			a[i] = even;
			b[i] = !even;
		}
		std::vector<double> v;
		for (size_t i = 0; i < out.bs.n; i++) {
			if ((out.bs.row[i]%2) == 0) {
				v = a;
				v.insert(v.end(), b.begin(), b.end());
			} else {
				v = b;
				v.insert(v.end(), b.begin(), b.end());
			}
			recycle(v, out.bs.nrows[i] * nc);
			if (!out.writeBlock(v, i)) return out;
		}
		//source[0].range_min.resize(1, 0);
		//source[0].range_max.resize(1, 1);
		//source[0].hasRange.resize(1, true);
	}

	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::init(std::vector<double> values, SpatOptions &opt) {
	SpatRaster out = geometry();
 	if (!out.writeStart(opt, filenames())) { return out; }
	unsigned nc = ncol();
	unsigned nl = nlyr();
	if (values.size() == 1) {
		double val = values[0];
		std::vector<double> v;
		for (size_t i = 0; i < out.bs.n; i++) {
			v.resize(out.bs.nrows[i]*nc*nl, val);
			if (!out.writeBlock(v, i)) return out;
		}
		//source[0].range_min.resize(1, val);
		//source[0].range_max.resize(1, val);
		//source[0].hasRange.resize(1, true);

	} else {
		int over = 0;
		for (size_t i = 0; i < out.bs.n; i++) {
			if (over > 0) {
				std::vector<double> newv(values.begin()+over, values.end());
				newv.insert(newv.end(), values.begin(), values.begin()+over);
				values = newv;
			}
			std::vector<double> v = values;
			recycle(v, out.bs.nrows[i]*nc);
			recycle(v, out.bs.nrows[i]*nc*nl);
			over = v.size() % values.size();
			if (!out.writeBlock(v, i)) return out;
		}
	}
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::rotate(bool left, SpatOptions &opt) {

	unsigned nc = ncol();
	unsigned nl = nlyr();
	unsigned hnc = (nc / 2);
	double addx = hnc * xres();
	if (left) {
		addx = -addx;
	}
	SpatRaster out = geometry(nlyr(), true, true, true);
	SpatExtent outext = out.getExtent();
	outext.xmin = outext.xmin + addx;
	outext.xmax = outext.xmax + addx;
	out.setExtent(outext, true, true, "");

	if (!hasValues()) return out;

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i=0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		std::vector<double> b;
		b.reserve(a.size());
		for (size_t j=0; j < nl; j++) {
			for (size_t r=0; r < out.bs.nrows[i]; r++) {
				unsigned s1 = j * out.bs.nrows[i] * nc + r * nc;
				unsigned e1 = s1 + hnc;
				b.insert(b.end(), a.begin()+e1, a.begin()+s1+nc);
				b.insert(b.end(), a.begin()+s1, a.begin()+e1);
			}
		}
		if (!out.writeBlock(b, i)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}



bool SpatRaster::shared_basegeom(SpatRaster &x, double tol, bool test_overlap) {
	if (!compare_origin(x.origin(), tol)) return false;
	if (!about_equal(xres(), x.xres(), xres() * tol)) return false;
	if (!about_equal(yres(), x.yres(), yres() * tol)) return false;
	if (test_overlap) {
		SpatExtent e = x.getExtent();
		e = e.intersect(getExtent());
		if (!e.valid()) return false;
	}
	return true;
}





SpatRaster SpatRaster::cover(SpatRaster x, std::vector<double> values, SpatOptions &opt) {

	unsigned nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl, true, true, true);

	bool rmatch = false;
	if (out.compare_geom(x, false, false, opt.get_tolerance(), true)) {
		rmatch = true;
	} else {
	//	if (!shared_basegeom(x, 0.1, true)) {
		out.setError("raster dimensions do not match");
		return(out);
	//	} else {
	//		out.msg.has_error = false;
	//		out.msg.error = "";
	//		SpatExtent e = getExtent();
	//		SpatExtent xe = x.getExtent();
	//		double prec = std::min(xres(), yres())/1000;
	//		if (!xe.compare(e, "<=", prec)) {
	//			SpatOptions xopt(opt);
	//			x = x.crop(e, "near", xopt);
	//		}
	//	}
	}


	if (!x.hasValues()) {
		return *this;
	}
	if (!hasValues()) {
		if (rmatch) {
			return x.deepCopy();
		} else {
			SpatExtent e = getExtent();
			return x.extend(e, "near", NAN, opt);
		}
	}

	std::vector<bool> cats = hasCategories();
	std::vector<bool> xcats = x.hasCategories();
	recycle(cats, xcats);
	for (size_t i =0; i<nl; i++) {
		if (cats[i]) {
			if (xcats[i]) {
				SpatCategories sc = getLayerCategories(i);
				SpatCategories xsc = x.getLayerCategories(i);
				if (sc.combine(xsc)) {
					out.source[0].cats[i] = sc;
				} else {
					std::string warn = "cannot merge categories of layer " + std::to_string(i+1);
					out.addWarning(warn);
					out.removeCategories(i);
				}
			}
		} else if (xcats[i]) {
			out.removeCategories(i);
		}
	}


	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		x.readStop();
		return out;
	}
	if (values.size() == 1) {
		double value=values[0];
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v, m;
			readValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol());
			x.readValues(m, out.bs.row[i], out.bs.nrows[i], 0, ncol());
			recycle(v, m);
			if (std::isnan(value)) {
				for (size_t j=0; j < v.size(); j++) {
					if (std::isnan(v[j])) {
						v[j] = m[j];
					}
				}
			} else {
				for (size_t j=0; j < v.size(); j++) {
					if (v[j] == value) {
						v[j] = m[j];
					}
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	} else {

		values = vunique(values);
		bool hasNA = false;
		for (int i = values.size()-1; i>=0; i--) {
			if (std::isnan(values[i])) {
				hasNA = true;
				values.erase(values.begin()+i);
			}
		}

		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v, m;
			readValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol());
			x.readValues(m, out.bs.row[i], out.bs.nrows[i], 0, ncol());
			recycle(v, m);
			for (size_t j=0; j < v.size(); j++) {
				if (hasNA) {
					if (std::isnan(v[j])) {
						v[j] = m[j];
						continue;
					}
				}
				for (size_t j=0; j<values.size(); j++) {
					if (v[j] == values[j]) {
						v[j] = m[j];
						continue;
					}
				}
			}
			if (!out.writeBlock(v, i)) return out;
		}
	}

	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}




SpatRaster SpatRaster::extend(SpatExtent e, std::string snap, double fill, SpatOptions &opt) {

	SpatRaster out = geometry_opt(nlyr(), true, true, true, true, opt);
	e = out.align(e, snap);
	SpatExtent extent = getExtent();
	e.unite(extent);

	out.setExtent(e, true, true, "");
	if (!hasValues() ) {
		if (!opt.get_filename().empty()) {
			out.addWarning("ignoring filename argument because there are no cell values");
		}
		return(out);
	}

	double tol = std::min(xres(), yres()) / 1000;
	if (extent.compare(e, "==", tol)) {
		// same extent
		if (opt.get_filename().empty()) {
			out = deepCopy();
		} else {
			out = writeRaster(opt);
		}
		return out;
	}


	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	out.fill(fill);
	BlockSize bs = getBlockSize(opt);
	for (size_t i=0; i<bs.n; i++) {
        std::vector<double> v;
		readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
        unsigned row1 = out.rowFromY(yFromRow(bs.row[i]));
        unsigned row2 = out.rowFromY(yFromRow(bs.row[i]+bs.nrows[i]-1));
        unsigned col1 = out.colFromX(xFromCol(0));
        unsigned col2 = out.colFromX(xFromCol(ncol()-1));
        if (!out.writeValuesRect(v, row1, row2-row1+1, col1, col2-col1+1)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::crop(SpatExtent e, std::string snap, bool expand, SpatOptions &opt) {

	SpatRaster out = geometry_opt(nlyr(), true, true, true, true, opt);

	if ( !e.valid() ) {
		out.setError("invalid extent");
		return out;
	}
	if ((e.xmin == e.xmax) && (e.ymin == e.ymax)) {
		out.setError("cannot crop a SpatRaster with an empty extent");
		return out;
	}
	SpatExtent ein = getExtent();
	SpatExtent fext = e;	
	e = e.intersect(ein);
	if ( !e.valid_notempty() ) {
		out.setError("extents do not overlap");
		return out;
	}

	SpatOptions ops;
	if (expand) {
		if ((fext.xmax <= ein.xmax)  && (fext.xmin >= ein.xmin) && (fext.ymax <= ein.ymax)  && (fext.ymin >= ein.ymin)) {
			expand = false;
		} else if ((fext.xmax >= ein.xmax)  && (fext.xmin <= ein.xmin) && (fext.ymax >= ein.ymax)  && (fext.ymin <= ein.ymin)) {
			return extend(fext, snap, NAN, opt);
		} else {
			ops = opt;
			opt = SpatOptions(opt);
		}
	}

	out.setExtent(e, true, false, snap);

	if (!hasValues() ) {
		if (expand) {
			if (!ops.get_filename().empty()) {
				out.addWarning("ignoring filename argument because there are no cell values");
			}
			out = out.extend(fext, snap, NAN, opt);
		} else {
			if (!opt.get_filename().empty()) {
				out.addWarning("ignoring filename argument because there are no cell values");
			}
		}
		return(out);
	}

	double hxr = xres() / 2;
	double hyr = yres() / 2;
	SpatExtent outext = out.getExtent();
	unsigned col1 = colFromX(outext.xmin + hxr);
	unsigned col2 = colFromX(outext.xmax - hxr);
	unsigned row1 = rowFromY(outext.ymax - hyr);
	unsigned row2 = rowFromY(outext.ymin + hyr);

	std::vector<bool> hw = hasWindow();
	bool haswin = hw[0];
	for (size_t i=1; i<nsrc(); i++) {
		haswin = (haswin | hw[i]);
	}

	if ((row1==0) && (row2==nrow()-1) && (col1==0) && (col2==ncol()-1) && (!haswin)) {
		// same extent
		if (opt.get_filename().empty()) {
			out = deepCopy();
		} else {
			out = writeRaster(opt);
		}
		return out;
	}

	std::vector<int> vt = getValueType(true);
	if (vt.size() == 1) {
		out.setValueType(vt[0]);
	}

	unsigned ncols = out.ncol();
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

//	opt.ncopies = 2;
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> v;
	for (size_t i = 0; i < out.bs.n; i++) {
		readValues(v, row1+out.bs.row[i], out.bs.nrows[i], col1, ncols);
		if (!out.writeBlock(v, i)) return out;
	}
	out.writeStop();
	readStop();

	if (expand) {
		out = out.extend(fext, snap, NAN, ops);
	}
	return(out);
}


SpatRaster SpatRaster::cropmask(SpatVector &v, std::string snap, bool touches, bool extend, SpatOptions &opt) {
	if (v.nrow() == 0) {
		SpatRaster out;
		out.setError("cannot crop a SpatRaster with an empty SpatVector");
		return out;
	}

	SpatOptions copt(opt);
	SpatRaster out = crop(v.extent, snap, extend, copt);
	if (out.hasError()) return out;
	// transfer warnings?
	return out.mask(v, false, NAN, touches, opt);
}



SpatRaster SpatRaster::flip(bool vertical, SpatOptions &opt) {

	SpatRaster out = geometry_opt(nlyr(), true, true, true, true, opt);
	if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	unsigned nc = ncol();
	unsigned nl = nlyr();

	if (vertical) {
		for (size_t i=0; i < out.bs.n; i++) {
			std::vector<double> a, b;
			size_t ii = out.bs.n - 1 - i;
			readBlock(a, out.bs, ii);
			b.reserve(a.size());
			for (size_t j=0; j < out.nlyr(); j++) {
				size_t offset = j * out.bs.nrows[ii] * nc;
				for (size_t k=0; k < out.bs.nrows[ii]; k++) {
					unsigned start = offset + (out.bs.nrows[ii] - 1 - k) * nc;
					b.insert(b.end(), a.begin()+start, a.begin()+start+nc);
				}
			}
			if (!out.writeBlock(b, i)) return out;
		}
	} else {
		for (size_t i=0; i < out.bs.n; i++) {
			std::vector<double> a, b;
			readBlock(a, out.bs, i);
			b.reserve(a.size());
			unsigned lyrrows = nl * out.bs.nrows[i];
			for (size_t j=0; j < lyrrows; j++) {
				unsigned start = j * nc;
				unsigned end = start + nc;
				std::vector<double> v(a.begin()+start, a.begin()+end);
				std::reverse(v.begin(), v.end());
				b.insert(b.end(), v.begin(), v.end());
			}
			if (!out.writeBlock(b, i)) return out;
		}
	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::reverse(SpatOptions &opt) {

	SpatRaster out = geometry_opt(nlyr(), true, true, true, true, opt);
	if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	unsigned nc = ncol();
	unsigned nl = nlyr();

	for (size_t i=0; i < out.bs.n; i++) {
		size_t ii = out.bs.n - 1 - i;
		std::vector<double> a, b;
		readBlock(a, out.bs, ii);
		b.reserve(a.size());
		unsigned lyrrows = nl * out.bs.nrows[ii];
		for (size_t j=0; j < lyrrows; j++) {
			unsigned start = (lyrrows - 1 - j) * nc;
			unsigned end = start + nc;
			std::vector<double> v(a.begin()+start, a.begin()+end);
			std::reverse(v.begin(), v.end());
			b.insert(b.end(), v.begin(), v.end());
		}
		if (!out.writeBlock(b, i)) return out;
	}

	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::shift(double x, double y, SpatOptions &opt) {
	SpatRaster out = deepCopy();
	SpatExtent outext = out.getExtent();
	outext.xmin = outext.xmin + x;
	outext.xmax = outext.xmax + x;
	outext.ymin = outext.ymin + y;
	outext.ymax = outext.ymax + y;
	out.setExtent(outext, true, true, "");
	return out;
}

bool SpatRaster::compare_origin(std::vector<double> x, double tol) {
	std::vector<double> y = origin();
	if (!about_equal(x[0], y[0], xres() * tol)) return false;
	if (!about_equal(x[1], y[1], yres() * tol)) return false;
	return true;
}




bool write_part(SpatRaster& out, SpatRaster& r, const double& hxr, unsigned& nl, bool notfirstlyr, bool warn, SpatOptions &opt) {
	BlockSize bs = r.getBlockSize(opt);
	if (!r.readStart()) {
	out.setError(r.getError());
		return false;
	}
	SpatExtent re = r.getExtent();
	if (!r.shared_basegeom(out, 0.1, true)) {
		SpatRaster temp = out.crop(re, "near", false, opt);
		std::vector<bool> hascats = r.hasCategories();
		std::string method = hascats[0] ? "near" : "bilinear";
		r = r.warper(temp, "", method, false, false, true, opt);
		if (r.hasError()) {
			out.setError(r.getError());
			return false;
		}
		warn = true;
		re = r.getExtent();
	}

	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v, vout;
		r.readBlock(v, bs, i);
		unsigned row1  = out.rowFromY(r.yFromRow(bs.row[i])); 
		unsigned row2  = out.rowFromY(r.yFromRow(bs.row[i]+bs.nrows[i]-1));
		unsigned col1  = out.colFromX(re.xmin + hxr);
		unsigned col2  = out.colFromX(re.xmax - hxr);
		unsigned ncols = col2-col1+1;
		unsigned nrows = row2-row1+1;
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


SpatRaster SpatRasterCollection::merge(bool first, bool narm, SpatOptions &opt) {

	SpatRaster out;
	unsigned n = size();
	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	if (n == 1) {
		out = ds[0].deepCopy();
		return(out);
	}

	std::vector<bool> hvals(n);
	hvals[0] = ds[0].hasValues();
	SpatExtent e = ds[0].getExtent();
	unsigned nl = ds[0].nlyr();
	for (size_t i=1; i<n; i++) {
		//  lyrs, crs, warncrs, ext, rowcol, res
		if (!ds[0].compare_geom(ds[i], false, false, opt.get_tolerance(), false, false, false, true)) {
			out.setError(ds[0].msg.error);
			return(out);
		}
		e.unite(ds[i].getExtent());
		hvals[i] = ds[i].hasValues();
		nl = std::max(nl, ds[i].nlyr());
	}

	out = ds[0].geometry(nl, true);
	//out = ds[0].geometry(1, false);
	out.setExtent(e, true, true, "");

	bool anyvals = false;
	for (size_t i=0; i<n; i++) {
		if (hvals[i]) {
			anyvals = true;
			break;
		}
	}
	if (!anyvals) {
		return out;
	}

	//out = out.geometry(nl, true);
	double hxr = out.xres()/2;

	std::vector<int> vt = getValueType(true);
		if (vt.size() == 1) {
		out.setValueType(vt[0]);
	}

 	if (!out.writeStart(opt, filenames())) { return out; }

	std::vector<size_t> seq(n);
	if (first) {
		std::iota(seq.rbegin(), seq.rend(), 0);
	} else {
		std::iota(seq.begin(), seq.end(), 0);
	}

	SpatOptions topt(opt);
	bool warn = false;
	bool notfirst = false;
	for (size_t i=0; i<n; i++) {
		if (!ds[seq[i]].hasValues()) continue;
		if (narm) {
			notfirst = i > 0;
		}
		if (!write_part(out, ds[seq[i]], hxr, nl, notfirst, warn, topt)) {
			return out;
		}
	}
	out.writeStop();
	if (warn) out.addWarning("rasters did not align and were resampled");

	return(out);
}



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


SpatRaster SpatRasterCollection::mosaic(std::string fun, SpatOptions &opt) {

	SpatRaster out;
	std::vector<std::string> f {"first", "last", "sum", "mean", "median", "min", "max"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("argument 'fun' is not a valid function name");
		return out;
	}
	if (fun == "first") {
		return merge(true, true, opt);
	}
	if (fun == "last") {
		return merge(false, true, opt);
	}
	unsigned n = size();

	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	if (n == 1) {
		out = ds[0].deepCopy();
		return(out);
	}

	std::vector<bool> hvals(n);
	hvals[0] = ds[0].hasValues();
	SpatExtent e = ds[0].getExtent();
	unsigned nl = ds[0].nlyr();
//std::vector<bool> resample(n, false);
		
	for (size_t i=1; i<n; i++) {
		SpatExtent ee = ds[i].getExtent();
									//  lyrs, crs, warncrs, ext, rowcol, res
		if (!ds[0].compare_geom(ds[i], false, false, opt.get_tolerance(), false, false, false, true)) {
			out.setError(ds[0].msg.error);
			return(out);
		}
		e.unite(ee);
		hvals[i] = ds[i].hasValues();
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
		return out;
	}

//	if (!overlaps(r1, r2, c1, c2)) {
//		return merge(true, true, opt);
//	}

	double ncl = 1000;
	if (n > 50) ncl = 500;
	if (n > 100) ncl = 250;
	double     ar = std::ceil(out.nrow() / ncl);
	unsigned arow = std::ceil(out.nrow() / ar);
	double     ac = std::ceil(out.ncol() / ncl);
	unsigned acol = std::ceil(out.ncol() / ac);

	SpatOptions sopt(opt);
	SpatRaster aout = out.aggregate({arow, acol}, "", true, sopt);
	SpatVector ve = aout.as_polygons(false, false, false, false, false, 0, sopt);

	SpatVector vcrp(out.getExtent(), "");
	ve = ve.intersect(vcrp, false);
	n = ve.nrow();
	bool warn = false;

 	if (!out.writeStart(opt, filenames())) { return out; }
	sopt.progressbar = false;

	std::vector<unsigned> use; 
	SpatRasterStack s;
	for (size_t i=0; i<n; i++) {
		SpatVector vi = ve.subset_rows(i);
		SpatExtent ce = vi.getExtent();
		SpatRasterCollection x = crop(ce, "near", true, use, sopt);
		if (x.empty()) {
			continue;
		} 
		s.ds = x.ds;
			//r = s.summary(fun, true, sopt);
// see #1159
//			if (i == 57 || i == 79 | i == 269) { // && (rcnt[i] == 6)) {
		SpatRaster r = s.collapse();
		r = r.summary(fun, true, sopt);	
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


/*
SpatRaster SpatRasterCollection::mosaic(std::string fun, SpatOptions &opt) {

	SpatRaster out;
	std::vector<std::string> f {"first", "last", "sum", "mean", "median", "min", "max"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("argument 'fun' is not a valid function name");
		return out;
	}
	if (fun == "first") {
		return merge(true, true, opt);
	}
	if (fun == "last") {
		return merge(false, true, opt);
	}
	unsigned n = size();

	if (n == 0) {
		out.setError("empty collection");
		return(out);
	}
	if (n == 1) {
		out = ds[0].deepCopy();
		return(out);
	}

	std::vector<bool> hvals(n);
	hvals[0] = ds[0].hasValues();
	SpatExtent e = ds[0].getExtent();
	unsigned nl = ds[0].nlyr();
	std::vector<bool> resample(n, false);
	for (size_t i=1; i<n; i++) {
									//  lyrs, crs, warncrs, ext, rowcol, res
		if (!ds[0].compare_geom(ds[i], false, false, opt.get_tolerance(), false, false, false, true)) {
			out.setError(ds[0].msg.error);
			return(out);
		}
		e.unite(ds[i].getExtent());
		hvals[i] = ds[i].hasValues();
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
		return out;
	}

//	SpatExtent eout = out.getExtent();
	double hxr = out.xres()/2;
	double hyr = out.yres()/2;
	std::vector<unsigned> r1, r2, c1, c2;
	r1.reserve(n); r2.reserve(n);
	c1.reserve(n); c2.reserve(n);
	SpatVector ve;
	ve.reserve(n);
	for (size_t i=0; i<n; i++) {
		SpatExtent ee = ds[i].getExtent();
		r1.push_back(out.rowFromY(ee.ymax - hyr));
		r2.push_back(out.rowFromY(ee.ymin + hyr));
		c1.push_back(out.colFromX(ee.xmin + hxr));
		c2.push_back(out.colFromX(ee.xmax - hxr));
		SpatVector v(ee, "");
		ve.geoms.push_back(v.geoms[0]);
	}

	if (!overlaps(r1, r2, c1, c2)) {
		return merge(true, true, opt);
	}

	ve = ve.unite();
	ve = ve.disaggregate(false);
	n = ve.nrow();
	std::vector<std::vector<unsigned>> rsti(n);
	for (size_t i=0; i<ve.ncol(); i++) {
		for (size_t j=0; j<n; j++) {
			if (ve.df.iv[i][j] == 1) {
				rsti[j].push_back(i);
			}
		}
	}
	std::vector<size_t> rcnt(n);
	for (size_t i=0; i<n; i++) {
		rcnt[i] = rsti[i].size();
	}
	std::vector<std::size_t> ord = sort_order_a(rcnt);
	permute(rcnt, ord);
	permute(rsti, ord);

	bool warn = false;
 	if (!out.writeStart(opt, filenames())) { return out; }
	SpatOptions sopt(opt);
	sopt.progressbar = false;

	for (size_t i=0; i<n; i++) {
		SpatRaster r;

		if (rcnt[i] == 1) {
			r = ds[rsti[i][0]];
		} else if (rcnt[i] > 1) {
			SpatVector vi = ve.subset_rows(ord[i]);

			
			SpatRasterCollection x = crop(vi.extent, "near", true, rsti[i], sopt);
			if (x.empty()) {
				continue;
			} 
			SpatRasterStack s;
			s.ds = x.ds;
			//r = s.summary(fun, true, sopt);
// see #1159
//			if (i == 57 || i == 79 | i == 269) { // && (rcnt[i] == 6)) {
			r = s.collapse();
			r = r.summary(fun, true, sopt);
			
			if (r.hasError()) {
				return r;
			}
		}
		if (!write_part(out, r, hxr, nl, false, warn, sopt)) {
			return out;
		}
	}
	out.writeStop();

	if (warn) out.addWarning("rasters did not align and were resampled");
	return out;
}
*/

SpatRaster SpatRasterCollection::morph(SpatRaster &x, SpatOptions &opt) {

	SpatRaster out;
	unsigned n = size();
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


void notisnan(const std::vector<double> &x, double &n) {
	for (size_t i=0; i<x.size(); i++) {
		n += !std::isnan(x[i]);
	}
}



void do_stat(std::vector<double> &v, std::string fun, bool narm, double &stat, double &stat2,double &n, size_t step) {
	if (v.empty()) return;
	if (fun == "sum") {
		double s = vsum(v, narm);
		if (step > 0) {
			std::vector<double> ss = {stat, s};
			stat = vsum(ss, narm);
		} else {
			stat = s;
		}
	} else if (fun == "mean") {
		double s = vsum(v, narm);
		if (step > 0) {
			std::vector<double> ss = {stat, s};
			stat = vsum(ss, narm);
		} else {
			stat = s;
		}
		if (narm) {
			notisnan(v, n);
		} else {
			n += v.size();
		}
	} else if (fun == "prod") {
		double p = vprod(v, narm);
		if (step > 0) {
			std::vector<double> pp = {stat, p};
			stat = vprod(pp, narm);
		} else {
			stat = p;
		}
	} else if (fun == "rms") {
		if (narm) {
			notisnan(v, n);
		} else {
			n += v.size();
		}
		double s = vsum2(v, narm);
		if (step > 0) {
			std::vector<double> ss = {stat, s};
			stat = vsum(ss, narm);
		} else {
			stat = s;
		}
	} else if (fun == "min") {
		double s = vmin(v, narm);
		if (step > 0) {
			std::vector<double> ss = {stat, s};
			stat = vmin(ss, narm);
		} else {
			stat = s;
		}
	} else if (fun == "max") {
		double s = vmax(v, narm);
		if (step > 0) {
			std::vector<double> ss = {stat, s};
			stat = vmax(ss, narm);
		} else {
			stat = s;
		}
	} else if (fun == "range") {
		double sn = vmin(v, narm);
		double sx = vmax(v, narm);
		if (step > 0) {
			std::vector<double> ss1 = {stat, sn};
			stat = vmin(ss1, narm);
			std::vector<double> ss2 = {stat2, sx};
			stat2 = vmax(ss2, narm);
		} else {
			stat = sn;
			stat2 = sx;
		}
	} else if (fun == "sd") {
		if (narm) {
			notisnan(v, n);
		} else {
			n += v.size();
		}
		double s1 = vsum(v, narm);
		double s2 = vsum2(v, narm);
		if (step > 0) {
			std::vector<double> ss1 = {stat, s1};
			stat = vsum(ss1, narm);
			std::vector<double> ss2 = {stat2, s2};
			stat2 = vsum(ss2, narm);
		} else {
			stat = s1;
			stat2 = s2;
		}
	} else if (fun == "notNA" || fun == "isNA") {
		notisnan(v, n);
	}
}



void do_mstats(std::vector<double> &v, size_t start, size_t end, std::vector<std::string> funs, bool narm, std::vector<double> &stat, std::vector<double> &stat2, double &n, bool first, bool last) {
	
	size_t nstat = funs.size();
	
	if (first) {
		stat.resize(0);
		stat.resize(nstat);
		stat2.resize(0);
		stat2.resize(nstat);
		n = 0;
	}
	
	if (v.empty()) return;

	double sum = 0;

	if (is_in_vector("sum", funs) || is_in_vector("mean", funs) || 
				is_in_vector("sd", funs) || is_in_vector("std", funs)) {
		if (narm) {
			sum = sum_se_rm(v, start, end);
		} else {
			sum = sum_se(v, start, end);
		}
	}
	size_t notna = 0;
	if (is_in_vector("mean", funs) || is_in_vector("rms", funs) || is_in_vector("sd", funs) || 
			is_in_vector("std", funs) || is_in_vector("notNA", funs) || is_in_vector("isNA", funs)) {
		if (narm) {
			notna = isnotna_se(v, start, end);
			n += notna;
		} else {
			n += (end - start);
		}
	}
	
	for (size_t i=0; i<nstat; i++) {
		std::string fun = funs[i];
		if (fun == "sum") {
			if (first) {
				stat[i] = sum;
			} else {
				std::vector<double> ss = {stat[i], sum};
				stat[i] = vsum(ss, narm);
			}
		} else if (fun == "mean") {
			if (first) {
				stat[i] = sum;
			} else {
				std::vector<double> ss = {stat[i], sum};
				stat[i] = vsum(ss, narm);
			}
			if (last) {
				if (n > 0) {
					stat[i] = stat[i] / n;
				} else {
					stat[i] = NAN;
				}
			}

		} else if (fun == "prod") {
			double p;
			if (narm) {
				 p = prod_se_rm(v, start, end);
			} else {
				p = prod_se(v, start, end);				
			}
			if (first) {
				stat[i] = p;
			} else {
				std::vector<double> pp = {stat[i], p};
				stat[i] = vprod(pp, narm);
			}
		} else if (fun == "rms") {
			double s;
			if (narm) {
				s = sum2_se_rm(v, start, end);
			} else {
				s = sum2_se(v, start, end);
			}
			if (first) {
				stat[i] = s;
			} else {
				std::vector<double> ss = {stat[i], s};
				stat[i] = vsum(ss, narm);
			}

			if (last) {
				// rms = sqrt(sum(x^2)/(n-1))
				if (n > 0) {
					stat[i] = sqrt(stat[i] / (n-1));
				} else {
					stat[i] = NAN;
				}
			}
		} else if (fun == "min") {
			double s;
			if (narm) {
				s = min_se_rm(v, start, end);
			} else {
				s = min_se(v, start, end);
			}
			if (first) {
				stat[i] = s;
			} else {
				std::vector<double> ss = {stat[i], s};
				stat[i] = vmin(ss, narm);
			}
		} else if (fun == "max") {
			double s;
			if (narm) {
				s = max_se_rm(v, start, end);
			} else {
				s = max_se(v, start, end);
			}
			if (first) {
				stat[i] = s;
			} else {
				std::vector<double> ss = {stat[i], s};
				stat[i] = vmax(ss, narm);
			}
		} else if ((fun == "sd") || (fun == "std")) {
			double s2;
			if (narm) {
				s2 = sum2_se_rm(v, start, end);
			} else {
				s2 = sum2_se(v, start, end);
			}
			if (first) {
				stat[i] = sum;
				stat2[i] = s2;
			} else {
				std::vector<double> ss1 = {stat[i], sum};
				stat[i] = vsum(ss1, narm);
				std::vector<double> ss2 = {stat2[i], s2};
				stat2[i] = vsum(ss2, narm);
			}

			if (last) {
				if (n > 0) {
					double mn = stat[i] / n;
					double mnsq = mn * mn;
					double mnsumsq = stat2[i] / n;
					if (fun == "std") {
						stat[i] = sqrt(mnsumsq - mnsq);
					} else {
						stat[i] = sqrt((mnsumsq - mnsq) * n/(n-1));
					}
				} else {
					stat[i] = NAN;
				}
			}
		} else if (fun == "notNA") {
			if (narm) {
				// if (last) {
				stat[i] = n;
			} else {
				stat[i] += isnotna_se(v, start, end);
			}
		} else if (fun == "isNA") {
			if (narm) {
				stat[i] += end - start - notna;
			} else {
				stat[i] += end - start - isnotna_se(v, start, end);
			}
		}
	}
}



SpatDataFrame SpatRaster::mglobal(std::vector<std::string> funs, bool narm, SpatOptions &opt) {

	SpatDataFrame out;
	std::vector<std::string> f {"sum", "mean", "min", "max", "prod", "rms", "sd", "std", "isNA", "notNA"};

	size_t nf = funs.size();
	for (size_t i=0; i<nf; i++) {
		if (std::find(f.begin(), f.end(), funs[i]) == f.end()) {
			out.setError(funs[i] + " is not a valid function");
			return(out);
		}
	}

	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return(out);
	}

	size_t nl = nlyr();
	
	std::vector<std::vector<double>> stats(nl, std::vector<double>(nf));
	std::vector<std::vector<double>> stats2(nl, std::vector<double>(nf));

	std::vector<double> n(nl);
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
	BlockSize bs = getBlockSize(opt);

	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v;
		readBlock(v, bs, i);
		unsigned off = bs.nrows[i] * ncol() ;
		for (size_t lyr=0; lyr<nl; lyr++) {
			unsigned offset = lyr * off;
			//std::vector<double> vv = { v.begin()+offset, v.begin()+offset+off };
			do_mstats(v, offset, (offset+off), funs, narm, stats[lyr], stats2[lyr], n[lyr], i==0, i==(bs.n-1));
		}
	}
	readStop();

	// transpose
	std::vector<std::vector<double>> tstat(nf, std::vector<double>(nl));
	std::vector<std::vector<double>> tstat2(nf, std::vector<double>(nl));
	for (size_t i=0; i<nl; i++) {
		for (size_t j=0; j<nf; j++) {
			tstat[j][i] = stats[i][j]; 
			tstat2[j][i] = stats2[i][j]; 
		}
	}

	for (size_t i=0; i<nf; i++) {
		if (funs[i]=="range") {
			out.add_column(tstat[i], "min");
			out.add_column(tstat2[i], "max");
		} else {
			out.add_column(tstat[i], funs[i]);
		}
	}
	
	return(out);
}



std::vector<std::vector<double>> SpatRaster::layerCor(std::string fun, bool narm, bool asSample, SpatOptions &opt) {

	std::vector<std::vector<double>> out(2);
	
	if (!hasValues()) {
		setError("SpatRaster has no values");
		return(out);
	}

	size_t nl = nlyr();

	if (fun == "pearson") {
		std::vector<double> means(nl*nl, NAN);
		std::vector<double> cor(nl*nl, 1);
		SpatOptions topt(opt);

		BlockSize bs = getBlockSize(topt);
		
		std::vector<std::string> gfuns = {"mean", "sd"};
		for (unsigned i=0; i<(nl-1); i++) {
			for (unsigned j=(i+1); j<nl; j++) {
				SpatRaster xi = subset({i}, topt);
				SpatRaster xj = subset({j}, topt);
				if (!xi.readStart()) {
					setError(getError());
					return(out);
				}
				if (!xj.readStart()) {
					setError(getError());
					return(out);
				}
				std::vector<std::vector<double>> stats(nl);
				std::vector<std::vector<double>> stats2(nl);
				std::vector<double> n(nl);
				std::vector<double> vi, vj;				
				for (size_t k=0; k<bs.n; k++) {
					xi.readBlock(vi, bs, k);
					xj.readBlock(vj, bs, k);
					if (narm) {
						for (size_t m=0; m<vi.size(); m++) {
							if (std::isnan(vi[m]) || std::isnan(vj[m])) {
								vi[m] = NAN;
								vj[m] = NAN;
							} 
						}
					}
					do_mstats(vi, 0, vi.size(), gfuns, narm, stats[0], stats2[0], n[0], k==0, k==(bs.n-1));
					do_mstats(vj, 0, vj.size(), gfuns, narm, stats[1], stats2[1], n[1], k==0, k==(bs.n-1));
				}
				double value = 0;
				for (long kk=bs.n; kk>0; kk--) {
					size_t k = kk-1;
					if (k < (bs.n-1)) {
						xi.readBlock(vi, bs, k);
						xj.readBlock(vj, bs, k);
					}
					if (narm) {
						for (size_t m=0; m<vi.size(); m++) {
							if (!std::isnan(vi[m])) {
								value += (vi[m] - stats[0][0]) * (vj[m]  - stats[1][0]);
							}
						}
					} else {
						for (size_t m=0; m<vi.size(); m++) {
							value += (vi[m] - stats[0][0]) * (vj[m]  - stats[1][0]);
						}
					}
				}
				value /= (n[0] - asSample) * (stats[0][1] * stats[1][1]);
				means[i*nl+j] = stats[0][0];
				means[j*nl+i] = stats[1][0];
				cor[i*nl+j] = value;
				cor[j*nl+i] = value;
				xi.readStop();
				xj.readStop();
			}
		}
		out[0] = cor;	
		out[1] = means;	
	}
	return(out);	
}





SpatDataFrame SpatRaster::global(std::string fun, bool narm, SpatOptions &opt) {

	SpatDataFrame out;
	std::vector<std::string> f {"sum", "mean", "min", "max", "range", "prod", "rms", "sd", "std", "stdpop", "isNA", "notNA"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("not a valid function");
		return(out);
	}

	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return(out);
	}

	std::string sdfun = fun;
	if ((fun == "std") || (fun == "sdpop")) {
		sdfun = "std";
		fun = "sd";
	}
	size_t nl = nlyr();
	std::vector<double> stats(nl);
	std::vector<double> stats2(nl);

	std::vector<double> n(nl);
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	BlockSize bs = getBlockSize(opt);
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v;
		readBlock(v, bs, i);
		unsigned off = bs.nrows[i] * ncol() ;
		for (size_t lyr=0; lyr<nl; lyr++) {
			unsigned offset = lyr * off;
			std::vector<double> vv = { v.begin()+offset, v.begin()+offset+off };
			do_stat(vv, fun, narm, stats[lyr], stats2[lyr], n[lyr], i);
		}
	}
	readStop();


	if (fun=="mean") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			if (n[lyr] > 0) {
				stats[lyr] = stats[lyr] / n[lyr];
			} else {
				stats[lyr] = NAN;
			}
		}
	} else if (fun=="rms") {
		// rms = sqrt(sum(x^2)/(n-1))
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			if (n[lyr] > 0) {
				stats[lyr] = sqrt(stats[lyr] / (n[lyr]-1));
			} else {
				stats[lyr] = NAN;
			}
		}
	} else if (fun == "sd") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			if (n[lyr] > 0) {
				double mn = stats[lyr] / n[lyr];
				double mnsq = mn * mn;
				double mnsumsq = stats2[lyr] / n[lyr];
				if (sdfun == "std") {
					stats[lyr] = sqrt(mnsumsq - mnsq);
				} else {
					stats[lyr] = sqrt((mnsumsq - mnsq) * n[lyr]/(n[lyr]-1));
				}
			} else {
				stats[lyr] = NAN;
			}
		}
	} else if ((fun == "notNA") || (fun == "isNA")) {
		for (size_t lyr=0; lyr<nl; lyr++) {
			stats[lyr] = n[lyr];
		}
	} else if (fun == "isNA") {
		double nc = ncell();
		for (size_t lyr=0; lyr<nl; lyr++) {
			stats[lyr] = nc - n[lyr];
		}
	}
	out.add_column(stats, fun);
	if (fun=="range") {
		out.add_column(stats2, "max");
	}
	return(out);
}



SpatDataFrame SpatRaster::global_weighted_mean(SpatRaster &weights, std::string fun, bool narm, SpatOptions &opt) {

	SpatDataFrame out;

	std::vector<std::string> f {"sum", "mean"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("not a valid function");
		return(out);
	}

	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return(out);
	}

	if (weights.nlyr() != 1) {
		out.setError("The weights raster must have 1 layer");
		return(out);
	}
	if (!compare_geom(weights, false, false, opt.get_tolerance(), true)) {
		out.setError( msg.getError() );
		return(out);
	}

	std::vector<double> stats(nlyr());
	double stats2 = 0;
	std::vector<double> n(nlyr());
	std::vector<double> w(nlyr());
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!weights.readStart()) {
		out.setError(weights.getError());
		return(out);
	}

	BlockSize bs = getBlockSize(opt);
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v, wv;
		readValues(v, bs.row[i], bs.nrows[i], 0, ncol());
		weights.readValues(wv, bs.row[i], bs.nrows[i], 0, ncol());

		unsigned off = bs.nrows[i] * ncol() ;
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			double wsum = 0;
			unsigned offset = lyr * off;
			std::vector<double> vv(v.begin()+offset,  v.begin()+offset+off);
			for (size_t j=0; j<vv.size(); j++) {
				if (!std::isnan(vv[j]) && !std::isnan(wv[j])) {
					vv[j] *= wv[j];
					wsum += wv[j];
				} else {
					vv[j] = NAN;
				}
			}
			do_stat(vv, fun, narm, stats[lyr], stats2, n[lyr], i);
			w[lyr] += wsum;
		}
	}
	readStop();
	weights.readStop();

	if (fun=="mean") {
		for (size_t lyr=0; lyr<nlyr(); lyr++) {
			if (n[lyr] > 0 && w[lyr] != 0) {
				stats[lyr] /= w[lyr];
			} else {
				stats[lyr] = NAN;
			}
		}
		out.add_column(stats, "weighted_mean");
	} else {
		out.add_column(stats, "weighted_sum");
	}

	return(out);
}


SpatRaster SpatRaster::scale(std::vector<double> center, bool docenter, std::vector<double> scale, bool doscale, SpatOptions &opt) {
	SpatRaster out;
	SpatOptions opts(opt);
	SpatDataFrame df;
	if (docenter) {
		if (center.empty()) {
			df = global("mean", true, opts);
			center = df.getD(0);
		}
		if (doscale) {
			out = arith(center, "-", false, false, opts);
		} else {
			out = arith(center, "-", false, false, opt);
		}
	}
	if (doscale) {
		if (scale.empty()) {
			// divide by sd if centered, and the root mean square otherwise.
			// rms = sqrt(sum(x^2)/(n-1)); if centered rms == sd
			if (docenter) {
				df = out.global("rms", true, opts);
			} else {
				df = global("rms", true, opts);
			}
			scale = df.getD(0);
		}
		if (docenter) {
			out = out.arith(scale, "/", false, false, opt);
		} else {
			out = arith(scale, "/", false, false, opt);
		}
	}
	return out;
}


/*
bool can_use_replace(const std::vector<double> &from, const std::vector<double> &to) {
	// test if any "to" later occurs in "from"
	size_t n = from.size();
	for (size_t i = 0; i < (n-1); i++) {
		for (size_t j = (i+1); j < n; j++) {
			if (to[i] == from[j]) {
				return false;
			}
		}
	}
	return true;
}
*/

SpatRaster SpatRaster::replaceValues(std::vector<double> from, std::vector<double> to, long nl, bool setothers, double others, bool keepcats, SpatOptions &opt) {

	SpatRaster out;
	if (from.empty()) {
		out.setError("argument 'from' cannot be empty");
		return out;
	}
	if (to.empty()) {
		out.setError("argument 'to' cannot be empty");
		return out;
	}

	bool mout = false;
	bool min = false;
	if (nl > 1) {
		if (nlyr() > 1) {
			out.setError("cannot create layer-varying replacement with multi-layer input");
			return out;
		}
		mout = true;
	} else if (nl < -1) {
		nl = abs(nl);
		if (nlyr() != (size_t) nl) {
			out.setError("nlyr() does not match ncol(from)");
			return out;
		}
		min = true;
	}

	if (min) {
		out = geometry(1);
		if (keepcats) {
			out.source[0].hasCategories[0] = source[0].hasCategories[0];
			out.source[0].cats[0] = source[0].cats[0];
			out.source[0].hasColors = source[0].hasColors;
			out.source[0].cols = source[0].cols;
			
		}
	} else {
		if (nl == 0) {
			out = geometry(nlyr());
			out.source[0].hasCategories = hasCategories();
			out.source[0].cats = getCategories();
			out.source[0].hasColors = hasColors();
			out.source[0].cols = getColors();
			
		} else {
			out = geometry(nl);
			if (keepcats) {
				for (long i=0; i<nl; i++) {
					out.source[0].hasCategories[i] = source[0].hasCategories[0];
					out.source[0].cats[i] = source[0].cats[0];
					out.source[0].hasColors[i] = source[0].hasColors[0];
					out.source[0].cols[i] = source[0].cols[0];
				}
			}
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	if (mout) {
		size_t tosz = to.size() / nl;
		size_t nlyr = out.nlyr();
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			size_t vs = v.size();
			v.reserve(vs * nlyr);
			for (size_t lyr = 1; lyr < nlyr; lyr++) {
				v.insert(v.end(), v.begin(), v.begin()+vs);
			}
			std::vector<double> vv;
			if (setothers) {
				vv.resize(v.size(), others);
			} else {
				vv = v;
			}
			for (size_t lyr = 0; lyr < nlyr; lyr++) {
				std::vector<double> tolyr(to.begin()+lyr*tosz, to.begin()+(lyr+1)*tosz);
				recycle(tolyr, from);
				size_t offset = lyr*vs;
				for (size_t j=0; j< from.size(); j++) {
					if (std::isnan(from[j])) {
						for (size_t k=offset; k<(offset+vs); k++) {
							vv[k] = std::isnan(v[k]) ? tolyr[j] : v[k];
						}
					} else {
						for (size_t k=offset; k<(offset+vs); k++) {
							if (v[k] == from[j]) {
								vv[k] = tolyr[j];
							}
						}
					}
				}
			}
			if (!out.writeBlock(vv, i)) return out;
		}
	} else if (min) {
		size_t n = from.size()/nl;
		size_t nlr = nl;
		recycle(to, n);
		std::vector<std::vector<double>> fro(n);
		for (size_t i=0; i<n; i++) {
			fro[i].reserve(nlr);
			for (size_t j=0; j<nlr; j++) {
				fro[i].push_back(from[i*nlr+j]);
			}
		}

		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			size_t nc = v.size() / nlr;
			std::vector<double> vv(nc, others);
			for (size_t j=0; j<nc; j++) {
				for (size_t m=0; m<n; m++) {
					bool match = true;
					for (size_t k=0; k<nlr; k++) {
						if (std::isnan(fro[m][k])) {
							if (!std::isnan(v[nc*k+j])) {
								match = false;
								break;
							}
						} else if (v[nc*k+j] != fro[m][k]) {
							match = false;
							break;
						}
					}
					if (match) {
						vv[j] = to[m];
						break;
					}
				}
			}
			if (!out.writeBlock(vv, i)) return out;
		}
	} else {
		recycle(to, from);
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			std::vector<double> vv;
			if (setothers) {
				vv.resize(v.size(), others);
			} else {
				vv = v;
			}
			for (size_t j=0; j< from.size(); j++) {
				if (std::isnan(from[j])) {
					for (size_t k=0; k<v.size(); k++) {
						if (std::isnan(v[k])) {
							vv[k] = to[j];
						}
					}
				} else {
					for (size_t k=0; k<v.size(); k++) {
						if (v[k] == from[j]) {
							vv[k] = to[j];
						}
					}
				}
			}
			if (!out.writeBlock(vv, i)) return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}


void reclass_vector(std::vector<double> &v, std::vector<std::vector<double>> rcl, bool right_closed, bool left_right_closed, bool lowest, bool others, double othersValue) {


	size_t nc = rcl.size(); // should be 2 or 3

	double NAval = NAN;

	size_t n = v.size();
	unsigned nr = rcl[0].size();

	if (nc == 1) {
		std::vector<double> rc = rcl[0];
		std::sort(rc.begin(), rc.end());
		if (right_closed) {
			if (lowest)	{
				for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						v[i] = NAval;
					} else if ((v[i] < rc[0]) || (v[i] > rc[nr-1])) {
						v[i] = NAval;
					} else {
						for (size_t j=1; j<nr; j++) {
							if (v[i] <= rc[j]) {
								v[i] = j-1;
								break;
							}
						}
					}
				}
			} else { // !lowest
				for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						v[i] = NAval;
					} else if ((v[i] <= rc[0]) || (v[i] > rc[nr-1])) {
						v[i] = NAval;
					} else {
						for (size_t j=1; j<nr; j++) {
							if (v[i] <= rc[j]) {
								v[i] = j-1;
								break;
							}
						}
					}
				}
			}
		} else { // left_closed
			if (lowest)	{ // which means highest in this context
				for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						v[i] = NAval;
					} else if ((v[i] < rc[0]) || (v[i] > rc[nr-1])) {
						v[i] = NAval;
					} else if (v[i] == rc[nr-1]) {
						v[i] = nr-2; // safe because there must be at least 2 classes
					} else {
						for (size_t j=1; j<nr; j++) {
							if (v[i] < rc[j]) {
								v[i] = j-1;
								break;
							}
						}
					}
				}
			} else { // not highest
				for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						v[i] = NAval;
					} else if ((v[i] < rc[0]) || (v[i] >= rc[nr-1])) {
						v[i] = NAval;
					} else {
						for (size_t j=1; j<nr; j++) {
							if (v[i] < rc[j]) {
								v[i] = j-1;
								break;
							}
						}
					}
				}
			}
		}

	// "is - becomes"
	} else if (nc == 2) {

		bool hasNAN = false;
		double replaceNAN = NAval;
		for (size_t j=0; j<nr; j++) {
			if (std::isnan(rcl[0][j])) {
				hasNAN = true;
				replaceNAN = rcl[1][j];
			}
		}
		for (size_t i=0; i<n; i++) {

			if (std::isnan(v[i])) {
				if (hasNAN) {
					v[i] = replaceNAN;
				} else {
					v[i] = NAval;
				}
			} else {
				bool found = false;
				for (size_t j=0; j<nr; j++) {
					if (v[i] == rcl[0][j]) {
						v[i] = rcl[1][j];
						found = true;
						break;
					}
				}
				if ((!found) && others) {
					v[i] = othersValue;
				}
			}
		}

	// "from - to - becomes"
	} else {

		bool hasNAN = false;
		double replaceNAN = NAval;
		for (size_t j=0; j<nr; j++) {
			if (std::isnan(rcl[0][j]) || std::isnan(rcl[1][j])) {
				hasNAN = true;
				replaceNAN = rcl[2][j];
			}
		}


		if (left_right_closed) {   // interval closed at left and right

			for (size_t i=0; i<n; i++) {
				if (std::isnan(v[i])) {
					if (hasNAN) {
						v[i] = replaceNAN;
					} else {
						v[i] = NAval;
					}
				} else {
					bool found = false;
					for (size_t j=0; j<nr; j++) {
						if ((v[i] >= rcl[0][j]) && (v[i] <= rcl[1][j])) {
							v[i] = rcl[2][j];
							found = true;
							break;
						}
					}
					if ((!found) && others) {
						v[i] = othersValue;
					}
				}
			}
		} else if (right_closed) {
			if (lowest) {  // include lowest value (left) of interval

				double lowval = rcl[0][0];
				double lowres = rcl[2][0];
				for (size_t i=1; i<nr; i++) {
					if (rcl[0][i] < lowval) {
						lowval = rcl[0][i];
						lowres = rcl[2][i];
					}
				}

				for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						if (hasNAN) {
							v[i] = replaceNAN;
						} else {
							v[i] = NAval;
						}
					} else if (v[i] == lowval) {
						v[i] = lowres;
					} else {
						bool found = false;
						for (size_t j=0; j<nr; j++) {
							if ((v[i] > rcl[0][j]) && (v[i] <= rcl[1][j])) {
								v[i] = rcl[2][j];
								found = true;
								break;
							}
						}
						if ((!found) && others) {
							v[i] = othersValue;
						}
					}
				}

			} else { // !lowest
					for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						if (hasNAN) {
							v[i] = replaceNAN;
						} else {
							v[i] = NAval;
						}
					} else {
						bool found = false;
						for (size_t j=0; j<nr; j++) {
							if ((v[i] > rcl[0][j]) && (v[i] <= rcl[1][j])) {
								v[i] = rcl[2][j];
								found = true;
								break;
							}
						}
						if ((!found) && others) {
							v[i] = othersValue;
						}
					}
				}
			}

		} else { // left closed

			if (lowest) { // which here means highest because right=FALSE

				double lowval = rcl[1][0];
				double lowres = rcl[2][0];
				for (size_t i=0; i<nr; i++) {
					if (rcl[1][i] > lowval) {
						lowval = rcl[1][i];
						lowres = rcl[2][i];
					}
				}

				for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						if (hasNAN) {
							v[i] = replaceNAN;
						} else {
							v[i] = NAval;
						}
					} else if (v[i] == lowval) {
						v[i] = lowres;
					} else {
						bool found = false;
						for (size_t j=0; j<nr; j++) {
							if ((v[i] >= rcl[0][j]) && (v[i] < rcl[1][j])) {
								v[i] = rcl[2][j];
								found = true;
								break;
							}
						}
						if ((!found) && others) {
							v[i] = othersValue;
						}
					}
				}

			} else { //!dolowest

				for (size_t i=0; i<n; i++) {
					if (std::isnan(v[i])) {
						if (hasNAN) {
							v[i] = replaceNAN;
						} else {
							v[i] = NAval;
						}
					} else {
						bool found = false;
						for (size_t j=0; j<nr; j++) {
							if ((v[i] >= rcl[0][j]) && (v[i] < rcl[1][j])) {
								v[i] = rcl[2][j];
								found = true;
								break;
							}
						}
						if ((!found) && others) {
							v[i] = othersValue;
						}
					}
				}
			}
		}
	}
}



SpatRaster SpatRaster::reclassify(std::vector<std::vector<double>> rcl, unsigned openclosed, bool lowest, bool others, double othersValue, bool bylayer, bool brackets, bool keepcats, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (keepcats) {
		out.source[0].hasCategories = hasCategories();
		out.source[0].cats = getCategories();
	}
	size_t nc = rcl.size();
	size_t nr = rcl[0].size();
	size_t nl = nlyr();
	if (nl == 1) bylayer = false;
	size_t maxnc = 3 + nl * bylayer;
	size_t rcldim = nc;

	if (bylayer) {
		if (((nc != maxnc) && (nc != (maxnc-1))) || nr < 1) {
			out.setError("reclass matrix is not correct. Should be nlyr(x) plus 1 or 2");
			return out;
		}
		rcldim = nc - (nl-1);
	} else {
		if (nc < 1 || nc > 3 || nr < 1) {
			out.setError("matrix must have 1, 2 or 3 columns, and at least one row");
			return out;
		}
	}


	//bool left = openclosed == 0;
	bool right = openclosed != 0 ;
	bool leftright = openclosed == 2;


	if (nc == 1) {
		if (nr == 1) {
			int breaks = rcl[0][0];
			if (breaks < 2) {
				out.setError("cannot classify with a single number that is smaller than 2");
				return out;
			}
			std::vector<bool> hr = hasRange();
			bool hasR = true;
			for (size_t i=0; i<hr.size(); i++) {
				if (!hr[i]) hasR = false;
			}

			if (!hasR) {
				SpatOptions xopt(opt);
				setRange(xopt, false);
			}
			std::vector<double> mn = range_min();
			std::vector<double> mx = range_max();
			double mnv = vmin(mn, true);
			double mxv = vmax(mx, true);
			rcl[0] = seq_steps(mnv, mxv, breaks);
		}

		size_t rn = rcl[0].size();
		if ((rn > 1) && (rn < 256)) {
			std::vector<std::string> s;
			if (brackets) {
				std::string bleft = ((!right) || lowest) ? "[" : "(";
				std::string bright = right ? "]" : ")";
				s.push_back(bleft+ double_to_string(rcl[0][0]) + "–" + double_to_string(rcl[0][1]) + bright);
				bleft = right ? "(" : "[";
				for (size_t i=2; i<(rn-1); i++) {
					s.push_back(bleft + double_to_string(rcl[0][i-1]) + "–" + double_to_string(rcl[0][i]) + bright);
				}
				bright = (right || lowest) ? "]" : ")";
				s.push_back(bleft + double_to_string(rcl[0][rn-2]) + "–" + double_to_string(rcl[0][rn-1]) + bright);
			} else {
				for (size_t i=1; i<rn; i++) {
					s.push_back(double_to_string(rcl[0][i-1]) + " – " + double_to_string(rcl[0][i]));
				}
			}
			std::vector<long> u(s.size());
			std::iota(u.begin(), u.end(), 0);
			std::vector<std::string> nms = getNames();
			for (size_t i=0; i<out.nlyr(); i++) {
				out.setLabels(i, u, s, nms[i]);
			}
		}
		nr = rcl[0].size();
	}
	for (size_t i=0; i<nc; i++) {
		if (rcl[i].size() != nr) {
			out.setError("matrix is not rectangular");
			return out;
		}
	}
	if (rcldim == 3) {
		for (size_t i=0; i<nr; i++) {
			if (rcl[0][i] > rcl[1][i]) {
				out.setError("'from' larger than 'to': (" + std::to_string(rcl[0][i]) + " - " + std::to_string(rcl[1][i]) +")");
				return out;
			}
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	if (bylayer) {
		std::vector<std::vector<double>> lyrrcl(rcldim+1);
		for (size_t i=0; i<rcldim; i++) {
			lyrrcl[i] = rcl[i];
		}
		for (size_t i = 0; i < out.bs.n; i++) {
			unsigned off = bs.nrows[i] * ncol() ;
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (size_t lyr = 0; lyr < nl; lyr++) {
				unsigned offset = lyr * off;
				lyrrcl[rcldim] = rcl[rcldim+lyr];
				std::vector<double> vx(v.begin()+offset, v.begin()+offset+off);
				reclass_vector(vx, lyrrcl, right, leftright, lowest, others, othersValue);
				std::copy(vx.begin(), vx.end(), v.begin()+offset);
			}
			if (!out.writeBlock(v, i)) return out;
		}
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			reclass_vector(v, rcl, right, leftright, lowest, others, othersValue);
			if (!out.writeBlock(v, i)) return out;
		}
	}

	readStop();
	out.writeStop();
	return(out);

}


SpatRaster SpatRaster::reclassify(std::vector<double> rcl, unsigned nc, unsigned openclosed, bool lowest, bool others, double othersValue, bool bylayer, bool brackets, bool keepcats, SpatOptions &opt) {

	SpatRaster out;
	if ((rcl.size() % nc) != 0) {
		out.setError("incorrect length of reclassify matrix");
		return(out);
	}
	size_t maxnc = 3 + bylayer * (nlyr() - 1);
	unsigned nr = rcl.size() / nc;
	if (nc > maxnc) {
		out.setError("incorrect number of columns in reclassify matrix");
		return(out);
	}
	std::vector< std::vector<double>> rc(nc);

	for (size_t i=0; i<nc; i++) {
		rc[i] = std::vector<double>(rcl.begin()+(i*nr), rcl.begin()+(i+1)*nr);
	}

	out = reclassify(rc, openclosed, lowest, others, othersValue, bylayer, brackets, keepcats, opt);
	return out;
}



std::vector<std::vector<double>> clump_getRCL(std::vector<std::vector<size_t>> rcl, size_t n) {
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


void clump_replace(std::vector<double> &v, size_t n, std::vector<double>& d, size_t cstart, std::vector<std::vector<size_t>>& rcl, size_t &ncps) {

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
	} else if (d[nd-1] == (ncps-1)) {  // or just "if" ?
		ncps--;
	}
}




void broom_clumps(std::vector<double> &v, std::vector<double>& above, const size_t &dirs, size_t &ncps, const size_t &nr, const size_t &nc, std::vector<std::vector<size_t>> &rcl, bool is_global) {

	size_t nstart = ncps;

	bool d4 = dirs == 4;
	size_t stopnc = nc-1;
	std::vector<double> d;

	//first row, no row above it, use "above"
	//first cell
	//Rcpp::Rcout << "r  x i v[i] nc v[i] nc" << std::endl;
	if ( !std::isnan(v[0]) ) {
		//Rcout << 0 << " ff " << 0 << " " << v[0] << " " << ncps << " " ;
		if (d4) {
			if (std::isnan(above[0])) {
				v[0] = ncps; // new patch
				ncps++;
			} else {
				v[0] = above[0]; // same as above
			}
		} else if (is_global) { //d8 global
			d = {above[0], above[1], above[stopnc]} ;
			clump_replace(v, 0, d, nstart, rcl, ncps);
		} else { //d8
			d = {above[0], above[1]} ;
			clump_replace(v, 0, d, nstart, rcl, ncps);
		}
		//Rcout << v[0] << " " << ncps << std::endl;
	}
	// other cells
	for (size_t i=1; i<stopnc; i++) {

		if (!std::isnan(v[i])) {
			//Rcout << 0 << " fm " << i << " " << v[i] << " " << ncps << " " ;
			if (d4) {
				d = {above[i], v[i-1]} ;
			} else {
				d = {above[i], above[i-1], above[i+1], v[i-1]} ;
			}
			clump_replace(v, i, d, nstart, rcl, ncps);
			//Rcout << v[i] << " " << ncps << std::endl;
		}
	}
	// last cell
	size_t i = stopnc;
	if (!std::isnan(v[i])) {
		//Rcout << 0 << " fl " << i << " " << v[i] << " " << ncps << " " ;

		if (is_global) {
			if (d4) {
				d = {above[i], v[i-1], v[0]} ;
			} else {
				d = {above[i], above[i-1], v[i-1], v[0], above[0]} ;
			}
		} else {
			if (d4) {
				d = {above[i], v[i-1]} ;
			} else {
				d = {above[i], above[i-1], v[i-1]} ;
			}
		}
		clump_replace(v, i, d, nstart, rcl, ncps);
		//Rcout << v[i] << " " << ncps << std::endl;
	}

	////////
	//other rows
	for (size_t r=1; r<nr; r++) {
		size_t start = r*nc;
		size_t i=start;
		// first cell
		if (!std::isnan(v[i])) {
			//Rcout << r << " f " << i << " " << v[i] << " " << ncps << " " ;
			if (is_global) {
				if (d4) {
					if (std::isnan(v[i-nc])) {
						v[i] = ncps;
						ncps++;
					} else {
						v[i] = v[i-nc];
					}
				} else {
					d = {v[i-1], v[i-nc], v[i-nc+1]}; //above and left-above
					clump_replace(v, i, d, nstart, rcl, ncps);
				}
			} else {
				if (d4) {
					if (std::isnan(v[i-nc])) {
						v[i] = ncps;
						ncps++;
					} else {
						v[i] = v[i-nc];
					}
				} else {
					d = {v[i-nc], v[i-nc+1]}; //above and left-above
					clump_replace(v, i, d, nstart, rcl, ncps);
				}
			}
			//Rcout << v[i] << " " << ncps << std::endl;
		}

		size_t stop = start + stopnc;

		// other cells
		for (size_t i=(start+1); i<stop; i++) {
			if (!std::isnan(v[i])) {
				//Rcout << r << " m " << i << " " << v[i] << " " << ncps << " " ;
				if (d4) {
					d = {v[i-nc], v[i-1]} ;
				} else {
					d = {v[i-nc], v[i-nc-1], v[i-nc+1], v[i-1]} ;
				}
				clump_replace(v, i, d, nstart, rcl, ncps);
				//Rcout << v[i] << " " << ncps << std::endl;
			}
		}

		// last cell
		i = stop;
		if (!std::isnan(v[i])) {
			//Rcout << r << " l " << i << " " << v[i] << " " << ncps << " " ;
			if (is_global) {
				if (d4) {
					d = {v[i-nc], v[i-1], v[start]} ;
				} else {
					d = {v[i-nc], v[i-nc-1], v[i-1], v[start], v[start-nc]} ;
				}
			} else if (!std::isnan(v[i])) {
				if (d4) {
					d = {v[i-nc], v[i-1]} ;
				} else {
					d = {v[i-nc], v[i-nc-1], v[i-1]} ;
				}
			}
			clump_replace(v, i, d, nstart, rcl, ncps);
			//Rcout << v[i] << " " << ncps << std::endl;
		}
	}
	size_t off = (nr-1) * nc;
	above = std::vector<double>(v.begin()+off, v.end());
}



SpatRaster SpatRaster::clumps(int directions, bool zeroAsNA, SpatOptions &opt) {

	SpatRaster out = geometry(1);

	if (nlyr() > 1) {
		SpatOptions ops(opt);
		std::vector<std::string> nms = getNames();
		if (ops.names.size() == nms.size()) {
			nms = opt.names;
		}
		for (size_t i=0; i<nlyr(); i++) {
			std::vector<unsigned> lyr = {(unsigned)i};
			ops.names = {nms[i]};
			SpatRaster x = subset(lyr, ops);
			x = x.clumps(directions, zeroAsNA, ops);
			out.addSource(x, false, ops);
		}
		if (!opt.get_filename().empty()) {
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
	if (!filename.empty()) {
		bool overwrite = opt.get_overwrite();
		std::string errmsg;
		if (!can_write({filename}, filenames(), overwrite, errmsg)) {
			out.setError(errmsg + " (" + filename +")");
			return(out);
		}
	}
	if (opt.names.empty()) {
		opt.names = {"patches"};
	}

	opt.set_filenames({""});
 	if (!out.writeStart(opt, filenames())) { return out; }
	size_t nc = ncol();
	size_t ncps = 1;
	std::vector<double> above(nc, NAN);
	std::vector<std::vector<size_t>> rcl(2);

	bool is_global = is_global_lonlat();
	for (size_t i = 0; i < out.bs.n; i++) {
		readBlock(v, out.bs, i);
		if (zeroAsNA) {
			std::replace(v.begin(), v.end(), 0.0, (double)NAN);
		}
		broom_clumps(v, above, directions, ncps, out.bs.nrows[i], nc, rcl, is_global);
		if (!out.writeBlock(v, i)) return out;
		// perhaps here keep track of unique values, so that gaps can be removed
	}
	out.writeStop();
	readStop();

	opt.set_filenames({filename});
	if (!rcl[0].empty()) {
		std::vector<std::vector<double>> rc = clump_getRCL(rcl, ncps);
		out = out.reclassify(rc, 3, true, false, 0.0, false, false, false, opt);
	} else if (!filename.empty()) {
		out = out.writeRaster(opt);
	}
	return out;
}



bool SpatRaster::replaceCellValues(std::vector<double> &cells, std::vector<double> &v, bool bylyr, SpatOptions &opt) {

	size_t cs = cells.size();
	double nce = ncell() - 1;
	for (size_t i=0; i<cs; i++) {
	if ((cells[i] < 0) || (cells[i] > nce)) {
			setError("cell number(s) out of range");
			return false;
		}
	}
	size_t vs = v.size();
	size_t nl = nlyr();
	if (vs == 1) {
		bylyr = false;
		recycle(v, cs);
	} else if (bylyr) {
		if (vs == nl) {
			rep_each(v, cs);
		} else if (vs != (cs*nl)) {
			setError("length of cells and values do not match");
			return false;
		}
	} else if (cs != vs) {
		if ((vs / nl) == cs) {
			bylyr = true;
		} else {
			setError("lengths of cells and values do not match");
			return false;
		}
	}
	size_t nc = ncell();
	size_t ns = nsrc();

	if (!hasValues()) {
		*this = init({NAN}, opt);
	}

	for (size_t i=0; i<ns; i++) {
		if (!source[i].memory) {
			if (!canProcessInMemory(opt)) {
				try {
					readAll();
				} catch(...) {
					setError("cannot process this raster in memory");
					return false;
				}
			} else {
				readAll();
			}
			break;
		}
	}
	if (bylyr) {
		size_t addlyr = 0;
		for (size_t i=0; i<ns; i++) {
			size_t nl = source[i].nlyr;
			for (size_t j=0; j<nl; j++) {
				size_t off = nc * j;
				size_t koff = cs * (j+addlyr);
				for (size_t k=0; k<cs; k++) {
					if (!std::isnan(cells[k])) {
						source[i].values[off + cells[k]] = v[koff + k];
					}
				}
			}
			source[i].setRange();
			addlyr += nl;
		}
	} else {
		//double minv = vmin(v, true);
		//double maxv = vmax(v, true);
		for (size_t i=0; i<ns; i++) {
			size_t nl = source[i].nlyr;
			for (size_t j=0; j<nl; j++) {
				size_t off = nc * j;
				for (size_t k=0; k<cs; k++) {
					if (!std::isnan(cells[k])) {
						source[i].values[off + cells[k]] = v[k];
					}
				}
			}
			source[i].setRange();
		}
	}
	return true;
}


bool SpatRaster::replaceCellValuesLayer(std::vector<size_t> layers, std::vector<double> &cells, std::vector<double> &v, bool bylyr, SpatOptions &opt) {


	size_t cs = cells.size();
	double nce = ncell() - 1;
	for (size_t i=0; i<cs; i++) {
		if ((cells[i] < 0) || (cells[i] > nce)) {
			setError("cell number(s) out of range");
			return false;
		}
	}

	size_t nl = layers.size();

	size_t maxnl = nlyr()-1;
	for (size_t i=0; i<nl; i++) {
		if (layers[i] > maxnl) {
			setError("invalid layer number");
			return(false);
		}
	}

	size_t vs = v.size();
	if (vs == 1) {
		bylyr = false;
		recycle(v, cs);
	} else if (bylyr) {
		if (vs != (cs*nl)) {
			setError("length of cells and values do not match");
			return false;
		}
	} else if (cs != vs) {
		if ((vs / nl) == cs) {
			bylyr = true;
		} else {
			setError("lengths of cells and values do not match");
			return false;
		}
	}
	size_t nc = ncell();


	if (!hasValues()) {
		*this = init({NAN}, opt);
	}

	std::vector<size_t> srcs;
	srcs.reserve(nl);
	for (size_t i=0; i<nl; i++) {
	    std::vector<unsigned> sl = findLyr(layers[i]);
		size_t src = sl[0];
		size_t lyr = sl[1];

		srcs.push_back(src);

		if (!source[src].memory) {
			// if sources is a temp file update the file?
			// or create a tmp file?
			try {
				readAll();
			} catch(...) {
				setError("cannot process this raster in memory");
				return false;
			}
		}

		size_t off = nc * lyr;
		if (bylyr) {
			size_t koff = cs * i;
			for (size_t k=0; k<cs; k++) {
				if (!std::isnan(cells[k])) {
					source[src].values[off + cells[k]] = v[koff + k];
				}
			}
		} else {
			for (size_t k=0; k<cs; k++) {
				if (!std::isnan(cells[k])) {
					source[src].values[off + cells[k]] = v[k];
				}
			}
		}
	}
	std::sort(srcs.begin(), srcs.end());
	srcs.erase(std::unique(srcs.begin(), srcs.end()), srcs.end());
	for (size_t i=0; i<srcs.size(); i++) {
		source[i].setRange();
	}
	return true;
}



SpatRaster SpatRaster::rgb2hsx(std::string type, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) {
		out.setError("no cell values");
		return out;
	}
	if ((!rgb)  || (rgblyrs.size() < 3)) {
		out.setError("no RGB channels");
		return out;
	}

	bool hsv=false;
	bool hsi=false;
	//, hsl;
	std::vector<std::string> nms;
	if (type == "hsv") {
		nms = {"hue", "saturation", "value"};
		hsv = true;
	} else if (type == "hsi") {
		nms = {"hue", "saturation", "intensity"};
		hsi = true;
	} else if (type == "hsl") {
		nms = {"hue", "saturation", "lightness"};
		//hsl = true;
	} else {
		out.setError("unknown type. Should be one of 'hsv', 'hsi' or 'hsl'");
		return out;
	}

	out.setNames(nms);
	out.rgbtype = type;

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
 	if (!out.writeStart(opt, filenames())) { return out; }

	size_t nc=ncol();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		size_t n = out.bs.nrows[i] * nc;
		size_t n2 = n * 2;
		size_t iR = rgblyrs[0] * n;
		size_t iG = rgblyrs[1] * n;
		size_t iB = rgblyrs[2] * n;
		for (size_t j = 0; j < n; j++) {
			double R = v[j + iR] / 255.;
			double G = v[j + iG] / 255.;
			double B = v[j + iB] / 255.;
			double m = std::min(std::min(R, G), B);
			double M = std::max(std::max(R, G), B);
			double C = (M - m);

			if ((M == 0) || (C == 0)) {
				v[j] = 0; // H (hue)
				v[j+n] = 0;  // S (saturation)
				if (hsv) {
					v[j+n2] = M; // V
				} else if (hsi) {
					v[j+n2] = (R + G + B) / 3; // I
				} else {
					v[j+n2] = (M + m) / 2; // L
				}
			} else {
				// S
				if (hsv) {
					v[j+n] = C / M;
					v[j+n2] = M; // value
				} else if (hsi) {
					v[j+n2] = (R + G + B) / 3; // I
					v[j+n] = 1 - m / v[j+n2];
				} else {
					double L = (M + m) / 2;
					v[j+n] = C / (1 - std::fabs(2 * L - 1));
					v[j+n2] = L;
				}
				// H
				if (hsi) {
					double H = ((R-G)+(R-B))/2.0;
					H = H/sqrt((R-G)*(R-G) + (R-B)*(G-B));
					H = acos(H);
					if (B > G) {
						H = 2 * M_PI - H;
					}
					v[j] = H/(2 * M_PI);
					//v[j] = acos( sqrt((((R-G) + (R-B)) / 2) /  pow((R - G),2) + (R-B)*(G-B)) );
				} else {
					if (M == R) {
						v[j] = 60 * (G - B) / C;
					} else if (M == G) {
						v[j] = 60 * ((B - R) / C) + 120;
					} else {
						v[j] = 60 * ((R - G) / C) + 240;
					}
					v[j] = v[j] < 0 ? (v[j] + 360) / 360 : v[j] / 360;
				}
			}
		}
		if (!out.writeBlock(v, i)) return out;
	}
	out.writeStop();
	readStop();
	return out;
}



SpatRaster SpatRaster::hsx2rgb(SpatOptions &opt) {
	SpatRaster out = geometry();
	if (nlyr() != 3) {
		out.setError("x must have three layers");
		return out;
	}
	if (!hasValues()) {
		out.setError("no cell values");
		return out;
	}
	bool hsv=false;
	bool hsl=false;

	if (rgbtype == "hsv") {
		hsv = true;
	} else if (rgbtype == "hsl") {
		hsl = true;
	} else if (rgbtype != "hsi") {
		out.setError("input color scheme should be one of 'hsv', 'hsi' or 'hsl'");
		return out;
	}

	std::vector<std::string> nms={"red", "green", "blue"};
	out.setNames(nms);
	out.rgb = true;
	out.rgblyrs = {0,1,2};
	out.rgbtype = "rgb";

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
 	if (!out.writeStart(opt, filenames())) { return out; }
	size_t nc=ncol();

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		size_t n = out.bs.nrows[i] * nc;
		size_t n2 = n * 2;
		for (size_t j = 0; j < n; j++) {
			if (std::isnan(v[j])) continue;

			double H = v[j] * 360;
			double S = v[j+n];
			double X, C, m;
			if (hsv) {
				double V = v[j + n2];
				C = V * S;
				m = V - C;
				X = C * (1 - std::fabs(std::fmod((H / 60.), 2) - 1));
			} else if (hsl) {
				double L = v[j + n2];
				C = (1 - std::fabs(2*L-1)) * S;
				m = L - C/2;
				X = C * (1 - std::fabs(std::fmod((H / 60.), 2) - 1));
			} else { // hsi
				double I = v[j + n2];
				double Z = 1 - std::fabs((std::fmod(H/60., 2.)) -1);
				C = (3 * I * S) / (1 + Z);
				X = C * Z;
				m = I * (1-S);
			}
			if (H < 60) { v[j]=C; v[j+n]=X; v[j+n2]=0; }
			else if (H < 120) { v[j]=X; v[j+n]=C; v[j+n2]=0; }
			else if (H < 180) { v[j]=0; v[j+n]=C; v[j+n2]=X; }
			else if (H < 240) { v[j]=0; v[j+n]=X; v[j+n2]=C; }
			else if (H < 300) { v[j]=X; v[j+n]=0; v[j+n2]=C; }
			else { v[j]=C; v[j+n]=0; v[j+n2]=X; }

			v[j] = (v[j] + m) * 255;
			v[j+n] = (v[j+n] + m) * 255;
			v[j+n2] = (v[j+n2] + m) * 255;

		}
		if (!out.writeBlock(v, i)) return out;
	}
	out.writeStop();
	readStop();
	return out;
}


SpatRaster SpatRaster::sort(bool decreasing, bool order, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) {
		return out;
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	unsigned nl = out.nlyr();
	std::vector<double> v(nl);
	unsigned nc;
	if (order) {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> a;
			readBlock(a, out.bs, i);
			nc = out.bs.nrows[i] * out.ncol();
			std::vector<size_t> knc;
			knc.reserve(nl);
			for (size_t k=0; k<nl; k++) {
				knc.push_back(k*nc);
			}
			std::vector<size_t> ord;
			for (size_t j=0; j<nc; j++) {
				for (size_t k=0; k<nl; k++) {
					v[k] = a[j+knc[k]];
				}
				if (decreasing) {
					ord = sort_order_d(v);
				} else {
					ord = sort_order_a(v);
				}
				for (size_t k=0; k<v.size(); k++) {
					a[j+knc[k]] = ord[k];
				}
			}
			if (!out.writeBlock(a, i)) return out;
		}
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> a;
			readBlock(a, out.bs, i);
			nc = out.bs.nrows[i] * out.ncol();
			for (size_t j=0; j<nc; j++) {
				for (size_t k=0; k<nl; k++) {
					v[k] = a[j+k*nc];
				}
				if (decreasing) {
					std::sort(v.rbegin(), v.rend());
				} else {
					std::sort(v.begin(), v.end());
				}
				for (size_t k=0; k<v.size(); k++) {
					a[j+k*nc] = v[k];
				}
			}
			if (!out.writeBlock(a, i)) return out;

		}
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::combineCats(SpatRaster x, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	unsigned nl = std::max(nlyr(), x.nlyr());
	if (nl > 1) {
		out.setError("can only do this for a single layer SpatRasters");
	}

	if (!out.compare_geom(x, false, false, opt.get_tolerance(), true)) {
		out.setError("raster dimensions do not match");
		return(out);
	}
	if (!x.hasValues() || !hasValues()) {
		out.setError("both SpatRasters must have cell values");
	}
	std::vector<bool> cats = hasCategories();
	std::vector<bool> xcats = x.hasCategories();
	if ((cats[0]) && (xcats[0])) {
		SpatCategories sc = getLayerCategories(0);
		SpatCategories xsc = x.getLayerCategories(0);
		if (sc.concatenate(xsc)) {
			SpatOptions topt(opt);
			x.addSource(*this, false, topt);
			std::vector<double> from, to;
			to = sc.d.as_double(0);
			for (size_t i=0; i<to.size(); i++) {
				from.push_back(sc.d.iv[2][i]);
				from.push_back(sc.d.iv[1][i]);
			}
			opt.names = {sc.d.names[sc.index]};
			std::vector<unsigned> cr = {0,1};
			sc.d = sc.d.subset_cols(cr);
			x.source[0].cats[0] = sc;
			x.source[0].hasCategories[0] = true;

			x = x.replaceValues(from, to, -2, false, NAN, true, opt);
			return x;
		} else {
			out.setError("cannot concatenate categories");
			return out;
		}
	} else {
		out.setError("both SpatRasters must be categorical");
		return out;
	}
	//return(out);
}


SpatRaster SpatRaster::intersect(SpatRaster &x, SpatOptions &opt) {
	
	size_t nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);
	out.setValueType(3);
	
	if (!hasValues()) return out;
	if (!x.hasValues()) return out;

	if (!out.compare_geom(x, false, false, opt.get_tolerance(), true)) {
		if (!shared_basegeom(x, 0.1, true)) {
			out.setError("rasters are not aligned");
			return(out);
		} else {
			out.msg.has_error = false;
			out.msg.error = "";
			SpatExtent e = getExtent();
			e = e.intersect(x.getExtent());
			if (e.empty()) {
				out.setError("rasters do not intersect");
				return(out);
			}
			SpatOptions xopt(opt);
			x = x.crop(e, "near", false, xopt);
			SpatRaster y = crop(e, "near", false, xopt);
			return y.intersect(x, opt);
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		x.readStop();
		return out;
	}

	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> a, b;
		readValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		x.readValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(a, b);
		std::vector<double> d(a.size());
		for (size_t j=0; j<a.size(); j++) {
			if (std::isnan(a[j]) && std::isnan(b[j])) {
				d[j] = NAN;
			} else if (std::isnan(a[j]) || std::isnan(b[j])) {
				d[j] = 0;				
			} else {
				d[j] = 1;
			}
		}
		if (!out.writeBlock(d, i)) return out;
	}

	out.writeStop();
	readStop();
	x.readStop();
	return(out);

}


SpatRaster SpatRaster::fill_range(long limit, bool circular, SpatOptions &opt) {
	
	size_t nl = limit;
	SpatRaster out = geometry(nl, false, false, false);
	
	if (limit < 3) {
		out.setError("limit must be larger than 3");
		return out;
	}
	if (nlyr() != 2) {
		out.setError("the input raster must have two layers");
		return out;		
	}
	if (!hasValues()) {
		out.setError("the input raster must have values");
		return out;		
	}
		
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	
	for (size_t i=0; i<out.bs.n; i++) {
		size_t nc = out.bs.nrows[i] * ncol();
		std::vector<double> v;
		readValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		std::vector<double> d((v.size() / 2) * nl);
		if (circular) {
			for (size_t j=0; j<nc; j++) {
				size_t jnc = j+nc;
				size_t start = v[j]-1;
				size_t end = v[jnc];
				if (std::isnan(v[j]) || std::isnan(v[jnc])) {
					for (size_t k=0; k<nl; k++) {
						d[k*nc+j] = NAN;
					}									
				} else {
					bool circ = false;
					if (start > end) {
						std::swap(start, end);
						circ = true;
					}
					if (end > nl) {
						for (size_t k=0; k<nl; k++) {
							d[k*nc+j] = NAN;
						}
					} else {
						if (circ) {
							for (size_t k=start; k<nl; k++) {
								d[k*nc+j] = 1;
							}
							for (size_t k=0; k<end; k++) {
								d[k*nc+j] = 1;
							}
						} else {
							for (size_t k=start; k<end; k++) {
								d[k*nc+j] = 1;
							}
						}
					}
				}
			}	
		} else {
			for (size_t j=0; j<nc; j++) {
				size_t jnc = j+nc;
				if (std::isnan(v[j]) || std::isnan(v[jnc]) || (v[j] < 1) || (v[jnc] > nl) || (v[jnc] < v[j])) {
					for (size_t k=0; k<nl; k++) {
						d[k*nc+j] = NAN;
					}				
				} else {
					for (size_t k=(v[j]-1); k<v[jnc]; k++) {
						d[k*nc+j] = 1;
					}
				}
			}
		}
		if (!out.writeBlock(d, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
	
}
