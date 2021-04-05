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

//#include <vector>
#include <limits>
#include <cmath>
#include "spatRaster.h"
#include "vecmathfun.h"

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


std::vector<double> compute_aggregates(std::vector<double> &in, size_t nr, size_t nc, size_t nl, std::vector<unsigned> dim, std::function<double(std::vector<double>&, bool)> fun, bool narm) {

// dim 0, 1, 2, are the aggregations factors dy, dx, dz
// and 3, 4, 5 are the new nrow, ncol, nlyr

	size_t dy = dim[0], dx = dim[1], dz = dim[2];
//	size_t bpC = dim[3];
// adjust for chunk
	size_t bpC = std::ceil(double(nr) / dim[0]);
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
	std::vector<double> out(nblocks, NAN);

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
	return(out);
}



SpatRaster SpatRaster::aggregate(std::vector<unsigned> fact, std::string fun, bool narm, SpatOptions &opt) {

	SpatRaster out;
	std::string message = "";
	bool success = get_aggregate_dims(fact, message);

// fact 0, 1, 2, are the aggregation factors dy, dx, dz
// and  3, 4, 5 are the new nrow, ncol, nlyr
	if (!success) {
		out.setError(message);
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

	unsigned outnc = out.ncol();

	//BlockSize bs = getBlockSize(4, opt.get_memfrac());
	BlockSize bs = getBlockSize(opt);
	//bs.n = floor(nrow() / fact[0]); # ambiguous on solaris
	bs.n = std::floor(static_cast <double> (nrow() / fact[0]));

	bs.nrows = std::vector<size_t>(bs.n, fact[0]);
	bs.row.resize(bs.n);
	for (size_t i =0; i<bs.n; i++) {
		bs.row[i] = i * fact[0];
	}
	size_t lastrow = bs.row[bs.n - 1] + bs.nrows[bs.n - 1] + 1;
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
	
	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}

	size_t nc = ncol();
	for (size_t i = 0; i < bs.n; i++) {
        std::vector<double> vin = readValues(bs.row[i], bs.nrows[i], 0, nc);
		std::vector<double> v  = compute_aggregates(vin, bs.nrows[i], nc, nlyr(), fact, agFun, narm);
		if (!out.writeValues(v, i, 1, 0, outnc)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}


