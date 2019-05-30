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


#include <vector>
#include <limits>
#include <cmath>
#include <functional>
#include "spatRaster.h"
#include "vecmath.h"


template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

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
		message = "all values in argument 'fact' are 1, nothing to aggregate";
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
	fact[0] = std::max(unsigned(1), std::min(fact[0], nrow()));
	fact[1] = std::max(unsigned(1), std::min(fact[1], ncol()));
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

std::vector<std::vector<double> > SpatRaster::get_aggregates(std::vector<double> &in, size_t nr, std::vector<unsigned> dim) {

// adjust for chunk
	//dim[3] = std::ceil(double(nr) / dim[0]);

	size_t dy = dim[0], dx = dim[1], dz = dim[2];
	size_t bpC = dim[3];
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
    size_t ncells = ncell();
    size_t nl = nlyr();
    size_t nc = ncol();
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


SpatRaster SpatRaster::aggregate(std::vector<unsigned> fact, std::string fun, bool narm, SpatOptions &opt) {

	std::string message = "";
	bool success = get_aggregate_dims(fact, message);
	if (!success) {
		SpatRaster er = geometry();
		er.setError(message);
		return er;
	}

	double xmax = extent.xmin + fact[4] * fact[1] * xres();
	double ymin = extent.ymax - fact[3] * fact[0] * yres();
	SpatExtent e = SpatExtent(extent.xmin, xmax, ymin, extent.ymax);
	SpatRaster out = SpatRaster(fact[3], fact[4], fact[5], e, crs);

	if (!source[0].hasValues) { return out; }
	if (!out.writeStart(opt)) { return out ;}


	std::vector<std::string> f {"sum", "mean", "min", "max", "median", "modal"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown function argument");
		return out;
	}

	BlockSize bs = getBlockSize(4);
	if (bs.n > 1) {
		unsigned nr = floor(bs.nrows[0] / fact[1]) * fact[1];
		for (size_t i =0; i<bs.n; i++) {
			bs.row[i] = i * nr;
			bs.nrows[i] = nr;
		}
		while (true) {
			unsigned lastrow = bs.row[bs.n] + bs.nrows[bs.n] + 1;
			if (lastrow < nrow()) {
				bs.row.push_back(lastrow);
				bs.nrows.push_back(std::min(bs.nrows[bs.n], nrow()-lastrow));
				bs.n += 1;
			}
		}
	}

	std::function<double(std::vector<double>&, bool)> agFun;
	if (fun == "mean") {
		agFun = vmean<double>;
	} else if (fun == "sum") {
		agFun = vsum<double>;
	} else if (fun == "min") {
		agFun = vmin<double>;
	} else if (fun == "max") {
		agFun = vmax<double>;
	} else if (fun == "median") {
		agFun = vmedian<double>;
	} else if (fun == "modal") {
		agFun = vmodal<double>;
	} else {
		agFun = vmean<double>;
	}

	size_t row, col, cell, lyrcell, nr, nc, ncells;
	nr = nrow();
	nc = ncol();
	ncells = nc * nr;
	readStart();
	for (size_t b = 0; b < out.bs.n; b++) {
		std::vector<double> in = readBlock(bs, b);
		std::vector<double > v(fact[3] * fact[4] * fact[5]);
		std::vector<std::vector< double > > a = get_aggregates(in, bs.nrows[b], fact);
		size_t nblocks = a.size();
		for (size_t i = 0; i < nblocks; i++) {
			row = (i / nc) % nr;
			col = i % nc;
			cell = row * nc + col;
			lyrcell = std::floor(i / (ncells)) * ncells + cell;
			v[lyrcell] = agFun(a[i], narm);
		}

		if (!out.writeValues(v, bs.row[b])) return out;

	}
	out.writeStop();
	return(out);
}


