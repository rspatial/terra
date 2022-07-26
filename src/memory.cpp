// Copyright (c) 2018-2022  Robert J. Hijmans
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
#include "ram.h"



bool SpatRaster::canProcessInMemory(SpatOptions &opt) {
	if (opt.get_todisk()) return false;
	double demand = size() * opt.ncopies;
	if (demand < opt.get_memmin()) {
		return true;
	}
	double supply;
	if (opt.get_memmax() > 0) {
		supply = opt.get_memmax() * opt.get_memfrac();
		//supply = std::min(supply, availableRAM());
	} else {
		supply = availableRAM() * opt.get_memfrac();
	}
	std::vector<double> v;
	double maxsup = v.max_size(); //for 32 bit systems
	supply = std::min(supply, maxsup);
	return (demand < supply);
}


size_t SpatRaster::chunkSize(SpatOptions &opt) {
	double n = opt.ncopies;
	double frac = opt.get_memfrac();

	double demand = size() * n;
	if (demand < opt.get_memmin()) {
		return nrow();
	}

	double cells_in_row = ncol() * nlyr() * n;
	double supply;

	if (opt.get_memmax() > 0) {
		supply = opt.get_memmax() * opt.get_memfrac();
		//supply = std::min(supply, availableRAM());
	} else {
		supply = availableRAM() * opt.get_memfrac();
	}
	double rows = supply * frac / cells_in_row;
	//double maxrows = 10000;
	//rows = std::min(rows, maxrows);
	size_t urows = floor(rows);
	urows = std::max(urows, (size_t)opt.minrows);
	if (urows < 1) return (1);
	if (urows < nrow()){
		return(urows);
	} else {
		return (nrow());
	}
}


std::vector<double> SpatRaster::mem_needs(SpatOptions &opt) {
	//returning bytes
	unsigned n = opt.ncopies;
	double memneed  = ncell() * (nlyr() * n);
	double memavail;
	if (opt.get_memmax() > 0) {
		memavail = opt.get_memmax();
		//memavail = std::min(memavail, availableRAM());
	} else {
		memavail = availableRAM();
	}
	double frac = opt.get_memfrac();
	double csize = chunkSize(opt);
	double inmem = canProcessInMemory(opt);
	std::vector<double> out = {memneed, memavail, frac, csize, inmem} ;
	return out;
}

//BlockSize SpatRaster::getBlockSize(unsigned n, double frac, unsigned steps) {
BlockSize SpatRaster::getBlockSize( SpatOptions &opt) {

	unsigned steps = opt.get_steps();

	BlockSize bs;
	size_t cs;

	if (steps > 0) {
		if (steps > nrow()) {
			steps = nrow();
		}
		bs.n = steps;
		cs = nrow() / steps;
	} else {
		cs = chunkSize(opt);
		bs.n = std::ceil(nrow() / double(cs));
	}
	bs.row = std::vector<size_t>(bs.n);
	bs.nrows = std::vector<size_t>(bs.n, cs);
	size_t r = 0;
	for (size_t i =0; i<bs.n; i++) {
		bs.row[i] = r;
		r += cs;
	}
	bs.nrows[bs.n-1] = cs - ((bs.n * cs) - nrow());

	return bs;
}

