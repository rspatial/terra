// Copyright (c) 2018-2020  Robert J. Hijmans
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
	double supply = (availableRAM()) * opt.get_memfrac();
	std::vector<double> v;
	double maxsup = v.max_size(); //for 32 bit systems
	supply = std::min(supply, maxsup);
	return (demand < supply);
}


uint_64 SpatRaster::chunkSize(unsigned n, double frac) {
	double cells_in_row = n * ncol() * nlyr();
	double rows = (availableRAM()) * frac / cells_in_row;
	//double maxrows = 10000;
	//rows = std::min(rows, maxrows);
	uint_64 urows = floor(rows);
	if (rows < 1) return (1);
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
	double memavail = availableRAM(); 
	double frac = opt.get_memfrac();
	double csize = chunkSize(n, frac);
	double inmem = canProcessInMemory(opt); 
	std::vector<double> out = {memneed, memavail, frac, csize, inmem} ;
	return out;
}

//BlockSize SpatRaster::getBlockSize(unsigned n, double frac, unsigned steps) {
BlockSize SpatRaster::getBlockSize( SpatOptions &opt) {

	unsigned n = opt.get_blocksizemp();
	double frac = opt.get_memfrac();
	unsigned steps = opt.get_steps();
	
	BlockSize bs;
	uint_64 cs;
	
	if (steps > 0) {
		if (steps > nrow()) {
			steps = nrow();
		}
		bs.n = steps;
		cs = nrow() / steps;
	} else {
		cs = chunkSize(n, frac);
		bs.n = std::ceil(nrow() / double(cs));
	}

	bs.row = std::vector<uint_64>(bs.n);
	bs.nrows = std::vector<uint_64>(bs.n, cs);
	uint_64 r = 0;
	for (size_t i =0; i<bs.n; i++) {
		bs.row[i] = r;
		r += cs;
	}
	bs.nrows[bs.n-1] = cs - ((bs.n * cs) - nrow());
	
	return bs;
}

