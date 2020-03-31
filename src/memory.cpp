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


bool SpatRaster::canProcessInMemory(unsigned n) {
	double f = 0.4;
	return (n * size()) < (availableRAM() * f);
}

unsigned SpatRaster::chunkSize(unsigned n) {
	double f = 0.2;
	unsigned cells_in_row = n * ncol() * nlyr();
	unsigned rows = availableRAM() * f / cells_in_row;
	unsigned maxrows = 1000;
	rows = std::min(rows, maxrows);
	return rows == 0 ? 1 : std::min(rows, nrow());
	
}

BlockSize SpatRaster::getBlockSize(unsigned n, unsigned steps) {
	BlockSize bs;
	unsigned cs;
	
	if (steps > 0) {
		steps = std::min(steps, nrow());
		bs.n = steps;
		cs = nrow() / steps;
	} else {
		cs = chunkSize(n);
		bs.n = std::ceil(nrow() / double(cs));
	}

	bs.row = std::vector<unsigned>(bs.n);
	bs.nrows = std::vector<unsigned>(bs.n, cs);
	unsigned r = 0;
	for (size_t i =0; i<bs.n; i++) {
		bs.row[i] = r;
		r += cs;
	}
	bs.nrows[bs.n-1] = cs - ((bs.n * cs) - nrow());
	return bs;
}
