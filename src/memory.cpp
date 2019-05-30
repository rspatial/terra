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

#include "spatRaster.h"
#include "ram.h"


bool SpatRaster::canProcessInMemory(unsigned n) {
	double f = 0.5;
	return (n * size()) < (availableRAM() * f);
}

unsigned SpatRaster::chunkSize(unsigned n) {
	double f = 0.25;
	unsigned cells_in_row = n * ncol() * nlyr();
	unsigned rows = availableRAM() * f / cells_in_row;
	return rows == 0 ? 1 : std::min(rows, nrow());
}

BlockSize SpatRaster::getBlockSize(unsigned n) {
	BlockSize bs;

//	if (source[0].filename == "") {
	// in memory
//		bs.row = {0};
//		bs.nrows = {nrow};
//		bs.n = 1;

//	} else {

		unsigned cs = chunkSize(n);
		unsigned chunks = ceil(nrow() / double(cs));
		bs.n = chunks;
		bs.row = std::vector<unsigned>(chunks);
		bs.nrows = std::vector<unsigned>(chunks, cs);

		unsigned r = 0;
		for (size_t i =0; i<chunks; i++) {
			bs.row[i] = r;
			r += cs;
		}
		bs.nrows[chunks-1] = cs - ((chunks * cs) - nrow());

//	}
	return bs;
}
