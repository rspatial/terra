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
#include "read_rst.h"


std::vector<double> SpatRaster::readValuesBinary(unsigned src, unsigned row, unsigned nrows, unsigned col, unsigned ncols) {
	unsigned nl = source[src].nlyr;
	std::string dtype = source[src].datatype;
	std::string filename = source[src].filename;
	std::string bndorder = source[src].bandorder;
	std::vector<double> out;
	if (ncols == ncol()) {
		if (nrows == nrow()) {
			out = readBinAll(filename, dtype, nrow(), ncol(), nl, bndorder);
		} else {
			out = readBinRows(filename, dtype, row, nrows, nrow(), ncol(), nl, bndorder);
		}
	} else {
		out = readBinBlock(filename, dtype, row, nrows, col, ncols, nrow(), ncol(), nl, bndorder);
	}
	return out;
}


