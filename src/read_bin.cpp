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
	unsigned nl = source[src].nlyrfile;
	std::string dtype = source[src].datatype;
	std::string filename = source[src].filename;
	std::string bndorder = source[src].bandorder;
	std::vector<unsigned> lyrs = source[src].layers;
	std::vector<double> out;
	if (ncols == ncol()) {
		if (nrows == nrow()) {
			out = readBinAll(filename, dtype, lyrs, nrow(), ncol(), nl, bndorder);
		} else {
			out = readBinRows(filename, dtype, row, nrows, lyrs, nrow(), ncol(), nl, bndorder);
		}
	} else {
		out = readBinBlock(filename, dtype, row, nrows, col, ncols, lyrs, nrow(), ncol(), nl, bndorder);
	}
	return out;
}


std::vector<double> SpatRaster::readSampleBinary(unsigned src, unsigned srows, unsigned scols) {

	unsigned nl = source[src].nlyrfile;
	std::string dtype = source[src].datatype;
	std::string filename = source[src].filename;
	std::string bndorder = source[src].bandorder;
	std::vector<unsigned> lyrs = source[src].layers;

	double rowstep = nrow() / (srows+1);
	std::vector<unsigned> steprows, stepcols;
	for (size_t i = 0; i < srows; i++) {
		steprows[i] = round((i+0.5) * rowstep);
	}
	double colstep = ncol() / (scols+1);
	for (size_t i = 0; i < scols; i++) {
		stepcols[i] = round((i+0.5) * colstep);
	}

	std::vector<double> cells = cellFromRowCol(steprows, stepcols);
	std::vector<std::vector<double>> v = readBinCell(filename, dtype, cells, lyrs, nrow(), ncol(), nl, bndorder);
	std::vector<double> out;
	for (size_t i = 0; i<lyrs.size(); i++) {
		out.insert(out.end(), v[i].begin(), v[i].end());
	}
	return out;
}


std::vector<std::vector<double>> SpatRaster::readCellsBinary(unsigned src, std::vector<double> cells) {
	
	unsigned nl = source[src].nlyrfile;
	std::string dtype = source[src].datatype;
	std::string filename = source[src].filename;
	std::string bndorder = source[src].bandorder;
	std::vector<unsigned> lyrs = source[src].layers;

	std::vector<std::vector<double>> out = readBinCell(filename, dtype, cells, lyrs, nrow(), ncol(), nl, bndorder);

	return out;
	
}


