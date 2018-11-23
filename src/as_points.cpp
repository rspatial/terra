// Copyright (c) 2018  Robert J. Hijmans
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
#include <cmath>
#include "spatRaster.h"

SpatVector SpatRaster::as_points(bool values, bool narm) {

// for now assuming one layer

	BlockSize bs = getBlockSize(4);
	std::vector<double> v, vout;
	vout.reserve(v.size());
	SpatVector pv;
	if (!values && !narm) {

	}

	readStart();
	unsigned nc = ncol();
	unsigned lyr;
	for (size_t i = 0; i < bs.n; i++) {
		v = readValues(bs.row[i], bs.nrows[i], 0, nc);
		unsigned ncells = bs.nrows[i] * nc;
		std::vector<double> cells;
		cells.reserve(ncells);
		if (narm) {
			unsigned off1 = bs.row[i] + nc;
			std::vector<double> boolcells(ncells, false);
			for (size_t j=0; j<v.size(); j++) {
				lyr = j / ncells;
				unsigned off2 = off1 - lyr*nc;
				if (!std::isnan(v[j])) {
					boolcells[off2 + j] = true;
				}
			}
			for (size_t j=0; j<cells.size(); j++) {
				if (boolcells[j]) {
					cells.push_back(j+off1);
				}
			}
		} else {
	   		std::iota(cells.begin(), cells.end(), 0);
		}

		std::vector<std::vector<double>> xy = xyFromCell(cells);
		std::vector<double> vals = extractCell(cells);
	}
	readStop();
	return(pv);
}

