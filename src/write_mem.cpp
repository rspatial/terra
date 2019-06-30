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


#include "SpatRaster.h"


bool SpatRaster::writeValuesMem(std::vector<double> &vals, unsigned startrow, unsigned nrows, unsigned startcol, unsigned ncols) {
	size_t sz = source[0].values.size();
	size_t start = startrow * ncol() * nlyr();
	if (sz == 0) { // first or all
		source[0].values = vals;
	} else if (sz == start) { // in chunks
		source[0].values.insert(source[0].values.end(), vals.begin(), vals.end());
	} else { // async
		if (start+vals.size() > sz) {
			source[0].values.resize(start+vals.size());
		}
		for (size_t i=0; i<vals.size(); i++) {
			source[0].values[start+i] = vals[i];
		}
	}
	return true;
}


