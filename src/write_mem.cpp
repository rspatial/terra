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
	if ((startcol==0) & (ncols==ncol())) {
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

	} else { // block writing
		size_t sz = source[0].values.size();
		if (sz == 0) {
			source[0].values.resize(ncol() * nrow() * nlyr(), NAN);
		}
		unsigned nc1 = nrows*ncols;
		unsigned nc2 = ncell();
		for (size_t i=0; i<nlyr(); i++) {
			unsigned off = i*nc1;
			for (size_t r=0; r<nrows; r++) {
				size_t start = r * ncols + off;
				std::vector<double> v(vals.begin()+start, vals.begin()+start+ncols);
				start = (startrow+r)*ncol() + i*nc2 + startcol;
				std::copy(v.begin(), v.end(), source[0].values.begin()+start);
			}
		}
	}
	return true;
}


