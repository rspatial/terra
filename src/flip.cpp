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
#include "spatRaster.h"

SpatRaster SpatRaster::flip(bool vertical, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return out;
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	std::vector<double> b;
	unsigned nc = ncol();
	unsigned nl = nlyr();
	
	if (vertical) {
		for (size_t i=0; i < out.bs.n; i++) {
			size_t ii = out.bs.n - 1 - i;
			std::vector<double> a = readBlock(out.bs, ii);
			unsigned lyrrows = nl * out.bs.nrows[ii];
			for (size_t j=0; j < lyrrows; j++) {
				unsigned start = (lyrrows - 1 - j) * nc;
				unsigned end = start + nc;
				b.insert(b.end(), a.begin()+start, a.begin()+end);
			}
			if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
			b.resize(0);
		}		
	} else {	
		for (size_t i=0; i < out.bs.n; i++) {
			std::vector<double> a = readBlock(out.bs, i);
			unsigned lyrrows = nl * out.bs.nrows[i];
			for (size_t j=0; j < lyrrows; j++) {
				unsigned start = j * nc;
				unsigned end = start + nc;
				b.insert(b.end(), a.begin()+start, a.begin()+end);
			}
			if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
			b.resize(0);
		}
	}	
	out.writeStop();
	readStop();
	return(out);
}




