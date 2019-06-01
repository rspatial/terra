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
#include "math_utils.h"


std::vector<std::vector<double>> SpatRaster::unique(unsigned decimals) {
	readStart();
	BlockSize bs = getBlockSize(4);
	std::vector<double> v, x;
	unsigned nc = ncol();
	unsigned nl = nlyr();
	std::vector<std::vector<double>> m(nl);
	std::vector<std::vector<double>> out(nl);
	for (size_t i = 0; i < bs.n; i++) {
		unsigned nr = bs.nrows[i];
		v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
		for (double &d : v) d = roundn(d, decimals);
		for (size_t j=0; j<nl; j++) {
			for (size_t r=0; r<nr; r++) {
				unsigned offset = nl*nc;
				m[j].insert(m[j].end(), v.begin()+offset, v.begin()+offset+nc);
			}
		}
		std::sort(m.begin(), m.end());
		m.erase(std::unique(m.begin(), m.end()), m.end());
		for (size_t j=0; j<nl; j++) {
			out[j].insert(out[j].end(), m[j].begin(), m[j].end());
		}
		m.resize(nl, std::vector<double>(0));
		std::sort(out.begin(), out.end());
		out.erase(std::unique(out.begin(), out.end()), out.end());	
	}
	readStop();
	return(out);
}

