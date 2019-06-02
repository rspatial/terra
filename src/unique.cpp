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


std::vector<std::vector<double>> SpatRaster::unique(bool bylayer) {
	readStart();
	BlockSize bs = getBlockSize(4);
	std::vector<double> v, x;
	unsigned nc = ncol();
	unsigned nl = nlyr();
	std::vector<std::vector<double>> out;
	
	if (bylayer) {
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<std::vector<double>> out(nl);
			v = readValues(bs.row[i], bs.nrows[i], 0, nc);
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*n;
				std::vector<double> m(v.begin()+off, v.begin()+off+n);
				std::sort(m.begin(), m.end());
				m.erase(std::unique(m.begin(), m.end()), m.end());
				out.insert(out.end(), m);
			}
		}
	} else {
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<std::vector<double>> m(n, std::vector<double>(nl));
			v = readValues(bs.row[i], bs.nrows[i], 0, nc);
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*n;
				for (size_t j=0; j<n; j++) {
					m[j][lyr] = v[off+j];
				}
			}
			std::sort(m.begin(), m.end());
			m.erase(std::unique(m.begin(), m.end()), m.end());
			for (size_t j=0; j<m.size(); j++) {
				out.insert(out.end(), m[j]);
			}
		}
	}
	readStop();
	std::sort(out.begin(), out.end());
	out.erase(std::unique(out.begin(), out.end()), out.end());	
	return(out);
}

