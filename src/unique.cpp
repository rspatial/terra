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
#include <limits>


std::vector<std::vector<double>> SpatRaster::unique(bool bylayer) {

	std::vector<std::vector<double>> out;
	if (!hasValues()) return out;

	constexpr double lowest_double = std::numeric_limits<double>::lowest();
	
	BlockSize bs = getBlockSize(4);
	unsigned nc = ncol();
	unsigned nl = nlyr();
	readStart();
	
	if (bylayer) {
		out.resize(nl);
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, nc);

			for (size_t j = 0; j < v.size(); j++) {
				if (std::isnan(v[j])) v[j] = lowest_double;
			}
			
			for (size_t lyr=0; lyr<nl; lyr++) {
				unsigned off = lyr*n;
				out[lyr] = std::vector<double>(v.begin()+off, v.begin()+off+n);
				std::sort(out[lyr].begin(), out[lyr].end());
				out[lyr].erase(std::unique(out[lyr].begin(), out[lyr].end()), out[lyr].end());
			}
		}
		for (size_t lyr=0; lyr<nl; lyr++) {
			std::sort(out[lyr].begin(), out[lyr].end());
			out[lyr].erase(std::unique(out[lyr].begin(), out[lyr].end()), out[lyr].end());
		}
	} else {
		for (size_t i = 0; i < bs.n; i++) {
			unsigned n = bs.nrows[i] * nc;
			std::vector<std::vector<double>> m(n, std::vector<double>(nl));
			std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, nc);
			for (size_t j = 0; j < v.size(); j++) {
				if (std::isnan(v[j])) v[j] = lowest_double;
			}

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
		std::sort(out.begin(), out.end());
		out.erase(std::unique(out.begin(), out.end()), out.end());
	}
	readStop();

	for (size_t i = 0; i < out.size(); i++) {
		for (size_t j = 0; j < out[0].size(); j++) {
			if (out[i][j] == lowest_double) {
				out[i][j] = NAN;
			}
		}
	}
	return(out);
}

