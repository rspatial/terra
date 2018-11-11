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
#include "spatraster.h"


SpatRaster SpatRaster::mask(SpatRaster x, std::string filename, std::string format, std::string datatype, bool overwrite) {

// check for size; need for recycling
	//SpatRaster out = *this;
	SpatRaster out = geometry();

//	out.source.resize(1);
 //   our.source[0].nlyr = 1;
//	out.values.resize(0);
  	out.writeStart(filename, format, datatype, overwrite);
	readStart();
	x.readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol, 0, nlyr());
		m = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol, 0, nlyr());
		for (size_t i=0; i < v.size(); i++) {
			if (std::isnan(m[i])) {
				v[i] = NAN;
			}
		}
		out.writeValues(v, out.bs.row[i]);
	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}


/*

std::vector<std::vector<double> > matrix(int nrow, int ncol) {
	std::vector<std::vector<double> > m (nrow, std::vector<double>(ncol));
	return(m);
}


int main() {
	std::vector<vector<double> > d = matrix(10, 2);
	std::vector<double> m (10);
	m[1] = 1;
	m[5] = 1;
	d = mask(d, m, 1, 9, false);
	for (int i=0; i < d.size(); i++) {
		for (int j=0; j < d[0].size(); j++) {
			std::cout << ' ' << d[i][j];
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}
*/
