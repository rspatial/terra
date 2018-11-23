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
#include "spatRaster.h"


SpatRaster SpatRaster::disaggregate(std::vector<unsigned> fact, SpatOptions opt) {

    SpatRaster out = geometry();
	std::string message = "";
	bool success = get_aggregate_dims(fact, message);
	if (!success) {
		out.setError(message);
		return out;
	}
	fact.resize(3);
    fact[3] = 1; // at least for now
    out.source[0].nrow = out.source[0].nrow  * fact[0];
    out.source[0].ncol = out.source[0].ncol  * fact[1];
    out.source[0].nlyr = out.source[0].nlyr  * fact[2];

    if (!hasValues()) {
        return out;
    }

	BlockSize bs = getBlockSize(4*fact[0]*fact[1]*fact[2]);
	std::vector<double> v, vout;
	double nc=ncol();
	std::vector<double> newrow(nc*fact[1]);
  	readStart();
  	out.writeStart(opt);
	for (size_t i = 0; i < bs.n; i++) {
		v = readValues(bs.row[i], bs.nrows[i], 0, nc);
		for (size_t row=0; row<bs.nrows[i]; row++) {
            double off = row*nc;
            for (size_t j=0; j<nc; j++) {
                for (size_t k=0; k<fact[1]; k++) {
                    newrow[j*fact[1]+k] = v[off+j];
                }
            }
            for (size_t col=0; col<fact[0]; col++) {
                vout.insert(vout.end(), newrow.begin(), newrow.end());
            }
		}
		out.writeValues(vout, bs.row[i]);
	}
	out.writeStop();
	readStop();
	return(out);
}


