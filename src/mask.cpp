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
#include "recycle.h"

SpatRaster SpatRaster::mask(SpatRaster x, bool inverse, double maskvalue, double updatevalue, SpatOptions &opt) {

	unsigned nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);

	if (!compare_geom(x, false, true, false, true, true, false)) {
		out.setError("dimensions and/or extent do not match");
		return(out);
	}

	readStart();
	x.readStart();
  	if (!out.writeStart(opt)) { return out; }
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		m = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(v, m);
		if (inverse) {
			if (std::isnan(maskvalue)) {
				for (size_t i=0; i < v.size(); i++) {
					if (!std::isnan(m[i])) {
						v[i] = updatevalue;
					}
				}
			} else {
				for (size_t i=0; i < v.size(); i++) {
					if (m[i] != maskvalue) {
						v[i] = updatevalue;
					}
				}
			}		
		} else {
			if (std::isnan(maskvalue)) {
				for (size_t i=0; i < v.size(); i++) {
					if (std::isnan(m[i])) {
						v[i] = updatevalue;
					}
				}
			} else {
				for (size_t i=0; i < v.size(); i++) {
					if (m[i] == maskvalue) {
						v[i] = updatevalue;
					}
				}
			}
		}
		if (!out.writeValues(v, out.bs.row[i])) return out;
		
	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}


SpatRaster SpatRaster::mask(SpatVector x, bool inverse, double maskvalue, double updatevalue, SpatOptions &opt) {
	std::string filename = opt.get_filename();
	opt.set_filename("");	
	SpatRaster m = rasterize(x, 0, opt);
	opt.set_filename(filename);
	SpatRaster out = mask(m, inverse, 0, updatevalue, opt);
	return(out);
}

	
