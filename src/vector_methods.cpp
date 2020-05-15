// Copyright (c) 2018-2020  Robert J. Hijmans
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

#include "SpatVector.h"
#include "string_utils.h"


SpatVector SpatVector::aggregate(std::string field, bool dissolve) {

	SpatVector out;
	
	int i = where_in_vector(field, get_names());
	if (i < 0) {
		out.setError("cannot find field");
		return out;		
	}
	SpatDataFrame uv;
	std::vector<int> idx = lyr.df.getIndex(i, uv);

	for (size_t i=0; i<uv.nrow(); i++) {
		SpatGeom g;
		g.gtype = lyr.geoms[0].gtype;
		for (size_t j=0; j<idx.size(); j++) {
			if (i == (size_t)idx[j]) {
				g.unite( getGeom(j) );
			}	
		}
		out.addGeom(g);
	}
	
	out.lyr.srs = lyr.srs;
	out.lyr.df  = uv; 

	if (dissolve) {
		out.addWarning("cannot dissolve yet");
	}
	
	return out;
}

