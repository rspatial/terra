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

#include "ogr_spatialref.h"
#include <vector>
#include "spatMessages.h"

SpatMessages transform_coordinates(std::vector<double> &x, std::vector<double> &y, std::string fromCRS, std::string toCRS) {
	
	SpatMessages m;
	
	OGRSpatialReference sourceCRS, targetCRS;
	OGRErr erro = sourceCRS.importFromProj4(&fromCRS[0]); 
	if (erro == 4) { 
		m.setError("crs is not valid");
		return(m);
	}
	erro = targetCRS.importFromProj4(&toCRS[0]); 
	if (erro == 4) { 
		m.setError("crs is not valid");
		return(m);
	}
	
	OGRCoordinateTransformation *poCT;
	poCT = OGRCreateCoordinateTransformation(&sourceCRS, &targetCRS );

	if( poCT == NULL )	{
		m.setError( "Transformation failed" );
		return (m);
	}
	
	unsigned failcount = 0;
	for (size_t i=0; i < x.size(); i++) {
		if( !poCT->Transform( 1, &x[i], &y[i] ) ) {
			x[i] = NAN;
			y[i] = NAN;
			failcount++;
		}
	}
	if (failcount > 0) {
		m.addWarning(std::to_string(failcount) + " failed transformations");
	}
	return m;
}

