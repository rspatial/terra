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


#include "ogr_spatialref.h"
#include <vector>
#include <string>
//#include "spatMessages.h"
#include "spatVector.h"
#include "gdalhelp.h"

SpatMessages transform_coordinates(std::vector<double> &x, std::vector<double> &y, std::string fromCRS, std::string toCRS) {

	SpatMessages m;
	OGRSpatialReference source, target;
	const char *pszDefFrom = fromCRS.c_str();
	OGRErr erro = source.SetFromUserInput(pszDefFrom);
	if (erro != OGRERR_NONE) {
		m.setError("input crs is not valid");
		return m;
	}
	const char *pszDefTo = toCRS.c_str();
	erro = target.SetFromUserInput(pszDefTo);
	if (erro != OGRERR_NONE) {
		m.setError("output crs is not valid");
		return m;
	}

	OGRCoordinateTransformation *poCT;
	poCT = OGRCreateCoordinateTransformation(&source, &target);

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



SpatVector SpatVector::project(std::string crs) {

	SpatVector s;

    #ifndef useGDAL
		s.setError("GDAL is not available");
		return(s);
	#else
	SpatDataFrame d = getGeometryDF();

	std::vector<double> x = d.dv[0];
	std::vector<double> y = d.dv[1];

	s.msg = transform_coordinates(x, y, getCRS(), crs);

	if (!s.msg.has_error) {
		unsigned n = d.iv[0].size();
		std::vector<unsigned> a, b, c;
		for (size_t i=0; i<n; i++) {
			a.push_back(d.iv[0][i]);
			b.push_back(d.iv[1][i]);
			c.push_back(d.iv[2][i]);
		}
		s.setGeometry(type(), a, b, x, y, c);
		std::vector<std::string> refs = srefs_from_string(crs);
		s.setCRS(refs[0]);
		s.setPRJ(refs[1]);
		s.lyr.df = lyr.df;
	}
	#endif
	return s;
}
