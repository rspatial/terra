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
#include <vector>
#include <string>
//#include "spatMessages.h"
#include "spatRaster.h"
#include "string_utils.h"


#ifndef useGDAL

bool SpatSRS::set(std::string txt, std::string &msg) {
	proj4 = txt;
	wkt = "";
	return true;
}
#else


#include "ogr_spatialref.h"
#include <gdal_priv.h> // GDALDriver


bool is_ogr_error(OGRErr err, std::string &msg) {
	if (err != OGRERR_NONE) {
		switch (err) {
			case OGRERR_NOT_ENOUGH_DATA:
				msg = "OGR: Not enough data";
			case OGRERR_UNSUPPORTED_GEOMETRY_TYPE:
				msg = "OGR: Unsupported geometry type";
			case OGRERR_CORRUPT_DATA:
				msg = "OGR: Corrupt data";
			case OGRERR_FAILURE:
				msg = "OGR: Invalid index";
			default:
				msg = "OGR: Error";
		}
		return true;
	}
	return false;
}



bool wkt_from_spatial_reference(const OGRSpatialReference *srs, std::string &wkt, std::string &msg) {
	char *cp;
#if GDAL_VERSION_MAJOR >= 3
	const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
	OGRErr err = srs->exportToWkt(&cp, options);
#else
	OGRErr err = srs->exportToWkt(&cp);
#endif
	if (is_ogr_error(err, msg)) {
		CPLFree(cp);
		return false;
	}
	wkt = std::string(cp);
	CPLFree(cp);
	return true;
}

bool prj_from_spatial_reference(const OGRSpatialReference *srs, std::string &prj, std::string &msg) {
	char *cp;
	OGRErr err = srs->exportToProj4(&cp);
	if (is_ogr_error(err, msg)) {
		CPLFree(cp);
		return false;
	}
	prj = std::string(cp);
	CPLFree(cp);
	return true;
}


bool string_from_spatial_reference(const OGRSpatialReference *srs, std::vector<std::string> &out, std::string &msg) {
	out = std::vector<std::string>(2, "");
	char *cp;
#if GDAL_VERSION_MAJOR >= 3
	const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
	OGRErr err = srs->exportToWkt(&cp, options);
#else
	OGRErr err = srs->exportToPrettyWkt(&cp);
#endif
	if (is_ogr_error(err, msg)) {
		CPLFree(cp);
		return false;
	}
	out[0] = std::string(cp);

	err = srs->exportToProj4(&cp);
	if (is_ogr_error(err, msg)) {
		CPLFree(cp);
		return false;
	}
	out[1] = std::string(cp);

	CPLFree(cp);
	return true;
}

/*
bool SpatSRS::set(OGRSpatialReference *poSRS, std::string &msg) {
	wkt="";
	proj4="";
	if (poSRS) {
		if (! wkt_from_spatial_reference(poSRS, wkt, msg)) {
			msg = "can't get wkt from srs";
			return false;
		};
		if (! prj_from_spatial_reference(poSRS, proj4, msg)) {
			msg = "can't get proj4 from srs";
			return false;
		};
	}
	return true;
}
*/

bool SpatSRS::set(std::string txt, std::string &msg) {
	wkt="";
	proj4="";
	lrtrim(txt);
	if (txt == "") {
		return true;
	} else {
		OGRSpatialReference *srs = new OGRSpatialReference;
		const char* s = txt.c_str();
		if (is_ogr_error(srs->SetFromUserInput(s), msg)) {
			delete srs;
			msg = "empty srs";
			return false;
		}
		if (! wkt_from_spatial_reference(srs, wkt, msg)) {
			delete srs;
			msg = "can't  get wkt from srs";
			return false;
		};
		if (! prj_from_spatial_reference(srs, proj4, msg)) {
			delete srs;
			msg = "can't  get proj4 from srs";
			return false;
		};
		delete srs;
		return true;
	}
	return false;
}

/*
bool SpatSRS::set(std::vector<std::string> txt, std::string &msg) {
	wkt="";
	proj4="";
	input="";
	if (txt.size() == 3) {
		proj4 == txt[0];
		wkt = txt[1];
		input = txt[2];
		return true;
	} else if (txt.size() == 2) {
		proj4 == txt[0];
		wkt = txt[1];
		return true;
	} else if (txt.size() == 1) { 
		input=txt[0];
		if (input != "") {
			OGRSpatialReference *srs = new OGRSpatialReference;
			const char* s = input.c_str();
			if (is_ogr_error(srs->SetFromUserInput(s), msg)) {
				delete srs;
				return false;
			}
			if (! wkt_from_spatial_reference(srs, wkt, msg)) {
				delete srs;
				return false;
			};
			if (! prj_from_spatial_reference(srs, proj4, msg)) {
				delete srs;
				return false;
			};
			delete srs;
			return true;
		}
	}
	return false;
}
*/


bool wkt_from_string(std::string input, std::string& wkt, std::string& msg) {
	lrtrim(input);
	wkt="";
	bool success = false;
	if (input != "") {
		OGRSpatialReference *srs = new OGRSpatialReference;
		const char* s = input.c_str();
		if (is_ogr_error(srs->SetFromUserInput(s), msg)) {
			delete srs;
			return false;
		}
		success = wkt_from_spatial_reference(srs, wkt, msg);
		delete srs;
	}
	return success;
}


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

	std::string srs = getSRS("wkt");
	std::string outwkt, msg;
	if (!wkt_from_string(crs, outwkt, msg)) {
		s.setError(msg);
		return s;
	}
	
	s.msg = transform_coordinates(x, y, srs, outwkt);

	if (!s.msg.has_error) {
		unsigned n = d.iv[0].size();
		std::vector<unsigned> a, b, c;
		for (size_t i=0; i<n; i++) {
			a.push_back(d.iv[0][i]);
			b.push_back(d.iv[1][i]);
			c.push_back(d.iv[2][i]);
		}
		s.setGeometry(type(), a, b, x, y, c);
		//std::vector<std::string> refs = srefs_from_string(crs);
		std::string msg;
		s.setSRS(crs);
		//s.setPRJ(refs[1]);
		s.df = df;
	}
	#endif
	return s;
}

#endif

