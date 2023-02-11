// Copyright (c) 2018-2023  Robert J. Hijmans
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


bool wkt_from_spatial_reference(const OGRSpatialReference srs, std::string &wkt, std::string &msg) {
	char *cp;
#if GDAL_VERSION_MAJOR >= 3
	const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
	OGRErr err = srs.exportToWkt(&cp, options);
#else
	OGRErr err = srs.exportToWkt(&cp);
#endif
	if (is_ogr_error(err, msg)) {
		CPLFree(cp);
		return false;
	}
	wkt = std::string(cp);
	CPLFree(cp);
	return true;
}



bool prj_from_spatial_reference(const OGRSpatialReference srs, std::string &prj, std::string &msg) {
	char *cp;
	OGRErr err = srs.exportToProj4(&cp);
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
	const char *options[3] = { "MULTILINE=NO", "FORMAT=WKT2", NULL };
	OGRErr err = srs->exportToWkt(&cp, options);
#else
	OGRErr err = srs->exportToWkt(&cp);
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


double SpatSRS::to_meter() {
	double out;
	OGRSpatialReference x;
	if (wkt.size() < 2) {
		return NAN;
	}
	OGRErr erro = x.SetFromUserInput(wkt.c_str());
	if (erro != OGRERR_NONE) {
		return NAN;
	}
	if (x.IsGeographic()) {
		return 0;
	}
	out = x.GetLinearUnits();
	return out;
}

bool SpatSRS::is_same(SpatSRS other, bool ignoreempty) {
	if (ignoreempty) {
		if (is_empty() || other.is_empty()) {
			return true;
		}
	}
	OGRSpatialReference x, y;
	OGRErr erro = x.SetFromUserInput(wkt.c_str());
	if (erro != OGRERR_NONE) {
		return false;
	}
	erro = y.SetFromUserInput(other.wkt.c_str());
	if (erro != OGRERR_NONE) {
		return false;
	}
	return x.IsSame(&y);
}


bool SpatSRS::is_same(std::string other, bool ignoreempty) {

	if (wkt.empty() && other.empty()) {
		return true;
	} else if (wkt.empty() || other.empty()) {
		return ignoreempty ? true : false;
	}

	OGRSpatialReference x, y;
	OGRErr erro = x.SetFromUserInput(wkt.c_str());
	if (erro != OGRERR_NONE) {
		return false;
	}
	erro = y.SetFromUserInput(other.c_str());
	if (erro != OGRERR_NONE) {
		return false;
	}
	return x.IsSame(&y);
}


bool SpatSRS::is_lonlat() {
	OGRSpatialReference x;
	if (wkt.size() < 2) {
		return false;
	}
	OGRErr erro = x.SetFromUserInput(wkt.c_str());
	if (erro != OGRERR_NONE) {
		return false;
	}
	return x.IsGeographic();
}


bool SpatSRS::set(std::string txt, std::string &msg) {
	wkt="";
	proj4="";
	lrtrim(txt);

	if (txt.empty()) {
		return true;
	} else {
		OGRSpatialReference srs;
		OGRErr e = srs.SetFromUserInput(txt.c_str());
		if (is_ogr_error(e, msg)) {
			msg = "empty srs";
			return false;
		}
		if (! wkt_from_spatial_reference(srs, wkt, msg)) {
			msg = "can't get wkt from srs";
			return false;
		};
		if (! prj_from_spatial_reference(srs, proj4, msg)) {
			msg = "";
			//msg = "can't get proj4 from srs";
			//return false;
		};
		return true;
	}
	return false;
}


bool wkt_from_string(std::string input, std::string& wkt, std::string& msg) {
	lrtrim(input);
	wkt="";
	bool success = false;
	if (!input.empty()) {
		OGRSpatialReference srs;
		OGRErr e = srs.SetFromUserInput(input.c_str());
		if (is_ogr_error(e, msg)) {
			return false;
		}
		success = wkt_from_spatial_reference(srs, wkt, msg);
	}
	return success;
}




bool can_transform(std::string fromCRS, std::string toCRS) {

	OGRSpatialReference source, target;
	const char *pszDefFrom = fromCRS.c_str();
	OGRErr erro = source.SetFromUserInput(pszDefFrom);
	if (erro != OGRERR_NONE) {
		return false;
	}
	const char *pszDefTo = toCRS.c_str();
	erro = target.SetFromUserInput(pszDefTo);
	if (erro != OGRERR_NONE) {
		return false;
	}

	OGRCoordinateTransformation *poCT;
	poCT = OGRCreateCoordinateTransformation(&source, &target);
	if( poCT == NULL )	{
		OCTDestroyCoordinateTransformation(poCT);
		return false;
	}
	OCTDestroyCoordinateTransformation(poCT);
	return true;
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
		m.setError( "Cannot do this coordinate transformation" );
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

	OCTDestroyCoordinateTransformation(poCT);
	if (failcount > 0) {
		m.addWarning(std::to_string(failcount) + " failed transformations");
	}
	return m;
}


std::vector<double> SpatVector::project_xy(std::vector<double> x, std::vector<double> y, std::string fromCRS, std::string toCRS) {

	msg = transform_coordinates(x, y, fromCRS, toCRS);
	x.insert(x.end(), y.begin(), y.end());
	return x;

}


SpatVector SpatVector::project(std::string crs) {

	SpatVector s;
	s.reserve(size());

    #ifndef useGDAL
		s.setError("GDAL is not available");
		return(s);
	#else

	OGRSpatialReference source, target;
	std::string vsrs = getSRS("wkt");
	const char *pszDefFrom = vsrs.c_str();
	OGRErr erro = source.SetFromUserInput(pszDefFrom);
	if (erro != OGRERR_NONE) {
		s.setError("input crs is not valid");
		return s;
	}
	const char *pszDefTo = crs.c_str();
	erro = target.SetFromUserInput(pszDefTo);
	if (erro != OGRERR_NONE) {
		s.setError("output crs is not valid");
		return s;
	}

	//CPLSetConfigOption("OGR_CT_FORCE_TRADITIONAL_GIS_ORDER", "YES");
	OGRCoordinateTransformation *poCT;
	poCT = OGRCreateCoordinateTransformation(&source, &target);

	if( poCT == NULL )	{
		s.setError( "Cannot do this transformation" );
		return(s);
	}

	s.setSRS(crs);
	s.df = df;
	std::vector<unsigned> keeprows;


	for (size_t i=0; i < size(); i++) {
		SpatGeom g = getGeom(i);
		SpatGeom gg;
		gg.gtype = g.gtype;
		for (size_t j=0; j < g.size(); j++) {
			SpatPart p = g.getPart(j);
			if (poCT->Transform(p.x.size(), &p.x[0], &p.y[0]) ) {
				SpatPart pp(p.x, p.y);
				if (p.hasHoles()) {
					for (size_t k=0; k < p.nHoles(); k++) {
						SpatHole h = p.getHole(k);
						if (poCT->Transform(h.x.size(), &h.x[0], &h.y[0])) {
							pp.addHole(h.x, h.y);
						}
					}
				}
				gg.addPart(pp);
			}
		}
		keeprows.push_back(i);
		s.addGeom(gg);
	}
	s.df = df.subset_rows(keeprows);
	OCTDestroyCoordinateTransformation(poCT);

	#endif
	return s;
}


#endif

