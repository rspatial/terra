// Copyright (c) 2018-2026  Robert J. Hijmans
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
#include <sstream>
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

#if GDAL_VERSION_MAJOR >= 3
#include "proj.h"
#endif

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


bool SpatSRS::m_dist(double &m, bool lonlat, std::string unit) {
	m = 1;
	if (!lonlat) {
		m = to_meter();
		m = std::isnan(m) ? 1 : m;
	}
	std::vector<std::string> ss {"m", "km"};
	if (std::find(ss.begin(), ss.end(), unit) == ss.end()) {
		return false;
	}
	if (unit == "km")	{
		m /= 1000;
	}
	return true;
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



bool SpatSRS::is_not_lonlat() {
	if (wkt.size() < 2) {
		return false;
	}
	OGRSpatialReference x;
	OGRErr erro = x.SetFromUserInput(wkt.c_str());
	if (erro != OGRERR_NONE) {
		return false;
	}
	return !x.IsGeographic();
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




void EmptyErrorHandler(CPLErr eErrClass, int errNo, const char *msg) {
    // do nothing
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

	OGRCoordinateTransformation *poCT = nullptr;
    CPLPushErrorHandler(EmptyErrorHandler);
	try {
		poCT = OGRCreateCoordinateTransformation(&source, &target);
	} catch(...) {
		OCTDestroyCoordinateTransformation(poCT);
		return false;
	}
    CPLPopErrorHandler();
	if (poCT == NULL) {
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


void transform_coordinates_partial(std::vector<double> &x, std::vector<double> &y, OGRCoordinateTransformation *poCT) {
	std::vector<double> X, Y;
	X.reserve(x.size());
	Y.reserve(y.size());
	std::vector<size_t> fails;
	for (size_t i=0; i < x.size(); i++) {
		if( poCT->Transform( 1, &x[i], &y[i] ) ) {
			X.push_back(x[i]);
			Y.push_back(y[i]);
		}
	}
	x = X;
	y = Y;
}


SpatVector SpatVector::project(std::string crs, bool partial, std::string pipeline,
		std::vector<double> AOI, double desired_accuracy, bool allow_ballpark) {

	bool remove_empty = false;

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

	bool use_pipeline = !pipeline.empty();
	if (!use_pipeline) {
		const char *pszDefTo = crs.c_str();
		erro = target.SetFromUserInput(pszDefTo);
		if (erro != OGRERR_NONE) {
			s.setError("output crs is not valid");
			return s;
		}
	}

	OGRCoordinateTransformation *poCT;
#if GDAL_VERSION_NUM >= 3000000
	OGRCoordinateTransformationOptions opts;
	bool use_opts = false;
	if (use_pipeline) {
		if (!opts.SetCoordinateOperation(pipeline.c_str(), false)) {
			s.setError("pipeline not accepted");
			return s;
		}
		use_opts = true;
	}
	if (AOI.size() == 4) {
		if (!opts.SetAreaOfInterest(AOI[0], AOI[1], AOI[2], AOI[3])) {
			s.setError("area of interest not accepted");
			return s;
		}
		use_opts = true;
	}
#if GDAL_VERSION_NUM >= 3030000
	if (desired_accuracy >= 0) {
		opts.SetDesiredAccuracy(desired_accuracy);
		use_opts = true;
	}
	if (!allow_ballpark) {
		opts.SetBallparkAllowed(false);
		use_opts = true;
	}
#endif
	OGRSpatialReference *pTarget = use_pipeline ? nullptr : &target;
	if (use_opts) {
		poCT = OGRCreateCoordinateTransformation(&source, pTarget, opts);
	} else {
		poCT = OGRCreateCoordinateTransformation(&source, &target);
	}
#else
	if (use_pipeline || AOI.size() == 4) {
		s.setError("pipeline and AOI require GDAL >= 3");
		return s;
	}
	poCT = OGRCreateCoordinateTransformation(&source, &target);
#endif

	if( poCT == NULL )	{
		s.setError( "Cannot do this transformation" );
		return(s);
	}
	
	s.setSRS(crs);
	s.df = df;
	std::vector<size_t> keeprows;
	
	if (partial) {
		#if GDAL_VERSION_MAJOR >= 2 && GDAL_VERSION_MINOR > 1
		poCT->SetEmitErrors(false);
		#endif
		std::string gt = type();
		size_t minpts = gt == "polygons" ? 3 : (gt == "lines" ? 2 : 1);
		for (size_t i=0; i < size(); i++) {
			SpatGeom g = getGeom(i);
			SpatGeom gg;
			gg.gtype = g.gtype;
			bool empty = true;
			for (size_t j=0; j < g.size(); j++) {
				SpatPart p = g.getPart(j);
				transform_coordinates_partial(p.x, p.y, poCT);	
				if (p.x.size() >= minpts) {
					SpatPart pp(p.x, p.y);
					if (p.hasHoles()) {
						for (size_t k=0; k < p.nHoles(); k++) {
							SpatHole h = p.getHole(k);
							transform_coordinates_partial(h.x, h.y, poCT);				
							if (h.x.size() >= 3) {
								pp.addHole(h.x, h.y);
							}
						}
					}
					gg.addPart(pp);
					empty = false;
				}
			}
			if (empty) {
				if (remove_empty) {
					keeprows.push_back(i);
				} else {
					s.addGeom(gg);
				}
			} else {	
				s.addGeom(gg);
			}
		}
		
		
	} else {

		for (size_t i=0; i < size(); i++) {
			SpatGeom g = getGeom(i);
			SpatGeom gg;
			gg.gtype = g.gtype;
			bool empty = true;
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
					empty = false;
				}
			}
			if (empty) {
				if (remove_empty) {
					keeprows.push_back(i);
				} else {
					s.addGeom(gg);
				}
			} else {
				s.addGeom(gg);
			}
		}
	}
	
	OCTDestroyCoordinateTransformation(poCT);

	if (size() > 0 && std::isnan(s.extent.xmin)) {
		s.setError("transformation failed (may need a grid file that could not be downloaded)");
		return s;
	}

	if (remove_empty) {
		s.df = df.subset_rows(keeprows);
	} else {
		s.df = df;		
	}
	#endif
	return s;
}


SpatDataFrame get_proj_pipelines(std::string source_crs, std::string target_crs,
		std::string authority, std::vector<double> AOI, std::string use, std::string grid_availability, 
		double desired_accuracy, bool strict_containment, bool axis_order_authority_compliant) {

	SpatDataFrame out;

#if GDAL_VERSION_MAJOR < 3
	out.setError("GDAL >= 3 and PROJ >= 7.1 required");
	return out;
#else
	#if PROJ_VERSION_MAJOR < 7 || (PROJ_VERSION_MAJOR == 7 && PROJ_VERSION_MINOR < 1)
	out.setError("PROJ >= 7.1 required");
	return out;
#endif

	PJ *src = proj_create(PJ_DEFAULT_CTX, source_crs.c_str());
	if (src == nullptr) {
		out.setError("invalid source crs");
		return out;
	}
	PJ *tgt = proj_create(PJ_DEFAULT_CTX, target_crs.c_str());
	if (tgt == nullptr) {
		proj_destroy(src);
		out.setError("invalid target crs");
		return out;
	}

	const char *auth = authority.empty() ? nullptr : authority.c_str();
	PJ_OPERATION_FACTORY_CONTEXT *ctx =
		proj_create_operation_factory_context(PJ_DEFAULT_CTX, auth);

	if (desired_accuracy >= 0.0) {
		proj_operation_factory_context_set_desired_accuracy(PJ_DEFAULT_CTX, ctx, desired_accuracy);
	}

	if (AOI.size() == 4) {
		proj_operation_factory_context_set_area_of_interest(PJ_DEFAULT_CTX, ctx,
			AOI[0], AOI[1], AOI[2], AOI[3]);
	} else {
		if (use == "BOTH") {
			proj_operation_factory_context_set_crs_extent_use(PJ_DEFAULT_CTX, ctx, PJ_CRS_EXTENT_BOTH);
		} else if (use == "INTERSECTION") {
			proj_operation_factory_context_set_crs_extent_use(PJ_DEFAULT_CTX, ctx, PJ_CRS_EXTENT_INTERSECTION);
		} else if (use == "SMALLEST") {
			proj_operation_factory_context_set_crs_extent_use(PJ_DEFAULT_CTX, ctx, PJ_CRS_EXTENT_SMALLEST);
		} else {
			proj_operation_factory_context_set_crs_extent_use(PJ_DEFAULT_CTX, ctx, PJ_CRS_EXTENT_NONE);
		}
	}

	if (strict_containment) {
		proj_operation_factory_context_set_spatial_criterion(PJ_DEFAULT_CTX, ctx,
			PROJ_SPATIAL_CRITERION_STRICT_CONTAINMENT);
	} else {
		proj_operation_factory_context_set_spatial_criterion(PJ_DEFAULT_CTX, ctx,
			PROJ_SPATIAL_CRITERION_PARTIAL_INTERSECTION);
	}

	if (grid_availability == "DISCARD") {
		proj_operation_factory_context_set_grid_availability_use(PJ_DEFAULT_CTX, ctx,
			PROJ_GRID_AVAILABILITY_DISCARD_OPERATION_IF_MISSING_GRID);
	} else if (grid_availability == "IGNORED") {
		proj_operation_factory_context_set_grid_availability_use(PJ_DEFAULT_CTX, ctx,
			PROJ_GRID_AVAILABILITY_IGNORED);
	} else if (grid_availability == "AVAILABLE") {
		proj_operation_factory_context_set_grid_availability_use(PJ_DEFAULT_CTX, ctx,
			PROJ_GRID_AVAILABILITY_KNOWN_AVAILABLE);
	} else {
		proj_operation_factory_context_set_grid_availability_use(PJ_DEFAULT_CTX, ctx,
			PROJ_GRID_AVAILABILITY_USED_FOR_SORTING);
	}

	PJ_OBJ_LIST *ops = proj_create_operations(PJ_DEFAULT_CTX, src, tgt, ctx);
	int n = proj_list_get_count(ops);

	std::vector<std::string> v_id(n);
	std::vector<std::string> v_description(n);
	std::vector<std::string> v_definition(n);
	std::vector<int8_t>      v_has_inverse(n);
	std::vector<double>      v_accuracy(n);
	std::vector<long>        v_grid_count(n);
	std::vector<int8_t>      v_instantiable(n);

	for (int i = 0; i < n; i++) {
		PJ *op = proj_list_get(PJ_DEFAULT_CTX, ops, i);

		PJ_PROJ_INFO info = proj_pj_info(op);
		v_id[i] = (info.id != nullptr) ? info.id : "";
		v_description[i] = (info.description != nullptr) ? info.description : "";
		v_has_inverse[i] = (info.has_inverse != 0) ? 1 : 0;
		v_accuracy[i] = info.accuracy;

		// Always extract the authority-compliant definition (before
		// any normalization) so that GDAL's OGR_CT_FORCE_TRADITIONAL_GIS_ORDER
		// can handle axis order without double-swapping.
		std::string def;
		if (info.definition != nullptr && info.definition[0] != '\0') {
			def = info.definition;
		} else {
			const char *ps = proj_as_proj_string(PJ_DEFAULT_CTX, op,
				PJ_PROJ_5, nullptr);
			if (ps != nullptr) def = ps;
		}
		// normalize: ensure every token starts with "+"
		{
			std::string norm;
			std::istringstream iss(def);
			std::string token;
			while (iss >> token) {
				if (!norm.empty()) norm += ' ';
				if (token[0] != '+') norm += '+';
				norm += token;
			}
			def = norm;
		}
		v_definition[i] = def;
		v_grid_count[i] = proj_coordoperation_get_grid_used_count(PJ_DEFAULT_CTX, op);
		v_instantiable[i] = proj_coordoperation_is_instantiable(PJ_DEFAULT_CTX, op) ? 1 : 0;

		proj_destroy(op);
	}

	proj_list_destroy(ops);
	proj_operation_factory_context_destroy(ctx);
	proj_destroy(tgt);
	proj_destroy(src);

	//out.add_column(v_id, "id");
	out.add_column(v_description, "description");
	out.add_column(v_definition, "definition");
	out.add_column(v_has_inverse, "has_inverse");
	out.add_column(v_accuracy, "accuracy");
	out.add_column(v_grid_count, "grid_count");
	out.add_column(v_instantiable, "instantiable");

	return out;
#endif
}


#endif

