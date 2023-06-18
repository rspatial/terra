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


#include "spatVector.h"
#include "file_utils.h"

#include "ogrsf_frmts.h"
#include "ogr_spatialref.h"

#include "crs.h"
#include "NA.h"

#include "string_utils.h"

std::string geomType(OGRLayer *poLayer) {
	std::string s = "";
    poLayer->ResetReading();
    OGRFeature *poFeature;
    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		const char* gname = poGeometry->getGeometryName();
		s = gname;
		break;
	}
	OGRFeature::DestroyFeature( poFeature );
	return s;
}


SpatDataFrame readAttributes(OGRLayer *poLayer, bool as_proxy) {
	SpatDataFrame df;

    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
	size_t nfields = poFDefn->GetFieldCount();
	if (nfields == 0) return df;

	OGRFieldType ft;
    poLayer->ResetReading();
    OGRFeature *poFeature;
	OGRFieldDefn *poFieldDefn;
	df.resize_cols(nfields);
	bool first = true;
	unsigned dtype;
	long longNA = NA<long>::value;

    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		if (first) {
			for (size_t i = 0; i < nfields; i++ ) {
				poFieldDefn = poFDefn->GetFieldDefn(i);
				std::string fname = poFieldDefn->GetNameRef();
				ft = poFieldDefn->GetType();
				if (ft == OFTReal) {
					dtype = 0;
				} else if ((ft == OFTInteger) | (ft == OFTInteger64)) {
					if (poFieldDefn->GetSubType() == OFSTBoolean) {
						dtype = 3;
					} else {
						dtype = 1;
					}
				} else {
					dtype = 2;
				}
				df.add_column(dtype, fname);
			}
			first = false;
		}

		for (size_t i = 0; i < nfields; i++ ) {
			poFieldDefn = poFDefn->GetFieldDefn( i );
			unsigned j = df.iplace[i];
			int not_null = poFeature->IsFieldSetAndNotNull(i);
			switch( poFieldDefn->GetType() ) {
				case OFTReal:
					if (not_null) {
						df.dv[j].push_back(poFeature->GetFieldAsDouble(i));
					} else {
						df.dv[j].push_back(NAN);
					}
					break;
				case OFTInteger:
					if (poFieldDefn->GetSubType() == OFSTBoolean) {
						if (not_null) {
							df.bv[j].push_back(poFeature->GetFieldAsInteger(i));
						} else {
							df.bv[j].push_back(2);
						}
					} else {
						if (not_null) {
							df.iv[j].push_back(poFeature->GetFieldAsInteger(i));
						} else {
							df.iv[j].push_back(longNA);
						}
					}
					break;
				case OFTInteger64:
					if (not_null) {
						df.iv[j].push_back(poFeature->GetFieldAsInteger64(i));
					} else {
						df.iv[j].push_back(longNA);
					}
					break;
	//          case OFTString:
				default:
					if (not_null) {
						df.sv[j].push_back(poFeature->GetFieldAsString(i));
					} else {
						df.sv[j].push_back(df.NAS);
					}
					break;
			}
		}
		OGRFeature::DestroyFeature(poFeature);
		if (as_proxy) break;
	}
	return df;
}


/*
std::string getDs_WKT(GDALDataset *poDataset) {
	std::string wkt = "";
	char *cp;
#if GDAL_VERSION_MAJOR >= 3
	const OGRSpatialReference *srs = poDataset->GetSpatialRef();
	const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
	OGRErr err = srs->exportToWkt(&cp, options);
	if (err == OGRERR_NONE) {
		wkt = std::string(cp);
		CPLFree(cp);
	}
#else
	const char *pszSrc = GDALGetProjectionRef( poDataset );
	if (pszSrc != NULL) {
		wkt = std::string(pszSrc);
	}

//	if (poDataset->GetProjectionRef() != NULL) {
//		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
//		OGRErr err = oSRS.exportToPrettyWkt(&cp);
//		if (err == OGRERR_NONE) {
//			wkt = std::string(cp);
//			CPLFree(cp);
//		}
//	}

#endif
	return wkt;
}

std::string getDs_PRJ(GDALDataset *poDataset) {
	std::string prj = "";
#if GDAL_VERSION_MAJOR >= 3
	char *cp;
	const OGRSpatialReference *srs = poDataset->GetSpatialRef();
	OGRErr err = srs->exportToProj4(&cp);
	if (err == OGRERR_NONE) {
		prj = std::string(cp);
		CPLFree(cp);
	}
#else
	if( poDataset->GetProjectionRef() != NULL ) {
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		char *pszPRJ = NULL;
		oSRS.exportToProj4(&pszPRJ);
		prj = pszPRJ;
	}
#endif
	return prj;
}
*/


SpatGeom getPointGeom(OGRGeometry *poGeometry) {
	SpatGeom g(points);
	if (poGeometry->IsEmpty()) {
		return g;
	}
	#if GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(2,3,0)
	OGRPoint *poPoint = poGeometry->toPoint();
	#else
	OGRPoint *poPoint = (OGRPoint *) poGeometry;
	#endif
	double x = poPoint->getX();
	double y = poPoint->getY();
	SpatPart p(x, y);
	g.addPart(p);
	return g;
}

SpatGeom getMultiPointGeom(OGRGeometry *poGeometry) {
	OGRMultiPoint *poMultipoint = ( OGRMultiPoint * )poGeometry;
	unsigned ng = poMultipoint->getNumGeometries();
	std::vector<double> X(ng);
	std::vector<double> Y(ng);
	for (size_t i=0; i<ng; i++) {
	   	OGRGeometry *poMpGeometry = poMultipoint->getGeometryRef(i);
		#if GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(2,3,0)
			OGRPoint *poPoint = poMpGeometry->toPoint();
		#else
			OGRPoint *poPoint = (OGRPoint *) poMpGeometry;
		#endif
		X[i] = poPoint->getX();
		Y[i] = poPoint->getY();
	}
	SpatPart p(X, Y);
	SpatGeom g(points);
	g.addPart(p);
	return g;
}


SpatGeom getLinesGeom(OGRGeometry *poGeometry) {

	OGRLineString *poGeom = (OGRLineString *) poGeometry;
	unsigned np = poGeom->getNumPoints();
	std::vector<double> X(np);
	std::vector<double> Y(np);
	OGRPoint ogrPt;
	for (size_t i=0; i<np; i++) {
		poGeom->getPoint(i, &ogrPt);
		X[i] = ogrPt.getX();
		Y[i] = ogrPt.getY();
	}
	SpatPart p(X, Y);
	SpatGeom g(lines);
	g.addPart(p);
	return g;
}

SpatGeom getMultiLinesGeom(OGRGeometry *poGeometry) {
	SpatGeom g(lines);
	OGRMultiLineString *poGeom = ( OGRMultiLineString * )poGeometry;
	unsigned ng = poGeom->getNumGeometries();
	OGRPoint ogrPt;
	for (size_t i=0; i<ng; i++) {
		OGRGeometry *poLineGeometry = poGeom->getGeometryRef(i);
		OGRLineString *poLine = ( OGRLineString * )poLineGeometry;
		unsigned np = poLine->getNumPoints();
		std::vector<double> X(np);
		std::vector<double> Y(np);
		for (size_t j=0; j<np; j++ ) {
			poLine->getPoint(j, &ogrPt);
			X[j] = ogrPt.getX();
			Y[j] = ogrPt.getY();
		}
		SpatPart p(X, Y);
		g.addPart(p);
	}
	return g;
}

//#include "Rcpp.h"

SpatGeom getPolygonsGeom(OGRGeometry *poGeometry) {
	SpatGeom g(polygons);
	OGRPoint ogrPt;
//	OGRwkbGeometryType geomtype = poGeometry->getGeometryType();
//	if ( geomtype == wkbPolygon ) {
		OGRPolygon *poGeom = ( OGRPolygon * )poGeometry;
		OGRLinearRing *poRing = poGeom->getExteriorRing();
		unsigned np = poRing->getNumPoints();
		std::vector<double> X(np);
		std::vector<double> Y(np);
		for (size_t i=0; i<np; i++) {
			poRing->getPoint(i, &ogrPt);
			X[i] = ogrPt.getX();
			Y[i] = ogrPt.getY();
		}
		SpatPart p(X, Y);
		unsigned nh = poGeom->getNumInteriorRings();
		for (size_t i=0; i<nh; i++) {
			OGRLinearRing *poHole = poGeom->getInteriorRing(i);
			unsigned np = poHole->getNumPoints();
			std::vector<double> X(np);
			std::vector<double> Y(np);
			for (size_t j=0; j<np; j++) {
				poHole->getPoint(j, &ogrPt);
				X[j] = ogrPt.getX();
				Y[j] = ogrPt.getY();
			}
			p.addHole(X, Y);
		}
		g.addPart(p);
//	}
	return g;
}


SpatGeom getMultiPolygonsGeom(OGRGeometry *poGeometry) {
	OGRMultiPolygon *poGeom = ( OGRMultiPolygon * )poGeometry;
	OGRPoint ogrPt;
	unsigned ng = poGeom->getNumGeometries();
	SpatGeom g(polygons);
	for (size_t i=0; i<ng; i++) {
		OGRGeometry *poPolygonGeometry = poGeom->getGeometryRef(i);
		OGRPolygon *poPolygon = ( OGRPolygon * )poPolygonGeometry;
		OGRLinearRing *poRing = poPolygon->getExteriorRing();
		unsigned np = poRing->getNumPoints();
		std::vector<double> X(np);
		std::vector<double> Y(np);
		for (size_t j=0; j<np; j++ ) {
			poRing->getPoint(j, &ogrPt);
			X[j] = ogrPt.getX();
			Y[j] = ogrPt.getY();
		}
		SpatPart p(X, Y);
		unsigned nh = poPolygon->getNumInteriorRings();
		for (size_t j=0; j<nh; j++) {
			OGRLinearRing *poHole = poPolygon->getInteriorRing(j);
			np = poHole->getNumPoints();
			std::vector<double> X(np);
			std::vector<double> Y(np);
			for (size_t k = 0; k < np; k++ ) {
				poHole->getPoint(k, &ogrPt);
				X[k] = ogrPt.getX();
				Y[k] = ogrPt.getY();
			}
			p.addHole(X, Y);
		}
		g.addPart(p);
	}
	return g;
}

std::vector<std::string> SpatVector::layer_names(std::string filename) {

	std::vector<std::string> out;

	if (filename.empty()) {
		setError("empty filename");
		return out;
	}
	// a gdb is a folder...
	//if (!file_exists(filename)) {
	//	setError("file does not exist");
	//	return out;
	//}

    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR,
				NULL, NULL, NULL ));

    if( poDS == NULL ) {
        setError("Cannot open this dataset" );
		return out;
    }

	size_t n = poDS->GetLayerCount();
	out.reserve(n);
	for (size_t i=0; i<n; i++) {
		OGRLayer *poLayer = poDS->GetLayer(i);
		if (poLayer == NULL) {
			out.push_back("");
		} else {
			out.push_back((std::string)poLayer->GetName());
		}
	}

	GDALClose(poDS);
	return out;
}

SpatGeom emptyGeom() {
	SpatGeom g;
	g.gtype = null;
	g.extent.xmin=NAN;
	g.extent.xmax=NAN;
	g.extent.ymin=NAN;
	g.extent.ymax=NAN;
	return g;
}


bool layerQueryFilter(GDALDataset *&poDS, OGRLayer *&poLayer, std::string &layer, std::string &query, std::vector<double> &extent, SpatVector &filter, std::string &errmsg, std::vector<std::string> &wrms) {

	if (query.empty()) {
		if (layer.empty()) {
			#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
				// do nothing
			#else
			std::vector<std::string> lyrnms;
			for ( auto&& poLayer: poDS->GetLayers() ) {
				lyrnms.push_back((std::string)poLayer->GetName());
			}
			if (lyrnms.size() > 1) {
				std::string lyrsel = lyrnms[0];
				lyrnms.erase(lyrnms.begin());
				std::string ccat = concatenate(lyrnms, ", ");
				wrms.push_back("Reading layer: " + lyrsel + "\nOther layers: " + ccat);
			}
			#endif

			poLayer = poDS->GetLayer(0);
			if (poLayer == NULL) {
				errmsg = "dataset has no layers";
				return false;
			}
		} else {
			poLayer = poDS->GetLayerByName(layer.c_str());
			if (poLayer == NULL) {
				errmsg = layer + " is not a valid layer name";
			#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
				// do nothing
			#else
				errmsg += "\nChoose one of: ";
				for ( auto&& poLayer: poDS->GetLayers() ) {
					errmsg += (std::string)poLayer->GetName() + ", ";
				}
				errmsg = errmsg.substr(0, errmsg.size()-2);
			#endif
				return false;
			}
		}
	} else {
		poLayer = poDS->ExecuteSQL(query.c_str(), NULL, NULL);
		if (poLayer == NULL) {
			errmsg = "Query failed";
			return false;
		}
	}
	
	if (filter.nrow() > 0) {
		if (filter.type() != "polygons") {
			filter = filter.hull("convex");
		} else if (filter.nrow() > 1) {
			filter = filter.aggregate(true);
		}
		GDALDataset *filterDS = filter.write_ogr("", "lyr", "Memory", false, true, std::vector<std::string>());
		if (filter.hasError()) {
			//setError(filter.getError());
			GDALClose(filterDS);
			errmsg = "filter has error";
			return false;
		}
		OGRLayer *fLayer = filterDS->GetLayer(0);
		fLayer->ResetReading();
		OGRFeature *fFeature = fLayer->GetNextFeature();
		if (fFeature != NULL ) {
			OGRGeometry *fGeometry = fFeature->StealGeometry();
			poLayer->SetSpatialFilter(fGeometry);
			OGRGeometryFactory::destroyGeometry(fGeometry);
		}
		OGRFeature::DestroyFeature( fFeature );
		GDALClose(filterDS);
	} else if (!extent.empty()) {
		poLayer->SetSpatialFilterRect(extent[0], extent[2], extent[1], extent[3]);
	}

	return true;
}


bool SpatVector::read_ogr(GDALDataset *&poDS, std::string layer, std::string query, std::vector<double> extent, SpatVector filter, bool as_proxy, std::string what) {

	if (poDS == NULL) {
		setError("dataset is empty");
		return false;
	}

	std::string crs = "";
	OGRLayer *poLayer;
	poLayer = poDS->GetLayer(0);

	read_query = query;
	read_extent = extent;
	std::string errmsg;
	std::vector<std::string> wrnmsg;

	if (!layerQueryFilter(poDS, poLayer, layer, query, extent, filter, errmsg, wrnmsg)) {
		setError(errmsg);
		return false;
	} else if (!wrnmsg.empty()) {
		for (size_t i=0; i < wrnmsg.size(); i++) addWarning(wrnmsg[i]);
	}

	OGRSpatialReference *poSRS = poLayer->GetSpatialRef();
	if (poSRS) {
		char *psz = NULL;
	#if GDAL_VERSION_MAJOR >= 3
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = poSRS->exportToWkt(&psz, options);
	#else
		OGRErr err = poSRS->exportToWkt(&psz);
	#endif
		if (err == OGRERR_NONE) {
			crs = psz;
		}
		setSRS(crs);
		CPLFree(psz);
	}

	if (what != "geoms") { 
		df = readAttributes(poLayer, as_proxy);
	}
	if (what == "attributes") {
		if (!query.empty()) {
			poDS->ReleaseResultSet(poLayer);
		}
		return true;
	}

	//const char* lname = poLayer->GetName();
	OGRwkbGeometryType wkbgeom = wkbFlatten(poLayer->GetGeomType());
	OGRFeature *poFeature;

	poLayer->ResetReading();
	poFeature = poLayer->GetNextFeature();
	if (poFeature != NULL) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		if (poGeometry != NULL) {
			if (poGeometry->Is3D()) {
				addWarning("Z coordinates ignored");
			}
			if (poGeometry->IsMeasured()) {
				addWarning("M coordinates ignored");
			}
		}
	}
	source_layer = poLayer->GetName();

	if (as_proxy) {
		SpatGeom g;
		if ((wkbgeom == wkbPoint) | (wkbgeom == wkbMultiPoint)) {
			//SpatPart p(0,0);
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			if (poGeometry != NULL) {
				if ( wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ) {
					g = getPointGeom(poGeometry);
				} else {
					g = getMultiPointGeom(poGeometry);
				}
			} else {
				g = emptyGeom();
			}
			addGeom(g);
			OGRFeature::DestroyFeature( poFeature );
		} else if (wkbgeom == wkbLineString || wkbgeom == wkbMultiLineString) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			if (poGeometry != NULL) {
				if (wkbFlatten ( poGeometry ->getGeometryType() ) == wkbLineString) {
					g = getLinesGeom(poGeometry);
				} else {
					g = getMultiLinesGeom(poGeometry);
				}
			} else {
				g = emptyGeom();
			}
			addGeom(g);
			OGRFeature::DestroyFeature( poFeature );
		} else if ( wkbgeom == wkbPolygon || wkbgeom == wkbMultiPolygon) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			if (poGeometry != NULL) {
				wkbgeom = wkbFlatten(poGeometry->getGeometryType());
				if (wkbgeom == wkbPolygon) {
					g = getPolygonsGeom(poGeometry);
				} else if (wkbgeom == wkbMultiPolygon ) {
					g = getMultiPolygonsGeom(poGeometry);
				} // else ?
			} else {
				g = emptyGeom();
			}
			addGeom(g);
			OGRFeature::DestroyFeature( poFeature );
		} else if (wkbgeom != wkbNone) {
			const char *geomtypechar = OGRGeometryTypeToName(wkbgeom);
			std::string strgeomtype = geomtypechar;
			std::string s = "cannot read this geometry type: "+ strgeomtype;
			setError(s);
			return false;
		}
		geom_count = poLayer->GetFeatureCount();
		is_proxy = true;
		return true;
	}

	OGRFeature::DestroyFeature( poFeature );
	poLayer->ResetReading();
	SpatGeom g;

	if ((wkbgeom == wkbPoint) | (wkbgeom == wkbMultiPoint)) {
		//SpatPart p(0,0);
		while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			if (poGeometry != NULL) {
				if ( wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ) {
					g = getPointGeom(poGeometry);
				} else {
					g = getMultiPointGeom(poGeometry);
				}
			} else {
				g = emptyGeom();
			}
			addGeom(g);
			OGRFeature::DestroyFeature( poFeature );
		}
	} else if (wkbgeom == wkbLineString || wkbgeom == wkbMultiLineString) {
		while ( (poFeature = poLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			if (poGeometry != NULL) {
				if (wkbFlatten ( poGeometry ->getGeometryType() ) == wkbLineString) {
					g = getLinesGeom(poGeometry);
				} else {
					g = getMultiLinesGeom(poGeometry);
				}
			} else {
				g = emptyGeom();
			}
			addGeom(g);
			OGRFeature::DestroyFeature( poFeature );
		}
	} else if ( wkbgeom == wkbPolygon || wkbgeom == wkbMultiPolygon) {
		while ( (poFeature = poLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			if (poGeometry != NULL) {
				wkbgeom = wkbFlatten(poGeometry->getGeometryType());
				if (wkbgeom == wkbPolygon) {
					g = getPolygonsGeom(poGeometry);
				} else if (wkbgeom == wkbMultiPolygon ) {
					g = getMultiPolygonsGeom(poGeometry);
				}
			} else {
				g = emptyGeom();
			}
			addGeom(g);
			OGRFeature::DestroyFeature( poFeature );
		}
	} else if (wkbgeom == wkbUnknown) {
		SpatVectorCollection sv;
		std::vector<double> dempty;
		SpatVector filter2;
		sv.read_ogr(poDS, "", "", dempty, filter2); 

		if (sv.size() > 0) {
			*this = sv.v[0];
			if (sv.v.size() > 1) {
				std::string gt = type();
				addWarning("returning " + gt + " ignoring additional geometry types. Use 'svc' to get all geometries");
			}
			return true;
		}
		if (sv.hasError()) {
			setError(sv.getError());
		}
		return false;
	} else if (wkbgeom != wkbNone) {
		const char *geomtypechar = OGRGeometryTypeToName(wkbgeom);
		std::string strgeomtype = geomtypechar;
		std::string s = "cannot read this geometry type: "+ strgeomtype;
		setError(s);
		return false;
	}


	if (!query.empty()) {
		poDS->ReleaseResultSet(poLayer);
	}

 	return true;
}


bool SpatVector::read(std::string fname, std::string layer, std::string query, std::vector<double> extent, SpatVector filter, bool as_proxy, std::string what) {
    //OGRRegisterAll();
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL ));
    if( poDS == NULL ) {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else {
			setError("Cannot open this file as a SpatVector: " + fname);
		}
		return false;
    }
	bool success = read_ogr(poDS, layer, query, extent, filter, as_proxy, what);
	if (poDS != NULL) GDALClose( poDS );
	source = fname;
	return success;
}

SpatVector SpatVector::fromDS(GDALDataset *poDS) {
	SpatVector out, fvct;
	std::vector<double> fext;
	out.read_ogr(poDS, "", "", fext, fvct, false, "");
	return out;
}


SpatVector::SpatVector(std::vector<std::string> wkt) {

	OGRGeometryFactory ogr;

	SpatGeom g;
	bool haveGeomt = false;
	SpatGeomType geomt = null;
	for (size_t i=0; i<wkt.size(); i++) {
		if (wkt[i] == "EMPTY") {
			g = emptyGeom();
			addGeom(g);
			continue;
		}

		OGRGeometry *poGeometry;

#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
		char *cstring = &wkt[i][0];
		std::vector<char*> cstr = { cstring };
		OGRErr err = ogr.createFromWkt(&cstr[0], NULL, &poGeometry );
#else
		const char* pszWKT = wkt[i].c_str();
		OGRErr err = ogr.createFromWkt( pszWKT, NULL, &poGeometry );
#endif

		if (err == OGRERR_NONE) {
			//const char* gname = poGeometry->getGeometryName();
			if (poGeometry != NULL) {
				OGRwkbGeometryType gtype = wkbFlatten(poGeometry->getGeometryType());
				if ( gtype == wkbPoint ) {
					g = getPointGeom(poGeometry);
				} else if ( gtype == wkbMultiPoint ) {
					g = getMultiPointGeom(poGeometry);
				} else if (gtype == wkbLineString) {
					g = getLinesGeom(poGeometry);
				} else if (gtype == wkbMultiLineString) {
					g = getMultiLinesGeom(poGeometry);
				} else if (gtype == wkbPolygon) {
					g = getPolygonsGeom(poGeometry);
				} else if (gtype == wkbMultiPolygon ) {
					g = getMultiPolygonsGeom(poGeometry);
				} else {
					const char *geomtypechar = OGRGeometryTypeToName(gtype);
					std::string strgeomtype = geomtypechar;
					std::string s = "cannot read geometry type: "+ strgeomtype;
					setError(s);
					return;
				}
				if (!haveGeomt) {
					haveGeomt = true;
					geomt = g.gtype;
				} else if (geomt != g.gtype) {
					setError("a SpatVector can only have a single geometry type");
					return;			
				}
				addGeom(g);
				OGRGeometryFactory::destroyGeometry(poGeometry);
			}
		} else {
			setError("not WKT");
			return;
		}
	}
}

bool SpatVectorCollection::read_ogr(GDALDataset *&poDS, std::string layer, std::string query, std::vector<double> extent, SpatVector filter) {
	
	OGRLayer *poLayer;
	poLayer = poDS->GetLayer(0);
	
	std::string errmsg;
	std::vector<std::string> wrnmsg;

	if (!layerQueryFilter(poDS, poLayer, layer, query, extent, filter, errmsg, wrnmsg)) {
		setError(errmsg);
		return false;
	} else if (!wrnmsg.empty()) {
		for (size_t i=0; i < wrnmsg.size(); i++) addWarning(wrnmsg[i]);
	}
	std::string crs = "";
	OGRSpatialReference *poSRS = poLayer->GetSpatialRef();
	if (poSRS) {
		char *psz = NULL;
	#if GDAL_VERSION_MAJOR >= 3
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = poSRS->exportToWkt(&psz, options);
	#else
		OGRErr err = poSRS->exportToWkt(&psz);
	#endif
		if (err == OGRERR_NONE) {
			crs = psz;
		}
		CPLFree(psz);
	}

	//const char* lname = poLayer->GetName();
//	OGRwkbGeometryType wkbgeom = wkbFlatten(poLayer->GetGeomType());
	OGRFeature *poFeature;

	poLayer->ResetReading();
	poFeature = poLayer->GetNextFeature();
	if (poFeature != NULL) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		if (poGeometry != NULL) {
			if (poGeometry->Is3D()) {
				addWarning("Z coordinates ignored");
			}
			if (poGeometry->IsMeasured()) {
				addWarning("M coordinates ignored");
			}
		}
	}

	std::string source_layer = poLayer->GetName();

	OGRFeature::DestroyFeature( poFeature );
	SpatDataFrame df = readAttributes(poLayer, false);
	poLayer->ResetReading();

	SpatVector points, lines, polygons;
	std::vector<unsigned> pnt, lin, pol;
	SpatGeom g;
	size_t i = 0;
	while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		if (poGeometry != NULL) {
			OGRwkbGeometryType wkb = wkbFlatten(poGeometry->getGeometryType());
			if (wkb  == wkbPoint ) {
				g = getPointGeom(poGeometry);
				points.addGeom(g);
				pnt.push_back(i);
			} else if ((wkb  == wkbMultiPoint) || (wkb  == wkbMultiPointZM) || (wkb  == wkbMultiPointM)) {
				g = getMultiPointGeom(poGeometry);
				points.addGeom(g);
				pnt.push_back(i);
			} else if (wkb == wkbLineString) {
				g = getLinesGeom(poGeometry);
				lines.addGeom(g);
				lin.push_back(i);
			} else if ((wkb == wkbMultiLineString) || (wkb == wkbMultiLineStringZM) || (wkb == wkbMultiLineStringM)) {
				g = getMultiLinesGeom(poGeometry);
				lines.addGeom(g);
				lin.push_back(i);
			} else if (wkb == wkbPolygon) {
				g = getPolygonsGeom(poGeometry);
				polygons.addGeom(g);
				pol.push_back(i);
			} else if ((wkb == wkbMultiPolygon) || (wkb == wkbMultiPolygonZM) || (wkb == wkbMultiPolygonM)) {
				g = getMultiPolygonsGeom(poGeometry);
				polygons.addGeom(g);
				pol.push_back(i);
			} else {
				// g = emptyGeom();
			}
			OGRFeature::DestroyFeature( poFeature );
			i++;
		}
	}
	if (!query.empty()) {
		poDS->ReleaseResultSet(poLayer);
	}

	if (polygons.size() > 0) {
		polygons.setSRS(crs);
		polygons.read_query = query;
		polygons.read_extent= extent;
		polygons.source_layer = source_layer;
		if (df.ncol() > 0) polygons.df = df.subset_rows(pol);
		v.push_back(polygons);
	}
	if (lines.size() > 0) {
		lines.setSRS(crs);
		lines.read_query = query;
		lines.read_extent= extent;
		lines.source_layer = source_layer;
		if (df.ncol() > 0) lines.df = df.subset_rows(lin);
		v.push_back(lines);
	}
	if (points.size() > 0) {
		points.setSRS(crs);
		points.read_query = query;
		points.read_extent= extent;
		points.source_layer = source_layer;
		if (df.ncol() > 0) points.df = df.subset_rows(pnt);
		v.push_back(points);
	}
	return true;
}


bool SpatVectorCollection::read(std::string fname, std::string layer, std::string query, std::vector<double> extent, SpatVector filter) {
    //OGRRegisterAll();
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL ));
    if( poDS == NULL ) {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else {
			setError("Cannot open this file as a SpatVector: " + fname);
		}
		return false;
    }
	bool success = read_ogr(poDS, layer, query, extent, filter);
	if (poDS != NULL) GDALClose( poDS );
	return success;
}


SpatVectorCollection::SpatVectorCollection(std::string filename, std::string layer, std::string query, std::vector<double> extent, SpatVector filter) {
	read(filename, layer, query, extent, filter);
}
