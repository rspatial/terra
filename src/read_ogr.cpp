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
	OGRFieldDefn *poFieldDefn;
	df.resize_cols(nfields);
	unsigned dtype;
	long longNA = NA<long>::value;
	SpatTime_t timeNA = NA<SpatTime_t>::value;

	for (size_t i = 0; i < nfields; i++ ) {
		poFieldDefn = poFDefn->GetFieldDefn(i);
		std::string fname = poFieldDefn->GetNameRef();
		ft = poFieldDefn->GetType();
		// OFTInteger64 may be too large 
		if ((ft == OFTReal) || (ft == OFTInteger64)) {
			dtype = 0;
		} else if (ft == OFTInteger) {
			if (poFieldDefn->GetSubType() == OFSTBoolean) {
				dtype = 3;
			} else {
				dtype = 1;
			}
		} else if ((ft == OFTDate) || (ft == OFTDateTime)) {
			dtype = 4;
		} else {
			dtype = 2;
		}
		df.add_column(dtype, fname);
	}

    OGRFeature *poFeature;
    poLayer->ResetReading();
    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {

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
						df.dv[j].push_back(poFeature->GetFieldAsInteger64(i));
					} else {
						df.dv[j].push_back(NAN);
					}
					break;
				case OFTDate:  
					if (i == 0) {
						df.tv[j].step = "days";
					}
					if (not_null) {
						int pnYear, pnMonth, pnDay, pnHour, pnMinute, pnTZFlag;
						float pfSecond;
						poFeature->GetFieldAsDateTime(i, &pnYear, &pnMonth, &pnDay, &pnHour, &pnMinute, &pfSecond, &pnTZFlag);
						SpatTime_t d = get_time(pnYear, pnMonth, pnDay, 0, 0, 0);
						df.tv[j].x.push_back(d);
					} else {
						df.tv[j].x.push_back(timeNA);
					}
					break;
				case OFTDateTime:
					if (i == 0) {
						df.tv[j].step = "seconds";
					}
					if (not_null) {
						int pnYear, pnMonth, pnDay, pnHour, pnMinute, pnTZFlag;
						float pfSecond;
						poFeature->GetFieldAsDateTime(i, &pnYear, &pnMonth, &pnDay, &pnHour, &pnMinute, &pfSecond, &pnTZFlag);
						SpatTime_t d = get_time(pnYear, pnMonth, pnDay, pnHour, pnMinute, (int)pfSecond);
						df.tv[j].x.push_back(d);
					} else {
						df.tv[j].x.push_back(timeNA);
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

SpatGeom emptyGeom() {
	SpatGeom g;
	g.gtype = null;
	g.extent.xmin=NAN;
	g.extent.xmax=NAN;
	g.extent.ymin=NAN;
	g.extent.ymax=NAN;
	return g;
}


SpatGeom getPointGeom(OGRGeometry *poGeometry) {
	SpatGeom g(points);
	if (poGeometry->IsEmpty()) {
		//SpatPart p(NAN, NAN);
		//g.addPart(p);
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
	SpatGeom g(points);
	if (poGeometry->IsEmpty()) {
		return g;
	}
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
		double x = poPoint->getX();
		double y = poPoint->getY();
		SpatPart p(x, y);
		g.addPart(p);
	}
	return g;
}


SpatGeom getLinesGeom(OGRGeometry *poGeometry) {
	SpatGeom g(lines);
	if (poGeometry->IsEmpty()) {
		return g;
	}
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
	g.addPart(p);
	return g;
}

SpatGeom getMultiLinesGeom(OGRGeometry *poGeometry) {
	SpatGeom g(lines);
	if (poGeometry->IsEmpty()) {
		return g;
	}
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


SpatGeom getPolygonsGeom(OGRGeometry *poGeometry) {
	SpatGeom g(polygons);
	if (poGeometry->IsEmpty()) {
		return g;
	}
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
	SpatGeom g(polygons);
	if (poGeometry->IsEmpty()) {
		return g;
	}
	OGRMultiPolygon *poGeom = ( OGRMultiPolygon * )poGeometry;
	OGRPoint ogrPt;
	unsigned ng = poGeom->getNumGeometries();
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


static OGRGeometry* linearize_geom(OGRGeometry *poGeometry) {
	if (poGeometry == nullptr || poGeometry->IsEmpty()) return nullptr;

	OGRwkbGeometryType gt = wkbFlatten(poGeometry->getGeometryType());

	if (gt == wkbPoint || gt == wkbMultiPoint ||
		gt == wkbLineString || gt == wkbMultiLineString ||
		gt == wkbPolygon || gt == wkbMultiPolygon) {
		return poGeometry->clone();
	}

	OGRGeometry *lin = poGeometry->getLinearGeometry();
	if (lin != nullptr) {
		OGRwkbGeometryType lt = wkbFlatten(lin->getGeometryType());
		if (lt == wkbPoint || lt == wkbMultiPoint ||
			lt == wkbLineString || lt == wkbMultiLineString ||
			lt == wkbPolygon || lt == wkbMultiPolygon) {
			return lin;
		}
		OGRGeometryFactory::destroyGeometry(lin);
	}

	if (gt == wkbGeometryCollection) {
		OGRGeometry *forced;
		forced = OGRGeometryFactory::forceTo(poGeometry->clone(), wkbMultiPolygon);
		if (forced && !forced->IsEmpty()) return forced;
		if (forced) OGRGeometryFactory::destroyGeometry(forced);

		forced = OGRGeometryFactory::forceTo(poGeometry->clone(), wkbMultiLineString);
		if (forced && !forced->IsEmpty()) return forced;
		if (forced) OGRGeometryFactory::destroyGeometry(forced);

		forced = OGRGeometryFactory::forceTo(poGeometry->clone(), wkbMultiPoint);
		if (forced && !forced->IsEmpty()) return forced;
		if (forced) OGRGeometryFactory::destroyGeometry(forced);
	}

	return nullptr;
}


static SpatGeom linearize_to_spatgeom(OGRGeometry *poGeometry, bool &ok) {
	ok = false;
	SpatGeom g;
	if (poGeometry == nullptr || poGeometry->IsEmpty()) {
		ok = true;
		return g;
	}
	OGRGeometry *lin = linearize_geom(poGeometry);
	if (lin == nullptr) return g;

	OGRwkbGeometryType gt = wkbFlatten(lin->getGeometryType());
	if (gt == wkbPoint) {
		g = getPointGeom(lin);
		ok = true;
	} else if (gt == wkbMultiPoint) {
		g = getMultiPointGeom(lin);
		ok = true;
	} else if (gt == wkbLineString) {
		g = getLinesGeom(lin);
		ok = true;
	} else if (gt == wkbMultiLineString) {
		g = getMultiLinesGeom(lin);
		ok = true;
	} else if (gt == wkbPolygon) {
		g = getPolygonsGeom(lin);
		ok = true;
	} else if (gt == wkbMultiPolygon) {
		g = getMultiPolygonsGeom(lin);
		ok = true;
	}
	OGRGeometryFactory::destroyGeometry(lin);
	return g;
}


void addOGRgeometry(SpatVector &x, OGRGeometry *poGeometry) {
	SpatGeom g;
		
	OGRwkbGeometryType gtype = wkbFlatten(poGeometry->getGeometryType());
	if (gtype == wkbPoint) {
		g = getPointGeom(poGeometry);
	} else if (gtype == wkbMultiPoint) {
		g = getMultiPointGeom(poGeometry);
	} else if (gtype == wkbLineString) {
		g = getLinesGeom(poGeometry);
	} else if (gtype == wkbMultiLineString) {
		g = getMultiLinesGeom(poGeometry);
	} else if (gtype == wkbPolygon) {
		g = getPolygonsGeom(poGeometry);
	} else if (gtype == wkbMultiPolygon) {
		g = getMultiPolygonsGeom(poGeometry);
	} else {
		const char *geomtypechar = OGRGeometryTypeToName(gtype);
		std::string strgeomtype = geomtypechar;
		x.setError("cannot read geometry type: " + strgeomtype);
		return;
	}
	if ((x.size() > 1) && (x.geoms[0].gtype != g.gtype)) {
		if (x.geoms[0].gtype != null) {
			x.setError("a SpatVector can only have a single geometry type");
			return;
		}
	}
	x.addGeom(g);
	OGRGeometryFactory::destroyGeometry(poGeometry);
}


bool SpatVector::addRawGeoms(std::vector<unsigned char*> wkbs, std::vector<size_t> sizes) {
	
	if (wkbs.size() == 0) {
		SpatGeom g = emptyGeom();
		addGeom(g);
		return true;
	}
	
	for (size_t i=0; i<wkbs.size(); i++) {
		OGRGeometry *poGeometry;
		OGRErr err = OGRGeometryFactory::createFromWkb(wkbs[i], NULL, &poGeometry, sizes[i]);
		if (err == OGRERR_NONE) {
			if (poGeometry != NULL) {				
				addOGRgeometry(*this, poGeometry);
			}
		} else {
			setError("not valid WKB");
			return false;
		}
	}
	
	return true;
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


bool layerQueryFilter(GDALDataset *&poDS, OGRLayer *&poLayer, std::string &layer, std::string &query, std::string &dialect, std::vector<double> &ext, SpatVector &filter, std::string &errmsg, std::vector<std::string> &wrms) {

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
		poLayer = poDS->ExecuteSQL(query.c_str(), NULL, dialect.c_str());
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
	} else if (!ext.empty()) {
		poLayer->SetSpatialFilterRect(ext[0], ext[2], ext[1], ext[3]);
	}

	return true;
}



bool SpatVector::read_ogr(GDALDataset *&poDS, std::string layer, std::string query, std::vector<double> ext, SpatVector filter, bool as_proxy, std::string what, std::string dialect) {

	if (poDS == NULL) {
		setError("dataset is empty");
		return false;
	}

	std::string crs = "";
	OGRLayer *poLayer;
	poLayer = poDS->GetLayer(0);

	read_query = query;
	read_extent = ext;
	std::string errmsg;
	std::vector<std::string> wrnmsg;

	if (!layerQueryFilter(poDS, poLayer, layer, query, dialect, ext, filter, errmsg, wrnmsg)) {
		setError(errmsg);
		return false;
	} else if (!wrnmsg.empty()) {
		for (size_t i=0; i < wrnmsg.size(); i++) addWarning(wrnmsg[i]);
	}

	const OGRSpatialReference *poSRS = poLayer->GetSpatialRef();
	if (poSRS) {
		char *psz = NULL;
	#if GDAL_VERSION_MAJOR >= 3
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = poSRS->exportToWkt(&psz, options);
	#else
		OGRErr err = poSRS->exportToWkt(&psz);
	#endif
		if (err == OGRERR_NONE && psz != NULL) {
			crs = psz;
		}
		CPLFree(psz);
		if (!crs.empty()) {
			setSRS(crs);
		}
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
		if (poFeature == NULL) {
			g = emptyGeom();
			addGeom(g);
		} else {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			bool ok = false;
			g = linearize_to_spatgeom(poGeometry, ok);
			if (!ok) {
				const char *gtn = poGeometry ?
					OGRGeometryTypeToName(wkbFlatten(poGeometry->getGeometryType())) : "NULL";
				setError("cannot read this geometry type: " + std::string(gtn));
				OGRFeature::DestroyFeature(poFeature);
				return false;
			}
			addGeom(g);
			OGRFeature::DestroyFeature(poFeature);
		}
		geom_count = poLayer->GetFeatureCount();

		OGREnvelope oExt;
		if (poLayer->GetExtent(&oExt, FALSE) == OGRERR_NONE) {
			extent.xmin = oExt.MinX;
			extent.xmax = oExt.MaxX;
			extent.ymin = oExt.MinY;
			extent.ymax = oExt.MaxY;
		} else {
			extent.xmin = NAN;
			extent.xmax = NAN;
			extent.ymin = NAN;
			extent.ymax = NAN;
		}

		is_proxy = true;
		return true;
	}

	OGRFeature::DestroyFeature(poFeature);
	poLayer->ResetReading();

	SpatVector points, lines, polygons;
	std::vector<size_t> pnt, lin, pol;
	size_t skip_count = 0;
	size_t i = 0;

	while ((poFeature = poLayer->GetNextFeature()) != NULL) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		if (poGeometry == NULL || poGeometry->IsEmpty()) {
			SpatGeom g = emptyGeom();
			polygons.addGeom(g);
			pol.push_back(i);
		} else {
			OGRGeometry *lin_geom = linearize_geom(poGeometry);
			if (lin_geom == NULL) {
				skip_count++;
			} else {
				OGRwkbGeometryType gt = wkbFlatten(lin_geom->getGeometryType());
				SpatGeom g;
				if (gt == wkbPoint) {
					g = getPointGeom(lin_geom);
					points.addGeom(g);
					pnt.push_back(i);
				} else if (gt == wkbMultiPoint) {
					g = getMultiPointGeom(lin_geom);
					points.addGeom(g);
					pnt.push_back(i);
				} else if (gt == wkbLineString) {
					g = getLinesGeom(lin_geom);
					lines.addGeom(g);
					lin.push_back(i);
				} else if (gt == wkbMultiLineString) {
					g = getMultiLinesGeom(lin_geom);
					lines.addGeom(g);
					lin.push_back(i);
				} else if (gt == wkbPolygon) {
					g = getPolygonsGeom(lin_geom);
					polygons.addGeom(g);
					pol.push_back(i);
				} else if (gt == wkbMultiPolygon) {
					g = getMultiPolygonsGeom(lin_geom);
					polygons.addGeom(g);
					pol.push_back(i);
				} else {
					skip_count++;
				}
				OGRGeometryFactory::destroyGeometry(lin_geom);
			}
		}
		OGRFeature::DestroyFeature(poFeature);
		i++;
	}

	if (skip_count > 0) {
		addWarning(std::to_string(skip_count) + " geometries could not be converted and were skipped");
	}

	size_t npol = polygons.size();
	size_t nlin = lines.size();
	size_t npnt = points.size();

	if (npol == 0 && nlin == 0 && npnt == 0) {
		if (!query.empty()) {
			poDS->ReleaseResultSet(poLayer);
		}
		return true;
	}

	size_t ntypes = (npol > 0) + (nlin > 0) + (npnt > 0);

	if (ntypes == 1) {
		if (npol > 0) {
			geoms = polygons.geoms;
			extent = polygons.extent;
			if (df.ncol() > 0 && pol.size() != i) df = df.subset_rows(pol);
		} else if (nlin > 0) {
			geoms = lines.geoms;
			extent = lines.extent;
			if (df.ncol() > 0 && lin.size() != i) df = df.subset_rows(lin);
		} else {
			geoms = points.geoms;
			extent = points.extent;
			if (df.ncol() > 0 && pnt.size() != i) df = df.subset_rows(pnt);
		}
	} else {
		if (npol > 0) {
			geoms = polygons.geoms;
			extent = polygons.extent;
			if (df.ncol() > 0) df = df.subset_rows(pol);
			std::string msg = "returning polygons. Ignoring ";
			if (nlin > 0) msg += std::to_string(nlin) + " line";
			if (nlin > 0 && npnt > 0) msg += " and ";
			if (npnt > 0) msg += std::to_string(npnt) + " point";
			msg += " geometries. Use 'svc' to get all geometries";
			addWarning(msg);
		} else if (nlin > 0) {
			geoms = lines.geoms;
			extent = lines.extent;
			if (df.ncol() > 0) df = df.subset_rows(lin);
			addWarning("returning lines. Ignoring " + std::to_string(npnt) + " point geometries. Use 'svc' to get all geometries");
		}
	}

	if (!query.empty()) {
		poDS->ReleaseResultSet(poLayer);
	}

	return true;
}


bool SpatVector::read(std::string fname, std::string layer, std::string query, std::vector<double> ext, SpatVector filter, bool as_proxy, std::string what, std::string dialect, std::vector<std::string> options) {

	char ** openops = NULL;
	for (size_t i=0; i<options.size(); i++) {
		std::vector<std::string> opt = strsplit(options[i], "=");
		if (opt.size() == 2) {
			openops = CSLSetNameValue(openops, opt[0].c_str(), opt[1].c_str());
		}
	}
		
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_VECTOR, NULL, openops, NULL ));
    if( poDS == NULL ) {
		if (!file_exists(fname)) {
			setError("file does not exist: " + fname);
		} else {
			setError("Cannot open this file as a SpatVector: " + fname);
		}
		return false;
    }
	bool success = read_ogr(poDS, layer, query, ext, filter, as_proxy, what, dialect);
	if (poDS != NULL) GDALClose( poDS );
	if (fname.substr(0, 1) != "{") { // not json
		source = fname;
	}
	return success;
}

SpatVector SpatVector::fromDS(GDALDataset *poDS) {
	SpatVector out, fvct;
	std::vector<double> fext;
	out.read_ogr(poDS, "", "", fext, fvct, false, "", "");
	return out;
}


SpatVector::SpatVector(std::vector<std::string> wkt) {

	OGRGeometryFactory ogr;

	extent.xmin = extent.xmax = extent.ymin = extent.ymax = NAN;
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

bool SpatVectorCollection::read_ogr(GDALDataset *&poDS, std::string layer, std::string query, std::string dialect, std::vector<double> extent, SpatVector filter) {
	
	OGRLayer *poLayer;
	poLayer = poDS->GetLayer(0);
	
	std::string errmsg;
	std::vector<std::string> wrnmsg;

	if (!layerQueryFilter(poDS, poLayer, layer, query, dialect, extent, filter, errmsg, wrnmsg)) {
		setError(errmsg);
		return false;
	} else if (!wrnmsg.empty()) {
		for (size_t i=0; i < wrnmsg.size(); i++) addWarning(wrnmsg[i]);
	}
	std::string crs = "";
	const OGRSpatialReference *poSRS = poLayer->GetSpatialRef();
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
	std::vector<size_t> pnt, lin, pol;
	SpatGeom g;
	size_t i = 0;
	size_t skip_count = 0;
	while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		if (poGeometry != NULL && !poGeometry->IsEmpty()) {
			OGRGeometry *lin_geom = linearize_geom(poGeometry);
			if (lin_geom == NULL) {
				skip_count++;
			} else {
				OGRwkbGeometryType wkb = wkbFlatten(lin_geom->getGeometryType());
				if (wkb == wkbPoint) {
					g = getPointGeom(lin_geom);
					points.addGeom(g);
					pnt.push_back(i);
				} else if (wkb == wkbMultiPoint) {
					g = getMultiPointGeom(lin_geom);
					points.addGeom(g);
					pnt.push_back(i);
				} else if (wkb == wkbLineString) {
					g = getLinesGeom(lin_geom);
					lines.addGeom(g);
					lin.push_back(i);
				} else if (wkb == wkbMultiLineString) {
					g = getMultiLinesGeom(lin_geom);
					lines.addGeom(g);
					lin.push_back(i);
				} else if (wkb == wkbPolygon) {
					g = getPolygonsGeom(lin_geom);
					polygons.addGeom(g);
					pol.push_back(i);
				} else if (wkb == wkbMultiPolygon) {
					g = getMultiPolygonsGeom(lin_geom);
					polygons.addGeom(g);
					pol.push_back(i);
				} else {
					skip_count++;
				}
				OGRGeometryFactory::destroyGeometry(lin_geom);
			}
		} else {
			g = emptyGeom();
			polygons.addGeom(g);
			pol.push_back(i);
		}
		OGRFeature::DestroyFeature( poFeature );
		i++;
	}
	if (skip_count > 0) {
		addWarning(std::to_string(skip_count) + " geometries could not be converted and were skipped");
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


bool SpatVectorCollection::read(std::string fname, std::string layer, std::string query, std::string dialect, std::vector<double> extent, SpatVector filter) {
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
	bool success = read_ogr(poDS, layer, query, dialect, extent, filter);
	if (poDS != NULL) GDALClose( poDS );
	return success;
}


SpatVectorCollection::SpatVectorCollection(std::string filename, std::string layer, std::string query, std::string dialect, std::vector<double> extent, SpatVector filter) {
	read(filename, layer, query, dialect, extent, filter);
}
