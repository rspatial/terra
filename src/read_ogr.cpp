// Copyright (c) 2018-2021  Robert J. Hijmans
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


SpatDataFrame readAttributes(OGRLayer *poLayer) {
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
    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		if (first) {
			for (size_t i = 0; i < nfields; i++ ) {
				poFieldDefn = poFDefn->GetFieldDefn(i);
				std::string fname = poFieldDefn->GetNameRef();
				ft = poFieldDefn->GetType();
				if (ft == OFTReal) {
					dtype = 0;
				} else if ((ft == OFTInteger) | (ft == OFTInteger64)) {
					dtype = 1;
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
			switch( poFieldDefn->GetType() ) {
				case OFTReal:
					df.dv[j].push_back(poFeature->GetFieldAsDouble(i));
					break;
				case OFTInteger:
					df.iv[j].push_back(poFeature->GetFieldAsInteger( i ));
					break;
				case OFTInteger64:
					df.iv[j].push_back(poFeature->GetFieldAsInteger64( i ));
					break;
	//          case OFTString:
				default:
					df.sv[j].push_back(poFeature->GetFieldAsString( i ));
					break;
			}
		}
		OGRFeature::DestroyFeature( poFeature );
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


bool SpatVector::read_ogr(GDALDataset *poDS, std::string layer,  std::string query, std::vector<double> extent, SpatVector filter) {

	std::string crs = "";

	OGRSpatialReference *poSRS = poDS->GetLayer(0)->GetSpatialRef();
	if (poSRS) {
		char *psz = NULL;
		OGRErr err = poSRS->exportToWkt(&psz);
		if (err == OGRERR_NONE) {
			crs = psz;
		}
		setSRS(crs);
		CPLFree(psz);
	}

	OGRLayer *poLayer;

	if (query != "") {
		poLayer = poDS->ExecuteSQL(query.c_str(), NULL, NULL);
		if (poLayer == NULL) {
			setError("Query failed");
			return false;
		}
	} else {
		if (layer == "") {
			poLayer = poDS->GetLayer(0);
			if (poLayer == NULL) {
				setError("dataset has no layers");
				return false;
			}
		} else {
			poLayer = poDS->GetLayerByName(layer.c_str());
			if (poLayer == NULL) {
				std::string msg = layer + " is not a valid layer name";
			#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
				// do nothing
			#else
				msg += "\nChoose one of: ";
				for ( auto&& poLayer: poDS->GetLayers() ) {
					msg += (std::string)poLayer->GetName() + ", ";
				}
				msg = msg.substr(0, msg.size()-2);
			#endif
				setError(msg);
				return false;
			}
		}
	}

	if (filter.nrow() > 0) {
		if (filter.type() != "polygons") {
			filter = filter.hull("convex");
		} else {
			if (filter.nrow() > 1) {
				filter = filter.aggregate(true);
			}
		}
		GDALDataset *filterDS = filter.write_ogr("", "lyr", "Memory", true, std::vector<std::string>());
		if (filter.hasError()) {
			setError(filter.getError());
			GDALClose(filterDS);
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
	} else if (extent.size() > 0) {
		poLayer->SetSpatialFilterRect(extent[0], extent[2], extent[1], extent[3]);
	}

	//const char* lname = poLayer->GetName();

	df = readAttributes(poLayer);
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
				SpatPart p;
				g = SpatGeom();
				g.addPart(p);
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
				SpatPart p;
				g = SpatGeom();
				g.addPart(p);
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
				g = SpatGeom();
			}
			addGeom(g);
			OGRFeature::DestroyFeature( poFeature );
		}
	} else if (wkbgeom != wkbNone) {
		const char *geomtypechar = OGRGeometryTypeToName(wkbgeom);
		std::string strgeomtype = geomtypechar;
		std::string s = "cannot read this geometry type: "+ strgeomtype;
		setError(s);
		return false;
	}

	if (query != "") {
		poDS->ReleaseResultSet(poLayer);
	}

 	return true;
}


bool SpatVector::read(std::string fname, std::string layer, std::string query, std::vector<double> extent, SpatVector filter) {
    //OGRRegisterAll();
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL ));
    if( poDS == NULL ) {
        setError("Cannot open this file as a SpatVector");
		return false;
    }
	bool success = read_ogr(poDS, layer, query, extent, filter);
	if (poDS != NULL) GDALClose( poDS );
	return success;
}

SpatVector SpatVector::fromDS(GDALDataset *poDS) {
	SpatVector out, fvct;
	std::vector<double> fext;
	out.read_ogr(poDS, "", "", fext, fvct);
	return out;
}


SpatVector::SpatVector(std::vector<std::string> wkt) {

	OGRGeometryFactory ogr;

	SpatGeom g;
	for (size_t i=0; i<wkt.size(); i++) {

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
				addGeom(g);
				OGRGeometryFactory::destroyGeometry(poGeometry);

			}
		} else {
			setError("not WKT");
			return;
		}
	}
}

