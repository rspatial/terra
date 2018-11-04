#include "spatVector.h"
#include "util.h"

#include "ogrsf_frmts.h"
#include "ogr_spatialref.h"

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
	OGRFieldType ft;
    poLayer->ResetReading();
    OGRFeature *poFeature;
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
	unsigned nfields = poFDefn->GetFieldCount();
	OGRFieldDefn *poFieldDefn;
	SpatDataFrame df;
	df.itype.resize(nfields);
	df.iplace.resize(nfields);
	unsigned dcnt = 0;
	unsigned icnt = 0; 
	unsigned scnt = 0;
	bool first = true;
    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		if (first) {
			for (size_t i = 0; i < nfields; i++ ) {
				poFieldDefn = poFDefn->GetFieldDefn(i);
				df.names.push_back(poFieldDefn->GetNameRef());
				ft = poFieldDefn->GetType();
				if (ft == OFTReal) {
					df.itype[i] = 0;
					df.iplace[i] = dcnt;
					dcnt++;
				} else if ((ft == OFTInteger) | (ft == OFTInteger64)) {
					df.itype[i] = 1;
					df.iplace[i] = icnt;
					icnt++;
				} else {
					df.itype[i] = 2;
					df.iplace[i] = scnt;
					scnt++;
				}
			}
			df.dv.resize(dcnt);
			df.iv.resize(icnt);
			df.sv.resize(scnt);	
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
	}
    OGRFeature::DestroyFeature( poFeature );
	return df;
}	




bool SpatLayer::read(std::string fname) {

    GDALAllRegister();
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL ));
    if( poDS == NULL ) {
        error_message = "Cannot open file";
		error = true;
		return false;
    }
	std::string crs = "";
	OGRSpatialReference *poSRS = poDS->GetLayer(0)->GetSpatialRef();
	if (poSRS) {
		char *pszPRJ = NULL;
		poSRS->exportToProj4(&pszPRJ);
		crs = pszPRJ;
	}
	
	OGRLayer *poLayer = poDS->GetLayerByName( basename(fname).c_str() );
	df = readAttributes(poLayer);

//	std::string geomtype = geomType(poLayer);
	OGRwkbGeometryType wkbgeom = wkbFlatten( poLayer ->GetGeomType());
	OGRFeature *poFeature;
	OGRPoint ogrPt;
	poLayer->ResetReading();
	unsigned np, nh, ng;
	
	if (wkbgeom == wkbPoint) {
		gtype = POINTS;
		SpatPart p(0,0);
		while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();     
			if( poGeometry != NULL) { // && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ) {
			#if GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(2,3,0)
				OGRPoint *poPoint = poGeometry->toPoint();
			#else
				OGRPoint *poPoint = (OGRPoint *) poGeometry;
			#endif
				p.x[0] = poPoint->getX();
				p.y[0] = poPoint->getY();
			} else {
				p.x[0] = NAN;
				p.y[0] = NAN;
			}
			//SpatPart p(x, y);
			SpatGeom g(p);
			addGeom(g);
		}

	} else if (wkbgeom == wkbLineString) {
		gtype = LINES;
		while ( (poFeature = poLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();     
			SpatGeom g;
			if (poGeometry != NULL) { 
				OGRLineString *poGeom = (OGRLineString *) poGeometry;
		        np = poGeom->getNumPoints();
				std::vector<double> X(np);
				std::vector<double> Y(np);
				for (size_t i=0; i<np; i++) {
					poGeom->getPoint(i, &ogrPt);
					X[i] = ogrPt.getX();
					Y[i] = ogrPt.getY();
				}
				SpatPart p(X, Y);
				g.addPart(p);
			}
			addGeom(g);
		}
		OGRFeature::DestroyFeature( poFeature );

	} else if (wkbgeom == wkbPolygon) {
		gtype = POLYGONS;
		while ( (poFeature = poLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = poFeature ->GetGeometryRef();
			SpatGeom g;
			if (poGeometry != NULL) { // && wkbFlatten ( poGeometry ->getGeometryType() ) == wkbPolygon )
				OGRPolygon *poGeom = ( OGRPolygon * )poGeometry;
				OGRLinearRing *poRing = poGeom->getExteriorRing();
				np = poRing->getNumPoints();				
				std::vector<double> X(np);
				std::vector<double> Y(np);
				for (size_t i=0; i<np; i++) {
					poRing->getPoint(i, &ogrPt);
					X[i] = ogrPt.getX();
					Y[i] = ogrPt.getY();
				}
				SpatPart p(X, Y);
				
				nh = poGeom->getNumInteriorRings();
				for (size_t i=0; i<nh; i++) {
					OGRLinearRing *poHole = poGeom->getInteriorRing(i);
					np = poHole->getNumPoints();
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
			}
			addGeom(g);
		}
	} else if (wkbgeom == wkbMultiPolygon) {
		gtype = POLYGONS;
		while ( (poFeature = poLayer->GetNextFeature()) != NULL ) {
			OGRGeometry *poGeometry = poFeature->GetGeometryRef();
			SpatGeom g;
			if ( poGeometry != NULL ) { // && wkbFlatten ( poGeometry ->getGeometryType() ) == wkbPolygon )
				OGRMultiPolygon *poGeom = ( OGRMultiPolygon * )poGeometry;
				ng = poGeom->getNumGeometries();
				for (size_t i=0; i<ng; i++) {
					OGRGeometry *poPolygonGeometry = poGeom->getGeometryRef(i);
					OGRPolygon *poPolygon = ( OGRPolygon * )poPolygonGeometry;
					OGRLinearRing *poRing = poPolygon->getExteriorRing();
					np = poRing->getNumPoints();
					std::vector<double> X(np);
					std::vector<double> Y(np);
					for (size_t j=0; j<np; j++ ) {
						poRing->getPoint(j, &ogrPt);
						X[j] = ogrPt.getX();
						Y[j] = ogrPt.getY();
					}
					SpatPart p(X, Y);

					nh = poPolygon->getNumInteriorRings();
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
			}
			addGeom(g);
		}
	} else {
		std::string gt = geomType(poLayer);		
		printf("unknown geomtype: %s \n", gt.c_str());		
	}
	
	OGRFeature::DestroyFeature( poFeature );
    GDALClose( poDS );
	setCRS(crs);
	return true;
}


 // get the extent
 //       OGREnvelope oExt;
 //       if( poLayer->GetExtent( &oExt, TRUE ) == OGRERR_NONE ){
 //           cout << "Extent: (" << oExt.MinX << ", " << oExt.MinY << ") - (" << oExt.MaxX << ", " << oExt.MaxY << ")" << endl;
