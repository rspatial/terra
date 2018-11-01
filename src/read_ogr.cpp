#include "spatvector.h"
#include "util.h"

#include "ogrsf_frmts.h"
#include "ogr_spatialref.h"
using namespace std;


string geomType(OGRLayer *poLayer) {
	string s = "";
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



void readPoints (OGRLayer *poLayer, std::vector<double> &x, std::vector<double> &y) {
	OGRFeature *poFeature;
    poLayer->ResetReading();
	size_t i=0;
    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
        OGRGeometry *poGeometry = poFeature->GetGeometryRef();     
        if( poGeometry != NULL) { // && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ) {
			OGRPoint *poPoint = (OGRPoint *) poGeometry;
            x[i] = poPoint->getX();
			y[i] = poPoint->getY();
        } else {
            x[i] = NAN;
			y[i] = NAN;
        }
		i++;
    }
    OGRFeature::DestroyFeature( poFeature );
}



bool SpatVector::read(std::string fname) {

    GDALAllRegister();
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL ));
    if( poDS == NULL ) {
        error_message = "Cannot open file";
		error = true;
		return false;
    }
	string crs = "";
	OGRSpatialReference *poSRS = poDS->GetLayer(0)->GetSpatialRef();
	if (poSRS) {
		char *pszPRJ = NULL;
		poSRS->exportToProj4(&pszPRJ);
		crs = pszPRJ;
	}
	
	OGRLayer  *poLayer = poDS->GetLayerByName( basename(fname).c_str() );
	df = readAttributes(poLayer);

	string geomtype = geomType(poLayer);
	if (geomtype == "POINT") {
		gtype = 1;
		std::vector<double> X(df.nrow());
		std::vector<double> Y(df.nrow());
		readPoints (poLayer, X, Y);
		pts.set(X, Y);
	} else {
		//printf("unknown geomtype");		
	}
	
    GDALClose( poDS );
	setCRS(crs);
	return true;
}

