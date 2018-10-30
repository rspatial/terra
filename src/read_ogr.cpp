using namespace std;
#include "spatvector.h"
#include "util.h"

#include "ogrsf_frmts.h"
#include "ogr_spatialref.h"


bool SpatVector::read(std::string fname) {

    GDALAllRegister();
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx( fname.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL ));
    if( poDS == NULL ) {
        error_message = "Cannot open file";
		error = true;
		return false;
    }
    OGRLayer  *poLayer = poDS->GetLayerByName( basename(fname).c_str() );
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    poLayer->ResetReading();
    OGRFeature *poFeature;
	unsigned nfields = poFDefn->GetFieldCount();
	df.itype.resize(nfields);
	df.iplace.resize(nfields);
	unsigned dcnt = 0;
	unsigned icnt = 0; 
	unsigned scnt = 0;
	OGRFieldType ft;
	std::vector<double> X, Y;
	bool first = true;

    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		if (first) {
			for (int i = 0; i < nfields; i++ ) {
				OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(i);
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
			first = false;
			
			df.dv.resize(dcnt);
			df.iv.resize(icnt);
			df.sv.resize(scnt);	
		} 
			
			
		for (int i = 0; i < nfields; i++ ) {
            OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( i );
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
//                case OFTString:
                default:
					df.sv[j].push_back(poFeature->GetFieldAsString( i ));
                    break;
            }
        }

        OGRGeometry *poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ) {
			OGRPoint *poPoint = (OGRPoint *) poGeometry;
            X.push_back( poPoint->getX() );
			Y.push_back( poPoint->getY() );
        } else {
            //printf( "no point geometry\n" );
        }
        OGRFeature::DestroyFeature( poFeature );
    }
	

    GDALClose( poDS );
	pts.set(X, Y);
	return true;
}

