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
    while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
        for(int iField = 0; iField < poFDefn->GetFieldCount(); iField++ ) {
            OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
            switch( poFieldDefn->GetType() ) {
                case OFTInteger:
                    printf( "%d,", poFeature->GetFieldAsInteger( iField ) );
                    break;
                case OFTInteger64:
                    printf( CPL_FRMT_GIB ",", poFeature->GetFieldAsInteger64( iField ) );
                    break;
                case OFTReal:
                    printf( "%.3f,", poFeature->GetFieldAsDouble(iField) );
                    break;
                case OFTString:
                    printf( "%s,", poFeature->GetFieldAsString(iField) );
                    break;
                default:
                    printf( "%s,", poFeature->GetFieldAsString(iField) );
                    break;
            }
        }

        OGRGeometry *poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ) {
            OGRPoint *poPoint = (OGRPoint *) poGeometry;
            printf( "%.3f, %.3f\n", poPoint->getX(), poPoint->getY() );
        } else {
            printf( "no point geometry\n" );
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    GDALClose( poDS );
	return true;
}

