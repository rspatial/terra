#include "spatraster.h"
#include "util.h"
#include "ogrsf_frmts.h"
using namespace std;


bool SpatLayer::write(std::string filename, bool overwrite) {

    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALAllRegister();
    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ) {
        error_message = "driver not available";
        error = true;
		return false; 
    }
    GDALDataset *poDS;
    poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        error_message = "cannot write file";
        error = true;
		return false; 
    }
    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer(basename(filename).c_str(), NULL, wkbPoint, NULL );
    if( poLayer == NULL ) {
        error_message = "Layer creation failed";
        error = true;
		return false; 
    }
    OGRFieldDefn oField( "Name", OFTString );
    oField.SetWidth(32);
    if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
        error_message = "Field creation failed";
        error = true;
		return false; 
    }
    double x, y;
    char szName[33];
    while( !feof(stdin) && fscanf( stdin, "%lf,%lf,%32s", &x, &y, szName ) == 3 )  {
        OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
        poFeature->SetField( "Name", szName );
        OGRPoint pt;
        pt.setX( x );
        pt.setY( y );
        poFeature->SetGeometry( &pt );
        if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ) {
			error_message = "Failed to create feature in shapefile";
			error = true;
			return false; 			
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    GDALClose( poDS );
	return true;
}


