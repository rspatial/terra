#include "spatraster.h"
#include "util.h"
#include "ogrsf_frmts.h"


bool SpatLayer::write(std::string filename, bool overwrite) {

	msg.success = true;

    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;
    GDALAllRegister();
    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ) {
        setError("driver not available");
		return false; 
    }
    GDALDataset *poDS;
    poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        setError("cannot write file");
		return false; 
    }
    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer(basename(filename).c_str(), NULL, wkbPoint, NULL );
    if( poLayer == NULL ) {
        setError("Layer creation failed");
		return false; 
    }
    OGRFieldDefn oField( "Name", OFTString );
    oField.SetWidth(32);
    if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
        setError("Field creation failed");
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
			setError("Failed to create feature in shapefile");
			return false; 			
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    GDALClose( poDS );
	return true;
}


