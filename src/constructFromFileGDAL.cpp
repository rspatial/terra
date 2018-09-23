using namespace std;
#include <string>
#include <vector>
#include "spat.h"


#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ogr_spatialref.h"

bool SpatRaster::constructFromFileGDAL(std::string fname) {


    GDALDataset  *poDataset;
    GDALAllRegister();
	const char* pszFilename =  fname.c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	
    if( poDataset == NULL )  {
		return false;
	}	
		
	ncol = poDataset->GetRasterXSize();
	nrow = poDataset->GetRasterYSize();
	unsigned nlyrs = poDataset->GetRasterCount();
	source.nlayers = { nlyrs };	

	double adfGeoTransform[6];
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None ) {
		double xmin = adfGeoTransform[0]; /* top left x */
		double xmax = xmin + adfGeoTransform[1] * ncol; /* w-e pixel resolution */
		double ymax = adfGeoTransform[3]; /* top left y */
		double ymin = ymax + nrow * adfGeoTransform[5]; /* n-s pixel resolution (negative value) */
		SpatExtent e(xmin, xmax, ymin, ymax);
		setExtent(e, false);
	}
		
	source.memory.push_back(false);
	source.filename.push_back( fname );
			
	source.driver = {"gdal"};
	

	
	string crs;
	if( poDataset->GetProjectionRef()  != NULL ) {
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		char *pszPRJ = NULL;
		oSRS.exportToProj4(&pszPRJ);
		crs = pszPRJ;
	} else {
		crs = "";
	}
	setCRS(crs);
	setnlyr();	
	hasValues = true; 
	
	return true;
 
}


	
