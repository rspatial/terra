using namespace std;
#include "spat.h"
#include "util.h"

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
		// the rounding below is to address a design flaw in GDAL
		// GDAL provides the coordinates of one corner and the resolution, instead of the coordinates of all (two opposite) corners.
		// computation of the opposite corder coordinates is only approximate for large rasters with a high resolution.
		double xmin = adfGeoTransform[0]; /* top left x */
		xmin = roundn(xmin, 9);
		double xmax = xmin + adfGeoTransform[1] * ncol; /* w-e pixel resolution */
		xmax = roundn(xmax, 9);
		double ymax = adfGeoTransform[3]; /* top left y */
		ymax = roundn(ymax, 9);
		double ymin = ymax + nrow * adfGeoTransform[5]; /* n-s pixel resolution (negative value) */
		ymin = roundn(ymin, 9);
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
	
	GDALClose( (GDALDatasetH) poDataset );
	
	setCRS(crs);
	setnlyr();	
	// for now
	hasValues = false;
	//hasValues = true; 
	
	return true;
 
}


	
