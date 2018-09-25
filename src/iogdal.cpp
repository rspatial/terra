using namespace std;
#include "spat.h"
#include "util.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ogr_spatialref.h"





bool SpatRaster::constructFromFileGDAL(std::string fname) {

    GDALDataset  *poDataset;
    GDALAllRegister();
	const char* pszFilename = fname.c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	
    if( poDataset == NULL )  {	return false;}	
		
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
	setCRS(crs);
	setnlyr();	

	GDALRasterBand  *poBand;
	//int nBlockXSize, nBlockYSize;
	double adfMinMax[2];
	int bGotMin, bGotMax;
	
	// need to loop over bands here
	poBand = poDataset->GetRasterBand(1);
	//poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
	//GDALGetColorInterpretationName( poBand->GetColorInterpretation()) );
	std::string dtype = GDALGetDataTypeName(poBand->GetRasterDataType());
	
	adfMinMax[0] = poBand->GetMinimum( &bGotMin );
	adfMinMax[1] = poBand->GetMaximum( &bGotMax );
	if( (bGotMin && bGotMax) ) {
		hasRange.push_back(true);
		range_min.push_back( adfMinMax[0] );
		range_max.push_back( adfMinMax[1] );
	}
    //	else GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
	//if( poBand->GetOverviewCount() > 0 ) printf( "Band has %d overviews.\n", poBand->GetOverviewCount() );
	//if( poBand->GetColorTable() != NULL )	printf( "Band has a color table with %d entries.\n", 
    //         poBand->GetColorTable()->GetColorEntryCount() );	


	GDALClose( (GDALDatasetH) poDataset );

	hasValues = true; 
	return true;
 
}


std::vector<double> SpatRaster::readGDALvalues(unsigned row, unsigned nrows, unsigned col, unsigned ncols) {
    GDALDataset  *poDataset;
	GDALRasterBand  *poBand;
	
    GDALAllRegister();
	const char* pszFilename = source.filename[0].c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	poBand = poDataset->GetRasterBand(1);
	int nXSize = 10;
	int nYSize = 10;

	double *pafScanline;
	pafScanline = (double *) CPLMalloc(sizeof(double)*nXSize);
	poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, pafScanline, nXSize, nYSize, GDT_Float32, 0, 0 );
	
	std::vector<double> out;
	out.insert(out.end(), &pafScanline[0], &pafScanline[nXSize * nYSize]);

	GDALClose( (GDALDatasetH) poDataset );
	//CPLFree(poBand);
	return(out);
}

