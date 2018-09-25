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
	if( poDataset->GetProjectionRef() != NULL ) {
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
	
	source.layers.resize(1);
	for (size_t i = 0; i < nlyrs; i++) {
		poBand = poDataset->GetRasterBand(i+1);
		source.layers[0].push_back(i+1);
		//poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
		//GDALGetColorInterpretationName( poBand->GetColorInterpretation()) );
		
		std::string dtype = GDALGetDataTypeName(poBand->GetRasterDataType());
		
		adfMinMax[0] = poBand->GetMinimum( &bGotMin );
		adfMinMax[1] = poBand->GetMaximum( &bGotMax );
		if( (bGotMin && bGotMax) ) {
			hasRange.push_back(true);
			range_min.push_back( adfMinMax[0] );
			range_max.push_back( adfMinMax[1] );
		} else {
			hasRange.push_back(false);
			range_min.push_back( NAN );
			range_max.push_back( NAN );		
			// GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
		}
			
		//	
		//if( poBand->GetOverviewCount() > 0 ) printf( "Band has %d overviews.\n", poBand->GetOverviewCount() );
		//if( poBand->GetColorTable() != NULL )	printf( "Band has a color table with %d entries.\n", 
		//         poBand->GetColorTable()->GetColorEntryCount() );	
		names.push_back( "lyr" + to_string(i+1) ) ;
	}

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
	std::vector<double> out;
	unsigned ncell = ncols*nrows;
	double *pafScanline;
	pafScanline = (double *) CPLMalloc(sizeof(double)*ncell);
	
	for (size_t i=0; i < source.layers[0].size(); i++) {
	
		poBand = poDataset->GetRasterBand(i+1);
		CPLErr err = poBand->RasterIO( GF_Read, row, col, ncols, nrows, pafScanline, ncols, nrows, GDT_Float64, 0, 0 );	
		if (err == 4) {
			std::vector<double> errout;
			return  errout;
		}
		out.insert(out.end(), &pafScanline[0], &pafScanline[ncell]);
	}
	
	CPLFree(pafScanline);
	GDALClose( (GDALDatasetH) poDataset );
	return(out);
}

