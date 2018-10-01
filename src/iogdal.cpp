using namespace std;
#include "spat.h"
#include "util.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"


// void SpatRaster::setMinMaxGDAL() {
// GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
// }

bool SpatRaster::constructFromFileGDAL(std::string fname) {

    GDALDataset  *poDataset;
    GDALAllRegister();
	const char* pszFilename = fname.c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	
    if( poDataset == NULL )  {	return false;}	
		
	ncol = poDataset->GetRasterXSize();
	nrow = poDataset->GetRasterYSize();
	unsigned nlyrs = poDataset->GetRasterCount();

	RasterSource s;
	s.nlyr = nlyrs;	
	
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
		
	s.memory = false;
	s.filename = fname;	
	s.driver = "gdal";
	
	string crs;
	if( poDataset->GetProjectionRef() != NULL ) {
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		char *pszPRJ = NULL;
		oSRS.exportToProj4(&pszPRJ);
		crs = pszPRJ;
	} else {
		crs = "";
	}
	s.crs = crs;
	setnlyr();	

	GDALRasterBand  *poBand;
	//int nBlockXSize, nBlockYSize;
	double adfMinMax[2];
	int bGotMin, bGotMax;
	
//	s.layers.resize(1);
	for (size_t i = 0; i < nlyrs; i++) {
		poBand = poDataset->GetRasterBand(i+1);
//		source.layers[0].push_back(i+1);
		//poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
		
		std::string dtype = GDALGetDataTypeName(poBand->GetRasterDataType());
		
		adfMinMax[0] = poBand->GetMinimum( &bGotMin );
		adfMinMax[1] = poBand->GetMaximum( &bGotMax );
		if( (bGotMin && bGotMax) ) {
			s.hasRange.push_back(true);
			s.range_min.push_back( adfMinMax[0] );
			s.range_max.push_back( adfMinMax[1] );
		} else {
			s.hasRange.push_back(false);
			s.range_min.push_back( NAN );
			s.range_max.push_back( NAN );		
		}
			
		//if( poBand->GetOverviewCount() > 0 ) printf( "Band has %d overviews.\n", poBand->GetOverviewCount() );
		
		//GDALGetColorInterpretationName( poBand->GetColorInterpretation()) );
		
		GDALColorTable *ct;
		ct = poBand->GetColorTable(); 
		if( ct != NULL )	{
			s.hasCT.push_back(true);
		} else {
			s.hasCT.push_back(false);
		}
		
		GDALRasterAttributeTable *rat = poBand->GetDefaultRAT(); 	
		if( rat != NULL )	{  // does not appear to work
			s.hasRAT.push_back(true);
		} else {
			s.hasRAT.push_back(false);
		}
		
		s.names.push_back( "lyr" + to_string(i+1) ) ;
	}

	GDALClose( (GDALDatasetH) poDataset );

	s.hasValues = true;
	source = { s };
	return true;
 
}

std::vector<double> SpatRaster::readValuesGDAL(unsigned row, unsigned nrows, unsigned col, unsigned ncols) {
    GDALDataset  *poDataset;
	GDALRasterBand  *poBand;
    GDALAllRegister();
	
	const char* pszFilename = source[0].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	std::vector<double> out;
	unsigned ncell = ncols*nrows;
	double *pafScanline;
	pafScanline = (double *) CPLMalloc(sizeof(double)*ncell);
	
	for (size_t i=0; i < source[0].nlyr; i++) {
	
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



bool SpatRaster::writeValuesGDAL(std::string filename, std::vector<double> values, std::string format) {

    const char *pszFormat = format.c_str();
	const char *pszDstFilename = filename.c_str();
    GDALDriver *poDriver;
    char **papszMetadata;
	GDALAllRegister();
	
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == NULL ) return (false);
    papszMetadata = poDriver->GetMetadata();
    if(! CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) )  return (false);
 
	GDALDataset *poDstDS;
	char **papszOptions = NULL;
	poDstDS = poDriver->Create( pszDstFilename, 512, 512, 1, GDT_Byte, papszOptions );

	double adfGeoTransform[6] = { 444720, 30, 0, 3751320, 0, -30 };
	OGRSpatialReference oSRS;
	char *pszSRS_WKT = NULL;
	GDALRasterBand *poBand;
	GByte abyRaster[512*512];
	poDstDS->SetGeoTransform( adfGeoTransform );
	oSRS.SetUTM( 11, TRUE );
	oSRS.SetWellKnownGeogCS( "NAD27" );
	oSRS.exportToWkt( &pszSRS_WKT );
	poDstDS->SetProjection( pszSRS_WKT );
	CPLFree( pszSRS_WKT );
	poBand = poDstDS->GetRasterBand(1);
	CPLErr err = poBand->RasterIO( GF_Write, 0, 0, 512, 512, abyRaster, 512, 512, GDT_Byte, 0, 0 );
	
	GDALClose( (GDALDatasetH) poDstDS );
	
	if (err == 4) {
		return false;
	}	
	return true;
}


