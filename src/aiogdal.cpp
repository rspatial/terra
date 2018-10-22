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
		
	RasterSource s;
	s.ncol = poDataset->GetRasterXSize();
	s.nrow = poDataset->GetRasterYSize();
	s.nlyr = poDataset->GetRasterCount();
	
	double adfGeoTransform[6];
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None ) {
		// the rounding below is to address a design flaw in GDAL
		// GDAL provides the coordinates of one corner and the resolution, 
		// instead of the coordinates of all (two opposite) corners.
		// This makes computation of the opposite corner coordinates only 
		// approximate for large rasters with a high resolution.
		double xmin = adfGeoTransform[0]; /* left x */
		double xmax = xmin + adfGeoTransform[1] * s.ncol; /* w-e pixel resolution */
		//xmax = roundn(xmax, 9);
		double ymax = adfGeoTransform[3]; /* top y */
		double ymin = ymax + s.nrow * adfGeoTransform[5]; /* n-s pixel resolution (negative value) */
		//ymin = roundn(ymin, 9);
		SpatExtent e(xmin, xmax, ymin, ymax);
		s.extent = e;
	}
		
	s.memory = false;
	s.filename = fname;	
	s.driver = "gdal";
	

	if( poDataset->GetProjectionRef() != NULL ) {
		OGRSpatialReference oSRS(poDataset->GetProjectionRef());
		char *pszPRJ = NULL;
		oSRS.exportToProj4(&pszPRJ);
		s.crs = pszPRJ;
	} else {
		s.crs = "";
	}


	GDALRasterBand  *poBand;
	//int nBlockXSize, nBlockYSize;
	double adfMinMax[2];
	int bGotMin, bGotMax;
	
//	s.layers.resize(1);
	for (size_t i = 0; i < s.nlyr; i++) {
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
	setSource(s);
	return true;
 
}

std::vector<double> SpatRaster::readValuesGDAL(unsigned row, unsigned nrows, unsigned col, unsigned ncols, unsigned lyr, unsigned nlyrs) {
	
    GDALDataset  *poDataset;
	GDALRasterBand  *poBand;
    GDALAllRegister();
	
	const char* pszFilename = source[0].filename.c_str();
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	std::vector<double> out;
	unsigned ncell = ncols*nrows;
//	double *pafScanline;
//	pafScanline = (double *) CPLMalloc(sizeof(double)*ncell);
	std::vector<double> scan(ncell);
	
	for (size_t i=0; i < nlyrs; i++) {
		
		poBand = poDataset->GetRasterBand(lyr + i + 1);
		CPLErr err = poBand->RasterIO( GF_Read, row, col, ncols, nrows, &scan[0], ncols, nrows, GDT_Float64, 0, 0 );	
		if (err == 4) {
			std::vector<double> errout;
			return errout;
		}
		out.insert(out.end(), scan.begin(), scan.end());
	}
	
	//CPLFree(pafScanline);
	GDALClose( (GDALDatasetH) poDataset );
	return(out);
}



bool SpatRaster::writeRasterGDAL(std::string filename, bool overwrite) {

	std::string format = "GTiff"; 
    
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
	poDstDS = poDriver->Create( pszDstFilename, ncol, nrow, nlyr(), GDT_Float64, papszOptions );

	std::vector<double> rs = resolution();
	double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
	poDstDS->SetGeoTransform( adfGeoTransform );
	

	string prj = getCRS();
	OGRSpatialReference oSRS;
	//std::vector<char> chars(prj.c_str(), prj.c_str() + prj.size() + 1u);
	//char *pszProj4 = &chars[0];
	//OGRErr erro = oSRS.importFromProj4( pszProj4); 
	OGRErr erro = oSRS.importFromProj4( &prj[0] ); 
	if (erro == 4) { return false ; }	// ??
	
	char *pszSRS_WKT = NULL;	
	oSRS.exportToWkt( &pszSRS_WKT );
	poDstDS->SetProjection( pszSRS_WKT );
	CPLFree( pszSRS_WKT );

	std::vector<double> rmin = range_min();
	std::vector<double> rmax = range_max();
	CPLErr err;
	bool result = true;
	std::vector<double> vals;
//	double* v = ( double* ) CPLMalloc( sizeof(double) * ncell() );
	GDALRasterBand *poBand;
	for (size_t i=0; i < nlyr(); i++) {
	
		vals = readValues(0,nrow,0,ncol,i,1);
	//	v = &vals[0];	
		poBand = poDstDS->GetRasterBand(i+1);
		err = poBand->RasterIO( GF_Write, 0, 0, ncol, nrow, &vals[0], ncol, nrow, GDT_Float64, 0, 0 );
		if (err == 4) {
			result = false;
			break;
		}
		poBand->SetStatistics(rmin[i], rmax[i], -9999., -9999.);
	}
	
	GDALClose( (GDALDatasetH) poDstDS );
	
	return result;
}



bool SpatRaster::writeStartGDAL(std::string filename, bool overwrite) {

	std::string format = "GTiff";     
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
	poDstDS = poDriver->Create( pszDstFilename, ncol, nrow, nlyr(), GDT_Float64, papszOptions );

	std::vector<double> rs = resolution();
	double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
	poDstDS->SetGeoTransform( adfGeoTransform );
	
	string prj = getCRS();
	OGRSpatialReference oSRS;
	OGRErr erro = oSRS.importFromProj4( &prj[0] ); 
	if (erro == 4) { return false ; }	// ??
	
	char *pszSRS_WKT = NULL;	
	oSRS.exportToWkt( &pszSRS_WKT );
	poDstDS->SetProjection( pszSRS_WKT );
	CPLFree( pszSRS_WKT );

	gdalconnection = poDstDS;

	source[0].range_min.resize(nlyr());
	source[0].range_max.resize(nlyr());
	for (size_t i =0; i<nlyr(); i++) {
		source[0].range_min[i] = std::numeric_limits<double>::max();
		source[0].range_max[i] = std::numeric_limits<double>::lowest();
	}

	return true;
}


bool SpatRaster::writeStopGDAL() {
	GDALRasterBand *poBand;
	for (size_t i=0; i < nlyr(); i++) {
		poBand = gdalconnection->GetRasterBand(i+1);
		poBand->SetStatistics(source[0].range_min[i], source[0].range_max[i], -9999., -9999.);
	}
	GDALClose( (GDALDatasetH) gdalconnection );
	return true;

}

bool SpatRaster::writeValuesGDAL(std::vector<double> vals, unsigned row){
	unsigned nrows = vals.size() / (nlyr() * ncol);
	unsigned start;
	CPLErr err;
	GDALRasterBand *poBand;
	double vmin, vmax;
	unsigned nc = nrows * ncol;
	for (size_t i=0; i < nlyr(); i++) {
		start = nc * i;
		poBand = gdalconnection->GetRasterBand(i+1);
		err = poBand->RasterIO( GF_Write, 0, row, ncol, nrows, &vals[start], ncol, nrows, GDT_Float64, 0, 0 );
		if (err == 4) break;
		minmax(vals.begin()+start, vals.begin()+start+nc, vmin, vmax);
		source[0].range_min[i] = std::min(source[0].range_min[i], vmin);
		source[0].range_max[i] = std::max(source[0].range_max[i], vmax); 
	}
	return true;
}

