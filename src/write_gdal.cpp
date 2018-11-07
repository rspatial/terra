#include "spatraster.h"
#include "math_utils.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"


bool SpatRaster::writeRasterGDAL(std::string filename, bool overwrite) {
	bool success;
	SpatRaster r = geometry();
	
	if (!hasValues()) {
		addWarning("none of the cells have values");
	}
	success = r.writeStartGDAL(filename, overwrite);
	if (!success) {
		setError("cannot open file");
		return false;
	}
	success = r.writeValuesGDAL(getValues(), 0);
	if (!success) {
		setError("cannot write values to file");
		return false;
	}
	success = r.writeStopGDAL();
	if (!success) {
		setError("cannot close file");
		return false;
	}
	return success;
}



bool SpatRaster::writeStartGDAL(std::string filename, bool overwrite) {

	std::string format = "GTiff";     
	const char *pszFormat = format.c_str();
	const char *pszDstFilename = filename.c_str();
    GDALDriver *poDriver;
    char **papszMetadata;
	GDALAllRegister();
	
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if(poDriver == NULL) return (false);
    papszMetadata = poDriver->GetMetadata();
    if(! CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE)) return (false);
 
	GDALDataset *poDstDS;
	char **papszOptions = NULL;
	poDstDS = poDriver->Create( pszDstFilename, ncol, nrow, nlyr(), GDT_Float64, papszOptions);

	std::vector<double> rs = resolution();
	double adfGeoTransform[6] = { extent.xmin, rs[0], 0, extent.ymax, 0, -1 * rs[1] };
	poDstDS->SetGeoTransform(adfGeoTransform);
	
	std::string prj = getCRS();
	OGRSpatialReference oSRS;
	OGRErr erro = oSRS.importFromProj4(&prj[0]); 
	if (erro == 4) { return false ; }	// ??
	
	char *pszSRS_WKT = NULL;	
	oSRS.exportToWkt(&pszSRS_WKT);
	poDstDS->SetProjection(pszSRS_WKT);
	CPLFree(pszSRS_WKT);

	gdalconnection = poDstDS;

	source[0].range_min.resize(nlyr());
	source[0].range_max.resize(nlyr());
	for (size_t i =0; i<nlyr(); i++) {
		source[0].range_min[i] = std::numeric_limits<double>::max();
		source[0].range_max[i] = std::numeric_limits<double>::lowest();
	}

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



bool SpatRaster::writeStopGDAL() {
	GDALRasterBand *poBand;
	source[0].hasRange.resize(nlyr());
	for (size_t i=0; i < nlyr(); i++) {
		poBand = gdalconnection->GetRasterBand(i+1);
		poBand->SetStatistics(source[0].range_min[i], source[0].range_max[i], -9999., -9999.);
		source[0].hasRange[i] = true;
	}
	GDALClose( (GDALDatasetH) gdalconnection );
	return true;
}


/*
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
	GDALRasterBand *poBand;
	for (size_t i=0; i < nlyr(); i++) {
	
		vals = readValues(0,nrow,0,ncol,i,1);
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
*/

