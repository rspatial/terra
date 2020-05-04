#include "gdalwarper.h"
#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "crs.h"

GDALDatasetH open_gdal(const SpatRaster& x) {
	std::string f = x.source[0].filename;
	GDALDatasetH hSrcDS = GDALOpen(f.c_str(), GA_ReadOnly );
	return hSrcDS;
}



bool gdal_ds_create(SpatRaster x, GDALDatasetH &hDS, bool fill, SpatOptions &opt, std::string &msg) {

	msg = "";
	
	char **papszOptions = NULL;
	for (size_t i=0; i<opt.gdal_options.size(); i++) {
		std::vector<std::string> wopt = strsplit(opt.gdal_options[i], "=");
		if (wopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, wopt[0].c_str(), wopt[1].c_str() );
		}
	}
	std::string filename = opt.filename;
	std::string format;
	if (filename == "") {
		format = "MEM";
	} else {
		format = "GTiff";
	}
		
	const char *pszFormat = format.c_str();
	GDALDriverH hDrv = GDALGetDriverByName(pszFormat);

	const char *pszFilename = filename.c_str();
	hDS = GDALCreate(hDrv, pszFilename, x.ncol(), x.nrow(), x.nlyr(), GDT_Float64, papszOptions );
	CSLDestroy( papszOptions );

	GDALRasterBandH hBand;
	std::vector<std::string> nms = x.getNames();
	for (size_t i=0; i < x.nlyr(); i++) {
		hBand = GDALGetRasterBand(hDS, i+1);
		GDALSetDescription(hBand, nms[i].c_str());
		GDALSetRasterNoDataValue(hBand, NAN);
		if (fill) GDALFillRaster(hBand, NAN, 0);
	}

	std::vector<double> rs = x.resolution();
	SpatExtent e = x.getExtent();
	double adfGeoTransform[6] = { e.xmin, rs[0], 0, e.ymax, 0, -1 * rs[1] };
	GDALSetGeoTransform( hDS, adfGeoTransform);

	OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
	std::vector<std::string> srs = x.getSRS();
	std::string prj = srs[1];
	OGRErr erro = OSRImportFromProj4(hSRS, &prj[0]);
	if (erro == 4) {
		msg = "CRS failure";
		OSRDestroySpatialReference( hSRS );
		return false;
	}
	
	char *pszSRS_WKT = NULL;
	OSRExportToWkt( hSRS, &pszSRS_WKT );
	GDALSetProjection( hDS, pszSRS_WKT );
	CPLFree(pszSRS_WKT);
	OSRDestroySpatialReference( hSRS );
	return true;
}




bool find_oputput_bounds(const GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, const std::string crs, std::string filename, std::string &msg) {

	msg = "";
	if ( hSrcDS == NULL ) {
		msg = "data source is NULL";
		return false;
	}

	// Create output with same datatype as first input band.
	GDALDataType eDT = GDALGetRasterDataType(GDALGetRasterBand(hSrcDS,1));

	// Get output driver (GeoTIFF format)

	// Get Source coordinate system.
	const char *pszSrcWKT = GDALGetProjectionRef( hSrcDS );
	if ( pszSrcWKT == NULL || strlen(pszSrcWKT) == 0 ) {
		msg = "data source has no WKT";
		return false;		
	}

	OGRSpatialReference* oSRS = new OGRSpatialReference;
	const char* s = crs.c_str();
	if (is_ogr_error(oSRS->SetFromUserInput(s), msg)) {
		return false;
	};

	//oSRS.SetWellKnownGeogCS( "WGS84" );
	
	char *pszDstWKT = NULL;
	oSRS->exportToWkt( &pszDstWKT );

	// Create a transformer that maps from source pixel/line coordinates
	// to destination georeferenced coordinates (not destination
	// pixel line).  We do that by omitting the destination dataset
	// handle (setting it to NULL).
	void *hTransformArg;
	hTransformArg =
		GDALCreateGenImgProjTransformer( hSrcDS, pszSrcWKT, NULL, pszDstWKT,
										 FALSE, 0, 1 );
	if (hTransformArg == NULL ) {
		msg = "cannot create TranformArg";
		return false;
	}

	// Get approximate output georeferenced bounds and resolution for file.
	double adfDstGeoTransform[6];
	int nPixels=0, nLines=0;
	CPLErr eErr = GDALSuggestedWarpOutput( hSrcDS,
									GDALGenImgProjTransform, hTransformArg,
									adfDstGeoTransform, &nPixels, &nLines );

	GDALDestroyGenImgProjTransformer( hTransformArg );
	if ( eErr != CE_None ) {
		msg = "cannot create warp output";
		return false;		
	}

	// Create the output DS.
	if (filename != "") {
		GDALDriverH hDriver = GDALGetDriverByName( "GTiff" );
		CPLAssert( hDriver != NULL );
		hDstDS = GDALCreate( hDriver, filename.c_str(), nPixels, nLines,
							GDALGetRasterCount(hSrcDS), eDT, NULL );
	} else {
		GDALDriverH hDriver = GDALGetDriverByName( "MEM" );
		CPLAssert( hDriver != NULL );

		hDstDS = GDALCreate(hDriver, NULL, nPixels, nLines,
							GDALGetRasterCount(hSrcDS), eDT, NULL );
	}
	
	if ( hDstDS == NULL ) {
		msg = "cannot create output dataset";
		return false;		
	}

	// Write out the projection definition.
	GDALSetProjection( hDstDS, pszDstWKT );
	GDALSetGeoTransform( hDstDS, adfDstGeoTransform );

	// Copy the color table, if required.
	GDALColorTableH hCT;
	hCT = GDALGetRasterColorTable( GDALGetRasterBand(hSrcDS,1) );
	if( hCT != NULL )
		GDALSetRasterColorTable( GDALGetRasterBand(hDstDS,1), hCT );
	
	return true;
}




bool gdal_warper(GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, std::string method, std::string msg) {

	msg="";

    // Setup warp options.
    GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
    psWarpOptions->hSrcDS = hSrcDS;
    psWarpOptions->hDstDS = hDstDS;
    psWarpOptions->nBandCount = 1;
    psWarpOptions->panSrcBands =
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
    psWarpOptions->panSrcBands[0] = 1;
    psWarpOptions->panDstBands =
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
    psWarpOptions->panDstBands[0] = 1;
    psWarpOptions->pfnProgress = GDALTermProgress;

    // Establish reprojection transformer.
    psWarpOptions->pTransformerArg =
        GDALCreateGenImgProjTransformer( hSrcDS,
                                        GDALGetProjectionRef(hSrcDS),
                                        hDstDS,
                                        GDALGetProjectionRef(hDstDS),
                                        FALSE, 0.0, 1 );
    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;

    // Initialize and execute the warp operation.
    GDALWarpOperation oOperation;
    oOperation.Initialize( psWarpOptions );
    oOperation.ChunkAndWarpImage( 0, 0,
                                GDALGetRasterXSize( hDstDS ),
                                GDALGetRasterYSize( hDstDS ) );
    GDALDestroyGenImgProjTransformer( psWarpOptions->pTransformerArg );
    GDALDestroyWarpOptions( psWarpOptions );
    GDALClose( hDstDS );
    GDALClose( hSrcDS );
    
	return true;
}	

bool is_valid_warp_method(const std::string &method) {
	std::vector<std::string> m { "near", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode", "max", "min", "med", "q1", "q3", "sum" };
	return (std::find(m.begin(), m.end(), method) != m.end());
}


SpatRaster SpatRaster::warp_crs(std::string crs, std::string method, SpatOptions &opt) {

	SpatRaster out;
	if (!is_valid_warp_method(method)) {
		out.setError("not a valid warp method");
		return out;
	}

	GDALDatasetH hSrcDS, hDstDS;
	hSrcDS = open_gdal(*this);

	std::string errmsg="";
	std::string filename = opt.get_filename();
	
	//test crs first 
	
	if (! find_oputput_bounds(hSrcDS, hDstDS, crs, filename, errmsg)) {
		out.setError(errmsg);
		return out;
	}
	
	std::string msg;
	bool success = gdal_warper(hSrcDS, hDstDS, method, msg);

	GDALClose( hSrcDS );
	GDALClose( hDstDS );

	if (!success) {
		out.setError(msg);
		return out;
	}
	return out;
}



SpatRaster SpatRaster::warp_rst(const SpatRaster &x, std::string method, SpatOptions &opt) {

	SpatRaster out;
	if (!is_valid_warp_method(method)) {
		out.setError("not a valid warp method");
		return out;
	}

	GDALDatasetH hSrcDS, hDstDS;
	hSrcDS = open_gdal(*this);

	std::string msg="";
	if (!gdal_ds_create(x, hDstDS, true, opt, msg)) {
		out.setError(msg);
		return(out);
	}

	bool success = gdal_warper(hSrcDS, hDstDS, method, msg);
	GDALClose( hSrcDS );
	GDALClose( hDstDS );

	if (!success) {
		out.setError(msg);
		return out;
	}
	return out;
}

