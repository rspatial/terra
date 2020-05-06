#include "gdalwarper.h"
#include "ogr_spatialref.h"

#include "spatRaster.h"
//#include "string_utils.h"
#include "crs.h"


bool find_oputput_bounds(const GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, const std::string crs, std::string filename, std::string driver, int nlyrs, std::string &msg) {

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
	if (is_ogr_error(oSRS->SetFromUserInput( crs.c_str() ), msg)) {
		return false;
	};

	char *pszDstWKT = NULL;
	oSRS->exportToWkt( &pszDstWKT );

	// Create a transformer that maps from source pixel/line coordinates
	// to destination georeferenced coordinates (not destination
	// pixel line).  We do that by omitting the destination dataset
	// handle (setting it to NULL).
	void *hTransformArg;
	hTransformArg =
		GDALCreateGenImgProjTransformer( hSrcDS, pszSrcWKT, NULL, pszDstWKT, FALSE, 0, 1 );
	if (hTransformArg == NULL ) {
		msg = "cannot create TranformArg";
		return false;
	}

	// Get approximate output georeferenced bounds and resolution for file.
	double adfDstGeoTransform[6];
	int nPixels=0, nLines=0;
	CPLErr eErr = GDALSuggestedWarpOutput( hSrcDS, GDALGenImgProjTransform, 
					hTransformArg, adfDstGeoTransform, &nPixels, &nLines );

	GDALDestroyGenImgProjTransformer( hTransformArg );
	if ( eErr != CE_None ) {
		msg = "cannot create warp output";
		return false;		
	}

	// Create the output DS.

	GDALDriverH hDriver = GDALGetDriverByName( driver.c_str() );
	if ( hDriver == NULL ) {
		msg = "empty driver";
		return false;
	}
	if (driver == "MEM") {
		hDstDS = GDALCreate(hDriver, "", nPixels, nLines, nlyrs, eDT, NULL );
	} else {
		hDstDS = GDALCreate( hDriver, filename.c_str(), nPixels, nLines, nlyrs, eDT, NULL );
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
    //psWarpOptions->pfnProgress = GDALTermProgress;

	psWarpOptions->papszWarpOptions =	
      CSLSetNameValue( psWarpOptions->papszWarpOptions, "INIT_DEST", "NO_DATA");

	//psWarpOptions->padfSrcNoDataReal =
	//   (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount );
    //psWarpOptions->padfSrcNoDataReal[0] = -3.4e+38;

	psWarpOptions->padfDstNoDataReal =
	    (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount );
    psWarpOptions->padfDstNoDataReal[0] = NAN;

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
    oOperation.ChunkAndWarpImage( 0, 0, GDALGetRasterXSize( hDstDS ), GDALGetRasterYSize( hDstDS ) );
    GDALDestroyGenImgProjTransformer( psWarpOptions->pTransformerArg );
    GDALDestroyWarpOptions( psWarpOptions );
    
	return true;
}	


bool is_valid_warp_method(const std::string &method) {
	std::vector<std::string> m { "near", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode", "max", "min", "med", "q1", "q3", "sum" };
	return (std::find(m.begin(), m.end(), method) != m.end());
}


/*
SpatRaster SpatRaster::warp_crs(std::string crs, std::string method, SpatOptions &opt) {

	SpatRaster out;
	if (!is_valid_warp_method(method)) {
		out.setError("not a valid warp method");
		return out;
	}

	GDALDatasetH hSrcDS, hDstDS;
	hSrcDS = open_gdal(*this);

	//test crs first 
	std::string filename = opt.filename;
	std::string driver = filename == "" ? "MEM" : "GTiff";
	std::string errmsg="";

	if (! find_oputput_bounds(hSrcDS, hDstDS, crs, filename, driver, nlyr(), errmsg)) {
		out.setError(errmsg);
		return out;
	}
	
	bool success = gdal_warper(hSrcDS, hDstDS, method, errmsg);
	if (!success) {
		out.setError(errmsg);
		GDALClose( hSrcDS );
		GDALClose( hDstDS );
		return out;
	}

	GDALClose( hSrcDS );

	if (driver == "MEM") {
		bool test = out.setValuesMEM(hDstDS, false); 
		GDALClose( hDstDS );
		if (!test) {
			out.setError("wat nu?");
			return out;
		}
	} else {
		for (size_t i=0; i < nlyr(); i++) {
			GDALRasterBandH hBand = GDALGetRasterBand(hDstDS, i+1);
			double adfMinMax[2];
			GDALComputeRasterMinMax(hBand, true, adfMinMax);
			GDALSetRasterStatistics(hBand, adfMinMax[0], adfMinMax[1], NAN, NAN);		
		}
		GDALClose( hDstDS );
		out = SpatRaster(filename);
	}
	
	// should not be needed (but it is)
	out.setSRS({crs});
	
	return out;
}
*/


SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, SpatOptions &opt) {
	
	SpatRaster out = x.geometry();
	if (!is_valid_warp_method(method)) {
		out.setError("not a valid warp method");
		return out;
	}

	GDALDatasetH hSrcDS, hDstDS;
	if (!open_gdal(hSrcDS)) {
		out.setError("cannot create dataset");
		return out;
	}

	std::string filename = opt.filename;
	std::string driver = filename == "" ? "MEM" : "GTiff";
	std::string errmsg="";
	bool use_crs = crs != "";  
	bool get_geom = false;

	if (use_crs) {  // use the crs, ignore argument "x"
		get_geom = true;
		//test crs first? 
		if (! find_oputput_bounds(hSrcDS, hDstDS, crs, filename, driver, nlyr(), errmsg)) {
			out.setError(errmsg);
			return out;
		}
	} else {
		if (!out.createDS(hDstDS, filename, driver, opt.gdal_options, true)) {
			return(out);
		}
	}
	
	bool success = gdal_warper(hSrcDS, hDstDS, method, errmsg);
	
	GDALClose( hSrcDS );
	if (!success) {
		GDALClose( hDstDS );
		out.setError(errmsg);
		return out;
	}

	if (driver == "MEM") {
		bool test = out.setValuesMEM(hDstDS, get_geom); 
		GDALClose( hDstDS );
		if (!test) {
			out.setError("wat nu?");
			return out;
		}
	} else {
		for (size_t i=0; i < nlyr(); i++) {
			GDALRasterBandH hBand = GDALGetRasterBand(hDstDS, i+1);
			double adfMinMax[2];
			GDALComputeRasterMinMax(hBand, true, adfMinMax);
			GDALSetRasterStatistics(hBand, adfMinMax[0], adfMinMax[1], NAN, NAN);		
		}
		GDALClose( hDstDS );
		out = SpatRaster(filename);
	}

	// should not be needed (but it is)
	if (use_crs) out.setSRS({crs});	
	return out;
}


/*
SpatRaster SpatRaster::tester(bool geom) {

	SpatRaster out = geometry();
	GDALDatasetH hDS;
	if (source[0].driver != "memory") {
		out.setError("mwaaa");
		return out;		
	}

	if (!open_gdal(hDS)) {
		out.setError("cannot create dataset");
		return out;
	}
	bool test = out.setValuesMEM(hDS, geom); 
	GDALClose( hDS );
	if (!test) {
		out.setError("wat nu?");
	}
	return out;
}
*/

