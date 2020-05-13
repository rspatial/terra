#include "gdalwarper.h"
#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "crs.h"
#include "file_utils.h"


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
		hDstDS = GDALCreate( hDriver, "", nPixels, nLines, nlyrs, eDT, NULL );
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


GDALResampleAlg getAlgo(std::string m) {
	GDALResampleAlg alg;
	if ( m == "near" ) { 
		alg = GRA_NearestNeighbour;
	} else if (m=="bilinear") {
		alg = GRA_Bilinear;
	} else if (m=="cubic") {		
		alg = GRA_Cubic;
	} else if (m=="cubicspline") {
		alg = GRA_CubicSpline; 	
	} else if (m=="lanczos") {
		alg = GRA_Lanczos;
	} else if (m=="mean") {
		alg = GRA_Average;	
//	} else if (m=="sum") {
//		alg = GRA_Sum;
	} else if (m=="mode") {
		alg = GRA_Mode; 	
	} else if (m=="max") {
		alg = GRA_Max;
	} else if (m=="min") {
		alg = GRA_Min;	
	} else if (m=="median") {
		alg = GRA_Med;	
	} else if (m=="q1") {
		alg = GRA_Q1;	
	} else if (m=="q3") {
		alg = GRA_Q3;
	} else { 
		alg = GRA_NearestNeighbour;
	}	
	return alg;
}


bool gdal_warper(GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, std::vector<unsigned> srcbands, std::vector<unsigned> dstbands, std::string method, std::string msg) {

	if (srcbands.size() != dstbands.size()) {
		msg = "number of source bands must match number of dest bands";
		return false;
	}
	int nbands = srcbands.size();

	GDALResampleAlg a = getAlgo(method);
	
    // Setup warp options.
    GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
    psWarpOptions->hSrcDS = hSrcDS;
    psWarpOptions->hDstDS = hDstDS;
	
	psWarpOptions->eResampleAlg = a;
	
    psWarpOptions->nBandCount = nbands;
    psWarpOptions->panSrcBands =
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
    psWarpOptions->panDstBands =
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
	psWarpOptions->padfDstNoDataReal =
	    (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount );
	
	for (int i=0; i<nbands; i++) {
		psWarpOptions->panSrcBands[i] = (int) srcbands[i]+1;
		psWarpOptions->panDstBands[i] = (int) dstbands[i]+1;
		//psWarpOptions->padfSrcNoDataReal[0] = -3.4e+38;
		psWarpOptions->padfDstNoDataReal[i] = NAN;
    }
	
	//psWarpOptions->pfnProgress = GDALTermProgress;

	psWarpOptions->papszWarpOptions =	
      CSLSetNameValue( psWarpOptions->papszWarpOptions, "INIT_DEST", "NO_DATA");

    // Establish reprojection transformer.
    psWarpOptions->pTransformerArg =
        GDALCreateGenImgProjTransformer( hSrcDS, GDALGetProjectionRef(hSrcDS),
                                        hDstDS, GDALGetProjectionRef(hDstDS),
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


SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, SpatOptions &opt) {

	SpatRaster out = x.geometry(nlyr());

	if (!is_valid_warp_method(method)) {
		out.setError("not a valid warp method");
		return out;
	}
	lrtrim(crs);
	std::string errmsg;
	std::string filename = opt.filename;

	bool use_crs = crs != "";  
	// should not be needed (need to fix)
	
	if ((!use_crs) & (!hasValues())) {
		if (filename != "") {
			addWarning("raster has no values, not writing to file");
		}
		return out;
	}

	if (filename == "") {
		if (!canProcessInMemory(4) || opt.get_todisk()) {
			filename = tempFile(opt.get_tempdir(), ".tif");
		} 
	} else {
		if (!can_write(filename, opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}
	}
	
	if (opt.names.size() == out.nlyr()) {
		out.setNames(opt.names);
	}
	
	if (!hasValues()) filename = ""; // for crs case
	std::string driver = filename == "" ? "MEM" : "GTiff";

	GDALDatasetH hSrcDS, hDstDS;
	size_t ns = nsrc();
	int bandstart = 0;
	
	for (size_t i=0; i<ns; i++) {
		
		if (!open_gdal(hSrcDS, i)) {
			out.setError("cannot create dataset from source");
			return out;
		}

		// create dest source, only once 
		if (i==0) {
			 // use the crs, ignore argument "x"
			if (use_crs) {
				if (! find_oputput_bounds(hSrcDS, hDstDS, crs, filename, driver, nlyr(), errmsg)) {
					out.setError(errmsg);
					GDALClose( hSrcDS );
					return out;
				}
				if (!hasValues()) {
					if (!out.from_gdalMEM(hDstDS, use_crs, false)) {
						out.setError("cannot get geometry from mem");
					} 
					GDALClose( hSrcDS );
					GDALClose( hDstDS );
					out.setSRS({crs});	 // fix the need for this
					return out;
				}
			} else {
				if (!out.create_gdalDS(hDstDS, filename, driver, false, NAN, opt.gdal_options)) {
					GDALClose( hSrcDS );
					//GDALClose( hDstDS );
					return out;
				}
			}
		}
		std::vector<unsigned> srcbands = source[i].layers;
		std::vector<unsigned> dstbands(srcbands.size()); 
		std::iota (dstbands.begin(), dstbands.end(), bandstart); 
		bandstart += dstbands.size();
		bool success = gdal_warper(hSrcDS, hDstDS, srcbands, dstbands, method, errmsg);
	
		GDALClose( hSrcDS );
		if (!success) {
			GDALClose( hDstDS );
			out.setError(errmsg);
			return out;
		}
	}
	
	if (driver == "MEM") {
		bool test = out.from_gdalMEM(hDstDS, use_crs, true); 
		GDALClose( hDstDS );
		if (!test) {
			out.setError("wat nu?");
			return out;
		}
	} else {
		for (size_t i=0; i < nlyr(); i++) {
			GDALRasterBandH hBand = GDALGetRasterBand(hDstDS, i+1);
			double adfMinMax[2];
			bool approx = ncell() > 10e+8;
			GDALComputeRasterMinMax(hBand, approx, adfMinMax);
			GDALSetRasterStatistics(hBand, adfMinMax[0], adfMinMax[1], NAN, NAN);		
		}
		GDALClose( hDstDS );
		out = SpatRaster(filename);
	}

	//if (use_crs) out.setSRS({crs});	// fix the need for this
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

