#include "gdalwarper.h"
#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "crs.h"


bool SpatRaster::open_gdal(GDALDatasetH &hDS) {
	// needs to loop over sources. thus should vector of GDALDatasetH
	// for now just doing the first
	if (!source[0].hasValues) {
		return false;
	
	} else if (source[0].driver == "gdal") {
		std::string f = source[0].filename;
		hDS = GDALOpen(f.c_str(), GA_ReadOnly);
		return(hDS != NULL);
		
	} else { // in memory
	
		size_t nl = nlyr();
		size_t ncls = nrow() * ncol();
		GDALDriverH hDrv = GDALGetDriverByName("MEM");

/*https://gis.stackexchange.com/questions/196048/how-to-reuse-memory-pointer-of-gdal-memory-driver
		char **papszOptions = NULL;
		hDS = GDALCreate(hDrv, "", ncol(), nrow(), 0, GDT_Float64, papszOptions);
		if (hDS == NULL) return false;
		std::vector<double> vals;
		for(size_t i=0; i<nl; i++)	{
			size_t off = ncls * i;
			vals = std::vector<double>(source[0].values.begin() +off, source[0].values.begin() +off+ncls);
			char szPtrValue[128] = { '\0' };
			int nRet = CPLPrintPointer( szPtrValue, reinterpret_cast<void*>(&vals[0]), sizeof(szPtrValue) );
			szPtrValue[nRet] = 0;
			papszOptions = CSLSetNameValue(papszOptions, "DATAPOINTER", szPtrValue);
			GDALAddBand(hDS, GDT_Float64, papszOptions);
		}
		CSLDestroy(papszOptions);
*/

		size_t nr = nrow();
		size_t nc = ncol();
		hDS = GDALCreate(hDrv, "", nc, nr, nl, GDT_Float64, NULL);
		if (hDS == NULL) return false;

		CPLErr err = CE_None;
		std::vector<double> vals;
		for (size_t i=0; i < nl; i++) {
			GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
			GDALSetRasterNoDataValue(hBand, NAN);
			size_t offset = ncls * i;
			vals = std::vector<double>(source[0].values.begin() + offset, source[0].values.begin() + offset + ncls);
			err = GDALRasterIO(hBand, GF_Write, 0, 0, nc, nr, &vals[0], nc, nr, GDT_Float64, 0, 0 );
			if (err != CE_None) {
				return false;
			}
		}
	}

	return true;
}


bool SpatRaster::setValuesMEM(GDALDatasetH hDS, bool set_geometry) {

	if (set_geometry) {
		RasterSource s;
		s.ncol = GDALGetRasterXSize( hDS );
		s.nrow = GDALGetRasterYSize( hDS );
		s.nlyr = GDALGetRasterCount( hDS );

		double adfGeoTransform[6];
		if( GDALGetGeoTransform( hDS, adfGeoTransform ) != CE_None ) {
			setError("Cannot get geotransform");
			return false;
		}
		double xmin = adfGeoTransform[0]; 
		double xmax = xmin + adfGeoTransform[1] * s.ncol; 
		double ymax = adfGeoTransform[3]; 
		double ymin = ymax + s.nrow * adfGeoTransform[5]; 
		s.extent = SpatExtent(xmin, xmax, ymin, ymax);

		s.driver = "memory";
		setSource(s);
		
		OGRSpatialReferenceH srs = GDALGetSpatialRef( hDS );
		char *cp;
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = OSRExportToWktEx(srs, &cp, options);
		std::string errmsg;
		if (is_ogr_error(err, errmsg)) {
			CPLFree(cp);
			return false;
		}
		std::string wkt = std::string(cp);
		CPLFree(cp);
		std::string msg;
		if (!s.srs.set({wkt}, msg)) {
			setError(msg);
			return false;
		}
	}
	
	
	source[0].values.reserve(ncell() * nlyr());
	CPLErr err = CE_None;
	int hasNA;
	for (size_t i=0; i < nlyr(); i++) {
		GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
		std::vector<double> lyrout( ncell() );
		err = GDALRasterIO(hBand, GF_Read, 0, 0, ncol(), nrow(), &lyrout[0], ncol(), nrow(), GDT_Float64, 0, 0);
		if (err != CE_None ) {
			setError("CE_None");
			return false;
		}
		
		//double naflag = -3.4e+38;
		double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
		if (hasNA) std::replace(lyrout.begin(), lyrout.end(), naflag, (double) NAN);
		source[0].values.insert(source[0].values.end(), lyrout.begin(), lyrout.end());
	}
	source[0].hasValues=TRUE;
	source[0].driver="memory";
	source[0].setRange();
	return true;
}


bool gdal_ds_create(SpatRaster x, GDALDatasetH &hDS, std::string filename, std::string driver, std::vector<std::string> foptions, bool fill, std::string &msg) {

	msg = "";
	
	char **papszOptions = NULL;
	for (size_t i=0; i < foptions.size(); i++) {
		std::vector<std::string> wopt = strsplit(foptions[i], "=");
		if (wopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, wopt[0].c_str(), wopt[1].c_str() );
		}
	}

	const char *pszFormat = driver.c_str();
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
		//GDALSetRasterNoDataValue(hBand, -3.4e+38);
		if (fill) GDALFillRaster(hBand, NAN, 0);
	}

	std::vector<double> rs = x.resolution();
	SpatExtent e = x.getExtent();
	double adfGeoTransform[6] = { e.xmin, rs[0], 0, e.ymax, 0, -1 * rs[1] };
	GDALSetGeoTransform( hDS, adfGeoTransform);

	OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
	std::vector<std::string> srs = x.getSRS();
	std::string wkt = srs[1];
	OGRErr erro = OSRSetFromUserInput(hSRS, wkt.c_str());
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
		if (!gdal_ds_create(out, hDstDS, filename, driver, opt.gdal_options, true, errmsg)) {
			out.setError(errmsg);
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

