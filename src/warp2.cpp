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




bool find_oputput_bounds(const GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, const std::string crs, std::string filename, std::string &errmsg) {


	CPLAssert( hSrcDS != NULL );

	// Create output with same datatype as first input band.
	GDALDataType eDT = GDALGetRasterDataType(GDALGetRasterBand(hSrcDS,1));

	// Get output driver (GeoTIFF format)
	GDALDriverH hDriver = GDALGetDriverByName( "GTiff" );
	CPLAssert( hDriver != NULL );

	// Get Source coordinate system.
	const char *pszSrcWKT = GDALGetProjectionRef( hSrcDS );
	CPLAssert( pszSrcWKT != NULL && strlen(pszSrcWKT) > 0 );

	errmsg = "";
	OGRSpatialReference* oSRS = new OGRSpatialReference;
	const char* s = crs.c_str();
	if (is_ogr_error(oSRS->SetFromUserInput(s), errmsg)) {
		out.setError(errmsg);
		return out;
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
	CPLAssert( hTransformArg != NULL );

	// Get approximate output georeferenced bounds and resolution for file.
	double adfDstGeoTransform[6];
	int nPixels=0, nLines=0;
	CPLErr eErr = GDALSuggestedWarpOutput( hSrcDS,
									GDALGenImgProjTransform, hTransformArg,
									adfDstGeoTransform, &nPixels, &nLines );
	CPLAssert( eErr == CE_None );
	GDALDestroyGenImgProjTransformer( hTransformArg );

	// Create the output DS.
	if (filename != "") {
		hDstDS = GDALCreate( hDriver, filename.c_str(), nPixels, nLines,
							GDALGetRasterCount(hSrcDS), eDT, NULL );
	} else {
		hDstDS = GDALCreate( "MEM", NULL, nPixels, nLines,
							GDALGetRasterCount(hSrcDS), eDT, NULL );
	}
	CPLAssert( hDstDS != NULL );

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



SpatRaster SpatRaster::warp2(const SpatRaster &x, std::string crs, std::string method, SpatOptions &opt) {

	SpatRaster out;

	GDALDatasetH hSrcDS, hDstDS;
	hSrcDS = open_gdal(*this);
	std::string errmsg="";
	std::string filename="test.tif";
	
	if (! find_oputput_bounds(hSrcDS, hDstDS, crs, filename, errmsg) {
		out.setError(msg);
		return out;
	}

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
    
	out = SpatRaster("out.tif");
	return out;
}	