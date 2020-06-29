// Copyright (c) 2018-2020  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.


#include "gdalwarper.h"
#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "file_utils.h"

#include "crs.h"

//#include <vector>
//#include "vecmath.h"


SpatVector SpatRaster::dense_extent() {
		
	std::vector<long> rows(nrow());
	std::iota(rows.begin(), rows.end(), 0);
	std::vector<long> cols(ncol());
	std::iota(cols.begin(), cols.end(), 0);

	std::vector<double> xcol = xFromCol(cols) ;
	std::vector<double> yrow = yFromRow(rows) ;

	std::vector<double> y0(ncol(), yFromRow(nrow()-1));
	std::vector<double> y1(ncol(), yFromRow(0));
	std::vector<double> x0(nrow(), xFromCol(0));
	std::vector<double> x1(nrow(), xFromCol(ncol()-1));

	std::vector<double> x = x0;
	std::vector<double> y = yrow;
	x.insert(x.end(), xcol.begin(), xcol.end());
	y.insert(y.end(), y0.begin(), y0.end());
	
	std::reverse(yrow.begin(), yrow.end());
	std::reverse(xcol.begin(), xcol.end());

	x.insert(x.end(), x1.begin(), x1.end());
	y.insert(y.end(), yrow.begin(), yrow.end() );
	x.insert(x.end(), xcol.begin(), xcol.end());
	y.insert(y.end(), y1.begin(), y1.end());
	
	SpatVector v(x, y, polygons, getSRS("wkt"));

	return v;
}


#if GDAL_VERSION_MAJOR >= 3

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
        (int *) CPLMalloc(sizeof(int) * nbands );
    psWarpOptions->panDstBands =
        (int *) CPLMalloc(sizeof(int) * nbands );
	psWarpOptions->padfSrcNoDataReal =
	    (double *) CPLMalloc(sizeof(double) * nbands );
	psWarpOptions->padfDstNoDataReal =
	    (double *) CPLMalloc(sizeof(double) * nbands );
	
	GDALRasterBandH hBand;
	int hasNA;
	for (int i=0; i<nbands; i++) {
		psWarpOptions->panSrcBands[i] = (int) srcbands[i]+1;
		psWarpOptions->panDstBands[i] = (int) dstbands[i]+1;

		hBand = GDALGetRasterBand(hSrcDS, srcbands[i]+1);
		double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
		if (hasNA) {
			psWarpOptions->padfSrcNoDataReal[i] = naflag;
		} else {
			psWarpOptions->padfSrcNoDataReal[i] = NAN;			
		}
		psWarpOptions->padfDstNoDataReal[i] = NAN;
    }
	
	//psWarpOptions->pfnProgress = GDALTermProgress;

	psWarpOptions->papszWarpOptions =	
      CSLSetNameValue( psWarpOptions->papszWarpOptions, "INIT_DEST", "NO_DATA");
	psWarpOptions->papszWarpOptions =	
      CSLSetNameValue( psWarpOptions->papszWarpOptions, "WRITE_FLUSH", "YES");

//GDALWarpInitSrcNoDataReal(GDALWarpOptions *psOptionsIn, double dNoDataReal)
//void GDALWarpInitDstNoDataReal(GDALWarpOptions *psOptionsIn, double dNoDataReal)

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


SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, bool mask, SpatOptions &opt) {

	SpatRaster out = x.geometry(nlyr());

	if (!is_valid_warp_method(method)) {
		out.setError("not a valid warp method");
		return out;
	}
	lrtrim(crs);
	std::string errmsg;
	SpatOptions mopt;
	if (mask) {
		mopt = opt;
		opt = SpatOptions(opt);
	}
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
		if (!canProcessInMemory(4, opt.get_memfrac()) || opt.get_todisk()) {
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
		out = SpatRaster(filename, -1, "");
	}
	
	if (mask) {
		SpatVector v = dense_extent();
		v = v.project(out.getSRS("wkt"));
		out = out.mask(v, false, NAN, mopt);
	}
	
	return out;
}


#else 
	

SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, bool mask, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);

	if (crs != "") {
		out.setError("You cannot project by specifying a crs with your version of GDAL");
		return out;
	}

	out.setNames(getNames());
	std::vector<std::string> f {"bilinear", "near"};
	if (std::find(f.begin(), f.end(), method) == f.end()) {
		out.setError("unknown warp method");
		return out;
	}
	if (!hasValues()) {
		return out;
	}

	std::string crsin = srs.wkt;
	std::string crsout = out.srs.wkt;
	bool do_prj = true;
	if ((crsin == crsout) || (crsin == "") || (crsout == "")) {
		do_prj = false;
	}

	if (!do_prj) {
		SpatExtent e = out.extent;
		e.intersect(extent);
		if (!e.valid()) {
			out.addWarning("No spatial overlap");
			return out;
		}
	}

	SpatRaster xx;
	if (do_prj) {
		xx = *this;
	} else {
		unsigned xq = x.xres() / xres();
		unsigned yq = x.yres() / yres();
		if (std::max(xq, yq) > 1) {
			xq = xq == 0 ? 1 : xq;
			yq = yq == 0 ? 1 : yq;
			std::vector<unsigned> agf = {yq, xq, 1};
			SpatOptions agopt;
			if (method == "bilinear") {
				xx = aggregate(agf, "mean", true, agopt);
			} else {
				xx = aggregate(agf, "modal", true, agopt);
			}
		} else {
			xx = *this;
		}
	}
	unsigned nc = out.ncol();

  	if (!out.writeStart(opt)) { return out; }
	for (size_t i = 0; i < out.bs.n; i++) {
        unsigned firstcell = out.cellFromRowCol(out.bs.row[i], 0);
		unsigned lastcell  = out.cellFromRowCol(out.bs.row[i]+out.bs.nrows[i]-1, nc-1);
		std::vector<double> cells(1+lastcell-firstcell);
		std::iota (std::begin(cells), std::end(cells), firstcell);
        std::vector<std::vector<double>> xy = out.xyFromCell(cells);
		if (do_prj) {
			#ifdef useGDAL
			out.msg = transform_coordinates(xy[0], xy[1], crsout, crsin);
			#else
			out.setError("GDAL is needed for crs transformation, but not available");
			return out;
			#endif
		}
		std::vector<std::vector<double>> v = xx.extractXY(xy[0], xy[1], method);
		if (!out.writeValues2(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
	}
	out.writeStop();
	return(out);
}



#endif
