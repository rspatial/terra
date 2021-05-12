// Copyright (c) 2018-2021  Robert J. Hijmans
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
#include "gdal_alg.h"
#include "ogrsf_frmts.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "file_utils.h"
#include "vecmath.h"

#include "crs.h"
#include "gdalio.h"
#include "recycle.h"


//#include <vector>
//#include "vecmath.h"


SpatVector SpatRaster::dense_extent() {

	std::vector<int_64> rows, cols;
	if (nrow() < 51) {
		rows.resize(nrow());
		std::iota(rows.begin(), rows.end(), 0);
	} else {
		rows = seq_steps((int_64) 0, (int_64) nrow(), 50);
		rows[rows.size()-1] = nrow()-1;
	} 
	if (ncol() < 20) {
		cols.resize(nrow());
		std::iota(cols.begin(), cols.end(), 0);
	} else {
		cols = seq_steps((int_64) 0, (int_64) ncol(), 50);
		cols[cols.size()-1] = ncol()-1;
	} 
	

	std::vector<double> xcol = xFromCol(cols) ;
	std::vector<double> yrow = yFromRow(rows) ;

	SpatExtent e = getExtent();
	std::vector<double> y0(cols.size(), e.ymin);
	std::vector<double> y1(cols.size(), e.ymax);
	std::vector<double> x0(rows.size(), e.xmin);
	std::vector<double> x1(rows.size(), e.xmax);

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

#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR < 2

SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, bool mask, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);
	out.setNames(getNames());

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

	std::string crsin = source[0].srs.wkt;
	std::string crsout = out.source[0].srs.wkt;
	bool do_prj = true;
	if ((crsin == crsout) || (crsin == "") || (crsout == "")) {
		do_prj = false;
	}

	if (!do_prj) {
		SpatExtent e = out.getExtent();
		e.intersect(getExtent());
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
		std::vector<std::vector<double>> v = xx.extractXY(xy[0], xy[1], method, false);
		if (!out.writeValues2(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
	}
	out.writeStop();
	return(out);
}




#else


bool find_output_bounds(const GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, std::string srccrs, const std::string dstcrs, std::string filename, std::string driver, int nlyrs, std::string datatype, std::string &msg) {

	msg = "";
	if ( hSrcDS == NULL ) {
		msg = "data source is NULL";
		return false;
	}

	// Create output with same datatype as first input band.
	//GDALDataType eDT = GDALGetRasterDataType(GDALGetRasterBand(hSrcDS,1));
	GDALDataType eDT;
	getGDALDataType(datatype, eDT);

	// Get output driver (GeoTIFF format)

	// Get Source coordinate system.
	// const char *pszSrcWKT = GDALGetProjectionRef( hSrcDS );
	const char *pszSrcWKT = srccrs.c_str();
	if ( pszSrcWKT == NULL || strlen(pszSrcWKT) == 0 ) {
		msg = "data source has no WKT";
		return false;
	}

	OGRSpatialReference* oSRS = new OGRSpatialReference;
	if (is_ogr_error(oSRS->SetFromUserInput( dstcrs.c_str() ), msg)) {
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

	CPLFree(pszDstWKT);
 	delete oSRS;

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


bool gdal_warper(GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, std::vector<unsigned> srcbands, std::vector<unsigned> dstbands, std::string method, std::string srccrs, std::string msg, bool verbose) {

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
    psWarpOptions->panSrcBands = (int *) CPLMalloc(sizeof(int) * nbands );
    psWarpOptions->panDstBands = (int *) CPLMalloc(sizeof(int) * nbands );
	psWarpOptions->padfSrcNoDataReal = (double *) CPLMalloc(sizeof(double) * nbands );
	psWarpOptions->padfDstNoDataReal = (double *) CPLMalloc(sizeof(double) * nbands );
	psWarpOptions->padfSrcNoDataImag = (double *) CPLMalloc(sizeof(double) * nbands );
	psWarpOptions->padfDstNoDataImag = (double *) CPLMalloc(sizeof(double) * nbands );

	GDALRasterBandH hBand;
	int hasNA;
	for (int i=0; i<nbands; i++) {
		psWarpOptions->panSrcBands[i] = (int) srcbands[i]+1;
		psWarpOptions->panDstBands[i] = (int) dstbands[i]+1;

		hBand = GDALGetRasterBand(hSrcDS, srcbands[i]+1);
		double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
		if (verbose && i == 0) {
#ifdef useRcpp
			std::string hna = hasNA ? "true" : "false";
			Rcpp::Rcout << "hasNA         : " << hna << std::endl;
			Rcpp::Rcout << "NA flag       : " << naflag << std::endl;
#endif
		}
		if (hasNA) {
			psWarpOptions->padfSrcNoDataReal[i] = naflag;
			psWarpOptions->padfDstNoDataReal[i] = naflag;
			hBand = GDALGetRasterBand(hDstDS, dstbands[i]+1);
			GDALSetRasterNoDataValue(hBand, naflag);
		} else {
			psWarpOptions->padfSrcNoDataReal[i] = NAN;
			psWarpOptions->padfDstNoDataReal[i] = NAN;
		}
		psWarpOptions->padfSrcNoDataImag[i] = 0;
		psWarpOptions->padfDstNoDataImag[i] = 0;
    }

	//psWarpOptions->pfnProgress = GDALTermProgress;

	psWarpOptions->papszWarpOptions =
     CSLSetNameValue( psWarpOptions->papszWarpOptions, "INIT_DEST", "NO_DATA");
	psWarpOptions->papszWarpOptions =
      CSLSetNameValue( psWarpOptions->papszWarpOptions, "WRITE_FLUSH", "YES");


    psWarpOptions->pTransformerArg =
        GDALCreateGenImgProjTransformer( hSrcDS, srccrs.c_str(),
                                        hDstDS, GDALGetProjectionRef(hDstDS),
                                        FALSE, 0.0, 1 );
    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;

    GDALWarpOperation oOperation;
    if (oOperation.Initialize( psWarpOptions ) != CE_None) {
		return false;
	}
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
	out.setNames(getNames());
	if (method == "near") {
		out.source[0].hasColors = hasColors();
		out.source[0].cols = getColors();
		out.source[0].hasCategories = hasCategories();
		out.source[0].cats = getCategories();
		out.rgb = rgb;
		out.rgblyrs = rgblyrs;		
	}

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
	std::string filename = opt.get_filename();

	std::string srccrs = getSRS("wkt");
	//if (opt.verbose) {
	//	Rcpp::Rcout << "wkt" << std::endl;
	//	Rcpp::Rcout << srccrs << std::endl;
	//}
	if (srccrs == "") {
		out.setError("input raster CRS not set");
		return out;	
	}

	bool use_crs = crs != "";  
	// should not be needed (need to fix)

	if ((!use_crs) & (!hasValues())) {
		if (filename != "") {
			addWarning("raster has no values, not writing to file");
		}
		return out;
	}

	if (filename == "") {
		if (!canProcessInMemory(opt) || !out.canProcessInMemory(opt)) {
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

		if (!open_gdal(hSrcDS, i, opt)) {
			out.setError("cannot create dataset from source");
			return out;
		}

		// create dest source, only once 
		if (i==0) {
			 // use the crs, ignore argument "x"
			if (use_crs) {
				if (! find_output_bounds(hSrcDS, hDstDS, srccrs, crs, filename, driver, nlyr(), opt.get_datatype(), errmsg)) {
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
				if (!out.create_gdalDS(hDstDS, filename, driver, false, NAN, opt)) {
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

		bool success = gdal_warper(hSrcDS, hDstDS, srcbands, dstbands, method, srccrs, errmsg, opt.get_verbose());

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
			out.setError("cannot do this transformation");
			return out;
		}
	} else {
		std::vector<std::string> nms = getNames();
		for (size_t i=0; i < nlyr(); i++) {
			GDALRasterBandH hBand = GDALGetRasterBand(hDstDS, i+1);
			double adfMinMax[2];
			bool approx = ncell() > 10e+8;
			GDALComputeRasterMinMax(hBand, approx, adfMinMax);
			GDALSetRasterStatistics(hBand, adfMinMax[0], adfMinMax[1], NAN, NAN);
			GDALSetDescription(hBand, nms[i].c_str());
		}
		GDALClose( hDstDS );
		out = SpatRaster(filename, {-1}, {""});
	}

	if (mask) {
		SpatVector v = dense_extent();
		v = v.project(out.getSRS("wkt"));
		out = out.mask(v, false, NAN, true, mopt);
	}

	return out;
}




#endif



/*

void SpatRaster::resample2(SpatRaster &out, const std::string &method, SpatOptions &opt) {

	unsigned nc = out.ncol();
  	if (!out.writeStart(opt)) { return; }
	for (size_t i = 0; i < out.bs.n; i++) {
        unsigned firstcell = out.cellFromRowCol(out.bs.row[i], 0);
		unsigned lastcell  = out.cellFromRowCol(out.bs.row[i]+out.bs.nrows[i]-1, nc-1);
		std::vector<double> cells(1+lastcell-firstcell);
		std::iota (std::begin(cells), std::end(cells), firstcell);
        std::vector<std::vector<double>> xy = out.xyFromCell(cells);
		std::vector<std::vector<double>> v = extractXY(xy[0], xy[1], method);
		if (!out.writeValues2(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return;
	}
	out.writeStop();

}


SpatRaster SpatRaster::resample1(SpatRaster &x, const std::string &method, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);
	out.setNames(getNames());
	std::vector<std::string> f {"bilinear", "ngb"};
	if (std::find(f.begin(), f.end(), method) == f.end()) {
		out.setError("unknown resample method");
		return out;
	}
	if (!hasValues()) {
		return out;
	}

	if ((!source[0].srs.is_empty()) && (!out.source[0].srs.is_empty())) {
		if (!source[0].srs.is_equal(out.source[0].srs)) {
			out.addWarning("Rasters have different crs");
		}
	}

	unsigned xq = x.xres() / xres();
	unsigned yq = x.yres() / yres();
	if (std::max(xq, yq) > 1) {
		SpatRaster xx;
		xq = xq == 0 ? 1 : xq;
		yq = yq == 0 ? 1 : yq;
		std::vector<unsigned> agf = {yq, xq, 1};
		SpatOptions agopt;
		if (method == "bilinear") {
			xx = aggregate(agf, "mean", true, agopt);
		} else {
			xx = aggregate(agf, "modal", true, agopt);
		}
		xx.resample2(out, method, opt);
	} else {
		resample2(out, method, opt);
	}
	return out;
}

*/





SpatRaster SpatRaster::rectify(std::string method, SpatRaster aoi, unsigned useaoi, bool snap, SpatOptions &opt) {
	SpatRaster out = geometry(0);

	if (nsrc() > 1) {
		out.setError("you can transform only one data source at a time");
		return(out);
	}
	if (!source[0].rotated) {
		out.setError("this source is not rotated");
		return(out);
	} 
	GDALDataset *poDataset;
	std::string fname = source[0].filename;
	poDataset = (GDALDataset *) GDALOpen(fname.c_str(), GA_ReadOnly );
	if( poDataset == NULL )  {
		setError("cannot read from " + fname );
		return out;
	}
	double gt[6];
	if( poDataset->GetGeoTransform(gt) != CE_None ) {
		out.setError("can't get geotransform");
		GDALClose( (GDALDatasetH) poDataset );
		return out;
	}
	GDALClose( (GDALDatasetH) poDataset );
	//SpatExtent e = getExtent();
	//std::vector<double> x = {e.xmin, e.xmin, e.xmax, e.xmax };
	//std::vector<double> y = {e.ymin, e.ymax, e.ymin, e.ymax };
	double nc = ncol();
	double nr = nrow();
	std::vector<double> x = {0, 0, nc, nc};
	std::vector<double> y = {0, nr, 0, nr};
	std::vector<double> xx(4);
	std::vector<double> yy(4);
	for (size_t i=0; i<4; i++) {
		xx[i] = gt[0] + x[i]*gt[1] + y[i]*gt[2];
		yy[i] = gt[3] + x[i]*gt[4] + y[i]*gt[5];
	}
	double xmin = vmin(xx, TRUE);
	double xmax = vmax(xx, TRUE);
	double ymin = vmin(yy, TRUE);
	double ymax = vmax(yy, TRUE);
	SpatExtent en(xmin, xmax, ymin, ymax);
	out = out.setResolution(gt[1], -gt[5]);
	out.setExtent(en, false, "out");

	if (useaoi == 1) { // use extent
		en = aoi.getExtent();
		if (snap) {
			en = out.align(en, "near");
			out.setExtent(en, false, "near");
		} else {
			out.setExtent(en, false, "");
		}
	} else if (useaoi == 2){  // extent and resolution
		out = aoi.geometry(0);
	} // else { // if (useaoi == 0) // no aoi

	
	out = warper(out, "", method, false, opt);

	return(out);
}



SpatVector SpatRaster::polygonize(bool trunc, bool values, bool narm, bool aggregate, SpatOptions &opt) {

	SpatVector out;
	SpatOptions topt(opt);

	if (nlyr() > 1) {
		out.addWarning("only the first layer is polygonized when 'dissolve=TRUE'");
	}
	SpatRaster tmp = subset({0}, topt);

	bool usemask = false;
	SpatRaster mask;
	if (narm) {
		usemask = true;
		mask = tmp.isfinite(topt);	
	} else if (trunc) {
		tmp = tmp.math("trunc", topt);
		trunc = false;
	} else if (tmp.sources_from_file()) {
		// for NAN and INT files. Should have a check for that
		//tmp = tmp.arith(0, "+", false, topt);
		// riskier  
		tmp.readAll();
	}
	
	
	GDALDatasetH rstDS;
	if (! tmp.sources_from_file() ) {
		if (!tmp.open_gdal(rstDS, 0, topt)) {
			out.setError("cannot open dataset");
			return out;
		}
	} else {
		std::string filename = tmp.source[0].filename;
		rstDS = GDALOpen( filename.c_str(), GA_ReadOnly);
		if (rstDS == NULL) {
			out.setError("cannot open dataset from file");
			return out;
		}
	}
    GDALDataset *srcDS=NULL;

#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
	srcDS = (GDALDataset *) rstDS;
#else
	srcDS = srcDS->FromHandle(rstDS);
#endif

	GDALDataset *maskDS=NULL;
	if (usemask) {
		GDALDatasetH rstMask;
		if (! mask.sources_from_file() ) {
			if (!mask.open_gdal(rstMask, 0, opt)) {
				out.setError("cannot open dataset");
				return out;
			}
		} else {
			std::string filename = mask.source[0].filename;
			rstMask = GDALOpen( filename.c_str(), GA_ReadOnly);
			if (rstMask == NULL) {
				out.setError("cannot open dataset from file");
				return out;		
			}
		}
#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
		maskDS = (GDALDataset *) rstMask;
#else
		maskDS = srcDS->FromHandle(rstMask);
#endif
	}

    GDALDataset *poDS = NULL;
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( "Memory" );
    if( poDriver == NULL )  {
        out.setError( "cannot create output dataset");
        return out;
    }
    poDS = poDriver->Create("", 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        out.setError("Creation of dataset failed" );
        return out;
    }
	std::vector<std::string> nms = getNames();
	std::string name = nms[0];

	OGRSpatialReference *SRS = NULL;
	std::string s = source[0].srs.wkt;
	if (s != "") {
		SRS = new OGRSpatialReference;
		OGRErr err = SRS->SetFromUserInput(s.c_str()); 
		if (err != OGRERR_NONE) {
			out.setError("crs error");
			delete SRS;
			return out;
		}
	}

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer(name.c_str(), SRS, wkbPolygon, NULL );
    if( poLayer == NULL ) {
        out.setError( "Layer creation failed" );
        return out;
    }
	if (SRS != NULL) SRS->Release();

	OGRFieldDefn oField(name.c_str(), trunc ?  OFTInteger : OFTReal);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
		out.setError( "Creating field failed");
		return out;
	}

	GDALRasterBand  *poBand;
	poBand = srcDS->GetRasterBand(1);

	//int hasNA=1;
	//double naflag = poBand->GetNoDataValue(&hasNA);

	CPLErr err;
	if (usemask) {
		GDALRasterBand  *maskBand;
		maskBand = maskDS->GetRasterBand(1);
		if (trunc) {
			err = GDALPolygonize(poBand, maskBand, poLayer, 0, NULL, NULL, NULL);
		} else {
			err = GDALFPolygonize(poBand, maskBand, poLayer, 0, NULL, NULL, NULL);
		}
		GDALClose(maskDS);
	} else {
		if (trunc) {
			err = GDALPolygonize(poBand, NULL, poLayer, 0, NULL, NULL, NULL);
		} else {
			err = GDALFPolygonize(poBand, NULL, poLayer, 0, NULL, NULL, NULL);
		}
	}
	if (err == 4) {
		out.setError("polygonize error");
		return out;
	}
	GDALClose(srcDS);

	out.read_ogr(poDS);
	GDALClose(poDS);

	if (aggregate && (out.nrow() > 0)) {
		out = out.aggregate(name, false);
	}

	if (!values) {
		out.df = SpatDataFrame();
	} 
	return out;
}

	
	
SpatRaster SpatRaster::rgb2col(size_t r,  size_t g, size_t b, SpatOptions &opt) {	
	SpatRaster out = geometry(1);
	if (nlyr() < 3) {
		out.setError("need at least three layers");
		return out;
	}
	size_t mxlyr = std::max(std::max(r, g), b);
	if (nlyr() < mxlyr) {
		out.setError("layer number for R, G, B, cannot exceed nlyr()");		
		return out;
	}
	
	std::string filename = opt.get_filename();
	opt.set_datatype("INT1U");
	std::string driver;
	if (filename == "") {
		if (canProcessInMemory(opt)) {
			driver = "MEM";
		} else {
			filename = tempFile(opt.get_tempdir(), ".tif");
			opt.set_filenames({filename});
			driver = "GTiff";
		} 
	} else {
		driver = opt.get_filetype();
		getGDALdriver(filename, driver);
		if (driver == "") {
			setError("cannot guess file type from filename");
			return out;
		}
		std::string errmsg;
		if (!can_write(filename, opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}	
	}
	
	std::vector<unsigned> lyrs = {(unsigned)r, (unsigned)g, (unsigned)b};
	SpatOptions ops(opt);
	*this = subset(lyrs, ops);
	*this = collapse_sources();	
	GDALDatasetH hSrcDS, hDstDS;
	if (!open_gdal(hSrcDS, 0, ops)) {
		out.setError("cannot create dataset from source");
		return out;
	}
	GDALRasterBandH R = GDALGetRasterBand(hSrcDS,1);
	GDALRasterBandH G = GDALGetRasterBand(hSrcDS,1);
	GDALRasterBandH B = GDALGetRasterBand(hSrcDS,1);

	GDALColorTableH hColorTable= GDALCreateColorTable(GPI_RGB);

	if (GDALComputeMedianCutPCT(R, G, B, NULL, 256, hColorTable, NULL, NULL) != CE_None) {
		out.setError("cannot create color table");
		GDALClose(hSrcDS);
		return out;
	}	
	
	if (!out.create_gdalDS(hDstDS, filename, driver, true, 0, opt)) {
		out.setError("cannot create new dataset");
		GDALClose(hSrcDS);
		return out;
	}
	
	GDALRasterBandH hTarget = GDALGetRasterBand(hDstDS, 1);
	GDALSetRasterColorInterpretation(hTarget, GCI_PaletteIndex);
	if (GDALDitherRGB2PCT(R, G, B, hTarget, hColorTable, NULL, NULL) != CE_None) {
		out.setError("cannot set color table");
		GDALClose(hSrcDS);
		GDALClose(hDstDS);
		return out;		
	}
	GDALClose(hSrcDS);

	if (driver == "MEM") {
		if (!out.from_gdalMEM(hDstDS, false, true)) {
			out.setError("conversion failed (mem)");
			GDALClose(hDstDS);
			return out;
		}
		SpatDataFrame cdf;
		cdf.add_column(1, "red");
		cdf.add_column(1, "green");
		cdf.add_column(1, "blue");
		cdf.add_column(1, "alpha");
		size_t nc = GDALGetColorEntryCount(hColorTable);
		cdf.reserve(nc);

		for (size_t i=0; i<nc; i++) {	
			const GDALColorEntry* col = GDALGetColorEntry(hColorTable, i);
			cdf.iv[0].push_back(col->c1);
			cdf.iv[1].push_back(col->c2);
			cdf.iv[2].push_back(col->c3);
			cdf.iv[3].push_back(col->c4);
		}
		out.source[0].hasColors.resize(1);
		out.source[0].hasColors[0] = true;
		out.source[0].cols.resize(1);
		out.source[0].cols[0] = cdf;		
	} else {
		if (GDALSetRasterColorTable(hTarget, hColorTable) != CE_None) {
			out.setError("cannot set color table");
			GDALClose(hDstDS);
			return out;
		}
	}
	GDALClose(hDstDS);
	if (driver != "MEM") {
		out = SpatRaster(filename, {-1}, {""});
	}
	return out;
}

	
SpatRaster SpatRaster::sievefilter(int threshold, int connections, SpatOptions &opt) {	
	SpatRaster out;
	std::string filename = opt.get_filename();
	std::string driver;
	if (filename == "") {
		if (canProcessInMemory(opt)) {
			driver = "MEM";
		} else {
			filename = tempFile(opt.get_tempdir(), ".tif");
			opt.set_filenames({filename});
			driver = "GTiff";
		} 
	} else {
		driver = opt.get_filetype();
		getGDALdriver(filename, driver);
		if (driver == "") {
			setError("cannot guess file type from filename");
			return out;
		}
		std::string errmsg;
		if (!can_write(filename, opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}	
	}

	SpatOptions ops(opt);
	GDALDatasetH hSrcDS, hDstDS;
	if (!open_gdal(hSrcDS, 0, ops)) {
		out.setError("cannot open input dataset");
		return out;
	}
	
	GDALDriverH hDriver = GDALGetDriverByName( driver.c_str() );
	if ( hDriver == NULL ) {
		out.setError("empty driver");
		return out;
	}
	if (!out.create_gdalDS(hDstDS, filename, driver, true, 0, opt)) {
		out.setError("cannot create new dataset");
		GDALClose(hSrcDS);
		return out;
	}

	GDALRasterBandH hSrcBand = GDALGetRasterBand(hSrcDS, 1);
	GDALRasterBandH hTargetBand = GDALGetRasterBand(hDstDS, 1);

	if (!GDALSieveFilter(hSrcBand, NULL, hTargetBand, threshold, connections, NULL, NULL, NULL)) {
		out.setError("sieve failed");
		GDALClose(hSrcDS);
		GDALClose(hDstDS);
		return out;
	}	
	
	GDALClose(hSrcDS);
	if (driver == "MEM") {
		if (!out.from_gdalMEM(hDstDS, false, true)) {
			out.setError("conversion failed (mem)");
			GDALClose(hDstDS);
			return out;
		}
	} else {
		out = SpatRaster(filename, {-1}, {""});
	}
	GDALClose(hDstDS);
	return out;
}


/*	
SpatRaster SpatRaster::fillna(int threshold, int connections, SpatOptions &opt) {	
CPLErr GDALFillNodata(GDALRasterBandH hTargetBand, GDALRasterBandH hMaskBand, doubledfMaxSearchDist, intbDeprecatedOption, intnSmoothingIterations, char**papszOptions, GDALProgressFuncpfnProgress, void*pProgressArg)
*/


