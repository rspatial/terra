// Copyright (c) 2018-2022  Robert J. Hijmans
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


SpatVector SpatRaster::dense_extent(bool inside, bool geobounds) {

	SpatExtent e = getExtent();
	if (geobounds && is_lonlat()) {		
		if ((e.ymin <= -90) || (e.ymax >= 90)) { 
			double fy = yres() / 10; // avoid Inf with Mercator
			SpatRaster g = geometry();
			e.ymin= std::max(e.ymin, -90.0+fy);
			e.ymax= std::min(e.ymax, 90.0-fy);
			g.source[0].extent = e;
			return g.dense_extent(inside, false);
		}
	}

	std::vector<int_64> rows, cols;
	if (nrow() < 51) {
		rows.resize(nrow());
		std::iota(rows.begin(), rows.end(), 0);
	} else {
		rows = seq_steps((int_64) 0, (int_64) nrow()-1, 50);
	}
	if (ncol() < 51) {
		cols.resize(nrow());
		std::iota(cols.begin(), cols.end(), 0);
	} else {
		cols = seq_steps((int_64) 0, (int_64) ncol()-1, 50);
	}


	std::vector<double> xcol = xFromCol(cols) ;
	std::vector<double> yrow = yFromRow(rows) ;

	double yr = yres() / 4;
	if (inside) {
		yrow.insert(yrow.begin(), e.ymax - yr);
		yrow.push_back(e.ymin + yr);
		std::vector<double> y0(xcol.size(), e.ymin+yr);
		std::vector<double> y1(xcol.size(), e.ymax-yr);
	} else {
		yrow.insert(yrow.begin(), e.ymax);
		yrow.push_back(e.ymin);
		std::vector<double> y0(xcol.size(), e.ymin);
		std::vector<double> y1(xcol.size(), e.ymax);
	}

	std::vector<double> y0(xcol.size(), e.ymin);
	std::vector<double> y1(xcol.size(), e.ymax);
	std::vector<double> x0(yrow.size(), e.xmin);
	std::vector<double> x1(yrow.size(), e.xmax);
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

	x.push_back(x[0]);
	y.push_back(y[0]);

	SpatVector v(x, y, polygons, getSRS("wkt"));

	return v;
}

#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR < 2

SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, bool mask, bool align, SpatOptions &opt) {
	SpatRaster out;
	out.setError("Not supported for this old version of GDAL");
	return(out);
}


#else



bool get_output_bounds(const GDALDatasetH &hSrcDS, std::string srccrs, const std::string dstcrs, SpatRaster &r) {

	if ( hSrcDS == NULL ) {
		r.setError("data source is NULL");
		return false;
	}

	// Get Source coordinate system.
	// const char *pszSrcWKT = GDALGetProjectionRef( hSrcDS );
	const char *pszSrcWKT = srccrs.c_str();
	if ( pszSrcWKT == NULL || strlen(pszSrcWKT) == 0 ) {
		r.setError("data source has no WKT");
		return false;
	}

	OGRSpatialReference* oSRS = new OGRSpatialReference;
	std::string msg = "";
	if (is_ogr_error(oSRS->SetFromUserInput( dstcrs.c_str() ), msg)) {
		r.setError(msg);
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
		r.setError("cannot create TranformArg");
		return false;
	}
	CPLFree(pszDstWKT);
 	delete oSRS;

	double adfDstGeoTransform[6];
	int nPixels=0, nLines=0;
	CPLErr eErr = GDALSuggestedWarpOutput( hSrcDS, GDALGenImgProjTransform,
					hTransformArg, adfDstGeoTransform, &nPixels, &nLines );

	GDALDestroyGenImgProjTransformer( hTransformArg );
	if ( eErr != CE_None ) {
		r.setError("cannot create warp output");
		return false;
	}


	r.source[0].ncol = nPixels;
	r.source[0].nrow = nLines;

	r.source[0].extent.xmin = adfDstGeoTransform[0]; /* left x */
	/* w-e pixel resolution */
	r.source[0].extent.xmax = r.source[0].extent.xmin + adfDstGeoTransform[1] * nPixels;
	r.source[0].extent.ymax = adfDstGeoTransform[3]; // top y
	r.source[0].extent.ymin = r.source[0].extent.ymax + nLines * adfDstGeoTransform[5];

	r.setSRS({dstcrs});

	return true;
}

/*
	// Create output with same datatype as first input band.
	GDALDataType eDT = GDALGetRasterDataType(GDALGetRasterBand(hSrcDS,1));
	GDALDataType eDT;
	getGDALDataType(datatype, eDT);

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
*/

bool getAlgo(GDALResampleAlg &alg, std::string m) {

	if (m=="sum") {
#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1
		alg = GRA_Sum;
		return true;
	}
#else
		return false;
	}
#endif

	if (m=="rms") {
#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 3
		alg = GRA_RMS;
		return true;
	}
#else
		return false;
	}
#endif


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
	} else if (m=="average") {
		alg = GRA_Average;
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
		return false;
	}
	return true;
}


bool is_valid_warp_method(const std::string &method) {
	std::vector<std::string> m { "near", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode", "max", "min", "med", "q1", "q3", "sum", "rms"};
	return (std::find(m.begin(), m.end(), method) != m.end());
}


bool set_warp_options(GDALWarpOptions *psWarpOptions, GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, std::vector<unsigned> srcbands, std::vector<unsigned> dstbands, std::string method, std::string srccrs, std::string msg, bool verbose, bool threads) {

	if (srcbands.size() != dstbands.size()) {
		msg = "number of source bands must match number of dest bands";
		return false;
	}
	int nbands = srcbands.size();

	GDALResampleAlg a;
	if (!getAlgo(a, method)) {
		if (method=="sum" || method=="rms") {
			msg = method + " not available in your version of GDAL";
		} else {
			msg = "unknown resampling algorithm";
		}
		return false;
	}

    // Setup warp options.
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

	if (threads) {
		psWarpOptions->papszWarpOptions =
			CSLSetNameValue( psWarpOptions->papszWarpOptions, "NUM_THREADS", "ALL_CPUS");
	}

    psWarpOptions->pTransformerArg =
        GDALCreateGenImgProjTransformer( hSrcDS, srccrs.c_str(),
                                        hDstDS, GDALGetProjectionRef(hDstDS),
                                        FALSE, 0.0, 1 );
    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;

	return true;
}


bool gdal_warper(GDALWarpOptions *psWarpOptions, GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS) {
    GDALWarpOperation oOperation;
    if (oOperation.Initialize( psWarpOptions ) != CE_None) {
		return false;
	}
    if (oOperation.ChunkAndWarpImage(0, 0, GDALGetRasterXSize(hDstDS), GDALGetRasterYSize(hDstDS)) != CE_None) {
		return false;
	}
    GDALDestroyGenImgProjTransformer( psWarpOptions->pTransformerArg );
    GDALDestroyWarpOptions( psWarpOptions );
	return true;
}


SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, bool mask, bool align, bool resample, SpatOptions &opt) {

	SpatRaster out = x.geometry(nlyr(), false, false);
	if (!is_valid_warp_method(method)) {
		out.setError("not a valid warp method");
		return out;
	}
	std::string srccrs = getSRS("wkt");
	if (resample) {
		out.setSRS(srccrs);
	}

	out.setNames(getNames());
	if (method == "near") {
		out.source[0].hasColors = hasColors();
		out.source[0].cols = getColors();
		out.source[0].hasCategories = hasCategories();
		out.source[0].cats = getCategories();
		out.rgb = rgb;
		out.rgblyrs = rgblyrs;
		out.rgbtype = rgbtype;
	}
	if (hasTime()) {
		out.source[0].hasTime = true;
		out.source[0].timestep = getTimeStep();
		out.source[0].timezone = getTimeZone();
		out.source[0].time = getTime();
	}

	bool use_crs = crs != "";
	if (use_crs) {
		align = false;
		resample = false;
	} else if (!hasValues()) {
		std::string fname = opt.get_filename();
		if (fname != "") {
			out.addWarning("raster has no values, not writing to file");
		}
		return out;
	}
	if (align) {
		crs = out.getSRS("wkt");
	}

	if (!resample) {
		if (srccrs == "") {
			out.setError("input raster CRS not set");
			return out;
		}
	}

	lrtrim(crs);
	SpatOptions sopt(opt);
	if (use_crs || align) {
		GDALDatasetH hSrcDS;
		if (!open_gdal(hSrcDS, 0, false, sopt)) {
			out.setError("cannot create dataset from source");
			return out;
		}
		out.setSRS(crs);
		if (!get_output_bounds(hSrcDS, srccrs, crs, out)) {
			GDALClose( hSrcDS );
			out.setError("cannot get output boundaries");
			return out;
		}
		GDALClose( hSrcDS );
	} else if (!resample) {
		OGRSpatialReference source, target;
		const char *pszDefFrom = srccrs.c_str();
		OGRErr erro = source.SetFromUserInput(pszDefFrom);
		if (erro != OGRERR_NONE) {
			out.setError("input crs is not valid");
			return out;
		}
		std::string targetcrs = out.getSRS("wkt");
		const char *pszDefTo = targetcrs.c_str();
		erro = target.SetFromUserInput(pszDefTo);
		if (erro != OGRERR_NONE) {
			out.setError("output crs is not valid");
			return out;
		}
		OGRCoordinateTransformation *poCT;
		poCT = OGRCreateCoordinateTransformation(&source, &target);
		if( poCT == NULL )	{
			out.setError( "Cannot do this transformation" );
			return(out);
		}
		OCTDestroyCoordinateTransformation(poCT);
	}

	if (align) {
		SpatExtent e = out.getExtent();
		e = x.align(e, "out");
		out.setExtent(e, false, true, "");
		std::vector<double> res = x.resolution();
		out = out.setResolution(res[0], res[1]);
	}
	if (!hasValues()) {
		return out;
	}

	SpatOptions mopt;
	if (mask) {
		mopt = opt;
		opt = SpatOptions(opt);
	}

	opt.ncopies += 4;
	if (!out.writeStart(opt)) {
		return out;
	}

	std::string errmsg;
	size_t ns = nsrc();
	SpatExtent eout = out.getExtent();


	std::vector<bool> has_so = source[0].has_scale_offset;
	std::vector<double> scale = source[0].scale;
	std::vector<double> offset = source[0].offset;
	for (size_t i=1; i<ns; i++) {
		has_so.insert(has_so.end(), source[0].has_scale_offset.begin(), source[0].has_scale_offset.end());
		scale.insert(scale.end(), source[0].scale.begin(), source[0].scale.end());
		offset.insert(offset.end(), source[0].offset.begin(), source[0].offset.end());
	}

	double halfy = out.yres() / 2;
	for (size_t i = 0; i < out.bs.n; i++) {
		int bandstart = 0;
		eout.ymax = out.yFromRow(out.bs.row[i]) + halfy;
		eout.ymin = out.yFromRow(out.bs.row[i] + out.bs.nrows[i]-1) - halfy;
		SpatRaster crop_out = out.crop(eout, "near", sopt);
		GDALDatasetH hDstDS;

		if (!crop_out.create_gdalDS(hDstDS, "", "MEM", false, NAN, has_so, scale, offset, sopt)) {
			return crop_out;
		}

		for (size_t j=0; j<ns; j++) {

			GDALDatasetH hSrcDS;

			if (!open_gdal(hSrcDS, j, false, sopt)) {
				out.setError("cannot create dataset from source");
				if( hDstDS != NULL ) GDALClose( (GDALDatasetH) hDstDS );
				return out;
			}
			std::vector<unsigned> srcbands = source[j].layers;
			std::vector<unsigned> dstbands(srcbands.size());
			std::iota (dstbands.begin(), dstbands.end(), bandstart);
			bandstart += dstbands.size();


			GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
			bool ok = set_warp_options(psWarpOptions, hSrcDS, hDstDS, srcbands, dstbands, method, srccrs, errmsg, opt.get_verbose(), opt.threads);
			if (!ok) {
				if( hDstDS != NULL ) GDALClose( (GDALDatasetH) hDstDS );
				out.setError(errmsg);
				return out;
			}

			//bool success = gdal_warper(hSrcDS, hDstDS, srcbands, dstbands, method, srccrs, errmsg, opt.get_verbose(), opt.threads);
			ok = gdal_warper(psWarpOptions, hSrcDS, hDstDS);
			if( hSrcDS != NULL ) GDALClose( (GDALDatasetH) hSrcDS );
			if (!ok) {
				if( hDstDS != NULL ) GDALClose( (GDALDatasetH) hDstDS );
				out.setError("warp failure");
				return out;
			}
		}

		bool ok = crop_out.from_gdalMEM(hDstDS, false, true);
		if( hDstDS != NULL ) GDALClose( (GDALDatasetH) hDstDS );
		if (!ok) {
			out.setError("cannot do this transformation (warp)");
			return out;
		}
		std::vector<double> v = crop_out.getValues(-1, opt);
		if (!out.writeBlock(v, i)) return out;
	}
	out.writeStop();


	if (mask) {
		SpatVector v = dense_extent(true, true);
		v = v.project(out.getSRS("wkt"));
		if (v.nrow() > 0) {
			out = out.mask(v, false, NAN, true, mopt);
		} else {
			out.addWarning("masking failed");
		}
	}
	return out;
}



#endif



SpatRaster SpatRaster::resample(SpatRaster x, std::string method, bool mask, bool agg, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);
	out.setNames(getNames());

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
		e = e.intersect(getExtent());
		if (!e.valid()) {
			out.addWarning("No spatial overlap");
			return out;
		}
	}

	if (agg) {
		if (do_prj) {
			// compare changes in true cell areas
			// if (some) output cells are much larger than input, we could
			//    a) disaggregate "x", warp, and aggregate the results
			//    b) or give a warning?
		} else {
			unsigned xq = x.xres() / xres();
			unsigned yq = x.yres() / yres();
			if (std::max(xq, yq) > 1) {
				xq = xq == 0 ? 1 : xq;
				yq = yq == 0 ? 1 : yq;
				std::vector<unsigned> agf = {yq, xq, 1};
				SpatOptions agopt(opt);
				SpatRaster xx;
				if (method == "bilinear") {
					xx = aggregate(agf, "mean", true, agopt);
				} else {
					xx = aggregate(agf, "modal", true, agopt);
				}
				return xx.resample(x, method, mask, false, opt);
			}
		}
	}

	SpatOptions mopt;
	if (mask) {
		mopt = opt;
		opt = SpatOptions(opt);
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
		std::vector<std::vector<double>> e = extractXY(xy[0], xy[1], method, false);
		std::vector<double> v = flatten(e);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i])) return out;
	}
	out.writeStop();


	if (mask) {
		SpatVector v = dense_extent(true, true);
		v = v.project(out.getSRS("wkt"));
		if (v.nrow() > 0) {
			out = out.mask(v, false, NAN, true, mopt);
		} else {
			out.addWarning("masking failed");
		}
	}


	return(out);
}




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
	GDALDataset *poDataset = openGDAL(source[0].filename, GDAL_OF_RASTER | GDAL_OF_READONLY, source[0].open_drivers, source[0].open_ops);

	if( poDataset == NULL )  {
		setError("cannot read from " + source[0].filename);
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

	out.setExtent(en, false, true, "out");
	SpatExtent e = out.getExtent();

	if (useaoi == 1) { // use extent
		en = aoi.getExtent();
		if (snap) {
			en = out.align(en, "near");
			out.setExtent(en, false, true, "near");
		} else {
			out.setExtent(en, false, true, "");
		}
	} else if (useaoi == 2){  // extent and resolution
		out = aoi.geometry(0);
	} // else { // if (useaoi == 0) // no aoi

	e = out.getExtent();
	out = warper(out, "", method, false, false, true, opt);

	return(out);
}



SpatVector SpatRaster::polygonize(bool trunc, bool values, bool narm, bool aggregate, SpatOptions &opt) {

	SpatVector out;
	out.srs = source[0].srs;
	SpatOptions topt(opt);

	SpatRaster tmp;
	if (nlyr() > 1) {
		out.addWarning("only the first layer is polygonized when 'dissolve=TRUE'");
		tmp = subset({0}, topt);
	} else {
		tmp = *this;
	}

//	bool usemask = false;
	SpatRaster mask;
	if (narm) {
//		usemask = true;
		SpatOptions mopt(topt);
		mopt.set_datatype("INT1U");
		mask = tmp.isfinite(mopt);
	} else if (trunc) {
		tmp = tmp.math("trunc", topt);
		trunc = false;
	} else if (tmp.sources_from_file()) {
		// for NAN and INT files. Should have a check for that
		//tmp = tmp.arith(0, "+", false, topt);
		// riskier
		tmp.readAll();
	}
	if (tmp.source[0].extset) {
		tmp = tmp.hardCopy(topt);
	}

	GDALDatasetH rstDS;
	if (!tmp.open_gdal(rstDS, 0, false, topt)) {
		out.setError("cannot open dataset");
		return out;
	}

    GDALDataset *srcDS=NULL;

#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
	srcDS = (GDALDataset *) rstDS;
#else
	srcDS = srcDS->FromHandle(rstDS);
#endif

	GDALDataset *maskDS=NULL;
	GDALDatasetH rstMask;
	if (narm) {
		if (!mask.open_gdal(rstMask, 0, false, opt)) {
			out.setError("cannot open dataset");
			return out;
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
        out.setError( "cannot create output driver");
        return out;
    }
    poDS = poDriver->Create("", 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        out.setError("Creation of output dataset failed" );
        return out;
    }
	std::vector<std::string> nms = getNames();
	std::string name = nms[0];

	OGRSpatialReference *SRS = NULL;

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
	if (narm) {
		GDALRasterBand *maskBand;
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

	std::vector<double> fext;
	SpatVector fvct;
	out.read_ogr(poDS, "", "", fext, fvct, false);
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
			filename = tempFile(opt.get_tempdir(), opt.pid, ".tif");
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
	if (!open_gdal(hSrcDS, 0, false, ops)) {
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

	if (!out.create_gdalDS(hDstDS, filename, driver, true, 0, {false}, {0.0}, {1.0}, opt)) {
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
		out = SpatRaster(filename, {-1}, {""}, {}, {});
	}
	return out;
}


SpatRaster SpatRaster::sieveFilter(int threshold, int connections, SpatOptions &opt) {
	SpatRaster out = geometry(1, true, true, true);
	if ((connections != 4) && (connections != 8)) {
		out.setError("connections should be 4 or 8");
		return out;
	}
	if (threshold < 1) {
		out.setError("threshold should be > 0");
		return out;
	}

	std::string filename = opt.get_filename();
	std::string driver;
	if (filename == "") {
		if (canProcessInMemory(opt)) {
			driver = "MEM";
		} else {
			filename = tempFile(opt.get_tempdir(), opt.pid, ".tif");
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
	if (!open_gdal(hSrcDS, 0, false, ops)) {
		out.setError("cannot open input dataset");
		return out;
	}

	GDALDriverH hDriver = GDALGetDriverByName( driver.c_str() );
	if ( hDriver == NULL ) {
		out.setError("empty driver");
		return out;
	}
	if (!out.create_gdalDS(hDstDS, filename, driver, true, 0, source[0].has_scale_offset, source[0].scale, source[0].offset, opt)) {
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
		out = SpatRaster(filename, {-1}, {""}, {}, {});
	}
	GDALClose(hDstDS);
	return out;
}


/*
SpatRaster SpatRaster::fillna(int threshold, int connections, SpatOptions &opt) {
CPLErr GDALFillNodata(GDALRasterBandH hTargetBand, GDALRasterBandH hMaskBand, doubledfMaxSearchDist, intbDeprecatedOption, intnSmoothingIterations, char**papszOptions, GDALProgressFuncpfnProgress, void*pProgressArg)
*/

/*
#include <gdalpansharpen.h>
SpatRaster SpatRaster::panSharpen(SpatRaster pan, SpatOptions &opt) {
	SpatRaster out = geometry();
	return out;
}
*/




/*
bool old_gdal_warper(GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS, std::vector<unsigned> srcbands, std::vector<unsigned> dstbands, std::string method, std::string srccrs, std::string msg, bool verbose, bool threads) {

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

	if (threads) {
		psWarpOptions->papszWarpOptions =
			CSLSetNameValue( psWarpOptions->papszWarpOptions, "NUM_THREADS", "ALL_CPUS");
	}

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
*/




/*
SpatRaster SpatRaster::very_old_warper(SpatRaster x, std::string crs, std::string method, bool mask, SpatOptions &opt) {

	SpatRaster out = x.geometry(nlyr(), false, false);
	out.setNames(getNames());
	if (method == "near") {
		out.source[0].hasColors = hasColors();
		out.source[0].cols = getColors();
		out.source[0].hasCategories = hasCategories();
		out.source[0].cats = getCategories();
		out.rgb = rgb;
		out.rgblyrs = rgblyrs;
	}
	if (hasTime()) {
		out.source[0].hasTime = true;
		out.source[0].timestep = getTimeStep();
		out.source[0].time = getTime();
	}


	bool use_crs = crs != "";
	std::string filename = opt.get_filename();
	if ((!use_crs) & (!hasValues())) {
		if (filename != "") {
			out.addWarning("raster has no values, not writing to file");
		}
		return out;
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
	filename = opt.get_filename();

	std::string srccrs = getSRS("wkt");
	if (srccrs == "") {
		out.setError("input raster CRS not set");
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
			if (use_crs) {
				if (!get_output_bounds(hSrcDS, srccrs, crs, out)) {
					GDALClose( hSrcDS );
					return out;
				}
			}
			if (!hasValues()) {
				GDALClose( hSrcDS );
				return out;
			}
			if (!out.create_gdalDS(hDstDS, filename, driver, false, NAN, source[i].has_scale_offset, source[i].scale, source[i].offset, opt)) {
				GDALClose( hSrcDS );
				//GDALClose( hDstDS );
				return out;
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
			out.setError("cannot do this transformation (warp)");
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
*/

