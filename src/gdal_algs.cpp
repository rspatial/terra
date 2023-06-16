// Copyright (c) 2018-2023  Robert J. Hijmans
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
#include <sstream>

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
		cols.resize(ncol());
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
#if GDAL_VERSION_MAJOR >= 3
	const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
	oSRS->exportToWkt( &pszDstWKT, options);
#else
	oSRS->exportToWkt( &pszDstWKT );
#endif

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
	} else if (m=="med") {
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
		if ((method=="sum") || (method=="rms")) {
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

/*
bool gdal_warper(GDALWarpOptions *psWarpOptions, GDALDatasetH &hSrcDS, GDALDatasetH &hDstDS) {
    GDALWarpOperation oOperation;
    if (oOperation.Initialize( psWarpOptions ) != CE_None) {
		return false;
	}
    if (oOperation.ChunkAndWarpImage(0, 0, GDALGetRasterXSize(hDstDS), GDALGetRasterYSize(hDstDS)) != CE_None) {
		return false;
	}
	return true;
}
*/


SpatRaster SpatRaster::warper(SpatRaster x, std::string crs, std::string method, bool mask, bool align, bool resample, SpatOptions &opt) {

/*
	if (extset) {
		std::vector<bool> m = inMemory();
		if (!vall(m)) {
			std::string fname = tempFile(opt.get_tempdir(), opt.pid, "_temp_rasterize.tif");
			SpatOptions xopt(opt);
			xopt.set_filenames({fname});
			SpatRaster y = x.writeRaster(opt);
			return warper(y, crs, method, mask, align, resample, opt);
		}
	}
*/
	size_t ns = nsrc();
	bool fixext = false;
	for (size_t j=0; j<ns; j++) {
		if (source[j].extset && (!source[j].memory)) {
			fixext = true;
			break;
		}
	}
	if (fixext) {
		SpatRaster r = *this;
		for (size_t j=0; j<ns; j++) {
			if (r.source[j].extset && (!r.source[j].memory)) {
				SpatRaster tmp(source[j]);
				//if (tmp.canProcessInMemory(opt)) {
				//	tmp.readAll();
				//} else {
				tmp = tmp.writeTempRaster(opt);
				r.source[j] = tmp.source[0]; 
			}
		}
		return r.warper(x, crs, method, mask, align, resample, opt);
	}


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

	bool use_crs = !crs.empty();
	if (use_crs) {
		align = false;
		resample = false;
	} else if (!hasValues()) {
		std::string fname = opt.get_filename();
		if (!fname.empty()) {
			out.addWarning("raster has no values, not writing to file");
		}
		return out;
	}
	if (align) {
		crs = out.getSRS("wkt");
	}

	if (!resample) {
		if (srccrs.empty()) {
			out.setError("input raster CRS not set");
			return out;
		}
	}

	lrtrim(crs);
	SpatOptions sopt(opt);
	if (use_crs || align) {
		GDALDatasetH hSrcDS;
		SpatRaster g = geometry(1);
		if (!g.open_gdal(hSrcDS, 0, false, sopt)) {
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
	if (!out.writeStart(opt, filenames())) {
		return out;
	}

	std::string errmsg;
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
		SpatRaster crop_out = out.crop(eout, "near", false, sopt);
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
			//ok = gdal_warper(psWarpOptions, hSrcDS, hDstDS);
			GDALWarpOperation oOperation;
			if (oOperation.Initialize( psWarpOptions ) != CE_None) {
				ok = false;
			} else if (oOperation.ChunkAndWarpImage(0, 0, GDALGetRasterXSize(hDstDS), GDALGetRasterYSize(hDstDS)) != CE_None) {
				ok = false;
			}
			GDALDestroyGenImgProjTransformer( psWarpOptions->pTransformerArg );
			GDALDestroyWarpOptions( psWarpOptions );

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
		v = v.project(out.getSRS("wkt"), true);
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
	if ((crsin == crsout) || crsin.empty() || crsout.empty()) {
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
  	if (!out.writeStart(opt, filenames())) { return out; }
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
		v = v.project(out.getSRS("wkt"), true);
		if (v.nrow() > 0) {
			out = out.mask(v, false, NAN, true, mopt);
		} else {
			out.addWarning("masking failed");
		}
	}


	return(out);
}



bool GCP_geotrans(GDALDataset *poDataset, double* adfGeoTransform) {
	int n = poDataset->GetGCPCount();
	if (n == 0) return false;
	const GDAL_GCP *gcp;
	gcp	= poDataset->GetGCPs();
	return GDALGCPsToGeoTransform(n, gcp, adfGeoTransform, true);
}

//#include <filesystem>

SpatRaster SpatRaster::rectify(std::string method, SpatRaster aoi, unsigned useaoi, bool snap, SpatOptions &opt) {
	SpatRaster out = geometry(0);

	if (nsrc() > 1) {
		out.setError("you can rectify only one data source at a time");
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
	
		if (GCP_geotrans(poDataset, gt)) {
		
			GDALClose( (GDALDatasetH) poDataset );
			std::string tmpfile = tempFile(opt.get_tempdir(), opt.pid, "_rect.tif");
			//++17 
			//std::filesystem::copy_file(source[0].filename, tmpfile);	

			std::ifstream  src(source[0].filename, std::ios::binary);
			std::ofstream  dst(tmpfile,   std::ios::binary);
			dst << src.rdbuf();	
		
			GDALDataset *poDataset = openGDAL(tmpfile, GDAL_OF_RASTER | GDAL_OF_READONLY, source[0].open_drivers, source[0].open_ops);
			poDataset->SetGeoTransform(gt);
			GDALClose( (GDALDatasetH) poDataset );
			SpatRaster tmp(tmpfile,  {-1}, {""}, {}, {});
			return tmp.rectify(method, aoi, useaoi, snap, opt);

		} else {
			out.setError("can't get the geotransform");
			GDALClose( (GDALDatasetH) poDataset );
			return out;
		}
	}
	GDALClose( (GDALDatasetH) poDataset );

//	gt[1] = std::abs(gt[1]);
	
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
	out = out.setResolution(fabs(gt[1]), fabs(gt[5]));

	out.setExtent(en, true, true, "out");
	//SpatExtent e = out.getExtent();

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

	//e = out.getExtent();

	out = warper(out, "", method, false, false, true, opt);

	return(out);
}



SpatVector SpatRaster::polygonize(bool round, bool values, bool narm, bool aggregate, int digits, SpatOptions &opt) {

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
		mask = tmp.isfinite(false, mopt);
	} 

		
	if (round && (digits > 0)) {
		tmp = tmp.math2("round", digits, topt);
		round = false;
	}
/*
	} else if (tmp.sources_from_file()) {
		// for NAN and INT files. Should have a check for that
		//tmp = tmp.arith(0, "+", false, topt);
		// riskier
		tmp.readAll();
	}
*/
	
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

	OGRFieldDefn oField(name.c_str(), round ? OFTInteger : OFTReal);
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
		if (round) {
			err = GDALPolygonize(poBand, maskBand, poLayer, 0, NULL, NULL, NULL);
		} else {
			err = GDALFPolygonize(poBand, maskBand, poLayer, 0, NULL, NULL, NULL);
		}
		GDALClose(maskDS);
	} else {
		if (round) {
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
	out.read_ogr(poDS, "", "", fext, fvct, false, "");
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
	if (filename.empty()) {
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
		if (driver.empty()) {
			out.setError("cannot guess file type from filename");
			return out;
		}

		std::string errmsg;
		if (!can_write({filename}, filenames(), opt.get_overwrite(), errmsg)) {
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





#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

SpatRaster SpatRaster::viewshed(const std::vector<double> obs, const std::vector<double> vals, const double curvcoef, const int mode, const double maxdist, const int heightmode, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (could_be_lonlat()) {
		out.setError("the method does not support lon/lat data");
		return out;
	}

	if (!hasValues()) {
		out.setError("input raster has no values");
		return out;
	}

	GDALViewshedOutputType outmode;
	if (heightmode==1) {
		outmode = GVOT_NORMAL; //= heightmode;
	} else if (heightmode==2) {
		outmode = GVOT_MIN_TARGET_HEIGHT_FROM_DEM;
	} else if (heightmode==3) {
		outmode = GVOT_MIN_TARGET_HEIGHT_FROM_GROUND;
	} else {
		out.setError("invalid output type");
		return out;
	}


	GDALViewshedMode emode;
	if (mode==1) {
		emode = GVM_Diagonal;
	} else if (mode==2) {
		emode = GVM_Edge;
	} else if (mode==3) {
		emode = GVM_Max;
	} else if (mode==4) {
		emode = GVM_Min;
	} else {
		out.setError("invalid mode");
		return out;
	}

	double minval = -9999;
	if (source[0].hasRange[0]) {
		minval = source[0].range_min[0] - 9999;
	}
	SpatOptions topt(opt);

	SpatRaster x;
	if (nlyr() > 1) {
		out.addWarning("viewshed is only done for the first layer");
		x = subset({0}, topt);
		x = x.replaceValues({NAN}, {minval}, 0, false, NAN, false, topt);
	} else {
		x = replaceValues({NAN}, {minval}, 0, false, NAN, false, topt);
	}

	std::string fname = opt.get_filename();
	std::string driver;
	if (!fname.empty()) {
		driver = opt.get_filetype();
		getGDALdriver(fname, driver);
		if (driver.empty()) {
			setError("cannot guess file type from filename");
			return out;
		}
		std::string errmsg;
		if (!can_write({fname}, filenames(), opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}
	}

	std::string filename = tempFile(topt.get_tempdir(), topt.pid, ".tif");
	driver = "GTiff";

	GDALDatasetH hSrcDS;
	if (!x.open_gdal(hSrcDS, 0, false, topt)) {
		out.setError("cannot open input dataset");
		return out;
	}

	GDALDriverH hDriver = GDALGetDriverByName( driver.c_str() );
	if ( hDriver == NULL ) {
		out.setError("empty driver");
		return out;
	}

	GIntBig diskNeeded = ncell() * 4;
	char **papszOptions = set_GDAL_options(driver, diskNeeded, false, topt.gdal_options);

	GDALRasterBandH hSrcBand = GDALGetRasterBand(hSrcDS, 1);

	GDALDatasetH hDstDS = GDALViewshedGenerate(hSrcBand, driver.c_str(), filename.c_str(), papszOptions, obs[0], obs[1], obs[2], obs[3], vals[0], vals[1], vals[2], vals[3], curvcoef, emode, maxdist, NULL, NULL, outmode, NULL);

	if (hDstDS != NULL) {
		GDALClose(hDstDS);
		GDALClose(hSrcDS);
		CSLDestroy( papszOptions );
		out = SpatRaster(filename, {-1}, {""}, {}, {});
	} else {
		GDALClose(hSrcDS);
		CSLDestroy( papszOptions );
		out.setError("something went wrong");
	}
	if (heightmode==1) {
		out.setValueType(3);
		out.setNames({"viewshed"});
	} else if (heightmode==2) {
		out.setNames({"above_sea"});
	} else {
		out.setNames({"above_land"});
	}
	out = out.mask(*this, false, NAN, NAN, opt);
	return out;
}

#else


SpatRaster SpatRaster::viewshed(const std::vector<double> obs, const std::vector<double> vals, const double curvcoef, const int mode, const double maxdist, const int heightmode, SpatOptions &opt) {
	SpatRaster out;
	out.setError("viewshed is not available for your version of GDAL. Need 3.1 or higher");
	return out;
}


#endif

std::string doubleToAlmostChar(double value){
    std::stringstream ss ;
    ss << value;
	std::string out = ss.str();
	return out;
}

SpatRaster SpatRaster::proximity(double target, double exclude, bool keepNA, std::string unit, bool buffer, double maxdist, bool remove_zero, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (nlyr() > 1) {
		out.addWarning("only the first layer is processed");
	}

	if (!hasValues()) {
		out.setError("input raster has no values");
		return out;
	}

	std::string filename = opt.get_filename();
	std::string driver;
	if (filename.empty()) {
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
		if (driver.empty()) {
			setError("cannot guess file type from filename");
			return out;
		}
		std::string errmsg;
		if (!can_write({filename}, filenames(), opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}
	}

	// GDAL proximity algo fails with other drivers? See #1116
//	driver = "MEM";

	GDALDatasetH hSrcDS, hDstDS;

	GDALDriverH hDriver = GDALGetDriverByName( driver.c_str() );
	if ( hDriver == NULL ) {
		out.setError("empty driver");
		return out;
	}

	GIntBig diskNeeded = ncell() * 4;
	char **papszOptions = set_GDAL_options(driver, diskNeeded, false, opt.gdal_options);
	papszOptions = CSLSetNameValue(papszOptions, "DISTUNITS", "GEO");

	SpatOptions ops(opt);
	SpatRaster x;
	bool mask = false;
	std::vector<double> mvals;
	if (buffer) {
		if (remove_zero) {
			x = isnotnan(true, ops);
		}		
		papszOptions = CSLSetNameValue(papszOptions, "MAXDIST", doubleToAlmostChar(maxdist).c_str());
		papszOptions = CSLSetNameValue(papszOptions, "FIXED_BUF_VAL", doubleToAlmostChar(1.0).c_str());
	} else if (std::isnan(target)) {
		if (std::isnan(exclude)) { // no exclusions
			x = isnotnan(false, ops);
		} else { // exclusion becomes target and is masked later
			x = replaceValues({exclude, NAN}, {0, 0}, 1, true, 1, false, ops);
			mvals.push_back(exclude);
			mask = true;
		}
	} else {
		//option for keepNA does not work, perhaps because of int conversion
		//papszOptions = CSLSetNameValue(papszOptions, "USE_INPUT_NODATA", "YES");
		if (std::isnan(exclude)) {
			if (keepNA) {
				x = replaceValues({target, NAN}, {0, 0}, 1, true, 1, false, ops);
				mvals.push_back(exclude);
				mask = true;
			} else {
				x = replaceValues({target}, {0}, 1, true, 1, false, ops);
			}
		} else {
			x = replaceValues({exclude, target}, {0, 0}, 1, true, 1, false, ops);
			mvals.push_back(exclude);
			mask = true;
		}
	}
	if (x.hasValues()) {
		if (!x.open_gdal(hSrcDS, 0, false, ops)) {
			out.setError("cannot open input dataset");
			return out;
		}
	} else if (!open_gdal(hSrcDS, 0, false, ops)) {
		out.setError("cannot open input dataset");
		return out;
	}

	std::string tmpfile = tempFile(opt.get_tempdir(), opt.pid, ".tif");
	std::string fname = mask ? tmpfile : filename;
	if (!out.create_gdalDS(hDstDS, fname, driver, false, 0, {false}, {1}, {0}, ops)) {
		out.setError("cannot create new dataset");
		GDALClose(hSrcDS);
		return out;
	}

	GDALRasterBandH hSrcBand = GDALGetRasterBand(hSrcDS, 1);
	GDALRasterBandH hTargetBand = GDALGetRasterBand(hDstDS, 1);

	if (GDALComputeProximity(hSrcBand, hTargetBand, papszOptions, NULL, NULL) != CE_None) {
		out.setError("proximity algorithm failed");
		GDALClose(hSrcDS);
		GDALClose(hDstDS);
		CSLDestroy( papszOptions );
		return out;
	}

	GDALClose(hSrcDS);
	CSLDestroy( papszOptions );

	if (driver == "MEM") {
		if (!out.from_gdalMEM(hDstDS, false, true)) {
			out.setError("conversion failed (mem)");
			GDALClose(hDstDS);
			return out;
		}
		GDALClose(hDstDS);
	} else {
		if (!mask) {
			double adfMinMax[2];
			GDALComputeRasterMinMax(hTargetBand, true, adfMinMax);
			GDALSetRasterStatistics(hTargetBand, adfMinMax[0], adfMinMax[1], -9999, -9999);
		}
		GDALClose(hDstDS);
		out = SpatRaster(fname, {-1}, {""}, {}, {});
	}
	
	if (mask) {
		out = out.mask(*this, false, mvals, NAN, opt);
	} 
	return out;
}

SpatRaster SpatRaster::sieveFilter(int threshold, int connections, SpatOptions &opt) {

	SpatRaster out = geometry(1, true, true, true);

	if (!hasValues()) {
		out.setError("input raster has no values");
		return out;
	}
	if (!((connections == 4) || (connections == 8))) {
		out.setError("connections should be 4 or 8");
		return out;
	}
	if (threshold < 2) {
		out.setError("a threshold < 2 is not meaningful");
		return out;
	}

	std::string filename = opt.get_filename();
	std::string driver;
	if (filename.empty()) {
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
		if (driver.empty()) {
			setError("cannot guess file type from filename");
			return out;
		}
		std::string errmsg;
		if (!can_write({filename}, filenames(), opt.get_overwrite(), errmsg)) {
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

	//opt.datatype = "INT4S";
	if (!out.create_gdalDS(hDstDS, filename, driver, true, 0, source[0].has_scale_offset, source[0].scale, source[0].offset, opt)) {
		out.setError("cannot create new dataset");
		GDALClose(hSrcDS);
		return out;
	}

	GDALRasterBandH hSrcBand = GDALGetRasterBand(hSrcDS, 1);
	GDALRasterBandH hTargetBand = GDALGetRasterBand(hDstDS, 1);

	if (GDALSieveFilter(hSrcBand, nullptr, hTargetBand, threshold, connections, nullptr, NULL, NULL) != CE_None) {
		out.setError("sieve failed");
		GDALClose(hSrcDS);
		GDALClose(hDstDS);
		return out;
	}

	GDALClose(hSrcDS);
	if (driver == "MEM") {
		if (!out.from_gdalMEM(hDstDS, false, true)) {
			out.setError("conversion failed (mem)");
		}
		GDALClose(hDstDS);
		return out;
	}
	
	double adfMinMax[2];
	GDALComputeRasterMinMax(hTargetBand, true, adfMinMax);
	GDALSetRasterStatistics(hTargetBand, adfMinMax[0], adfMinMax[1], -9999, -9999);
	GDALClose(hDstDS);
	return SpatRaster(filename, {-1}, {""}, {}, {});
}



bool getGridderAlgo(std::string algo, GDALGridAlgorithm &a) {
	if (algo == "nearest") {
		a = GGA_NearestNeighbor;
	} else if (algo == "invdistpow") {
		a = GGA_InverseDistanceToAPower;
	} else if (algo == "invdistpownear") {
		a = GGA_InverseDistanceToAPowerNearestNeighbor;
	} else if (algo == "mean") {
		a = GGA_MovingAverage;
	} else if (algo == "min") {
		a = GGA_MetricMinimum;
	} else if (algo == "max") {
		a = GGA_MetricMaximum;
	} else if (algo == "range") {
		a = GGA_MetricRange;
	} else if (algo == "count") {
		a = GGA_MetricCount;
	} else if (algo == "distto") {
		a = GGA_MetricAverageDistance;
	} else if (algo == "distbetween") {
		a = GGA_MetricAverageDistancePts;
	} else if (algo == "linear") {
		a = GGA_Linear;
	} else {
		return false;
	}
	return true;
}


void *metricOptions(std::vector<double> op) {
	GDALGridDataMetricsOptions *poOptions = static_cast<GDALGridDataMetricsOptions *>(
		CPLCalloc(sizeof(GDALGridDataMetricsOptions), 1));

	#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 6
	#else
	poOptions->nSizeOfStructure = sizeof(GDALGridDataMetricsOptions);
	#endif 
		
	poOptions->dfRadius1 = op[0];
	poOptions->dfRadius2 = op[1];
	poOptions->dfAngle = op[2];
	poOptions->nMinPoints = std::max(0.0, op[3]);
	poOptions->dfNoDataValue = op[4];
	return poOptions;
}

void *invDistPowerOps(std::vector<double> op) {
	GDALGridInverseDistanceToAPowerOptions *poOptions = static_cast<GDALGridInverseDistanceToAPowerOptions *>(
		CPLCalloc(sizeof(GDALGridInverseDistanceToAPowerOptions), 1));

	#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 6
	#else
	poOptions->nSizeOfStructure = sizeof(GDALGridInverseDistanceToAPowerOptions);
	#endif
	
	poOptions->dfPower = op[0];
	poOptions->dfSmoothing = op[1];
	poOptions->dfRadius1 = op[2];
	poOptions->dfRadius2 = op[3];
	poOptions->dfAngle = op[4];
	poOptions->nMaxPoints = std::max(0.0, op[5]);
	poOptions->nMinPoints = std::max(0.0, op[6]);
	poOptions->dfNoDataValue = op[7];

	poOptions->dfAnisotropyRatio = 1;
	poOptions->dfAnisotropyAngle = 0;

	return poOptions;
}

void *invDistPowerNNOps(std::vector<double> op) {
	GDALGridInverseDistanceToAPowerNearestNeighborOptions *poOptions = static_cast<GDALGridInverseDistanceToAPowerNearestNeighborOptions *>(
		CPLCalloc(sizeof(GDALGridInverseDistanceToAPowerNearestNeighborOptions), 1));

	#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 6
	#else
	poOptions->nSizeOfStructure = sizeof(GDALGridInverseDistanceToAPowerNearestNeighborOptions);
	//poOptions->nMaxPointsPerQuadrant = 
	//poOptions->nMinPointsPerQuadrant = 
	#endif
	poOptions->dfPower = op[0];
	poOptions->dfSmoothing = op[1];
	poOptions->dfRadius = op[2];
	poOptions->nMaxPoints = std::max(0.0, op[3]);
	poOptions->nMinPoints = std::max(0.0, op[4]);
	poOptions->dfNoDataValue = op[5];
	return poOptions;
}

void *moveAvgOps(std::vector<double> op) {
	GDALGridMovingAverageOptions *poOptions = static_cast<GDALGridMovingAverageOptions *>(
		CPLCalloc(sizeof(GDALGridMovingAverageOptions), 1));

	#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 6
	#else
	poOptions->nSizeOfStructure = sizeof(GDALGridMovingAverageOptions);
	#endif


	poOptions->dfRadius1 = op[0];
	poOptions->dfRadius2 = op[1];
	poOptions->dfAngle = op[2];
	poOptions->nMinPoints = std::max(op[3], 0.0);
	poOptions->dfNoDataValue = op[4];
	return poOptions;
}


void *nearngbOps(std::vector<double> op) {
	GDALGridNearestNeighborOptions *poOptions = static_cast<GDALGridNearestNeighborOptions *>(
		CPLCalloc(sizeof(GDALGridNearestNeighborOptions), 1));

	#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 6
	#else
	poOptions->nSizeOfStructure = sizeof(GDALGridNearestNeighborOptions);
	#endif

	poOptions->dfRadius1 = op[0];
	poOptions->dfRadius2 = op[1];
	poOptions->dfAngle = op[2];
	poOptions->dfNoDataValue = op[3];
	return poOptions;
}

void *LinearOps(std::vector<double> op) {
	GDALGridLinearOptions *poOptions = static_cast<GDALGridLinearOptions *>(
		CPLCalloc(sizeof(GDALGridLinearOptions), 1));

	#if GDAL_VERSION_MAJOR <= 3 && GDAL_VERSION_MINOR < 6
	#else
	poOptions->nSizeOfStructure = sizeof(GDALGridLinearOptions);
	#endif

	poOptions->dfRadius = op[0];
	poOptions->dfNoDataValue = op[1];
	return poOptions;
}

SpatRaster SpatRaster::rasterizeWindow(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::string algo, std::vector<double> algops, SpatOptions &opt) {

	SpatRaster out = geometry(1);

	GDALGridAlgorithm eAlg;
		if (!getGridderAlgo(algo, eAlg)) {
		out.setError("unknown algorithm");
		return out;
	}
	void *poOptions;
	if (is_in_vector(algo, {"min", "max", "range", "count", "distto", "distbetween"})) {
		if (algops.size() != 5) {
			out.setError("incorrect algorithm options");
			return out;
		}
		poOptions = metricOptions(algops) ;
	} else if (algo == "mean") {
		if (algops.size() != 5) {
			out.setError("incorrect algorithm options");
			return out;
		}
		poOptions = moveAvgOps(algops);
	} else if (algo == "invdistpow") {
		if (algops.size() != 8) {
			out.setError("incorrect algorithm options");
			return out;
		}
		poOptions = invDistPowerOps(algops);
	} else if (algo == "invdistpownear") {
		if (algops.size() != 6) {
			out.setError("incorrect algorithm options");
			return out;
		}
		poOptions = invDistPowerNNOps(algops);
	} else if (algo == "nearest") {
		if (algops.size() != 4) {
			out.setError("incorrect algorithm options");
			return out;
		}
		poOptions = nearngbOps(algops);
	} else if (algo == "linear") {
		if (algops.size() != 2) {
			out.setError("incorrect algorithm options");
			return out;
		}
		poOptions = LinearOps(algops);
	} else {
		out.setError("unknown algorithm");
		return out;
	}

	SpatExtent e = out.getExtent();
	if (!out.writeStart(opt, out.filenames())) {
		return out;
	}

	GUInt32 npts = x.size();
	GDALGridContext *ctxt = GDALGridContextCreate(eAlg, poOptions, npts, &x[0], &y[0], &z[0], true);
	CPLFree( poOptions );

	double rsy = out.yres() / 2;
	size_t ncs = out.ncol();
	BlockSize bs = out.getBlockSize(opt);
	std::vector<double> v;
	for (size_t i=0; i < bs.n; i++) {
		double ymax = yFromRow(bs.row[i]) + rsy;
		double ymin = yFromRow(bs.row[i] + bs.nrows[i] - 1) - rsy;
		v.resize(bs.nrows[i] * ncs);

		CPLErr eErr = GDALGridContextProcess(ctxt, e.xmin, e.xmax, ymin, ymax, ncs, bs.nrows[i], GDT_Float64, &v[0], NULL, NULL);

		if ( eErr != CE_None ) {
			out.setError("something went wrong");
			GDALGridContextFree(ctxt);
			return out;
		}

		std::vector<double> f;
		f.reserve(v.size());
		for (size_t j=0; j < bs.nrows[i]; j++) {
			unsigned start = (bs.nrows[i] - 1 - j) * ncs;
			f.insert(f.end(), v.begin()+start, v.begin()+start+ncs);
		}
		if (!out.writeBlock(f, i)) {
			GDALGridContextFree(ctxt);
			return out;
		}
	}
	GDALGridContextFree(ctxt);
	out.writeStop();
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



