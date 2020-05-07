#include <algorithm>
#include "spatRaster.h"
#include "string_utils.h"
#include "file_utils.h"

#include "crs.h"

#ifdef useGDAL

	#include "cpl_port.h"
	#include "cpl_conv.h" // CPLFree()
	#include "gdal_version.h"
	#include "gdalwarper.h"
	#include "ogr_srs_api.h"

	#if (!(GDAL_VERSION_MAJOR == 2 && GDAL_VERSION_MINOR < 1))
		#include "gdal_utils.h"
		#define GDALutils
	#endif

	#include "gdal_errors.h"


std::vector<double> getValuesMEM(GDALDatasetH hDS, unsigned ncol, unsigned nrow, unsigned nlyr) {

	unsigned ncell = ncol * nrow;
	std::vector<double> out;
	out.reserve(ncell * nlyr);
	CPLErr err = CE_None;
	int hasNA;
	for (size_t i=0; i < nlyr; i++) {
		GDALRasterBandH hBand = GDALGetRasterBand(hDS, i+1);
		std::vector<double> lyrout(ncell);
		err = GDALRasterIO(hBand, GF_Read, 0, 0, ncol, nrow, &lyrout[0], ncol, nrow, GDT_Float64, 0, 0);
		if (err != CE_None ) {
			return out;
			//setError("CE_None");
			break;
		}
		double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
		if (hasNA) std::replace(lyrout.begin(), lyrout.end(), naflag, (double) NAN);
		out.insert(out.end(), lyrout.begin(), lyrout.end());
	}
	return out;
}


GDALDatasetH ds_create(SpatRaster &x, std::string format, bool fill, SpatOptions &opt) {

	GDALDatasetH hDS;

	char **papszOptions = NULL;
	for (size_t i=0; i<opt.gdal_options.size(); i++) {
		std::vector<std::string> wopt = strsplit(opt.gdal_options[i], "=");
		if (wopt.size() == 2) {
			papszOptions = CSLSetNameValue( papszOptions, wopt[0].c_str(), wopt[1].c_str() );
		}
	}
	const char *pszFormat = format.c_str();
	GDALDriverH hDrv = GDALGetDriverByName(pszFormat);

	std::string filename = opt.filename;
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
		x.setError("CRS failure");
		//return false ;
	}
	char *pszSRS_WKT = NULL;
	OSRExportToWkt( hSRS, &pszSRS_WKT );
	GDALSetProjection( hDS, pszSRS_WKT );
	CPLFree(pszSRS_WKT);

	OSRDestroySpatialReference( hSRS );
	return hDS;
}


bool valid_warp_method(std::string method) {
	std::vector<std::string> m { "near", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode", "max", "min", "med", "q1", "q3", "sum" };
	return (std::find(m.begin(), m.end(), method) != m.end());
}


SpatRaster SpatRaster::warp(SpatRaster x, const std::string &method, SpatOptions &opt) {

	SpatRaster out = x.geometry(nlyr());

	if (!valid_warp_method(method)) {
		out.setError("invalid warp method");
		return out;
	}

	std::string errmsg;
	std::string filename = opt.filename;
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


#ifdef GDALutils
	std::vector<std::string> tmpfs;
	SpatOptions topt(opt);
	SpatRaster inp = sources_to_disk(tmpfs, true, topt);
	size_t nsrc = inp.source.size();
	std::vector<GDALDatasetH> src(nsrc);
	for (size_t i = 0; i < nsrc; i++) {
		src[i] = GDALOpen((const char *) inp.source[i].filename.c_str(), GA_ReadOnly);
	}

	std::vector<std::string> sops;
	if (filename != "") {
		std::vector<double> e = out.extent.asVector();
		std::vector<std::string> srs = out.getSRS();
		sops = {"-t_srs", srs[1], //"-dstnodata", "NAN",
			"-ts", std::to_string(out.ncol()), std::to_string(out.nrow()),
			"-te", std::to_string(e[0]), std::to_string(e[2]), std::to_string(e[1]), std::to_string(e[3]),
			"-r", method, "-ovr", "NONE", "-r", method, "-overwrite"};
	} else {
		sops = {"-r", method, "-overwrite"};
	}
	//if (format != "") {
	//	sops.push_back("-of");
	//	sops.push_back(format);
	//}

	std::vector <char *> options_char = string_to_charpnt(sops);
	GDALWarpAppOptions* gopts = GDALWarpAppOptionsNew(options_char.data(), nullptr);

	int err = 0;
	GDALDatasetH result;
	//std::string format = opt.get_filetype();
	if (filename == "") {
		GDALDatasetH hDstDS = ds_create(out, "MEM", true, opt);
		result = GDALWarp(NULL, hDstDS, nsrc, src.data(), gopts, &err);
		if (result != NULL) {
			out.setValues(getValuesMEM(result, out.ncol(), out.nrow(), out.nlyr()));
			GDALClose(result);
		}
	} else {
		const char *pszFilename = filename.c_str();
		result = GDALWarp(pszFilename, NULL, nsrc, src.data(), gopts, &err);
		if (result != NULL) {
			GDALClose(result);
			out = SpatRaster(filename);
		}
	}
	for (size_t i=0; i<tmpfs.size(); i++) remove(tmpfs[i].c_str());

	GDALWarpAppOptionsFree(gopts);
	for (size_t i = 0; i < nsrc; i++) {
		if (src[i] != NULL) GDALClose(src[i]);
	}

	if (err) {
		out.setError("an error occurred");
	} else if (result == NULL) {
		out.setError("no output");
	}

#else
	// we could use native project and resample
	out.setError("not implemented for GDAL < 2.1 -- let us know if that is a problem");
#endif
	return out;
}


SpatRaster SpatRaster::warpcrs(std::string x, const std::string &method, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr());
	// always to file

	std::string msg, wkt;
	if (!wkt_from_string(x, wkt, msg)) {
		out.setError(msg);
		return out;
	}
	std::string errmsg;
	std::string filename = opt.filename;
	if (filename == "") {
		filename = tempFile(opt.get_tempdir(), ".tif");
	} else {
		if (!can_write(filename, opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}
	}

#ifdef GDALutils
	if (!valid_warp_method(method)) {
		out.setError("invalid warp method");
		return out;
	}

	std::vector<std::string> tmpfs;
	SpatOptions topt(opt);
	out = sources_to_disk(tmpfs, true, topt);


	//std::string format = opt.filetype();

	//GDALDatasetH src = GDALOpen((const char *) source[0].filename.c_str(), GA_ReadOnly);
	size_t nsrc = out.source.size();
	std::vector<GDALDatasetH> src(nsrc);
	for (size_t i = 0; i < nsrc; i++) {
		src[i] = GDALOpen((const char *) out.source[i].filename.c_str(), GA_ReadOnly);
	}

	//GDALformat(filename, format);
	//"-dstnodata", "NAN",
	std::vector<std::string> sops = {"-t_srs", wkt, "-r", method, "-ovr", "NONE", "-dstnodata", "NAN"};
	if (opt.overwrite) {
		sops.push_back("-overwrite");
	}
	//if (format != "") {
	//	sops.push_back("-of");
	//	sops.push_back(format);
	//}

	std::vector <char *> options_char = string_to_charpnt(sops);
	GDALWarpAppOptions* gopts = GDALWarpAppOptionsNew(options_char.data(), NULL);

	GDALDatasetH result;
	int err = 0;
	const char *pszFilename = filename.c_str();
	result = GDALWarp(pszFilename, NULL, nsrc, &src[0], gopts, &err);
	if (result != NULL) {
		GDALClose(result);
		out = SpatRaster(filename);
	}
	GDALWarpAppOptionsFree(gopts);
	for (size_t i = 0; i < nsrc; i++) {
		if (src[i] != NULL) GDALClose(src[i]);
	}

	if (err) {
		out.setError("an error occurred");
	} else if (result == NULL) {
		out.setError("no output");
	}
	for (size_t i=0; i<tmpfs.size(); i++) remove(tmpfs[i].c_str());
	return out;

#else
	// we could use native project and resample
	out.setError("not implemented for GDAL < 2.1 -- let us know if that is a problem");
#endif

	return out;
}


bool gdalwarp(std::string src, std::string dst,	std::vector<std::string> options, std::vector<std::string> oo, std::vector<std::string> doo) {

	int err = 0;

	std::vector <char *> oo_char = string_to_charpnt(oo); // input open options
	std::vector <char *> doo_char = string_to_charpnt(doo); // output open options

// to be replaced by SpatRaster::sources
	std::vector<GDALDatasetH> src_pt(src.size());
	for (size_t i = 0; i < src.size(); i++) {
		src_pt[i] = GDALOpenEx((const char *) src.c_str(), GA_ReadOnly, NULL, oo_char.data(), NULL);
	}
	GDALDatasetH dst_ds = GDALOpenEx((const char *) dst.c_str(), GDAL_OF_RASTER | GA_Update, NULL, doo_char.data(), NULL);

	std::vector <char *> options_char = string_to_charpnt(options);
	GDALWarpAppOptions* opt = GDALWarpAppOptionsNew(options_char.data(), NULL);


	GDALDatasetH result = GDALWarp(dst_ds == NULL ? (const char *) dst.c_str() : NULL, dst_ds, src.size(), src_pt.data(), opt, &err);
	GDALWarpAppOptionsFree(opt);
	for (size_t i = 0; i < src.size(); i++) {
		if (src_pt[i] != NULL)	GDALClose(src_pt[i]);
	}
	if (result != NULL)	GDALClose(result);
	return (result != NULL) || (!err) ;
}

#endif




/*
#include <vector>
#include "spatRaster.h"
#include "vecmath.h"

#ifdef useGDAL
	#include "crs.h"
#endif


SpatRaster SpatRaster::warp(SpatRaster x, std::string method, SpatOptions &opt) {

	unsigned nl = nlyr();
	SpatRaster out = x.geometry(nl);
	out.setNames(getNames());
	std::vector<std::string> f {"bilinear", "ngb"};
	if (std::find(f.begin(), f.end(), method) == f.end()) {
		out.setError("unknown warp method");
		return out;
	}
	if (!hasValues()) {
		return out;
	}

	std::string crsin = getCRS();
	std::string crsout = out.getCRS();
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


*/