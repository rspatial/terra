
#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "file_utils.h"

#include "gdal_alg.h"
#include "ogrsf_frmts.h"

#include "spatFactor.h"
#include "recycle.h"
#include "gdalio.h"


SpatRaster SpatRaster::rasterizePoints(SpatVector x, std::string fun, std::vector<double> values, double background, SpatOptions &opt) {
	
	SpatRaster out = geometry(1, false, false, false);
	if (!out.writeStart(opt)) {
		return out;
	}
	if (fun != "count" && (values.size() != x.size())) {
		out.setError("values do not match geometries");
		return out;
	}	

	size_t nc = ncol();
	std::vector<std::vector<double>> pxy = x.coordinates();
	std::vector<double> cells = cellFromXY(pxy[0], pxy[1]);
	if (fun == "count") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), 0);
			for (size_t j=0; j<cells.size(); j++) {
				if (cells[j] >= cmin && cells[j] <= cmax) {
					size_t k = cells[j] - cmin; 
					v[k]++; 
				}
			}
			if (background != 0) {
				for (size_t j=0; j<v.size(); j++) {
					if (v[j] == 0) {
						v[j] = background;
					}
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "sum") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), 0);
			for (size_t j=0; j<cells.size(); j++) {
				if (cells[j] >= cmin && cells[j] <= cmax) {
					size_t k = cells[j] - cmin; 
					v[k] += values[j]; 
				}
			}
			if (background != 0) {
				for (size_t j=0; j<v.size(); j++) {
					if (v[j] == 0) {
						v[j] = background;
					}
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else if (fun == "mean") {
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), 0);
			std::vector<double> cnt = v;
			for (size_t j=0; j<cells.size(); j++) {
				if (cells[j] >= cmin && cells[j] <= cmax) {
					size_t k = cells[j] - cmin; 
					v[k] += values[j]; 
					cnt[k]++;
				}
			}
			for (size_t j=0; j<cnt.size(); j++) {
				if (cnt[j] > 0) {
					v[j] /= cnt[j];
				}
			}
			if (background != 0) {
				for (size_t j=0; j<v.size(); j++) {
					if (v[j] == 0) {
						v[j] = background;
					}
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	} else { // "last"
		for (size_t i=0; i < out.bs.n; i++) {
			double cmin = out.bs.row[i] * nc;
			double cmax = (out.bs.row[i]+out.bs.nrows[i]) * nc - 1;
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
			for (size_t j=0; j<cells.size(); j++) {
				if (cells[j] >= cmin && cells[j] <= cmax) {
					size_t k = cells[j] - cmin; 
					v[k] = values[j]; 
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
	}
	out.writeStop();
	return out;
}


SpatRaster SpatRaster::rasterizeGeom(SpatVector x, std::string unit, std::string fun, SpatOptions &opt) {

	if (x.type() != "points") {

		SpatRaster out = geometry(1, false, false, false);
		SpatOptions ops(opt);

		std::vector<std::string> ss {"m", "km"};
		if (std::find(ss.begin(), ss.end(), unit) == ss.end()) {
			out.setError("invalid unit (not 'm' or 'km')");
			return out;
		}
		if ((x.type() == "lines")) {
			ss = {"count", "length", "crosses"};
			if (std::find(ss.begin(), ss.end(), fun) == ss.end()) {
				out.setError("invalid value for 'fun' (not 'count', 'crosses', or 'length')");
				return out;
			}
		} else {
			ss = {"area", "count"};
			if (std::find(ss.begin(), ss.end(), fun) == ss.end()) {
				out.setError("invalid value for 'fun' (not 'area' or 'count')");
				return out;
			}
		}

		SpatRaster empty = out.geometry();
		SpatExtent e = out.getExtent();
		double rsy = out.yres() / 2;

		double m = unit == "m" ? 1 : 1000;
		if (!x.is_lonlat()) {
			double tom = x.srs.to_meter();
			tom = std::isnan(tom) ? 1 : tom;
			m *= tom;
		}
		if (x.type() == "lines") {
			out.setNames({"length"});
		} else {
			out.setNames({"area"});	
			m *= m;
		}
		opt.ncopies = std::max(opt.ncopies, (unsigned)4) * 8;
		if (!out.writeStart(opt)) {
			return out;
		}
		for (size_t i=0; i < out.bs.n; i++) {
			e.ymax = yFromRow(out.bs.row[i]) + rsy;
			e.ymin = yFromRow(out.bs.row[i] + out.bs.nrows[i] - 1) - rsy;
			SpatRaster tmp = empty.crop(e, "near", ops);

			SpatVector p = tmp.as_polygons(true, false, false, false, false, ops);
			std::vector<double> v(out.bs.nrows[i] * out.ncol(), 0);

			if (fun == "crosses") {
				std::vector<int> r = p.relate(x, "crosses");
				size_t nx = x.size();
				for (size_t j=0; j< r.size(); j++) {
					size_t k= j / nx;
					v[k] += r[j];
				}
			} else {
				std::vector<long> cell(p.size());
				std::iota(cell.begin(), cell.end(), 0);
				p.df.add_column(cell, "cell");
				p = p.intersect(x);
				std::vector<double> stat;
				if (x.type() == "lines") {
					stat = p.length();
				} else {
					stat = p.area("m", false, {});
				}
				if (fun == "count") {
					for (size_t j=0; j<stat.size(); j++) {
						size_t k = p.df.iv[0][j]; 
						v[k]++; 
					}
				} else {
					for (size_t j=0; j<stat.size(); j++) {
						size_t k = p.df.iv[0][j]; 
						v[k] += (stat[j] / m); 
					}				
				}
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i]))  return out;
		}
		out.writeStop();
		return(out);

	} else {
		return rasterizePoints(x, "count", {}, 0.0, opt);
	}
}


SpatRaster SpatRaster::hardCopy(SpatOptions &opt) {
	SpatRaster out = geometry(-1, true, true);
	if (!hasValues()) {
		out.addWarning("raster has no values");
		return out;
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
 	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v; 
		readBlock(v, out.bs, i);
		if (!out.writeBlock(v, i)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}

/* gdalcopy
			GDALDatasetH hSrcDS = GDALOpenEx(out.source[0].filename.c_str(), GDAL_OF_RASTER | GDAL_OF_UPDATE, NULL, NULL, NULL);
			if(hSrcDS == NULL) {
				out.setError("cannot open dataset";)
				return false;
			}
			GDALDriverH hDriver = GDALGetDatasetDriver(hSrcDS);
			GDALDatasetH hDstDS = GDALCreateCopy( hDriver, filename.c_str(), hSrcDS, FALSE, NULL, NULL, NULL );
			GDALClose(hSrcDS);
			if(hDstDS == NULL) {
				out.setError("cannot create dataset";
				return false;
			}
			GDALClose(hDstDS);
*/

bool SpatRaster::getDSh(GDALDatasetH &rstDS, SpatRaster &out, std::string &filename, std::string &driver, double &naval, bool update, double background, SpatOptions &opt) {

	filename = opt.get_filename();
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
			out.setError("cannot guess file type from filename");
			return false;
		}
		std::string msg;
		if (!can_write(filename, opt.get_overwrite(), msg)) {
			out.setError(msg);
			return false;
		}
	}

	if (opt.names.size() == nlyr()) {
		out.setNames(opt.names);
	}


	if (update) {
		out = hardCopy(opt);
		//size_t ns = source.size();
		if (!out.open_gdal(rstDS, 0, true, opt)) {
			return false;
		}
	} else if (!out.create_gdalDS(rstDS, filename, driver, true, background, source[0].has_scale_offset, source[0].scale, source[0].offset, opt)) {
		out.setError("cannot create dataset");
		return false;
	}

	GDALRasterBandH hBand = GDALGetRasterBand(rstDS, 1);
	GDALDataType gdt = GDALGetRasterDataType(hBand);
	getNAvalue(gdt, naval);
	int hasNA;
	double naflag = GDALGetRasterNoDataValue(hBand, &hasNA);
	naval = hasNA ? naflag : naval;
	return true;
}





SpatRaster SpatRaster::rasterizeLyr(SpatVector x, double value, double background, bool touches, bool update, SpatOptions &opt) {

// not working well in some cases. See #552
	std::string gtype = x.type();
	SpatRaster out;
	out.setNames({"ID"});

	if ( !hasValues() ) update = false;
	if (update) { // all lyrs
		out = geometry();
	} else {
		out = geometry(1);
	}

	GDALDataset *vecDS = x.write_ogr("", "lyr", "Memory", false, true, std::vector<std::string>());
	if (x.hasError()) {
		out.setError(x.getError());
		return out;
	}

	OGRLayer *poLayer = vecDS->GetLayer(0);
#if GDAL_VERSION_MAJOR <= 2 && GDAL_VERSION_MINOR <= 2
	OGRLayerH hLyr = poLayer;
#else
	OGRLayerH hLyr = poLayer->ToHandle(poLayer);
#endif
    std::vector<OGRLayerH> ahLayers;
	ahLayers.push_back( hLyr );

	std::string driver, filename;
	GDALDatasetH rstDS;
	double naval;
	if (!getDSh(rstDS, out, filename, driver, naval, update, background, opt)) {
		return out;
	}
	if (std::isnan(value)) {
		// passing NULL instead may also work.
		value = naval;
	}

	std::vector<int> bands(out.nlyr());
	std::iota(bands.begin(), bands.end(), 1);
	std::vector<double> values(out.nlyr(), value);

	char** papszOptions = NULL;
	CPLErr err;
	if (touches) {
		papszOptions = CSLSetNameValue(papszOptions, "ALL_TOUCHED", "TRUE"); 
	}
	err = GDALRasterizeLayers(rstDS, static_cast<int>(bands.size()), &(bands[0]),
			1, &(ahLayers[0]), NULL, NULL, &(values[0]), papszOptions, NULL, NULL);

	CSLDestroy(papszOptions);

//	for (size_t i=0; i<ahGeometries.size(); i++) {
//		OGR_G_DestroyGeometry(ahGeometries[i]);
//	}
	GDALClose(vecDS);

	if ( err != CE_None ) {
		out.setError("rasterization failed");
		GDALClose(rstDS);
		return out;
	}

	if (driver == "MEM") {
		if (!out.from_gdalMEM(rstDS, false, true)) {
			out.setError("rasterization failed (mem)");
		}
	}

	GDALRasterBandH band = GDALGetRasterBand(rstDS, 1);
	double adfMinMax[2];
	GDALComputeRasterMinMax(band, false, adfMinMax);
	GDALSetRasterStatistics(band, adfMinMax[0], adfMinMax[1], -9999, -9999);

	GDALClose(rstDS);
	if (driver != "MEM") {
		out = SpatRaster(filename, {-1}, {""}, {});
	}
	return out;
}


SpatRaster SpatRaster::rasterize(SpatVector x, std::string field, std::vector<double> values, 
	double background, bool touches, bool add, bool weights, bool update, bool minmax, SpatOptions &opt) {

	std::string gtype = x.type();
	bool ispol = gtype == "polygons";
	if (weights) update = false;

	if (weights && ispol) {
		SpatOptions sopts(opt);
		SpatRaster wout = geometry(1);
		unsigned agx = 1000 / ncol();
		agx = std::max((unsigned)10, agx); 
		unsigned agy = 1000 / nrow();
		agy = std::max((unsigned)10, agy);
		wout = wout.disaggregate({agx, agy}, sopts);
		field = "";
		double f = agx * agy;
		wout = wout.rasterize(x, field, {1/f}, background, touches, add, false, false, false, sopts);
		wout = wout.aggregate({agx, agy}, "sum", true, opt);
		return wout;
	}

	SpatRaster out;
	if ( !hasValues() ) update = false;
	if (update) {
		out = hardCopy(opt);
	} else {
		out = geometry(1);
		if (field == "") {
			out.setNames({"layer"});
		} else {
			out.setNames({field});
		}
	}

	size_t nGeoms = x.size();
	if (nGeoms == 0) {
		if (!update) {
			out = out.init({background}, opt);
		}
		return out;
	}

	if (ispol && touches && add) {
		add = false;
		out.addWarning("you cannot use 'sum' and 'touches' at the same time");
	}

	if (field != "") {
		int i = x.df.get_fieldindex(field);
		if (i < 0) {
			out.setError("field " + field + " not found");
			return out;
		}
		std::string dt = x.df.get_datatype(field);
		if (dt == "string") {
			//std::vector<std::string> ss = ;
			SpatFactor f;
			f.set_values(x.df.getS(i));
			values.resize(f.v.size());
			for (size_t i=0; i<values.size(); i++) {
				values[i] = f.v[i];
			}
			if (!add && !update) {
				std::vector<long> u(f.labels.size());
				std::iota(u.begin(), u.end(), 0);
				std::vector<std::string> nms = getNames();
				out.setLabels(0, u, f.labels, field);
			}
			if (add) {
				add = false;
				addWarning("cannot add factors");
			}
		} else if (dt == "double") {
			values = x.df.getD(i);
		} else {
			std::vector<long> v = x.df.getI(i);
			values.resize(v.size());
			for (size_t i=0; i<values.size(); i++) {
				values[i] = v[i];
			}
		}
	}
	
	if (values.size() != nGeoms) {
		recycle(values, nGeoms);
	}

	GDALDataset *vecDS = x.write_ogr("", "lyr", "Memory", false, true, std::vector<std::string>());
	if (x.hasError()) {
		out.setError(x.getError());
		return out;
	}
    std::vector<OGRGeometry *> ogrGeoms;
	ogrGeoms.reserve(nGeoms);
	
	OGRLayer *poLayer = vecDS->GetLayer(0);
	poLayer->ResetReading();
	
	OGRFeature *poFeature;
	while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->GetGeometryRef();
		if (poGeometry != NULL) {
			OGRGeometry *copyGeom = poGeometry->clone();
			ogrGeoms.push_back( copyGeom );
		}
		OGRFeature::DestroyFeature( poFeature );
	}
	GDALClose(vecDS);

	std::string errmsg, driver, filename;
	GDALDatasetH rstDS;
	double naval;
	if (add) {	background = 0;	}

	if (!out.getDSh(rstDS, out, filename, driver, naval, update, background, opt)) {
		return out;
	}
	for (double &d : values) d = std::isnan(d) ? naval : d;
		// passing NULL instead may also work.

	std::vector<int> bands(out.nlyr());
	std::iota(bands.begin(), bands.end(), 1);
	rep_each(values, out.nlyr());

	char** papszOptions = NULL;
	CPLErr err;
	if (ispol && touches && (nGeoms > 1)) {
		// first to get the touches
		papszOptions = CSLSetNameValue(papszOptions, "ALL_TOUCHED", "TRUE"); 
		err = GDALRasterizeGeometries(rstDS, 
				static_cast<int>(bands.size()), &(bands[0]),
				static_cast<int>(ogrGeoms.size()), 
				(OGRGeometryH *) &(ogrGeoms[0]),
				NULL, NULL, &(values[0]), papszOptions, NULL, NULL);
		CSLDestroy(papszOptions);

		if ( err != CE_None ) {
			out.setError("rasterization failed");
			GDALClose(rstDS);
			for (size_t i=0; i<ogrGeoms.size(); i++) {
				OGR_G_DestroyGeometry(ogrGeoms[i]);
			}
			return out;
		}
		//GDALFlushCache(rstDS);
		// second time to fix the internal area
		err = GDALRasterizeGeometries(rstDS, 
				static_cast<int>(bands.size()), &(bands[0]),
				static_cast<int>(ogrGeoms.size()), 
				(OGRGeometryH *) &(ogrGeoms[0]),
				NULL, NULL, &(values[0]), NULL, NULL, NULL);

	} else {
		if (touches) {
			papszOptions = CSLSetNameValue(papszOptions, "ALL_TOUCHED", "TRUE"); 
		} else if (add) {
			papszOptions = CSLSetNameValue(papszOptions, "MERGE_ALG", "ADD"); 
		}
		err = GDALRasterizeGeometries(rstDS, 
				static_cast<int>(bands.size()), &(bands[0]),
				static_cast<int>(ogrGeoms.size()), 
				(OGRGeometryH *) &(ogrGeoms[0]),
				NULL, NULL, &(values[0]), papszOptions, NULL, NULL);

		CSLDestroy(papszOptions);
	}

	for (size_t i=0; i<ogrGeoms.size(); i++) {
		OGR_G_DestroyGeometry(ogrGeoms[i]);
	}

	if ( err != CE_None ) {
		out.setError("rasterization failed");
		GDALClose(rstDS);
		return out;
	}

	if (driver == "MEM") {
		if (!out.from_gdalMEM(rstDS, false, true)) {
			out.setError("rasterization failed (mem)");
		}
	}

	GDALRasterBandH band = GDALGetRasterBand(rstDS, 1);

	if (minmax) {
		double adfMinMax[2];
		GDALComputeRasterMinMax(band, false, adfMinMax);
		GDALSetRasterStatistics(band, adfMinMax[0], adfMinMax[1], -9999, -9999);
	}

	GDALClose(rstDS);
	if (driver != "MEM") {
		out = SpatRaster(filename, {-1}, {""}, {});
	} else {
		std::string fname = opt.get_filename();
		if ((fname != "") && (!update)) {
			out = out.writeRaster(opt);
		}
	}
	return out;
}


std::vector<double> SpatRaster::rasterizeCells(SpatVector &v, bool touches, SpatOptions &opt) { 
// note that this is only for lines and polygons
    SpatOptions ropt(opt);
	SpatRaster r = geometry(1);
	SpatExtent e = getExtent();
	e = e.intersect(v.getExtent());
	if ( !e.valid() ) {
		std::vector<double> out(1, NAN);
		return out;
	}

	SpatRaster rc = r.crop(e, "out", ropt);
	std::vector<double> feats(1, 1) ;
    SpatRaster rcr = rc.rasterize(v, "", feats, NAN, touches, false, false, false, false, ropt); 
	SpatVector pts = rcr.as_points(false, true, false, ropt);
	std::vector<double> cells;
	if (pts.size() == 0) {
		pts = v.as_points(false, true);
		SpatDataFrame vd = pts.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		cells = r.cellFromXY(x, y);
		cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
		if (cells.size() == 0) {
			cells.resize(1, NAN);
		}
	} else {
		SpatDataFrame vd = pts.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		cells = r.cellFromXY(x, y);
	}
	return cells;
}

void SpatRaster::rasterizeCellsWeights(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt) { 
// note that this is only for polygons
    SpatOptions ropt(opt);
	opt.progress = nrow()+1;
	SpatRaster rr = geometry(1);
	std::vector<unsigned> fact = {10, 10};
	SpatExtent e = getExtent();
	SpatExtent ve = v.getExtent();
	e = e.intersect(ve);
	if ( !e.valid() ) {
		return;
	}
	SpatRaster r = rr.crop(v.extent, "out", ropt);
	r = r.disaggregate(fact, ropt);
	std::vector<double> feats(1, 1) ;
	r = r.rasterize(v, "", feats, NAN, true, false, false, false, false, ropt); 
	r = r.arith(100.0, "/", false, ropt);
	r = r.aggregate(fact, "sum", true, ropt);
	SpatVector pts = r.as_points(true, true, false, ropt);
	if (pts.size() == 0) {
		weights.resize(1);
		weights[0] = NAN;
		cells.resize(1);
		cells[0] = NAN;
	} else {
		SpatDataFrame vd = pts.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		cells = rr.cellFromXY(x, y);
		weights = pts.df.dv[0];
	}
	return;
}

void SpatRaster::rasterizeCellsExact(std::vector<double> &cells, std::vector<double> &weights, SpatVector &v, SpatOptions &opt) { 
    
	SpatOptions ropt(opt);
	opt.progress = nrow()+1;
	SpatRaster r = geometry(1);
	r = r.crop(v.extent, "out", ropt);

//	if (r.ncell() < 1000) {
		std::vector<double> feats(1, 1) ;
		r = r.rasterize(v, "", feats, NAN, true, false, false, false, false, ropt); 

		SpatVector pts = r.as_points(true, true, false, ropt);
		if (pts.size() == 0) {
			weights.resize(1);
			weights[0] = NAN;
			cells.resize(1);
			cells[0] = NAN;
		} else {
			SpatDataFrame vd = pts.getGeometryDF();
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			cells = cellFromXY(x, y);

			SpatVector rv = r.as_polygons(false, false, false, true, false, ropt);
			std::vector<double> csize = rv.area("m", true, {});
			rv.df.add_column(csize, "area");
			rv.df.add_column(cells, "cells");
			rv = rv.crop(v);
			weights = rv.area("m", true, {});
			for (size_t i=0; i<weights.size(); i++) {
				weights[i] /= rv.df.dv[0][i];
			}
			cells = rv.df.dv[1];
		}
//	}

/*
// the below would need to remove all cells already included in the above
// because touches=false includes partly overlapped cells [#346]

	else {
		std::vector<double> feats(1, 1) ;
		SpatVector vv = v.as_lines();
		SpatRaster b = r.rasterize(vv, "", feats, NAN, true, false, false, false, false, opt); 
		SpatVector pts = b.as_points(true, true, opt);
		if (pts.nrow() > 0) {
			SpatDataFrame vd = pts.getGeometryDF();
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			cells = cellFromXY(x, y);

			SpatVector bv = b.as_polygons(false, false, false, true, opt);
			std::vector<double> csize = bv.area("m", true, {});
			bv.df.add_column(csize, "cellsize");
			bv.df.add_column(cells, "cellnr");
			bv = bv.crop(v);
			weights = bv.area("m", true, {});
			for (size_t i=0; i<weights.size(); i++) {
				weights[i] /= bv.df.dv[0][i];
			}
			cells = bv.df.dv[1];
		}
		// touches = false
		r = r.rasterize(v, "", feats, NAN, false, false, false, false, false, opt); 
		pts = r.as_points(true, true, opt);
		if (pts.nrow() > 0) {
			SpatDataFrame vd = pts.getGeometryDF();
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			std::vector<double> cells2 = cellFromXY(x, y);
			cells.insert(cells.end(), cells2.begin(), cells2.end());
			weights.resize(weights.size() + cells2.size(), 1);
		}

		if (cells.size() == 0) {
			weights.resize(1);
			weights[0] = NAN;
			cells.resize(1);
			cells[0] = NAN;
		}
	}
*/

}



