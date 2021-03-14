#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "file_utils.h"

#include "gdal_alg.h"
#include "ogrsf_frmts.h"

#include "spatFactor.h"
#include "string_utils.h"
#include "recycle.h"
#include "gdalio.h"

SpatRaster rasterizePoints(SpatVector p, SpatRaster r, std::vector<double> values, double background, SpatOptions &opt) {
	r.setError("not implemented yet");
	return(r);
}


SpatRaster SpatRaster::rasterize2(SpatVector x, std::string field, std::vector<double> values, 
	double background, bool touches, bool add, bool weights, SpatOptions &opt) {

	std::string gtype = x.type();

	if (weights && (gtype == "polygons")) {
		SpatOptions sopts(opt);
		SpatRaster wout = geometry(1);
		unsigned agx = 1000 / ncol();
		agx = std::max((unsigned)10, agx); 
		unsigned agy = 1000 / nrow();
		agy = std::max((unsigned)10, agy);
		wout = wout.disaggregate({agx, agy}, sopts);
		field = "";
		double f = agx * agy;
		wout = wout.rasterize2(x, field, {1/f}, background, touches, add, false, sopts);
		wout = wout.aggregate({agx, agy}, "sum", true, opt);
		return wout;
	}

	SpatRaster out;
//	if ( !hasValues() ) update = false;
//	if (update) {
//		out = geometry();
//	} else {
		out = geometry(1);
		out.setNames({field});
//	}

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
//			if (!update) {
				std::vector<double> levels(f.levels.size());
				for (size_t i=0; i<levels.size(); i++) {
					levels[i] = f.levels[i];
				}
				out.setCategories(0, levels, f.labels);
//			}
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
	size_t nGeoms = x.size();
	if (values.size() != nGeoms) {
		recycle(values, nGeoms);
	}

	GDALDataset *vecDS = x.write_ogr("", "lyr", "Memory", true);
	if (x.hasError()) {
		out.setError(x.getError());
		return out;
	}
    std::vector<OGRGeometryH> ahGeometries;
	OGRLayer *poLayer = vecDS->GetLayer(0);
	poLayer->ResetReading();
	OGRFeature *poFeature;
	while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
		OGRGeometry *poGeometry = poFeature->StealGeometry();
        OGRGeometryH hGeom = poGeometry;
        ahGeometries.push_back( hGeom );
	}
	GDALClose(vecDS);


	std::string errmsg;
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
		std::string driver = opt.get_filetype();
		getGDALdriver(filename, driver);
		if (driver == "") {
			setError("cannot guess file type from filename");
			return out;
		}
		if (!can_write(filename, opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}	
	}

	GDALDatasetH rstDS;
/*
	if (update) {
		size_t nsrc = source.size();
		if (driver == "MEM") {
			// force into single source
			SpatOptions svopt;
			std::vector<double> v = getValues();
			out.setValues(v, svopt);
			if (!out.open_gdal(rstDS, 0, opt)) {
				out.setError("cannot open dataset");
				return out;
			}
		} else {
			// make a copy first
			// including for the odd case that MEM is false but the source in memory
			if ( (nsrc > 1) || (!sources_from_file()) ) {
				SpatRaster out = writeRaster(opt);
			} else {
				// writeRaster should do the below? copyRaster?
				GDALDatasetH hSrcDS = GDALOpen(source[0].filename.c_str(), GA_ReadOnly );
				if(hSrcDS == NULL) {
					out.setError("cannot open source dataset");
					return out;
				}
				GDALDriverH hDriver = GDALGetDatasetDriver(hSrcDS);
				GDALDatasetH hDstDS = GDALCreateCopy( hDriver, filename.c_str(), hSrcDS, FALSE, NULL, NULL, NULL );
				GDALClose(hSrcDS);
				if(hDstDS == NULL) {
					out.setError("cannot create dataset");
					return out;
				}
				GDALClose(hDstDS);
			}
			rstDS = GDALOpen( filename.c_str(), GA_Update);	
		}
	
	} else {
*/		
		if (add) {
			background = 0;
		}
		if (!out.create_gdalDS(rstDS, filename, driver, true, background, opt)) {
			out.setError("cannot create dataset");
			return out;
		}
//	}

	std::vector<int> anBandList = {1};
//	papszOptions
//	std::vector <char *> options_char = string_to_charpnt(options);
//	GDALRasterizeOptions* ropts = GDALRasterizeOptionsNew(options_char.data(), NULL);
	char** papszOptions = NULL;
	if (touches) {
		papszOptions = CSLSetNameValue(papszOptions, "ALL_TOUCHED", "TRUE"); 
	}
	if (add) {
		papszOptions = CSLSetNameValue(papszOptions, "MERGE_ALG", "ADD"); 
	}
	
	
	CPLErr err = GDALRasterizeGeometries(rstDS, 
			static_cast<int>(anBandList.size()), &(anBandList[0]),
            static_cast<int>(ahGeometries.size()), &(ahGeometries[0]),
			NULL, NULL,
            &(values[0]), 
			papszOptions, 
			NULL, NULL);
			
	if ( err != CE_None ) {
		Rcpp::Rcout << err << std::endl;
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
		out = SpatRaster(filename, {-1}, {""});
	}
	return out;
}


std::vector<double> SpatRaster::rasterizeCells(SpatVector &v, bool touches) { 
// note that this is only for lines and polygons
    SpatOptions opt;
	SpatRaster r = geometry(1);
	SpatExtent e = getExtent();
	e.intersect(v.getExtent());
	if ( !e.valid() ) {
		std::vector<double> out;
		return out;
	}
	SpatRaster rc = r.crop(e, "out", opt);
#if GDAL_VERSION_MAJOR >= 3		
	std::vector<double> feats(1, 1) ;		
    SpatRaster rcr = rc.rasterize1(v, "", feats, {""}, NAN, false, touches, false, false, opt); 
#else
	std::vector<double> feats(v.size(), 1) ;		
    SpatRaster rcr = rc.rasterize(v, "", feats, {""}, NAN, false, touches, false, false, opt); 
#endif
	SpatVector pts = rcr.as_points(false, true, opt);
    SpatDataFrame vd = pts.getGeometryDF();
    std::vector<double> x = vd.getD(0);
    std::vector<double> y = vd.getD(1);
	std::vector<double> cells = r.cellFromXY(x, y);
	return cells;
}

std::vector<std::vector<double>> SpatRaster::rasterizeCellsWeights(SpatVector &v, bool touches) { 
// note that this is only for polygons
    SpatOptions opt;
	opt.progress = nrow()+1;
	SpatRaster rr = geometry(1);
	std::vector<unsigned> fact = {10, 10};
	SpatExtent e = getExtent();
	e.intersect(v.getExtent());
	std::vector<std::vector<double>> out(2);
	if ( !e.valid() ) {
		return out;
	}
	SpatRaster r = rr.crop(v.extent, "out", opt);
	r = r.disaggregate(fact, opt);
#if GDAL_VERSION_MAJOR >= 3
	std::vector<double> feats(1, 1) ;		
	r = r.rasterize1(v, "", feats, {""}, NAN, false, touches, false, false, opt); 
#else
	std::vector<double> feats(v.size(), 1) ;		
	r = r.rasterize(v, "", feats, {""}, NAN, false, touches, false, false, 	opt); 
#endif
	r = r.arith(100.0, "/", false, opt);
	r = r.aggregate(fact, "sum", true, opt);
	SpatVector pts = r.as_points(true, true, opt);
	SpatDataFrame vd = pts.getGeometryDF();
	std::vector<double> x = vd.getD(0);
	std::vector<double> y = vd.getD(1);
	out[0] = rr.cellFromXY(x, y);
	out[1] = pts.df.dv[0];
	return out;
}

