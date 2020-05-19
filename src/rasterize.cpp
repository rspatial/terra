#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "crs.h"
#include "file_utils.h"

#include "gdal_utils.h"


SpatRaster rasterizePoints(SpatVector p, SpatRaster r, std::vector<double> values, double background, SpatOptions &opt) {
	r.setError("not implemented yet");
	return(r);
}



#if GDAL_VERSION_MAJOR >= 3


SpatRaster SpatRaster::rasterize(SpatVector x, std::string field, std::vector<double> values, double background, bool update, bool touches, bool inverse, SpatOptions &opt) {

	SpatRaster out;
	if ( !hasValues() ) update = false;
	
	if (update) {
		out = geometry();
	} else {
		out = geometry(1);
		out.setNames({""});
	}

	std::string errmsg;
	std::string filename = opt.get_filename();
	std::string driver = filename == "" ? "MEM" : "GTiff";

	bool canRAM = canProcessInMemory(4);
	if (filename == "") {
		if (!canRAM || opt.get_todisk()) {
			filename = tempFile(opt.get_tempdir(), ".tif");
			opt.set_filename(filename);
			driver = "GTiff";
		} 
	} else {
		if (!can_write(filename, opt.get_overwrite(), errmsg)) {
			out.setError(errmsg);
			return out;
		}
		//if (canRAM) driver == "MEM";
	}


	std::vector<std::string> options; 
	if (inverse) options.push_back("-i");
	if (touches) options.push_back("-at");

	if (field != "") {
		std::vector<std::string> nms = x.get_names();
		if (!is_in_vector(field, nms)) {
			out.setError("field " + field + " not found");
			return out;
		}
		if (!update) out.setNames({field});
		options.push_back("-a");
		options.push_back(field);
	} else {
		if (values.size() == 1) {
			options.push_back("-burn");
			options.push_back(std::to_string(values[0]));
		} else if (values.size() == x.size()) {
			std::string burnvar = "rst_var";
			if (!x.lyr.df.add_column(values, burnvar)) {
				out.setError("this does not work??");
				return out;
			}
			options.push_back("-a");
			options.push_back(burnvar);
		} else {
			out.setError("the length of values must 1 or the number of features");
			return out;
		}
	}

	GDALDatasetH rstDS, vecDS;

	if (update) {
		size_t nsrc = source.size();
		if (driver == "MEM") {
			// force into single source
			out.setValues(getValues());
			if (!out.open_gdal(rstDS, 0)) {
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
		for (size_t i=0; i<nlyr(); i++) {
			options.push_back("-b");
			options.push_back(std::to_string(i+1));
		}
		
	} else {
		if (!out.create_gdalDS(rstDS, filename, driver, true, background, opt.gdal_options)) {
			out.setError("cannot create dataset");
			return out;
		}
	}
	
	
	GDALDataset *poDS = x.write_ogr("", "lyr", "Memory", true);
	vecDS = poDS->ToHandle(poDS);

	std::vector <char *> options_char = string_to_charpnt(options);
	GDALRasterizeOptions* ropts = GDALRasterizeOptionsNew(options_char.data(), NULL);

	int err = 0;
	GDALDatasetH hDst = GDALRasterize(NULL, rstDS, vecDS, ropts, &err);
	GDALRasterizeOptionsFree(ropts);
	if (err != 0) {
		setError("error "+ std::to_string(err));
	}
	
	if (driver == "MEM") {
		bool test = out.from_gdalMEM(hDst, false, true); 
		GDALClose( hDst );
		if (!test) {
			out.setError("wat nu?");
			return out;
		}
		if (update) {
			// seems to a bug that the input layers are returned as well
			out.source[0].values.erase(out.source[0].values.begin(), out.source[0].values.begin()+out.size());
		}
		if (filename != "") {
			writeRaster(opt);
		}
	} else {
		for (size_t i=0; i < nlyr(); i++) { //currently always one band
			GDALRasterBandH hBand = GDALGetRasterBand(hDst, i+1);
			double adfMinMax[2];
			bool approx = ncell() > 10e+8;
			GDALComputeRasterMinMax(hBand, approx, adfMinMax);
			GDALSetRasterStatistics(hBand, adfMinMax[0], adfMinMax[1], NAN, NAN);		
		}
		GDALClose( hDst );
		out = SpatRaster(filename);
	}

	return out;
}

#else 
	

std::vector<double> rasterize_polygon(std::vector<double> r, double value, const std::vector<double> &pX, const std::vector<double> &pY, const unsigned startrow, const unsigned nrows, const unsigned ncols, const double xmin, const double ymax, const double rx, const double ry) {

	unsigned n = pX.size();
	std::vector<unsigned> nCol(n);
	for (size_t row=0; row < (nrows); row++) {
		double y = ymax - (startrow+row+0.5) * ry;

		// find nodes.
		unsigned nodes = 0;
		size_t j = n-1;
		for (size_t i=0; i<n; i++) {
			if (((pY[i] < y) && (pY[j] >= y)) || ((pY[j] < y) && (pY[i] >= y))) {
			//	nCol[nodes++]=(int)  (((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx);
				double nds = ((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx;
				nds = nds < 0 ? 0 : nds;
		        nds = nds > ncols ? ncols : nds;
				nCol[nodes] = (unsigned) nds;
				nodes++;
			}
			j = i;
		}

		std::sort(nCol.begin(), nCol.begin()+nodes);
		unsigned ncell = ncols * row;

		//  Fill the cells between node pairs.
		for (size_t i=0; i < nodes; i+=2) {
			if (nCol[i+1] > 0 && nCol[i] < ncols) {
				for (size_t col = nCol[i]; col < nCol[i+1]; col++) {
					r[col + ncell] = value;
				}
			}
		}
	}
	return(r);
}




SpatRaster rasterizePolygons(SpatVector p, SpatRaster r, std::vector<double> value, double background, SpatOptions &opt) {

	SpatRaster out = r.geometry(1);

  	if (!out.writeStart(opt)) { return out; }
	double resx = out.xres();
	double resy = out.yres();
	SpatGeom poly;
	SpatPart part;
	SpatHole hole;
	unsigned n = p.size();
	unsigned nc = out.ncol();
	SpatExtent extent = out.getExtent();

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v(out.bs.nrows[i] * nc, background);

		for (size_t j = 0; j < n; j++) {
			poly = p.getGeom(j);
			unsigned np = poly.size();

			for (size_t k = 0; k < np; k++) {
				part = poly.getPart(k);
				if (part.hasHoles()) {
					std::vector<double> vv = rasterize_polygon(v, value[j], part.x, part.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
					for (size_t h=0; h < part.nHoles(); h++) {
						hole = part.getHole(h);
						vv = rasterize_polygon(vv, background, hole.x, hole.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
					}
					for (size_t q=0; q < vv.size(); q++) {
						if ((vv[q] != background) && (!std::isnan(vv[q]))) {
							v[q] = vv[q];
						}
					}
				} else {
					v = rasterize_polygon(v, value[j], part.x, part.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
				}
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;

	}
	out.writeStop();
	return(out);
}




std::vector<double> rasterize_line(std::vector<double> r, double value, const std::vector<double> &pX, const std::vector<double> &pY, const unsigned startrow, const unsigned nrows, const unsigned ncols, const double xmin, const double ymax, const double rx, const double ry) {
	unsigned n = pX.size();
	for (size_t row=0; row<nrows; row++) {
		double y = ymax - (startrow+row+0.5) * ry;
		unsigned ncell = ncols * row;
		for (size_t i=1; i<n; i++) {
            size_t j = i-1;
			if (((pY[i] < y) && (pY[j] >= y)) || ((pY[j] < y) && (pY[i] >= y))) {
				double col = ((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx;
				if ((col >= 0) & (col < ncols)) {
                    r[ncell + col] = value;
				}
			}
		}
	}
	return(r);
}



SpatRaster rasterizeLines(SpatVector p, SpatRaster r, std::vector<double> value, double background, SpatOptions &opt) {

	SpatRaster out = r.geometry(1);
  	if (!out.writeStart(opt)) { return out; }
	double resx = out.xres();
	double resy = out.yres();
	SpatGeom line;
	SpatPart part;
	unsigned n = p.size();
	SpatExtent extent = out.getExtent();

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
		for (size_t j = 0; j < n; j++) {
			line = p.getGeom(j);
			unsigned nln = line.size();
			for (size_t k = 0; k < nln; k++) {
				part = line.getPart(k);
				v = rasterize_line(v, value[j], part.x, part.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
	}
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::rasterize(SpatVector x, std::string field, std::vector<double> values, double background, bool update, bool touches, bool inverse, SpatOptions &opt) {

//SpatRaster SpatRaster::rasterize(SpatVector x, std::vector<double> values, double background, bool update, SpatOptions &opt) {


	std::string gtype = x.type();
	SpatRaster out = geometry(1);

	if (field != "") {
		std::vector<std::string> nms = x.get_names();
		if (!is_in_vector(field, nms)) {
			out.setError("field " + field + " not found");
			return out;
		}
		if (!update) out.setNames({field});
	} else {
		if (values.size() == 1) {
			values = std::vector<double>(x.size(), values[0]);
		} else if (values.size() != x.size()) {
			out.setError("the length of values must 1 or the number of features");
			return out;
		}
	}



	SpatOptions opts(opt);
	if (!update) {
		opts = opt;
	}
	if (gtype == "polygons") {
		out = rasterizePolygons(x, out, values, background, opts);
	} else if (gtype == "lines") {
		out = rasterizeLines(x, out, values, background, opts);
	}  else {
		out = rasterizePoints(x, out, values, background, opts);
	}
	if (update) out = cover(out, background, opt);

	if (touches) {
		out.addWarning("argument touches is not supported with your version of GDAL");	
	}
	if (inverse) {
		out.addWarning("argument inverse is not supported with your version of GDAL");	
	}

	return out;
}


#endif
