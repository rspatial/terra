#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "crs.h"
#include "file_utils.h"

#include "gdal_utils.h"


SpatRaster SpatRaster::grasterize(SpatVector x, std::string field, std::vector<double> values, double background, bool update, bool touches, bool inverse, SpatOptions &opt) {

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

