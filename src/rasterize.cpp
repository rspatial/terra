#include "gdalwarper.h"
#include "ogr_spatialref.h"

#include "spatRaster.h"
#include "string_utils.h"
#include "crs.h"
#include "file_utils.h"

#include "gdal_utils.h"


SpatRaster SpatRaster::grasterize(SpatVector x, std::string field, std::vector<double> values, bool touches, bool inverse, SpatOptions &opt) {

	SpatRaster out = geometry(1);

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


	std::vector<std::string> options; 
	if (inverse) options.push_back("-i");
	if (touches) options.push_back("-at");

	if (field != "") {
		std::vector<std::string> nms = x.get_names();
		if (!is_in_vector(field, nms)) {
			out.setError("field " + field + " not found");
			return out;
		}
		out.setNames({field});
		options.push_back("-a");
		options.push_back(field);
	} else {
		out.setNames({"rasterized"});
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

	std::string driver = filename == "" ? "MEM" : "GTiff";
	GDALDatasetH rstDS, vecDS;
	if (!out.create_gdalDS(rstDS, filename, driver, true, opt.gdal_options)) {
		return out;
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

