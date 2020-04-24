#include <vector>
#include <string>
#include "string_utils.h"

#include "cpl_port.h"
#include "cpl_conv.h" // CPLFree()
#include "gdal_version.h"

#include "gdal_errors.h" 

// code adapted from the 'sf' package by Edzer Pebesma et al

#if (!(GDAL_VERSION_MAJOR == 2 && GDAL_VERSION_MINOR < 1))
# include "gdal_utils.h" // requires >= 2.1



std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> openopts) {
	
	std::string ret = "";
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(openopts); // open options
	GDALInfoOptions* opt = GDALInfoOptionsNew(options_char.data(), NULL);
	GDALDatasetH ds = GDALOpenEx((const char *) filename.c_str(), GA_ReadOnly, NULL, oo_char.data(), NULL);
	if (ds == NULL) return ret; // #nocov
	char *ret_val = GDALInfo(ds, opt);
	ret = ret_val;
	CPLFree(ret_val);
	GDALInfoOptionsFree(opt);
	GDALClose(ds);
	return ret;
}

#else

std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> oo) {
	std::string out = "GDAL version >= 2.1 required for gdalinfo");
	return out;
}

#endif

