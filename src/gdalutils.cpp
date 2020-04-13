//#include "spatRaster.h"
#include "cpl_port.h"
#include "cpl_conv.h" // CPLFree()
#include "gdal_version.h"
#include "gdalhelp.h"
#include "ogr_srs_api.h"

#ifdef GDALutils
# include "gdal_utils.h" // requires >= 2.1

// code adapted from the 'sf' package by Edzer Pebesma et al
	
bool gdal_rasterize(std::string src, std::string dst, 
		std::vector<std::string> options, std::vector<std::string> oo, std::vector<std::string> doo,
		bool overwrite = false) {

	int err = 0;
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(oo); // open options
	GDALRasterizeOptions* opt =  GDALRasterizeOptionsNew(options_char.data(), NULL);

	GDALDatasetH src_pt = GDALOpenEx((const char *) src.c_str(), GDAL_OF_VECTOR | GA_ReadOnly, NULL, oo_char.data(), NULL);
	if (src_pt == NULL) {
		return false;
	}
	unset_error_handler();
	GDALDatasetH dst_pt = NULL;
	if (! overwrite) {
		std::vector <char *> doo_char = string_to_charpnt(doo); // open options
		dst_pt = GDALOpenEx((const char *) dst.c_str(), GDAL_OF_RASTER | GA_Update, NULL, doo_char.data(), NULL);
	}
	set_error_handler();
	GDALDatasetH result = GDALRasterize(dst_pt == NULL ? (const char *) dst.c_str() : NULL, dst_pt, src_pt, opt, &err);
	GDALRasterizeOptionsFree(opt);
	if (src_pt != NULL)	GDALClose(src_pt);
	if (result != NULL) GDALClose(result);
	return result == NULL || err;
}



bool CPL_gdalbuildvrt(std::string src, std::string dst,
		std::vector<std::string> options, std::vector<std::string> oo) {

	int err = 0;
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(oo); // open options
	GDALBuildVRTOptions* opt = GDALBuildVRTOptionsNew(options_char.data(), NULL);

	//std::vector<const char *> srcpt(src.size());
	//for (int i = 0; i < src.size(); i++)
	//	srcpt[i] = (const char *) src[i];

	std::vector<GDALDatasetH> srcpt(src.size());
	for (size_t i = 0; i < src.size(); i++) {
		srcpt[i] = GDALOpenEx((const char *) src.c_str(), GDAL_OF_RASTER | GA_ReadOnly, NULL, oo_char.data(), NULL);
	}
	GDALDatasetH result = GDALBuildVRT((const char *) dst.c_str(), src.size(), srcpt.data(), NULL, opt, &err);

	GDALBuildVRTOptionsFree(opt);
	for (size_t i = 0; i < src.size(); i++) {
		GDALClose(srcpt[i]);
	}
	if (result != NULL) GDALClose(result);
	return result == NULL || err;
}


bool gdaldemprocessing(std::string src, std::string dst,
		std::vector<std::string> options, std::vector<std::string> processing, std::vector<std::string> colorfilename, std::vector<std::string> oo) {

	int err = 0;
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(oo); // open options
	GDALDEMProcessingOptions* opt =  GDALDEMProcessingOptionsNew(options_char.data(), NULL);

	GDALDatasetH src_pt = GDALOpenEx((const char *) src.c_str(), GDAL_OF_RASTER | GA_ReadOnly, NULL, oo_char.data(), NULL);
	if (src_pt == NULL) {
		//Rcpp::stop("cannot open source dataset"); // #nocov
		return false;
	}
	GDALDatasetH result = GDALDEMProcessing((const char *) dst.c_str(), src_pt, 
		processing.size() == 0 ? NULL : (const char *) processing[0].c_str(), 
		colorfilename.size() == 0 ? NULL : (const char *) colorfilename[0].c_str(), 
		opt, &err);
	GDALDEMProcessingOptionsFree(opt);
	if (result != NULL)	GDALClose(result);
	if (src_pt != NULL)	GDALClose(src_pt);
	return result == NULL || err;
}

/*
bool CPL_gdalnearblack(std::vector<std::string> src, std::vector<std::string> dst,
		std::vector<std::string> options, std::vector<std::string> oo, std::vector<std::string> doo) {

	int err = 0;
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(oo); // open options
	std::vector <char *> doo_char = string_to_charpnt(doo); // open options
	GDALNearblackOptions* opt =  GDALNearblackOptionsNew(options_char.data(), NULL);

	// GDALDatasetH src_pt = GDALOpen((const char *) src.c_str(), GA_ReadOnly);
	GDALDatasetH src_pt = GDALOpenEx((const char *) src.c_str(), GDAL_OF_RASTER | GA_ReadOnly, NULL, oo_char.data(), NULL);
	GDALDatasetH dst_pt = GDALOpenEx((const char *) dst.c_str(), GDAL_OF_RASTER | GA_Update, NULL, doo_char.data(), NULL);
	GDALDatasetH result = GDALNearblack(dst_pt == NULL ? (const char *) dst.c_str() : NULL, dst_pt, src_pt, opt, &err);
	GDALNearblackOptionsFree(opt);
	if (src_pt != NULL) 
		GDALClose(src_pt);
	if (result != NULL)
		GDALClose(result);
	return result == NULL || err;
}

bool CPL_gdalgrid(std::vector<std::string> src, std::vector<std::string> dst,
		std::vector<std::string> options, std::vector<std::string> oo) {

	int err = 0;
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(oo); // open options
	GDALGridOptions* opt =  GDALGridOptionsNew(options_char.data(), NULL);

	GDALDatasetH src_pt = GDALOpenEx((const char *) src.c_str(), GDAL_OF_ALL | GA_ReadOnly, 
		NULL, oo_char.data(), NULL);
	GDALDatasetH result = GDALGrid((const char *) dst.c_str(), src_pt, opt, &err);
	GDALGridOptionsFree(opt);
	if (src_pt != NULL)
		GDALClose(src_pt);
	if (result != NULL)
		GDALClose(result);
	return result == NULL || err;
}



bool CPL_gdaltranslate(std::vector<std::string> src, std::vector<std::string> dst,
		std::vector<std::string> options, std::vector<std::string> oo) {

	int err = 0;
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(oo);
	GDALTranslateOptions* opt =  GDALTranslateOptionsNew(options_char.data(), NULL);

	GDALDatasetH src_pt = GDALOpenEx((const char *) src.c_str(), GDAL_OF_RASTER | GA_ReadOnly, 
		NULL, oo_char.data(), NULL);
	if (src_pt == NULL)
		return 1; // #nocov
	GDALDatasetH result = GDALTranslate((const char *) dst.c_str(), src_pt, opt, &err);
	if (src_pt != NULL)
		GDALClose(src_pt);
	GDALTranslateOptionsFree(opt);
	if (result != NULL)
		GDALClose(result);
	return result == NULL || err;
}

bool CPL_gdalvectortranslate(std::vector<std::string> src, std::vector<std::string> dst,
		std::vector<std::string> options, std::vector<std::string> oo, std::vector<std::string> doo) {

	int err = 0;
	std::vector <char *> options_char = string_to_charpnt(options);
	GDALVectorTranslateOptions* opt =  GDALVectorTranslateOptionsNew(options_char.data(), NULL);

	std::vector <char *> oo_char = string_to_charpnt(oo); // open options
	GDALDatasetH src_pt = GDALOpenEx((const char *) src.c_str(), GDAL_OF_VECTOR | GA_ReadOnly, NULL, 
		oo_char.data(), NULL);
	if (src_pt == NULL)
		return 1; // #nocov
	std::vector <char *> doo_char = string_to_charpnt(doo); // open options
	unset_error_handler();
	GDALDatasetH dst_pt = GDALOpenEx((const char *) dst.c_str(), GDAL_OF_VECTOR | GA_Update, NULL, doo_char.data(), NULL);
	set_error_handler();
	GDALDatasetH result = 
		GDALVectorTranslate(dst_pt == NULL ? (const char *) dst.c_str() : NULL, dst_pt, 1, &src_pt, opt, &err);
	GDALVectorTranslateOptionsFree(opt);
	GDALClose(src_pt);
	if (result != NULL)
		GDALClose(result);
	return result == NULL || err;
}

*/

#endif