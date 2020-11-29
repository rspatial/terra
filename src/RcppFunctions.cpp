#include <Rcpp.h>
//#include "spatRaster.h"
#include "spatRasterMultiple.h"

//#include <memory> //std::addressof
#include "gdal_priv.h"
# include "gdal_info.h"

#if GDAL_VERSION_MAJOR >= 3
#include "proj.h"
#endif



// [[Rcpp::export(name = ".geotransform")]]
std::vector<double> geotransform(std::string fname) {
	std::vector<double> out;
    GDALDataset *poDataset;
    poDataset = (GDALDataset *) GDALOpen(fname.c_str(), GA_ReadOnly );

    if( poDataset == NULL )  {
		Rcpp::Rcout << "cannot read from: " + fname << std::endl;
		return out;
	}

	double gt[6];	
	if( poDataset->GetGeoTransform( gt ) != CE_None ) {
		Rcpp::Rcout << "bad" << std::endl;
	}
	out = std::vector<double>(std::begin(gt), std::end(gt));
	GDALClose( (GDALDatasetH) poDataset );

	return out;
}


// [[Rcpp::export(name = ".gdalinfo")]]
std::string ginfo(std::string filename, std::vector<std::string> options, std::vector<std::string> oo) {
	std::string out = gdalinfo(filename, options, oo);
	return out;
}

// [[Rcpp::export(name = ".sdinfo")]]
std::vector<std::vector<std::string>> sd_info(std::string filename) {
	std::vector<std::vector<std::string>> sd = sdinfo(filename);
	return sd;
}

// [[Rcpp::export(name = ".gdalversion")]]
std::string gdal_version() {
	const char* what = "RELEASE_NAME";
	const char* x = GDALVersionInfo(what);
	std::string s = (std::string) x;
	return s;
}

// [[Rcpp::export(name = ".metadata")]]
std::vector<std::string> metatdata(std::string filename) {
	std::vector<std::string> m = get_metadata(filename);
	return m;
}

// [[Rcpp::export(name = ".sdsmetadata")]]
std::vector<std::string> sdsmetatdata(std::string filename) {
	std::vector<std::string> m = get_metadata_sds(filename);
	return m;
}

// [[Rcpp::export(name = ".parsedsdsmetadata")]]
std::vector<std::vector<std::string>> sdsmetatdataparsed(std::string filename) {
	std::vector<std::string> m = sdsmetatdata(filename);
	std::vector<std::vector<std::string>> s = parse_metadata_sds(m);
	return s;
}

template <typename... Args>
inline void warningNoCall(const char* fmt, Args&&... args ) {
    Rf_warningcall(R_NilValue, tfm::format(fmt, std::forward<Args>(args)... ).c_str());
}

template <typename... Args>
inline void NORET stopNoCall(const char* fmt, Args&&... args) {
    throw Rcpp::exception(tfm::format(fmt, std::forward<Args>(args)... ).c_str(), false);
}

static void __err_warning(CPLErr eErrClass, int err_no, const char *msg) {
	switch ( eErrClass ) {
        case 0:
            break; 
        case 1:
        case 2:
            warningNoCall("%s (GDAL %d)", msg, err_no); 
            break; 
        case 3:
            warningNoCall("%s (GDAL error %d)", msg, err_no); 
            break;
        case 4:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no); 
            break;
        default:
            warningNoCall("%s (GDAL error class %d, #%d)", msg, eErrClass, err_no); 
            break; 
    }
    return;	
}

static void __err_error(CPLErr eErrClass, int err_no, const char *msg) {
	switch ( eErrClass ) {
        case 0:
        case 1:
        case 2:
            break; 		
        case 3:
            warningNoCall("%s (GDAL error %d)", msg, err_no); 
            break;
        case 4:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no); 
            break;
        default:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no); 
            break; 
    }
    return;
}


static void __err_fatal(CPLErr eErrClass, int err_no, const char *msg) {
	switch ( eErrClass ) {
        case 0:
        case 1:
        case 2:
        case 3:
            break;
        case 4:
            stopNoCall("%s (GDAL unrecoverable error %d)", msg, err_no); 
            break;
        default:
            break; 
    }
    return;
}

static void __err_none(CPLErr eErrClass, int err_no, const char *msg) {
    return;
}



// [[Rcpp::export(name = ".set_gdal_warnings")]]
void set_gdal_warnings(int level) {
	if (level==4) {
		CPLSetErrorHandler((CPLErrorHandler)__err_none);
	} else if (level==1) {
		CPLSetErrorHandler((CPLErrorHandler)__err_warning);			
	} else if (level==2) {
		CPLSetErrorHandler((CPLErrorHandler)__err_error);			
	} else {
		CPLSetErrorHandler((CPLErrorHandler)__err_fatal);
	} 		
}


// [[Rcpp::export(name = ".gdalinit")]]
void gdal_init(std::string path) {
	set_gdal_warnings(3);
    GDALAllRegister();
    OGRRegisterAll(); 
	//GDALregistred = true;
#if GDAL_VERSION_MAJOR >= 3
	if (path != "") {
		const char *cp = path.c_str();
		proj_context_set_search_paths(PJ_DEFAULT_CTX, 1, &cp);
	}
#endif
}

