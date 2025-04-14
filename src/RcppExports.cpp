// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// have_TBB
bool have_TBB();
RcppExport SEXP _terra_have_TBB() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(have_TBB());
    return rcpp_result_gen;
END_RCPP
}
// proj_version
std::string proj_version();
RcppExport SEXP _terra_proj_version() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(proj_version());
    return rcpp_result_gen;
END_RCPP
}
// hex2rgb
std::vector<unsigned char> hex2rgb(std::string s);
RcppExport SEXP _terra_hex2rgb(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(hex2rgb(s));
    return rcpp_result_gen;
END_RCPP
}
// rgb2hex
std::string rgb2hex(std::vector<unsigned char> x);
RcppExport SEXP _terra_rgb2hex(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<unsigned char> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rgb2hex(x));
    return rcpp_result_gen;
END_RCPP
}
// sameSRS
bool sameSRS(std::string x, std::string y);
RcppExport SEXP _terra_sameSRS(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(sameSRS(x, y));
    return rcpp_result_gen;
END_RCPP
}
// getCRSname
std::vector<std::string> getCRSname(std::string s);
RcppExport SEXP _terra_getCRSname(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(getCRSname(s));
    return rcpp_result_gen;
END_RCPP
}
// getLinearUnits
double getLinearUnits(std::string s);
RcppExport SEXP _terra_getLinearUnits(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(getLinearUnits(s));
    return rcpp_result_gen;
END_RCPP
}
// geotransform
std::vector<double> geotransform(std::string fname);
RcppExport SEXP _terra_geotransform(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(geotransform(fname));
    return rcpp_result_gen;
END_RCPP
}
// gdal_setconfig
void gdal_setconfig(std::string option, std::string value);
RcppExport SEXP _terra_gdal_setconfig(SEXP optionSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type option(optionSEXP);
    Rcpp::traits::input_parameter< std::string >::type value(valueSEXP);
    gdal_setconfig(option, value);
    return R_NilValue;
END_RCPP
}
// gdal_getconfig
std::string gdal_getconfig(std::string option);
RcppExport SEXP _terra_gdal_getconfig(SEXP optionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type option(optionSEXP);
    rcpp_result_gen = Rcpp::wrap(gdal_getconfig(option));
    return rcpp_result_gen;
END_RCPP
}
// ginfo
std::string ginfo(std::string filename, std::vector<std::string> options, std::vector<std::string> oo);
RcppExport SEXP _terra_ginfo(SEXP filenameSEXP, SEXP optionsSEXP, SEXP ooSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type oo(ooSEXP);
    rcpp_result_gen = Rcpp::wrap(ginfo(filename, options, oo));
    return rcpp_result_gen;
END_RCPP
}
// sd_info
std::vector<std::vector<std::string>> sd_info(std::string filename);
RcppExport SEXP _terra_sd_info(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(sd_info(filename));
    return rcpp_result_gen;
END_RCPP
}
// gdal_version
std::string gdal_version();
RcppExport SEXP _terra_gdal_version() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(gdal_version());
    return rcpp_result_gen;
END_RCPP
}
// geos_version
std::string geos_version(bool runtime, bool capi);
RcppExport SEXP _terra_geos_version(SEXP runtimeSEXP, SEXP capiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type runtime(runtimeSEXP);
    Rcpp::traits::input_parameter< bool >::type capi(capiSEXP);
    rcpp_result_gen = Rcpp::wrap(geos_version(runtime, capi));
    return rcpp_result_gen;
END_RCPP
}
// metatdata
std::vector<std::string> metatdata(std::string filename);
RcppExport SEXP _terra_metatdata(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(metatdata(filename));
    return rcpp_result_gen;
END_RCPP
}
// sdsmetatdata
std::vector<std::string> sdsmetatdata(std::string filename);
RcppExport SEXP _terra_sdsmetatdata(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(sdsmetatdata(filename));
    return rcpp_result_gen;
END_RCPP
}
// sdsmetatdataparsed
std::vector<std::vector<std::string>> sdsmetatdataparsed(std::string filename);
RcppExport SEXP _terra_sdsmetatdataparsed(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(sdsmetatdataparsed(filename));
    return rcpp_result_gen;
END_RCPP
}
// gdal_drivers
std::vector<std::vector<std::string>> gdal_drivers();
RcppExport SEXP _terra_gdal_drivers() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(gdal_drivers());
    return rcpp_result_gen;
END_RCPP
}
// set_gdal_warnings
void set_gdal_warnings(int level);
RcppExport SEXP _terra_set_gdal_warnings(SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type level(levelSEXP);
    set_gdal_warnings(level);
    return R_NilValue;
END_RCPP
}
// seed_init
void seed_init(uint32_t seed_val);
RcppExport SEXP _terra_seed_init(SEXP seed_valSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint32_t >::type seed_val(seed_valSEXP);
    seed_init(seed_val);
    return R_NilValue;
END_RCPP
}
// gdal_init
void gdal_init(std::string projpath, std::string datapath);
RcppExport SEXP _terra_gdal_init(SEXP projpathSEXP, SEXP datapathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type projpath(projpathSEXP);
    Rcpp::traits::input_parameter< std::string >::type datapath(datapathSEXP);
    gdal_init(projpath, datapath);
    return R_NilValue;
END_RCPP
}
// percRank
std::vector<double> percRank(std::vector<double> x, std::vector<double> y, double minc, double maxc, int tail);
RcppExport SEXP _terra_percRank(SEXP xSEXP, SEXP ySEXP, SEXP mincSEXP, SEXP maxcSEXP, SEXP tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type minc(mincSEXP);
    Rcpp::traits::input_parameter< double >::type maxc(maxcSEXP);
    Rcpp::traits::input_parameter< int >::type tail(tailSEXP);
    rcpp_result_gen = Rcpp::wrap(percRank(x, y, minc, maxc, tail));
    return rcpp_result_gen;
END_RCPP
}
// clearVSIcache
void clearVSIcache(bool vsi);
RcppExport SEXP _terra_clearVSIcache(SEXP vsiSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type vsi(vsiSEXP);
    clearVSIcache(vsi);
    return R_NilValue;
END_RCPP
}
// setGDALCacheSizeMB
void setGDALCacheSizeMB(double x, bool vsi);
RcppExport SEXP _terra_setGDALCacheSizeMB(SEXP xSEXP, SEXP vsiSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type vsi(vsiSEXP);
    setGDALCacheSizeMB(x, vsi);
    return R_NilValue;
END_RCPP
}
// getGDALCacheSizeMB
double getGDALCacheSizeMB(bool vsi);
RcppExport SEXP _terra_getGDALCacheSizeMB(SEXP vsiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type vsi(vsiSEXP);
    rcpp_result_gen = Rcpp::wrap(getGDALCacheSizeMB(vsi));
    return rcpp_result_gen;
END_RCPP
}
// get_proj_search_paths
std::vector<std::string> get_proj_search_paths();
RcppExport SEXP _terra_get_proj_search_paths() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_proj_search_paths());
    return rcpp_result_gen;
END_RCPP
}
// set_proj_search_paths
bool set_proj_search_paths(std::vector<std::string> paths);
RcppExport SEXP _terra_set_proj_search_paths(SEXP pathsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type paths(pathsSEXP);
    rcpp_result_gen = Rcpp::wrap(set_proj_search_paths(paths));
    return rcpp_result_gen;
END_RCPP
}
// PROJ_network
std::string PROJ_network(bool enable, std::string url);
RcppExport SEXP _terra_PROJ_network(SEXP enableSEXP, SEXP urlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type enable(enableSEXP);
    Rcpp::traits::input_parameter< std::string >::type url(urlSEXP);
    rcpp_result_gen = Rcpp::wrap(PROJ_network(enable, url));
    return rcpp_result_gen;
END_RCPP
}
// pearson_cor
double pearson_cor(std::vector<double> x, std::vector<double> y, bool narm);
RcppExport SEXP _terra_pearson_cor(SEXP xSEXP, SEXP ySEXP, SEXP narmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type narm(narmSEXP);
    rcpp_result_gen = Rcpp::wrap(pearson_cor(x, y, narm));
    return rcpp_result_gen;
END_RCPP
}
// weighted_pearson_cor
double weighted_pearson_cor(std::vector<double> x, std::vector<double> y, std::vector<double> weights, bool narm);
RcppExport SEXP _terra_weighted_pearson_cor(SEXP xSEXP, SEXP ySEXP, SEXP weightsSEXP, SEXP narmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type narm(narmSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_pearson_cor(x, y, weights, narm));
    return rcpp_result_gen;
END_RCPP
}
// uniqueSymmetricRows
Rcpp::IntegerMatrix uniqueSymmetricRows(std::vector<size_t> x, std::vector<size_t> y);
RcppExport SEXP _terra_uniqueSymmetricRows(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<size_t> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<size_t> >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(uniqueSymmetricRows(x, y));
    return rcpp_result_gen;
END_RCPP
}
// dist2segmentPoint_geo
double dist2segmentPoint_geo(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double& ilon, double& ilat);
RcppExport SEXP _terra_dist2segmentPoint_geo(SEXP plonSEXP, SEXP platSEXP, SEXP lon1SEXP, SEXP lat1SEXP, SEXP lon2SEXP, SEXP lat2SEXP, SEXP ilonSEXP, SEXP ilatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type plon(plonSEXP);
    Rcpp::traits::input_parameter< double >::type plat(platSEXP);
    Rcpp::traits::input_parameter< double >::type lon1(lon1SEXP);
    Rcpp::traits::input_parameter< double >::type lat1(lat1SEXP);
    Rcpp::traits::input_parameter< double >::type lon2(lon2SEXP);
    Rcpp::traits::input_parameter< double >::type lat2(lat2SEXP);
    Rcpp::traits::input_parameter< double& >::type ilon(ilonSEXP);
    Rcpp::traits::input_parameter< double& >::type ilat(ilatSEXP);
    rcpp_result_gen = Rcpp::wrap(dist2segmentPoint_geo(plon, plat, lon1, lat1, lon2, lat2, ilon, ilat));
    return rcpp_result_gen;
END_RCPP
}
// intermediate
std::vector<std::vector<double>> intermediate(double lon1, double lat1, double lon2, double lat2, int n, double distance);
RcppExport SEXP _terra_intermediate(SEXP lon1SEXP, SEXP lat1SEXP, SEXP lon2SEXP, SEXP lat2SEXP, SEXP nSEXP, SEXP distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lon1(lon1SEXP);
    Rcpp::traits::input_parameter< double >::type lat1(lat1SEXP);
    Rcpp::traits::input_parameter< double >::type lon2(lon2SEXP);
    Rcpp::traits::input_parameter< double >::type lat2(lat2SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type distance(distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(intermediate(lon1, lat1, lon2, lat2, n, distance));
    return rcpp_result_gen;
END_RCPP
}
// ncdf_open
int ncdf_open(std::string filename, bool write);
RcppExport SEXP _terra_ncdf_open(SEXP filenameSEXP, SEXP writeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type write(writeSEXP);
    rcpp_result_gen = Rcpp::wrap(ncdf_open(filename, write));
    return rcpp_result_gen;
END_RCPP
}
// ncdf_close
bool ncdf_close(int ncid);
RcppExport SEXP _terra_ncdf_close(SEXP ncidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ncid(ncidSEXP);
    rcpp_result_gen = Rcpp::wrap(ncdf_close(ncid));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_spat();

static const R_CallMethodDef CallEntries[] = {
    {"_terra_have_TBB", (DL_FUNC) &_terra_have_TBB, 0},
    {"_terra_proj_version", (DL_FUNC) &_terra_proj_version, 0},
    {"_terra_hex2rgb", (DL_FUNC) &_terra_hex2rgb, 1},
    {"_terra_rgb2hex", (DL_FUNC) &_terra_rgb2hex, 1},
    {"_terra_sameSRS", (DL_FUNC) &_terra_sameSRS, 2},
    {"_terra_getCRSname", (DL_FUNC) &_terra_getCRSname, 1},
    {"_terra_getLinearUnits", (DL_FUNC) &_terra_getLinearUnits, 1},
    {"_terra_geotransform", (DL_FUNC) &_terra_geotransform, 1},
    {"_terra_gdal_setconfig", (DL_FUNC) &_terra_gdal_setconfig, 2},
    {"_terra_gdal_getconfig", (DL_FUNC) &_terra_gdal_getconfig, 1},
    {"_terra_ginfo", (DL_FUNC) &_terra_ginfo, 3},
    {"_terra_sd_info", (DL_FUNC) &_terra_sd_info, 1},
    {"_terra_gdal_version", (DL_FUNC) &_terra_gdal_version, 0},
    {"_terra_geos_version", (DL_FUNC) &_terra_geos_version, 2},
    {"_terra_metatdata", (DL_FUNC) &_terra_metatdata, 1},
    {"_terra_sdsmetatdata", (DL_FUNC) &_terra_sdsmetatdata, 1},
    {"_terra_sdsmetatdataparsed", (DL_FUNC) &_terra_sdsmetatdataparsed, 1},
    {"_terra_gdal_drivers", (DL_FUNC) &_terra_gdal_drivers, 0},
    {"_terra_set_gdal_warnings", (DL_FUNC) &_terra_set_gdal_warnings, 1},
    {"_terra_seed_init", (DL_FUNC) &_terra_seed_init, 1},
    {"_terra_gdal_init", (DL_FUNC) &_terra_gdal_init, 2},
    {"_terra_percRank", (DL_FUNC) &_terra_percRank, 5},
    {"_terra_clearVSIcache", (DL_FUNC) &_terra_clearVSIcache, 1},
    {"_terra_setGDALCacheSizeMB", (DL_FUNC) &_terra_setGDALCacheSizeMB, 2},
    {"_terra_getGDALCacheSizeMB", (DL_FUNC) &_terra_getGDALCacheSizeMB, 1},
    {"_terra_get_proj_search_paths", (DL_FUNC) &_terra_get_proj_search_paths, 0},
    {"_terra_set_proj_search_paths", (DL_FUNC) &_terra_set_proj_search_paths, 1},
    {"_terra_PROJ_network", (DL_FUNC) &_terra_PROJ_network, 2},
    {"_terra_pearson_cor", (DL_FUNC) &_terra_pearson_cor, 3},
    {"_terra_weighted_pearson_cor", (DL_FUNC) &_terra_weighted_pearson_cor, 4},
    {"_terra_uniqueSymmetricRows", (DL_FUNC) &_terra_uniqueSymmetricRows, 2},
    {"_terra_dist2segmentPoint_geo", (DL_FUNC) &_terra_dist2segmentPoint_geo, 8},
    {"_terra_intermediate", (DL_FUNC) &_terra_intermediate, 6},
    {"_terra_ncdf_open", (DL_FUNC) &_terra_ncdf_open, 2},
    {"_terra_ncdf_close", (DL_FUNC) &_terra_ncdf_close, 1},
    {"_rcpp_module_boot_spat", (DL_FUNC) &_rcpp_module_boot_spat, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_terra(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
