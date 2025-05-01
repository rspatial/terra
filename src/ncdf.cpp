
#include <string>
#include "Rcpp.h"


#if !defined(__APPLE__)
// [[Rcpp::export(name = ".ncdf_open")]]
int ncdf_open(std::string filename, bool write) {
	return(-1);
}

// [[Rcpp::export(name = ".ncdf_close")]]
bool ncdf_close(int ncid) {
	return(-1);
}

#else

#include <netcdf.h>

// [[Rcpp::export(name = ".ncdf_open")]]
int ncdf_open(std::string filename, bool write) {
	
	int nc_mode = NC_NOWRITE;
	if (write) {
		nc_mode = NC_WRITE;
	}	
	int ncid = -1;
	int status = NC_NOERR;
	status = nc_open(filename.c_str(), nc_mode, &ncid);
	if (status != NC_NOERR) {
		std::string s = nc_strerror(status);
		Rcpp::Rcout << s << std::endl;
		return false;
	}
	return ncid;
}

// [[Rcpp::export(name = ".ncdf_close")]]
bool ncdf_close(int ncid) {
	int err = nc_close(ncid);
	if (err != NC_NOERR) {
		std::string s = nc_strerror(err);
		Rcpp::Rcout << s << std::endl;
		return false;
	}
	return true;
}

#endif