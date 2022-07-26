/*

#include "spatRaster.h"
#include "netcdf.h"


// [[Rcpp::export(name = ".readncdf")]]
bool readncdf(std::string f) {
	int ncid;
	int err = nc_open(f.c_str(), 0, &ncid);
	Rcpp::Rcout << err << std::endl;
	if (err != 0) return false;
	err = nc_close(ncid);
	Rcpp::Rcout << err << std::endl;
	if (err != 0) return false;
	return true;
}

*/