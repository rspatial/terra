#include "spatraster.h"

SpatRaster SpatRaster::trim(unsigned padding, std::string filename, std::string format, std::string datatype, bool overwrite) {

	long nrl = nrow * nlyr();
	long ncl = ncol * nlyr();
	
	std::vector<double> v;
	unsigned r;
	for (r=0; r<nrow; r++) {
		v = readValues(r, 1, 0, ncol, 0, nlyr());
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
			break;
		}
	}
	
	if ( r == nrow) { //stop('only NA values found')
	}
	unsigned firstrow = std::min(std::max(r - padding, unsigned(0)), nrow);

	for (r=nrow-1; r>firstrow; r--) {
		v = readValues(r, 1, 0, ncol, 0, nlyr());
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < ncl) {
			break;
		}
	}
	
	unsigned lastrow = std::max(std::min(r+padding, nrow), unsigned(0));
	
	unsigned tmp;
	if (lastrow < firstrow) { 
		tmp = firstrow;
		firstrow = lastrow;
		lastrow = tmp;
	}
	unsigned c;
	for (c=0; c<ncol; c++) {
		v = readValues(0, nrow, c, 1, 0, nlyr());
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
			break;
		}
	}
	unsigned firstcol = std::min(std::max(c-padding, unsigned(0)), ncol);
	
	
	for (c=ncol-1; c>firstcol; c--) {
		v = readValues(0, nrow, c, 1, 0, nlyr());
		if (std::count_if( v.begin(), v.end(), [](double d) { return std::isnan(d); } ) < nrl) {
			break;
		}
	}
	unsigned lastcol = std::max(std::min(c+padding, ncol), unsigned(0));
	
	if (lastcol < firstcol) { 
		tmp = firstcol;
		firstcol = lastcol;
		lastcol = tmp;
	}
	
	std::vector<double> res = resolution();
	double xr = res[0];
	double yr = res[1];
	SpatExtent e = SpatExtent(xFromCol(firstcol)-0.5*xr, xFromCol(lastcol)+0.5*xr, yFromRow(lastrow)-0.5*yr, yFromRow(firstrow)+0.5*yr);
	
	return( crop(e, "near", filename, format, datatype, overwrite) ) ;
}

