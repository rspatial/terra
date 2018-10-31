using namespace std;
#include "spatraster.h"
#include "util.h"

void SpatRaster::setCRS(std::string _crs) {
	lrtrim(_crs);
	for (size_t i = 0; i < nsrc(); i++) { source[i].crs = _crs; }
	crs = _crs;
}

std::vector<double> SpatRaster::resolution() { 
	return std::vector<double> { (extent.xmax - extent.xmin) / ncol, (extent.ymax - extent.ymin) / nrow };
}

unsigned SpatRaster::nlyr() {
	unsigned x = 0;
	for (size_t i=0; i<source.size(); i++) { x += source[i].nlyr; }
	return(x);
}

std::vector<string> SpatRaster::filenames() {
	std::vector<string> x(source.size());
	for (size_t i=0; i<x.size(); i++) { x[i] = source[i].filename; }
	return(x);
}

std::vector<bool> SpatRaster::inMemory() {
	std::vector<bool> m(source.size());
	for (size_t i=0; i<m.size(); i++) { m[i] = source[i].memory; }
	return(m);
}

