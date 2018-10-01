using namespace std;
#include "spat.h"


SpatRaster::SpatRaster(std::string fname) {
	constructFromFile(fname);	
}


SpatRaster::SpatRaster(RasterSource s) {
source = { s };
	nrow = s.nrow;
	ncol = s.ncol;
	nlyr = s.nlyr;
	extent = s.extent;
	crs = s.crs;
}

		
SpatRaster::SpatRaster() {
	nrow=10; 
	ncol=10; 
	extent = SpatExtent();
	RasterSource s;
	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = 0;
	s.hasRange = { false };
	s.hasValues = false; 
	s.layers.resize(1,1);
	s.datatype = "";
	s.names = {"lyr.1"};
	source = { s };
	setnlyr();
}


SpatRaster::SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string _crs) {

	RasterSource s;
	s.nrow=rcl[0]; 
	s.ncol=rcl[1];
	s.extent.xmin = ext[0];
	s.extent.xmax = ext[1];
	s.extent.ymin = ext[2];
	s.extent.ymax = ext[3];
	s.hasValues = false; 
	s.hasRange = {false};

	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = rcl[2];
	s.layers.resize(1, 1);
	s.datatype = "";
	s.crs =_crs;
	for (unsigned i=0; i < rcl[2]; i++) { s.names.push_back("lyr." + std::to_string(i+1)) ; }
	source = { s };

	setnlyr();	
}


SpatRaster::SpatRaster(unsigned _nrow, unsigned _ncol, unsigned _nlyr, SpatExtent ext, std::string _crs) {
	nrow=_nrow; ncol=_ncol;
	extent = ext;
	hasValues = false; 
	RasterSource s;
	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = _nlyr;
	s.hasRange = { false };
	s.layers.resize(1, 1);
	s.datatype = "";
	s.crs=_crs;
	for (unsigned i=0; i < _nlyr; i++) {	s.names.push_back("lyr." + std::to_string(i+1)) ; }
	source = {s};
	setnlyr();	
}


