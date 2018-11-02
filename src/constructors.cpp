#include "spatraster.h"


SpatRaster::SpatRaster(std::string fname) {
	constructFromFile(fname);
}


void SpatRaster::setSources(std::vector<RasterSource> s) {
	source = s;
	nrow = s[0].nrow;
	ncol = s[0].ncol;
	extent = s[0].extent;
	crs = s[0].crs;
	open_write = false;
	open_read = false;
}


void SpatRaster::setSource(RasterSource s) {
	setSources({s});
}


SpatRaster::SpatRaster(RasterSource s) {
	setSources( {s} );
}


SpatRaster::SpatRaster() {

	RasterSource s;
	s.nrow = 10;
	s.ncol = 10;
	s.extent = SpatExtent();
	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = 1;
	s.hasRange = { false };
	s.hasValues = false;
	s.layers.resize(1,1);
	s.datatype = "";
	s.names = {"lyr.1"};
	s.crs = "+proj=longlat +datum=WGS84";

	setSource(s);
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

	setSource(s);
}


SpatRaster::SpatRaster(unsigned _nrow, unsigned _ncol, unsigned _nlyr, SpatExtent ext, std::string _crs) {

	RasterSource s;
	s.nrow = _nrow;
	s.ncol = _ncol;
	s.extent = ext;
	s.hasValues = false;
	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = _nlyr;
	s.hasRange = { false };
	s.layers.resize(1, 1);
	s.datatype = "";
	s.crs=_crs;
	for (unsigned i=0; i < _nlyr; i++) {	s.names.push_back("lyr." + std::to_string(i+1)) ; }
	setSource(s);
}


SpatRaster SpatRaster::geometry(long nlyrs) {
	RasterSource s;
	s.nrow = nrow;
	s.ncol = ncol;
	s.extent = extent;
	s.crs = crs;
	s.memory = true;
	if (nlyrs < 1) {
		s.nlyr = nlyr();
	} else {
		s.nlyr = nlyrs;
	}
	std::vector<std::string> nms(s.nlyr);
	for (size_t i=0; i < s.nlyr; i++) { nms[i] = "lyr" + std::to_string(i+1); }
	s.names = nms;
	SpatRaster out;
	out.setSource(s);
	out.error = false;
	out.warning = false;
	return out;
}


SpatRaster SpatRaster::deepCopy() {

	SpatRaster out = *this;
	out.error = false;
	out.warning = false;
	out.nrow = nrow;
	out.ncol = ncol;
	out.extent = extent;
	out.crs = crs;
	return out;
}

