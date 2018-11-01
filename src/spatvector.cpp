using namespace std;
#include "spatvector.h"
#include "util.h"

/*SpatVector::SpatVector() {
	df.names.push_back("test1");
	df.names.push_back("test2");
	std::vector<double> v = {10,11.5,12.5,13.5,15.5};
	std::vector<long> w = {1,2,3,4,5};
	df.dv.push_back(v);
	df.itype.push_back(0);
	df.iplace.push_back(0);
	df.iv.push_back(w);
	df.itype.push_back(1);
	df.iplace.push_back(0);
	df.dv.push_back(v);
	df.itype.push_back(0);
	df.iplace.push_back(1);
}
*/

std::vector<double> SpatVector::getDv(unsigned i) {
	unsigned j = df.iplace[i];
	return df.dv[j];
}

std::vector<long> SpatVector::getIv(unsigned i){
	unsigned j = df.iplace[i];
	return df.iv[j];
}

std::vector<string> SpatVector::getSv(unsigned i){
	unsigned j = df.iplace[i];
	return df.sv[j];
}

std::vector<unsigned> SpatVector::getItype(){
	return df.itype;
}

std::vector<unsigned> SpatVector::getIplace(){
	return df.iplace;
}

std::vector<string> SpatVector::names(){
	return df.names;
}

unsigned SpatVector::ncol() {
	return df.ncol();
}

unsigned SpatVector::nrow() {
// should use the geoms instead for when there are no atts
//	return size();
	return df.nrow();
}


SpatExtent SpatVector::extent(){
	if (gtype == 0) {
		return pos.extent;
	} else if (gtype == 1) {
		return pts.extent;
	} else if (gtype == 2) {
		return lns.extent;		
	} else {
		SpatExtent e;
		return e;
	}
}

std::string SpatVector::getCRS(){
	if (gtype == 0) {
		return pos.crs;
	} else if (gtype == 1) {
		return pts.crs;
	} else if (gtype == 2) {
		return lns.crs;		
	} else {
		return "?";
	}
}

void SpatVector::setCRS(std::string crs){
	lrtrim(crs);
	if (gtype == 0) {
		pos.crs = crs;
	} else if (gtype == 1) {
		pts.crs = crs;
	} else if (gtype == 2) {
		lns.crs = crs;
	} else {
		error=true;
		error_message = "unknown geometry type";
	}
}


string SpatVector::type(){
	if (gtype == 0) {
		return "polygons";
	} else if (gtype == 1) {
		return "points";
	} else if (gtype == 2) {
		return "lines";		
	} else {
		return("?");
	}
}

