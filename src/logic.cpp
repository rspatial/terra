#include <vector>
#include <algorithm>
#include <functional>
#include "spatraster.h"


template <typename T>
std::vector<T> operator&(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::logical_and<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
    return result;
}


template <typename T>
std::vector<T> operator|(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::logical_or<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
    return result;
}


SpatRaster SpatRaster::isnot(std::string filename, bool overwrite) {
	SpatRaster out = geometry();
  	out.writeStart(filename, overwrite);
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		for (size_t j=0; j<a.size(); j++) {
			a[i] = !a[i];
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();	
	return(out);
}



SpatRaster SpatRaster::logic(SpatRaster x, std::string oper, std::string filename, bool overwrite) {
	
	SpatRaster out = geometry();
	
	std::vector<std::string> f {"&", "|"}; 
	if (std::find(f.begin(), f.end(), oper) == f.end()) {
		out.error = true;
		out.error_message = "unknown logic function";
		return out;
	}

	if (!compare_geom(x, true, false)) {
		out.error = true;
		out.error_message = "dimensions and/or extent do not match";
		return(out);
	}
	
  	out.writeStart(filename, overwrite);
	readStart();
	x.readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		std::vector<double> b = x.readBlock(out.bs, i);
		if (oper == "&") {
			a = a & b; 
		} else if (oper == "|") {
			a = a | b; 
		} else {
			// stop
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();	
	x.readStop();	
	return(out);
}



SpatRaster SpatRaster::logic(bool x, std::string oper, std::string filename, bool overwrite) {

	SpatRaster out = geometry();
	
  	out.writeStart(filename, overwrite);
	readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		if (std::isnan(x)) {
			for(double& d : a)  d = NAN;
		} else if (oper == "&") {
//			for(double& d : a)  d & x;
		} else if (oper == "|") {
//			for(double& d : a)  d | x;
		} else {
			// stop
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();		
	return(out);
}

