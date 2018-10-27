using namespace std;
#include <vector>
#include "spat.h"

#include <algorithm>
#include <functional>

SpatRaster SpatRaster::math(std::string fun, std::string filename, bool overwrite) {

	SpatRaster out = geometry();
	
  	out.writeStart(filename, overwrite);
	readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		if (fun == "sqrt") {
			for(double& d : a) if (!std::isnan(d)) d = sqrt(d);
		} else if (fun == "abs") {
			for(double& d : a) if (!std::isnan(d)) d = abs(d);
		} else if (fun == "log") {
			for(double& d : a) if (!std::isnan(d)) d = log(d);
		} else if (fun == "log2") {
			for(double& d : a) if (!std::isnan(d)) d = log2(d);
		} else if (fun == "log10") {
			for(double& d : a) if (!std::isnan(d)) d = log10(d);
		} else if (fun == "exp") {
			for(double& d : a) if (!std::isnan(d)) d = exp(d);
//		} else if (fun == "sign") {
//			for(double& d : a) if (!std::isnan(d)) d = sign(d);

// not working, probably because double instead of int
		} else if (fun == "ceiling") {
			for(double& d : a) if (!std::isnan(d)) d = ceil(d);
		} else if (fun == "floor") {
			for(double& d : a) if (!std::isnan(d)) d = floor(d);
		} else if (fun == "trunc") {
			for(double& d : a) if (!std::isnan(d)) d = trunc(d);
		} else {
			// check is too late, should be done before opening files.
			out.error = true;
			out.error_message = "unknown mathematical function";
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();		
	return(out);
}


// "acos", "acosh",, "asinh", "atanh",  "cosh", "cospi",  "sinh", "sinpi",  "tanh", "tanpi", 

SpatRaster SpatRaster::trig(std::string fun, std::string filename, bool overwrite) {

	SpatRaster out = geometry();
	
  	out.writeStart(filename, overwrite);
	readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol, 0, nlyr());
		if (fun == "sin") {
			for(double& d : a) if (!std::isnan(d)) d = sin(d);
		} else if (fun == "cos") {
			for(double& d : a) if (!std::isnan(d)) d = cos(d);
		} else if (fun == "tan") {
			for(double& d : a) if (!std::isnan(d)) d = tan(d);
		} else if (fun == "asin") {
			for(double& d : a) if (!std::isnan(d)) d = asin(d);
		} else if (fun == "acos") {
			for(double& d : a) if (!std::isnan(d)) d = acos(d);
		} else if (fun == "atan") {
			for(double& d : a) if (!std::isnan(d)) d = atan(d);
		} else {
			// check is too late, should be done before opening files.
			out.error = true;
			out.error_message = "unknown mathematical function";
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();		
	return(out);
}

