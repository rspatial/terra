// Copyright (c) 2018-2019  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <algorithm>
#include <functional>
#include "spatRaster.h"


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


SpatRaster SpatRaster::isnot(SpatOptions &opt) {
	SpatRaster out = geometry();
  	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		for (size_t j=0; j<a.size(); j++) {
			a[i] = !a[i];
		}
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		
	}
	out.writeStop();
	readStop();	
	return(out);
}



SpatRaster SpatRaster::logic(SpatRaster x, std::string oper, SpatOptions &opt) {
	
	SpatRaster out = geometry();
	
	std::vector<std::string> f {"&", "|"}; 
	if (std::find(f.begin(), f.end(), oper) == f.end()) {
		out.setError("unknown logic function");
		return out;
	}

	if (!compare_geom(x, true, false)) {
		out.setError("dimensions and/or extent do not match");
		return(out);
	}
	
  	if (!out.writeStart(opt)) { return out; }
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
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		
	}
	out.writeStop();
	readStop();	
	x.readStop();	
	return(out);
}



SpatRaster SpatRaster::logic(bool x, std::string oper, SpatOptions &opt) {

	SpatRaster out = geometry();
	
  	if (!out.writeStart(opt)) { return out; }
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
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		
	}
	out.writeStop();
	readStop();		
	return(out);
}

