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

#include <functional>
#include "spatRaster.h"
#include "math_utils.h"
#include "recycle.h"


template <typename T> int sign(T value) {
    return (T(0) < value) - (value < T(0));
}


SpatRaster SpatRaster::math(std::string fun, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return out;

	std::vector<std::string> f {"abs", "sqrt", "ceiling", "floor", "trunc", "log", "log10", "log2", "log1p", "exp", "expm1", "sign"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown math function");
		return out;
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
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
		} else if (fun == "log1p") {
			for(double& d : a) if (!std::isnan(d)) d = log1p(d);
		} else if (fun == "exp") {
			for(double& d : a) if (!std::isnan(d)) d = exp(d);
		} else if (fun == "expm1") {
			for(double& d : a) if (!std::isnan(d)) d = expm1(d);
		} else if (fun == "sign") {
			for(double& d : a) if (!std::isnan(d)) d = sign(d);
		} else if (fun == "ceiling") {
			for(double& d : a) if (!std::isnan(d)) d = ceil(d);
		} else if (fun == "floor") {
			for(double& d : a) if (!std::isnan(d)) d = floor(d);
		} else if (fun == "trunc") {
			for(double& d : a) if (!std::isnan(d)) d = trunc(d);
		}
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::math2(std::string fun, unsigned digits, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return out;

	std::vector<std::string> f {"round", "signif"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown math2 function");
		return out;
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		if (fun == "round") {
			for(double& d : a) d = roundn(d, digits);
		} else if (fun == "signif") {
			for(double& d : a) if (!std::isnan(d)) d = signif(d, digits);
		} 
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::trig(std::string fun, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return out;

	std::vector<std::string> f {"acos", "asin", "atan", "cos", "sin", "tan", "acosh", "asinh", "atanh", "cosh", "cospi", "sinh", "sinpi", "tanh", "tanpi"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown trig function");
		return out;
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
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
		} else if (fun == "sinh") {
			for(double& d : a) if (!std::isnan(d)) d = sinh(d);
		} else if (fun == "cosh") {
			for(double& d : a) if (!std::isnan(d)) d = cosh(d);
		} else if (fun == "tanh") {
			for(double& d : a) if (!std::isnan(d)) d = tanh(d);
		} else if (fun == "asinh") {
			for(double& d : a) if (!std::isnan(d)) d = asinh(d);
		} else if (fun == "acosh") {
			for(double& d : a) if (!std::isnan(d)) d = acosh(d);
		} else if (fun == "atanh") {
			for(double& d : a) if (!std::isnan(d)) d = atanh(d);
		} else if (fun == "sinpi") {
			for(double& d : a) if (!std::isnan(d)) d = sin(d * M_PI);
		} else if (fun == "cospi") {
			for(double& d : a) if (!std::isnan(d)) d = cos(d * M_PI);
		} else if (fun == "tanpi") {
			for(double& d : a) if (!std::isnan(d)) d = tan(d * M_PI);
		}
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		
	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::atan_2(SpatRaster x, SpatOptions &opt) {
	SpatRaster out = geometry();
	if (!hasValues()) return out;
  	if (!out.writeStart(opt)) { return out; }
	readStart();
	x.readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		std::vector<double> b = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(a, b);
		std::vector<double> d(a.size());
		for (size_t i=0; i<a.size(); i++) {
			if (std::isnan(a[i]) || std::isnan(b[i])) {
				d[i] = NAN;
			} else {
				d[i] = atan2(a[i], b[i]);
			}
		}
		if (!out.writeValues(d, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}

