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


double dabs(double x) {
	return (x < 0 ? -1 * x : x);
}

SpatRaster SpatRaster::math(std::string fun, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return out;

	std::vector<std::string> f {"abs", "sqrt", "ceiling", "floor", "trunc", "log", "log10", "log2", "log1p", "exp", "expm1", "sign"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown math function");
		return out;
	}

	std::function<double(double)> mathFun;
	if (fun == "sqrt") {
		mathFun = static_cast<double(*)(double)>(sqrt);
	} else if (fun == "abs") {
		mathFun = dabs;
	} else if (fun == "log") {
		mathFun = static_cast<double(*)(double)>(log);
	} else if (fun == "log2") {
		mathFun = static_cast<double(*)(double)>(log2);
	} else if (fun == "log10") {
		mathFun = static_cast<double(*)(double)>(log10);
	} else if (fun == "log1p") {
		mathFun = static_cast<double(*)(double)>(log1p);
	} else if (fun == "exp") {
		mathFun = static_cast<double(*)(double)>(exp);
	} else if (fun == "expm1") {
		mathFun = static_cast<double(*)(double)>(expm1);
	} else if (fun == "sign") {
		mathFun = sign<double>;
	} else if (fun == "ceiling") {
		mathFun = static_cast<double(*)(double)>(ceil);
	} else if (fun == "floor") {
		mathFun = static_cast<double(*)(double)>(floor);
	} else if (fun == "trunc") {
		mathFun = static_cast<double(*)(double)>(trunc);
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		for(double& d : a) if (!std::isnan(d)) d = mathFun(d);
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



double sin_pi(double &x) {
	return sin(x * M_PI);
}

double cos_pi(double &x) {
	return sin(x * M_PI);
}

double tan_pi(double &x) {
	return sin(x * M_PI);
}



SpatRaster SpatRaster::trig(std::string fun, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) return out;

	std::vector<std::string> f {"acos", "asin", "atan", "cos", "sin", "tan", "acosh", "asinh", "atanh", "cosh", "cospi", "sinh", "sinpi", "tanh", "tanpi"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown trig function");
		return out;
	}

	std::function<double(double&)> trigFun;
	if (fun == "sin") {
		trigFun = static_cast<double(*)(double)>(sin);
	} else if (fun == "cos") {
		trigFun = static_cast<double(*)(double)>(cos);
	} else if (fun == "tan") {
		trigFun = static_cast<double(*)(double)>(tan);
	} else if (fun == "asin") {
		trigFun = static_cast<double(*)(double)>(asin);
	} else if (fun == "acos") {
		trigFun = static_cast<double(*)(double)>(acos);
	} else if (fun == "atan") {
		trigFun = static_cast<double(*)(double)>(atan);
	} else if (fun == "sinh") {
		trigFun = static_cast<double(*)(double)>(sinh);
	} else if (fun == "cosh") {
		trigFun = static_cast<double(*)(double)>(cosh);
	} else if (fun == "tanh") {
		trigFun = static_cast<double(*)(double)>(tanh);
	} else if (fun == "asinh") {
		trigFun = static_cast<double(*)(double)>(asinh);
	} else if (fun == "acosh") {
		trigFun = static_cast<double(*)(double)>(acosh);
	} else if (fun == "atanh") {
		trigFun = static_cast<double(*)(double)>(atanh);
	} else if (fun == "sinpi") {
		trigFun = sin_pi;
	} else if (fun == "cospi") {
		trigFun = cos_pi;
	} else if (fun == "tanpi") {
		trigFun = tan_pi;
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		for (double& d : a) if (!std::isnan(d)) d = trigFun(d);
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

