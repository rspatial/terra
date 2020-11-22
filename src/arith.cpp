// Copyright (c) 2018-2020  Robert J. Hijmans
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
#include "spatRasterMultiple.h"
#include "recycle.h"
#include "math_utils.h"
#include "vecmathfun.h"
#include "modal.h"


// need to take care of NAs here. OK for NAN, but not for int types
template <typename T>
void operator+(std::vector<T>& a, const std::vector<T>& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<T>());
}


template <typename T>
void operator-(std::vector<T>& a, const std::vector<T>& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<T>());
}


template <typename T>
void operator/(std::vector<T>& a, const std::vector<T>& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::divides<T>());
}

template <typename T>
void operator*(std::vector<T>& a, const std::vector<T>& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<T>());
}



template <typename T>
void operator%(std::vector<T>& a, const std::vector<T>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::modulus<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = std::fmod(a[i], b[i]);
		}
	}
}


template <typename T>
void operator==(std::vector<T>& a, const std::vector<T>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::equal_to<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] == b[i];
		}
	}
}

template <typename T>
void operator!=(std::vector<T>& a, const std::vector<T>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::not_equal_to<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] != b[i];
		}
	}
}

template <typename T>
void operator>=(std::vector<T>& a, const std::vector<T>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::greater_equal<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] >= b[i];
		}
	}
}

template <typename T>
void operator<=(std::vector<T>& a, const std::vector<T>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::less_equal<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] <= b[i];
		}
	}
 }


template <typename T>
void operator>(std::vector<T>& a, const std::vector<T>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::greater<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] > b[i];
		}
	}
}

template <typename T>
void operator<(std::vector<T>& a, const std::vector<T>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::less<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] < b[i];
		}
	}
}


template <typename T>
void power(std::vector<T>& a, const std::vector<T>& b) {
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = std::pow(a[i], b[i]);
		}
	}
}


bool smooth_operator(std::string oper) {
	std::vector<std::string> f {"+", "-", "*", "^", "/", "%", "==", "!=", ">", "<", ">=", "<="};
	return (std::find(f.begin(), f.end(), oper) != f.end());
}



SpatRaster SpatRaster::arith(SpatRaster x, std::string oper, SpatOptions &opt) {

	size_t nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);


	if (!smooth_operator(oper)) {
		out.setError("unknown arith function");
		return out;
	}

	if (!out.compare_geom(x, false, true)) {
		return(out);
	}
	
	if (!(hasValues() & x.hasValues())) {
		out.setError("raster has no values"); // or warn and treat as NA?
		return out;
	}
 	if (!out.writeStart(opt)) { return out; }
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		std::vector<double> b = x.readBlock(out.bs, i);
		recycle(a,b);
		if (oper == "+") {
			a + b;
		} else if (oper == "-") {
			a - b;
		} else if (oper == "*") {
			a * b;
		} else if (oper == "/") {
			a / b;
		} else if (oper == "^") {
			power(a, b);
		} else if (oper == "%") {
			 a % b;
		} else if (oper == "==") {
			a == b;
		} else if (oper == "!=") {
			a == b;
		} else if (oper == ">=") {
			a >= b;
		} else if (oper == "<=") {
			a <= b;
		} else if (oper == ">") {
			a > b;
		} else if (oper == "<") {
			a < b;
		}
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;

	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}



SpatRaster SpatRaster::arith(double x, std::string oper, bool reverse, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr());
	if (!smooth_operator(oper)) {
		out.setError("unknown arith function");
		return out;
	}	
	if (!hasValues()) {
		out.setError("raster has no values"); // or warn and treat as NA?
		return out;
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
  	if (!out.writeStart(opt)) { return out; }
	
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		if (std::isnan(x)) {
			for(double& d : a)  d = NAN;
		} else if (oper == "+") {
			for(double& d : a)  d += x;
		} else if (oper == "-") {
			if (reverse) {
				for(double& d : a)  d = x - d;
			} else {
				for(double& d : a)  d -= x;
			}
		} 
		else if (oper == "*") {
			for(double& d : a)  d *= x;
		} else if (oper == "/") {
			if (reverse) {
				for(double& d : a)  d = x / d;
			} else {
				for(double& d : a)  d /= x;
			}
		} else if (oper == "^") {
			if (reverse) {
				for(double& d : a)  d = std::pow(x, d);				
			} else {
				for(double& d : a)  d = std::pow(d, x);
			}
		} else if (oper == "%") {
			if (reverse) {
				for (size_t i=0; i<a.size(); i++) {
					a[i] = std::fmod(x, a[i]);
				}
			} else {
				for (size_t i=0; i<a.size(); i++) {
					a[i] = std::fmod(a[i], x);
				}
			}
		} else if (oper == "==") {
			for(double& d : a) if (!std::isnan(d)) d = d == x;
		} else if (oper == "!=") {
			for(double& d : a) if (!std::isnan(d)) d = d != x;
		} else if (oper == ">=") {
			if (reverse) {
				for(double& d : a) if (!std::isnan(d)) d = x >= d;
			} else {
				for(double& d : a) if (!std::isnan(d)) d = d >= x;
			}
		} else if (oper == "<=") {
			if (reverse) {
				for(double& d : a) if (!std::isnan(d)) d = x <= d;				
			} else {
				for(double& d : a) if (!std::isnan(d)) d = d <= x;
			}
		} else if (oper == ">") {
			if (reverse) {
				for(double& d : a) if (!std::isnan(d)) d = x > d;				
			} else {
				for(double& d : a) if (!std::isnan(d)) d = d > x;
			}
		} else if (oper == "<") {
			if (reverse) {
				for(double& d : a) if (!std::isnan(d)) d = x < d;				
			} else {
				for(double& d : a) if (!std::isnan(d)) d = d < x;
			}
		} else {
			// stop
		}
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::arith(std::vector<double> x, std::string oper, bool reverse, SpatOptions &opt) {

	if (x.size() == 1) {
		return(arith(x[0], oper, reverse, opt));	
	}
	
	SpatRaster out = geometry();
	if (!smooth_operator(oper)) {
		out.setError("unknown arith function");
		return out;
	}
	if (!hasValues()) {
		out.setError("raster has no values"); // or warn and treat as NA?
		return out;
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

		
  	if (!out.writeStart(opt)) { return out; }

	unsigned nl = nlyr();
	unsigned nc = ncol();
	recycle(x, nlyr());
	
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		std::vector<double> vv;
		unsigned off = out.bs.nrows[i] * nc;
		for (size_t j=0; j<nl; j++) {
			unsigned s = j * off;
			std::vector<double> a(v.begin()+s, v.begin()+s+off);
			if (std::isnan(x[j])) {
				for(double& d : a) d = NAN;
			} else if (oper == "+") {
				for(double& d : a) d += x[j];
			} else if (oper == "-") {
				if (reverse) {
					for(double& d : a) d = x[j] - d;
				} else {
					for(double& d : a) d -= x[j];
				}
			} else if (oper == "*") {
				for(double& d : a)  d *= x[j];
			} else if (oper == "/") {
				if (reverse) {
					for(double& d : a) d = x[j] / d;
				} else {
					for(double& d : a) d /= x[j];
				}
			} else if (oper == "^") {
				if (reverse) {
					for(double& d : a)  d = std::pow(x[j], d);				
				} else {
					for(double& d : a)  d = std::pow(d, x[j]);
				}
			} else if (oper == "%") {
				if (reverse) {
					for (size_t k=0; k<a.size(); k++) {
						a[k] = std::fmod(x[j], a[k]);
					}
				} else {
					for (size_t k=0; k<a.size(); k++) {
						a[k] = std::fmod(a[k], x[j]);
					}
				}
			} else if (oper == "==") {
				for(double& d : a) if (!std::isnan(d)) d = d == x[j];
			} else if (oper == "!=") {
				for(double& d : a) if (!std::isnan(d)) d = d != x[j];
			} else if (oper == ">=") {
				if (reverse) {
					for(double& d : a) if (!std::isnan(d)) d = x[j] >= d;				
				} else {
					for(double& d : a) if (!std::isnan(d)) d = d >= x[j];
				}
			} else if (oper == "<=") {
				if (reverse) {
					for(double& d : a) if (!std::isnan(d)) d = x[j] <= d;				
				} else {
					for(double& d : a) if (!std::isnan(d)) d = d <= x[j];
				}
			} else if (oper == ">") {
				if (reverse) {
					for(double& d : a) if (!std::isnan(d)) d = x[j] > d;				
				} else {
					for(double& d : a) if (!std::isnan(d)) d = d > x[j];
				}
			} else if (oper == "<") {
				if (reverse) {
					for(double& d : a) if (!std::isnan(d)) d = x[j] < d;				
				} else {
					for(double& d : a) if (!std::isnan(d)) d = d < x[j];
				}
			} else {
				// stop
			}
			vv.insert(vv.end(), a.begin(), a.end());
		}
		if (!out.writeValues(vv, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}





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

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
		
  	if (!out.writeStart(opt)) { return out; }
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

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
		
  	if (!out.writeStart(opt)) { return out; }
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

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
		
  	if (!out.writeStart(opt)) { return out; }
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
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
		
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}	
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
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}	
  	if (!out.writeStart(opt)) { return out; }
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

	if (!out.compare_geom(x, true, false)) {
		return(out);
	}
	
 	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}
	
		
 	if (!out.writeStart(opt)) { return out; }
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
	
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
		
  	if (!out.writeStart(opt)) { return out; }
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		if (std::isnan(x)) {
			for(double& d : a)  d = NAN;
		} else if (oper == "&") {
			for(double& d : a)  d = (d==1) & x;
		} else if (oper == "|") {
			for(double& d : a)  d = (d==1) | x;
		} else {
			out.setError("?");
			return out;
		}
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();		
	return(out);
}


SpatRaster SpatRaster::cum(std::string fun, bool narm, SpatOptions &opt) {

	SpatRaster out = geometry();

	std::vector<std::string> f {"sum", "prod", "min", "max"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown cum function");
		return out;
	}
	if (!hasValues()) {
	//	out.setError("raster has no values");
		return out;
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	
  	if (!out.writeStart(opt)) { return out; }
	unsigned nl = out.nlyr();
	std::vector<double> v(nl);
	unsigned nc;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		nc = out.bs.nrows[i] * out.ncol();
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}
			if (fun == "sum") {
				cumsum(v, narm);
			} else if (fun == "prod") {
				cumprod(v, narm);
			} else if (fun == "min") {
				cummin(v, narm);
			} else if (fun == "max") {
				cummax(v, narm);
			}
			for (size_t k=0; k<v.size(); k++) {
				a[j+k*nc] = v[k];
			}
		}
		if (!out.writeValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;

	}
	out.writeStop();
	readStop();
	return(out);
}


double vstdev(std::vector<double> v, bool narm) {
	double m = vmean(v, narm);
	for (double& d : v) d = pow(d - m, 2);
	m = vmean(v, narm);
	return sqrt(m);
}
	


SpatRaster SpatRaster::summary_numb(std::string fun, std::vector<double> add, bool narm, SpatOptions &opt) {

	SpatRaster out = geometry(1);

	std::vector<std::string> f {"sum", "mean", "median", "which.min", "which.max", "min", "max", "range", "prod", "any", "all", "stdev"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown summary function");
		return out;
	}

	if (fun == "range") {
		return range(add, narm, opt);
	} 
	out.source[0].names[0] = fun;
  	if (!hasValues()) { return out; }


	std::function<double(std::vector<double>&, bool)> sumFun;
	if (fun == "stdev") {
		sumFun = vstdev;
	} else {
		sumFun = getFun(fun);
	}
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	
	if (!out.writeStart(opt)) { return out; }
	unsigned nl = nlyr();
	std::vector<double> v(nl);
	if (add.size() > 0) v.insert( v.end(), add.begin(), add.end() );

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}
			b[j] = sumFun(v, narm);
		}
		if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;

	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::summary(std::string fun, bool narm, SpatOptions &opt) {
	std::vector<double> add;
	return summary_numb(fun, add, narm, opt);
}



SpatRaster SpatRaster::modal(std::vector<double> add, std::string ties, bool narm, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	out.source[0].names[0] = "modal" ;
  	if (!hasValues()) { return out; }


	std::vector<std::string> f {"lowest", "highest", "first", "random", "NA"};
	//std::vector<std::string>::iterator it; 
	auto it = std::find(f.begin(), f.end(), ties);
	if (it == f.end()) {
		out.setError("unknown ties choice");
		return out;
	} 
	size_t ities = std::distance(f.begin(), it);

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	
  	if (!out.writeStart(opt)) { return out; }

	uint32_t seed = 1;
	std::default_random_engine rgen(seed);
	std::uniform_real_distribution<double> dist (0.0,1.0);

	unsigned nl = nlyr();
	std::vector<double> v(nl);
	v.insert( v.end(), add.begin(), add.end() );

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}		
			b[j] = modal_value(v, ities, narm, rgen, dist);
		}
		if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::range(std::vector<double> add, bool narm, SpatOptions &opt) {
	SpatRaster out = geometry(2);
	out.source[0].names.resize(2);
	out.source[0].names[0] = "range_min" ;
	out.source[0].names[1] = "range_max" ;
  	if (!hasValues()) { return out; }

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	
  	if (!out.writeStart(opt)) { return out; }
	unsigned nl = nlyr();
	std::vector<double> v(nl);
	v.insert( v.end(), add.begin(), add.end() );

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc * 2);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}
			std::vector<double> rng = vrange(v, narm);
			b[j] = rng[0];
			b[j+nc] = rng[1];
		}
		if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;

	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRasterStack::summary_numb(std::string fun, std::vector<double> add, bool narm, SpatOptions &opt) {

	std::vector<unsigned> vnl = nlyr();
	unsigned nl = vmax(vnl, false);
	SpatRaster out = ds[0].geometry(nl);
	unsigned ns = nsds();

	std::vector<std::string> f {"sum", "mean", "median", "which.min", "which.max", "min", "max", "range", "prod", "any", "all", "stdev"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown summary function");
		return out;
	}

	if (fun == "range") {
		out.setError("parallel range not implemented, use min and max");
		return out;
	}
  	if (!ds[0].hasValues()) { return out; }

	std::function<double(std::vector<double>&, bool)> sumFun;
	if (fun == "stdev") {
		sumFun = vstdev;
	} else {
		sumFun = getFun(fun);
	}
	for (size_t i=0; i < ns; i++) {
		if (!ds[i].readStart()) {
			out.setError(ds[i].getError());
			return(out);
		}
	}

	
  	if (!out.writeStart(opt)) { return out; }
	std::vector<double> v(ns);
	if (add.size() > 0) v.insert( v.end(), add.begin(), add.end() );

	std::vector<std::vector<double>> a(ns); 
	for (size_t i=0; i < out.bs.n; i++) {
		unsigned nc = out.bs.nrows[i] * out.ncol() * nl;
		for (size_t j=0; j < ns; j++) {
			a[j] = ds[j].readBlock(out.bs, i);
			recycle(a[j], nc);
		}
		std::vector<double> b(nc);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<ns; k++) {
				v[k] = a[k][j];
			}
			b[j] = sumFun(v, narm);
		}
		if (!out.writeValues(b, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;

	}
	for (size_t i=0; i < ns; i++) {
		ds[i].readStop();
	}
	out.writeStop();
	return(out);
}



SpatRaster SpatRasterStack::summary(std::string fun, bool narm, SpatOptions &opt) {
	std::vector<double> add;
	return summary_numb(fun, add, narm, opt);
}

