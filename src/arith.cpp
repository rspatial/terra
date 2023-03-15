// Copyright (c) 2018-2023  Robert J. Hijmans
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
#include "vecmath.h"

//#include "modal.h"

/*
// need to take care of NAs here. OK for NAN, but not for int types
//template <typename T>
//void operator+(std::vector<T>& a, const std::vector<T>& b) {
void operator+(std::vector<double>& a, const std::vector<double>& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<double>());
}


//template <typename T>
void operator-(std::vector<double>& a, const std::vector<double>& b) {
   std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<double>());
}


//template <typename T>
void operator/(std::vector<double>& a, const std::vector<double>& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::divides<double>());
}

//template <typename T>
void operator*(std::vector<double>& a, const std::vector<double>& b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<double>());
}

*/

//template <typename T>
void operator%(std::vector<double>& a, const std::vector<double>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::modulus<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = std::fmod(a[i], b[i]);
		}
	}
}


//template <typename T>
void operator>=(std::vector<double>& a, const std::vector<double>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::greater_equal<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] >= b[i];
		}
	}
}

//template <typename T>
void operator<=(std::vector<double>& a, const std::vector<double>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::less_equal<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] <= b[i];
		}
	}
 }


//template <typename T>
void operator>(std::vector<double>& a, const std::vector<double>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::greater<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] > b[i];
		}
	}
}

//template <typename T>
void operator<(std::vector<double>& a, const std::vector<double>& b) {
//    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::less<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			a[i] = NAN;
		} else {
			a[i] = a[i] < b[i];
		}
	}
}



bool smooth_operator(std::string &oper, bool &logical, bool &reverse, bool &falseNA) {
	std::vector<std::string> f {"==", "!=", ">", "<", ">=", "<="};
	logical = std::find(f.begin(), f.end(), oper) != f.end();
	f = {"+", "-", "*", "^", "/", "%", "%%", "%/%"};
	bool ok = logical || (std::find(f.begin(), f.end(), oper) != f.end());
	if (ok) {
		if (oper == "%%") oper = "%";
		if (logical) {
			if (reverse) {
				if (oper == ">") {				
					oper = "<=";
				} else if (oper == "<") {				
					oper = ">=";
				} else if (oper == ">=") {				
					oper = "<";
				} else if (oper == "<=") {				
					oper = ">";
				}
				reverse = false;
			}
		} else {
			falseNA = false;
		}
	}
	return ok;
}




SpatRaster SpatRaster::arith(SpatRaster x, std::string oper, bool falseNA, SpatOptions &opt) {

	size_t nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);
	if (!(hasValues() && x.hasValues())) {
		out.setError("raster has no values"); // or warn and treat as NA?
		return out;
	}

	bool logical = false;
	bool reverse = false;
	if (!smooth_operator(oper, logical, reverse, falseNA)) {
		out.setError("unknown arith function");
		return out;
	}
	if (logical) {
		out.setValueType(3);
	} else if (oper != "/") {
		std::vector<int> v = getValueType(false);
		std::vector<int> vx = x.getValueType(false);
		bool is_int = true;
		for (size_t i = 0; i<v.size(); i++) {
			if ((v[i] != 1) && (v[i] != 3)) { 
				is_int = false;
				break;
			}
		}
		for (size_t i = 0; i<vx.size(); i++) {
			if (!is_int) break;
			if ((vx[i] != 1) && (vx[i] != 3)) { 
				is_int = false;
			}
		}
		if (is_int) {
			out.setValueType(1);
		}
	}

	if (!out.compare_geom(x, false, true, opt.get_tolerance(), true)) {
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
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		x.readStop();
		return out;
	}

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a, b;
		readBlock(a, out.bs, i);
		x.readBlock(b, out.bs, i);
		recycle(a,b);
		if (oper == "+") {
			std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<double>());
		} else if (oper == "-") {
			std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<double>());
		} else if (oper == "*") {
			std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<double>());
		} else if (oper == "/") {
			std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::divides<double>());
		} else if (oper == "^") {
			for (size_t i=0; i<a.size(); i++) {
				if (std::isnan(a[i]) || std::isnan(b[i])) {
					a[i] = NAN;
				} else {
					a[i] = std::pow(a[i], b[i]);
				}
			}
		} else if (oper == "%") {
			 a % b;
		} else if (oper == "%/%") {
			for (size_t i=0; i<a.size(); i++) {
				a[i] = trunc(a[i] / b[i]);
			}
		} else if (oper == "==") {
			for (size_t i=0; i<a.size(); i++) {
				if (std::isnan(a[i]) || std::isnan(b[i])) {
					a[i] = NAN;
				} else {
					a[i] = a[i] == b[i];
				}
			}
		} else if (oper == "!=") {
			for (size_t i=0; i<a.size(); i++) {
				if (std::isnan(a[i]) || std::isnan(b[i])) {
					a[i] = NAN;
				} else {
					a[i] = a[i] != b[i];
				}
			}
		} else if (oper == ">=") {
			a >= b;
		} else if (oper == "<=") {
			a <= b;
		} else if (oper == ">") {
			a > b;
		} else if (oper == "<") {
			a < b;
		}
		if (falseNA) {
			for (double& d : a) if (!d) d = NAN;
		}
		if (!out.writeBlock(a, i)) return out;

	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}


SpatRaster SpatRaster::arith(double x, std::string oper, bool reverse, bool falseNA, SpatOptions &opt) {

	SpatRaster out = geometry(nlyr());
	if (!hasValues()) {
		out.setError("raster has no values"); // or warn and treat as NA?
		return out;
	}

	bool logical;
	if (!smooth_operator(oper, logical, reverse, falseNA)) {
		out.setError("unknown arith function");
		return out;
	}
	if (logical) {
		out.setValueType(3);
	} else if (oper != "/") {
		std::vector<int> v = getValueType(false);
		bool is_int = true;
		for (size_t i = 0; i<v.size(); i++) {
			if ((v[i] != 1) && (v[i] != 3)) { 
				is_int = false;
				break;
			}
		}
		if (is_int && (x == std::round(x))) { 
			out.setValueType(1);
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		if (std::isnan(x)) {
			for(double& d : a)  d = NAN;
		} else if (oper == "+") {
			for (double& d : a)  d += x;
		} else if (oper == "-") {
			if (reverse) {
				for (double& d : a)  d = x - d;
			} else {
				for (double& d : a)  d -= x;
			}
		}
		else if (oper == "*") {
			for(double& d : a)  d *= x;
		} else if (oper == "/") {
			if (reverse) {
				for (double& d : a)  d = x / d;
			} else {
				for (double& d : a)  d /= x;
			}
		} else if (oper == "^") {
			if (reverse) {
				for (double& d : a)  d = std::pow(x, d);
			} else {
				for (double& d : a)  d = std::pow(d, x);
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
		} else if (oper == "%/%") {
			if (reverse) {
				for (double& d : a) d = trunc(x / d);
			} else {
				for (double& d : a) d = trunc(d / x);
			}
		} else if (oper == "==") {
			for (double& d : a) if (!std::isnan(d)) d = d == x;
		} else if (oper == "!=") {
			for(double& d : a) if (!std::isnan(d))  d = d != x;
		} else if (oper == ">=") {
			for (double& d : a) if (!std::isnan(d)) d = d >= x;
		} else if (oper == "<=") {
			for (double& d : a) if (!std::isnan(d)) d = d <= x;
		} else if (oper == ">") {
			for (double& d : a) if (!std::isnan(d)) d = d > x;
		} else if (oper == "<") {
			for (double& d : a) if (!std::isnan(d)) d = d < x;
		} else {
			out.setError("unknown arith function");
			return out;
		}
		if (falseNA) {
			for (double& d : a) if (!d) d = NAN;
		}

		if (!out.writeBlock(a, i)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}

SpatRaster SpatRaster::is_true(bool falseNA, SpatOptions &opt) {
	return arith(1, "==", false, false, opt);
}

SpatRaster SpatRaster::is_false(bool falseNA, SpatOptions &opt) {
	return arith(1, "!=", false, falseNA, opt);
}


SpatRaster SpatRaster::arith(std::vector<double> x, std::string oper, bool reverse, bool falseNA, SpatOptions &opt) {

	if (x.empty()) {
		SpatRaster out;
		out.setError("cannot compute with nothing");
		return out;
	}
	if (!hasValues()) {
		SpatRaster out;
		out.setError("raster has no values"); // or warn and treat as NA?
		return out;
	}

	if (x.size() == 1) {
		return(arith(x[0], oper, reverse, falseNA, opt));
	}

	unsigned innl = nlyr();
	unsigned outnl = innl;

	if (x.size() > innl) {
		outnl = x.size();
	}
	SpatRaster out = geometry(outnl);

	bool logical=false;
	if (!smooth_operator(oper, logical, reverse, falseNA)) {
		out.setError("unknown arith function");
		return out;
	}
	
	if (logical) {
		out.setValueType(3);
	} else if (oper != "/") {
		std::vector<int> v = getValueType(false);
		bool is_int = true;
		for (size_t i = 0; i<v.size(); i++) {
			if ((v[i] != 1) && (v[i] != 3)) { 
				is_int = false;
				break;
			}
		}
		for (size_t i = 0; i<x.size(); i++) {
			if (!is_int) break;
			is_int = x[i] == std::round(x[i]); 
		}
		if (is_int) { 
			out.setValueType(1);
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	unsigned nc = ncol();
	recycle(x, outnl);

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		if (outnl > innl) {
			recycle(v, outnl * out.bs.nrows[i] * nc);
		}
		//std::vector<double> vv;
		unsigned off = out.bs.nrows[i] * nc;
		for (size_t j=0; j<outnl; j++) {
			unsigned s = j * off;
			if (std::isnan(x[j])) {
				for (size_t k=0; k<off; k++) {
					v[s+k] = NAN;
				}
			} else if (oper == "+") {
				for (size_t k=0; k<off; k++) {
					v[s+k] += x[j];
				}
			} else if (oper == "-") {				
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = x[j] - v[s+k];
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] -= x[j];
					}
				}
			} else if (oper == "*") {
				for (size_t k=0; k<off; k++) {
					v[s+k] *= x[j];
				}
			} else if (oper == "/") {				
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = x[j] / v[s+k];
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] /= x[j];
					}
				}
			} else if (oper == "^") {				
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::pow(x[j], v[s+k]);
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::pow(v[s+k], x[j]);
					}
				}
			} else if (oper == "%") {
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::fmod(x[j], v[s+k]);
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::fmod(v[s+k], x[j]);
					}
				}
			} else if (oper == "==") {
				for (size_t k=0; k<off; k++) {
					v[s+k] = v[s+k] == x[j];
				}			
			} else if (oper == "!=") {
				for (size_t k=0; k<off; k++) {
					v[s+k] = v[s+k] != x[j];
				}			
			} else if (oper == ">=") {				
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] >= x[j];
				}			
			} else if (oper == "<=") {
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] <= x[j];
				}
			} else if (oper == ">") {
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] > x[j];
				}
			} else if (oper == "<") {
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] < x[j];
				}
			} else {
				out.setError("unknown arith function");
				return out;
			}
		}
		if (falseNA) {
			for (double& d : v) if (!d) d = NAN;
		}

		if (!out.writeBlock(v, i)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::arith_m(std::vector<double> x, std::string oper, std::vector<size_t> dim, bool reverse, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (!hasValues()) {
		out.setError("raster has no values"); // or warn and treat as NA?
		return out;
	}
	
	size_t nx = x.size();
	if (nx == 0) {
		out.setError("cannot compute with nothing");
		return out;
	}
	if (dim.size() != 2) {
		out.setError("incorrect dimensions"); 
		return out;
	}
	if ((dim[1] * dim[0]) != nx) {
		out.setError("incorrect matrix dimensions (dim(m) != length(x))"); 
		return out;
	}
	if (nx == 1) {
		return(arith(x[0], oper, reverse, false, opt));
	}
	// single cell
	if (dim[0] < 2) {
		return(arith(x, oper, reverse, false, opt));
	}
	if (dim[0] > ncell()) {
		out.setError("incorrect matrix dimensions (nrow > ncell(x))"); 
		return out;
	}

	size_t nl = nlyr();
	if (dim[1] > nl) {
		out.setError("incorrect matrix dimensions (ncol != nlyr(x))"); 
		return out;
	} else if (dim[1] < nl) {
		recycle(x, nl * dim[0]);
		dim[1] = nl;
	}
	
	bool logical;
	bool falseNA=false;
	if (!smooth_operator(oper, logical, reverse, falseNA)) {
		out.setError("unknown arith function");
		return out;
	}
	if (logical) {
		out.setValueType(3);
	} 
	
	if (logical) {
		out.setValueType(3);
	} else if (oper != "/") {
		std::vector<int> v = getValueType(false);
		bool is_int = true;
		for (size_t i = 0; i<v.size(); i++) {
			if ((v[i] != 1) && (v[i] != 3)) { 
				is_int = false;
				break;
			}
		}
		for (size_t i = 0; i<x.size(); i++) {
			if (!is_int) break;
			is_int = x[i] == std::round(x[i]); 
		}
		if (is_int) { 
			out.setValueType(1);
		}
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	unsigned nc = ncol();
	
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		readBlock(v, out.bs, i);
		unsigned off = out.bs.nrows[i] * nc;
		for (size_t j=0; j<nl; j++) {
			size_t s = j * off;
			size_t d = j * dim[0];
			std::vector<double> xj = {x.begin()+d, x.begin()+d+dim[0]};
			size_t start = (out.bs.row[i] * nc) % xj.size();
			if (start != 0) {
				std::rotate(xj.begin(), xj.begin()+start, xj.end());
			}
			recycle(xj, off);
			if (oper == "+") {
				for (size_t k=0; k<off; k++) {
					v[s+k] += xj[k];
				}
			} else if (oper == "-") {				
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = xj[k] - v[s+k];
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] -= xj[k];
					}
				}
			} else if (oper == "*") {
				for (size_t k=0; k<off; k++) {
					v[s+k] *= xj[k];
				}
			} else if (oper == "/") {				
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = xj[k] / v[s+k];
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] /= xj[k];
					}
				}
				
			} else if (oper == "^") {				
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::pow(xj[k], v[s+k]);
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::pow(v[s+k], xj[k]);
					}
				}
			} else if (oper == "%") {
				if (reverse) {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::fmod(xj[k], v[s+k]);
					}
				} else {
					for (size_t k=0; k<off; k++) {
						v[s+k] = std::fmod(v[s+k], xj[k]);
					}
				}
			} else if (oper == "==") {
				for (size_t k=0; k<off; k++) {
					v[s+k] = v[s+k] == xj[k];
				}			
			} else if (oper == "!=") {
				for (size_t k=0; k<off; k++) {
					v[s+k] = v[s+k] != xj[k];
				}			
			} else if (oper == ">=") {				
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] >= xj[k];
				}			
			} else if (oper == "<=") {
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] <= xj[k];
				}
			} else if (oper == ">") {
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] > xj[k];
				}
			} else if (oper == "<") {
				for (size_t k=0; k<off; k++) {
					if (!std::isnan(v[s+k])) v[s+k] = v[s+k] < xj[k];
				}
			} else {
				// stop
			}
		}
		if (!out.writeBlock(v, i)) return out;
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

	std::vector<std::string> f = {"ceiling", "floor", "trunc", "sign", "log", "log10", "log2", "log1p", "exp", "expm1", "abs", "sqrt"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown math function");
		return out;
	}
	f = {"ceiling", "floor", "trunc", "sign"};
	if (std::find(f.begin(), f.end(), fun) != f.end()) {
		out.setValueType(1);
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


  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		for(double& d : a) if (!std::isnan(d)) d = mathFun(d);
		if (!out.writeBlock(a, i)) return out;
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

	if (digits == 0) out.setValueType(1);

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}


  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		if (fun == "round") {
			for(double& d : a) d = roundn(d, digits);
		} else if (fun == "signif") {
			for(double& d : a) if (!std::isnan(d)) d = signif(d, digits);
		}
		if (!out.writeBlock(a, i)) return out;
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


  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		for (double& d : a) if (!std::isnan(d)) d = trigFun(d);
		if (!out.writeBlock(a, i)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::atan_2(SpatRaster x, SpatOptions &opt) {

	size_t nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);
	if ((!hasValues()) || (!x.hasValues())) {
		return out;
	} 
	
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!x.readStart()) {
		out.setError(x.getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		x.readStop();
		return out;
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a, b;
		readValues(a, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		x.readValues(b, out.bs.row[i], out.bs.nrows[i], 0, ncol());
		recycle(a, b);
		std::vector<double> d(a.size());
		for (size_t j=0; j<a.size(); j++) {
			if (std::isnan(a[j]) || std::isnan(b[j])) {
				d[j] = NAN;
			} else {
				d[j] = atan2(a[j], b[j]);
			}
		}
		if (!out.writeBlock(d, i)) return out;
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


SpatRaster SpatRaster::isnot(bool falseNA, SpatOptions &opt) {
	SpatRaster out = geometry();
	out.setValueType(3);
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	if (falseNA) {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> a;
			readBlock(a, out.bs, i);
			for (size_t j=0; j<a.size(); j++) {
				a[i] = a[i] ? NAN : 1;
			}
			if (!out.writeBlock(a, i)) return out;

		}
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> a;
			readBlock(a, out.bs, i);
			for (size_t j=0; j<a.size(); j++) {
				a[i] = !a[i];
			}
			if (!out.writeBlock(a, i)) return out;

		}
	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::logic(SpatRaster x, std::string oper, SpatOptions &opt) {

	size_t nl = std::max(nlyr(), x.nlyr());
	SpatRaster out = geometry(nl);

	out.setValueType(3);

	std::vector<std::string> f {"&", "|"};
	if (std::find(f.begin(), f.end(), oper) == f.end()) {
		out.setError("unknown logic function");
		return out;
	}

	if (!out.compare_geom(x, false, true, opt.get_tolerance(), true)) {
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

 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		std::vector<double> b;
		x.readBlock(b, out.bs, i);
		recycle(a, b);
		if (oper == "&") {
			a = a & b;
		} else if (oper == "|") {
			a = a | b;
		} else {
			// stop
		}
		if (!out.writeBlock(a, i)) return out;

	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}



SpatRaster SpatRaster::logic(bool x, std::string oper, SpatOptions &opt) {

	SpatRaster out = geometry();
	out.setValueType(3);

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		if (std::isnan(x)) {
			for(double& d : a)  d = NAN;
		} else if (oper == "&") {
			for(double& d : a)  d = (d==1) & x;
		} else if (oper == "|") {
			for(double& d : a)  d = (d==1) | x;
		} else if (oper == "istrue") {
			for(double& d : a)  d = d==1 ? 1 : 0;
		} else if (oper == "isfalse") {
			for(double& d : a)  d = d!=1 ? 1 : 0;
		} else {
			out.setError("unknown operator: " + oper);
			return out;
		}
		if (!out.writeBlock(a, i)) return out;
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


  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	unsigned nl = out.nlyr();
	std::vector<double> v(nl);
	unsigned nc;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
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
		if (!out.writeBlock(a, i)) return out;

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
  	if (!hasValues()) { return out; }

	if (fun == "range") {
		return range(add, narm, opt);
	}
	out.source[0].names[0] = fun;
	std::function<double(std::vector<double>&, bool)> sumFun;
	if (fun == "std") {
		sumFun = vstdev;
	} else {
		if (!haveFun(fun)) {
			out.setError("unknown function argument");
			return out;
		}
		sumFun = getFun(fun);
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}


	opt.ncopies = std::max(opt.ncopies, nlyr() * 2);
	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	unsigned nl = nlyr();
	std::vector<double> v(nl);
	if (!add.empty()) v.insert( v.end(), add.begin(), add.end() );

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}
			b[j] = sumFun(v, narm);
		}
		if (!out.writeBlock(b, i)) return out;

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

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	uint32_t seed = 1;
	std::default_random_engine rgen(seed);
	std::uniform_real_distribution<double> dist (0.0,1.0);

	unsigned nl = nlyr();
	std::vector<double> v(nl);
	v.insert( v.end(), add.begin(), add.end() );

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}
			b[j] = modal_value(v, ities, narm, rgen, dist);
		}
		if (!out.writeBlock(b, i)) return out;
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


  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	unsigned nl = nlyr();
	std::vector<double> v(nl);
	v.insert( v.end(), add.begin(), add.end() );

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a;
		readBlock(a, out.bs, i);
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
		if (!out.writeBlock(b, i)) return out;

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

	if (fun == "range") {
		out.setError("parallel range not implemented, use min and max");
		return out;
	}
  	if (!ds[0].hasValues()) { return out; }

	std::function<double(std::vector<double>&, bool)> sumFun;
	if (fun == "std") {
		sumFun = vstdev;
	} else {
		if (!haveFun(fun)) {
			out.setError("unknown function argument");
			return out;
		}
		sumFun = getFun(fun);
	}

	for (size_t i=0; i < ns; i++) {
		if (!ds[i].readStart()) {
			out.setError(ds[i].getError());
			return(out);
		}
	}


  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> v(ns);
	if (!add.empty()) v.insert( v.end(), add.begin(), add.end() );

	std::vector<std::vector<double>> a(ns);
	for (size_t i=0; i < out.bs.n; i++) {
		unsigned nc = out.bs.nrows[i] * out.ncol() * nl;
		for (size_t j=0; j < ns; j++) {
			ds[j].readBlock(a[j], out.bs, i);
			recycle(a[j], nc);
		}
		std::vector<double> b(nc);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<ns; k++) {
				v[k] = a[k][j];
			}
			b[j] = sumFun(v, narm);
		}
		if (!out.writeBlock(b, i)) return out;

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



SpatRaster SpatRaster::isnan(bool falseNA, SpatOptions &opt) {
	SpatRaster out = geometry();
	out.setValueType(3);

    if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	if (falseNA) {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = std::isnan(d) ? 1 : NAN;
			if (!out.writeBlock(v, i)) return out;
		}
	} else {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = std::isnan(d);
			if (!out.writeBlock(v, i)) return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::anynan(bool falseNA, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	out.setValueType(3);

    if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	double inval = falseNA ? NAN : 0;
	size_t nl = nlyr();
	size_t nc = ncol();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v, w;
		readBlock(v, out.bs, i);
		size_t off = out.bs.nrows[i] * nc;
		w.resize(off, inval);
		for (size_t j=0; j<off; j++) {
			for (size_t k=0; k<nl; k++) {
				size_t cell = j + k * off;
				if (std::isnan(v[cell])) {
					w[j] = 1;
					continue;
				}
			}
		}
		if (!out.writeBlock(w, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::nonan(bool falseNA, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	out.setValueType(3);

    if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	double inval = falseNA ? NAN : 0;
	size_t nl = nlyr();
	size_t nc = ncol();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v, w;
		readBlock(v, out.bs, i);
		size_t off = out.bs.nrows[i] * nc;
		w.resize(off, 1);
		for (size_t j=0; j<off; j++) {
			for (size_t k=0; k<nl; k++) {
				size_t cell = j + k * off;
				if (std::isnan(v[cell])) {
					w[j] = inval;
					continue;
				}
			}
		}
		if (!out.writeBlock(w, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::allnan(bool falseNA, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	out.setValueType(3);

    if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	double outval = falseNA ? NAN : 0;

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	size_t nl = nlyr();
	size_t nc = ncol();
	for (size_t i=0; i<out.bs.n; i++) {
		std::vector<double> v, w;
		readBlock(v, out.bs, i);
		size_t off = out.bs.nrows[i] * nc;
		w.resize(off, 1);
		for (size_t j=0; j<off; j++) {
			for (size_t k=0; k<nl; k++) {
				size_t cell = j + k * off;
				if (!std::isnan(v[cell])) {
					w[j] = outval;
					continue;
				}
			}
		}
		if (!out.writeBlock(w, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::isnotnan(bool falseNA, SpatOptions &opt) {
	SpatRaster out = geometry();
	out.setValueType(3);
    if (!hasValues()) return out;

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	if (falseNA) {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = std::isnan(d) ? NAN : 1;
			if (!out.writeBlock(v, i)) return out;
		}
	} else {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = !std::isnan(d);
			if (!out.writeBlock(v, i)) return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::countnan(long n, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (n > 0) {
		out.setValueType(3);
	}

    if (!hasValues()) return out;
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	size_t nl = nlyr();
	size_t nc = ncol();
	if (n > ((long) nlyr())) {
		out.addWarning("n > nlyr(x)");
		std::vector<double> w;
		for (size_t i=0; i<out.bs.n; i++) {
			size_t off = out.bs.nrows[i] * nc;
			w.resize(off, 0);
			if (!out.writeBlock(w, i)) return out;
		}
	} else if (n > 0) {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v, w;
			readBlock(v, out.bs, i);
			size_t off = out.bs.nrows[i] * nc;
			w.resize(off, 0);
			for (size_t j=0; j<off; j++) {
				long cnt = 0;
				for (size_t k=0; k<nl; k++) {
					size_t cell = j + k * off;
					if (std::isnan(v[cell])) {
						cnt++;
						if (cnt == n) {
							w[j] = 1;
						}
					}
				}
			}
			if (!out.writeBlock(w, i)) return out;
		}
	} else {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v, w;
			readBlock(v, out.bs, i);
			size_t off = out.bs.nrows[i] * nc;
			w.resize(off, 0);
			for (size_t j=0; j<off; j++) {
				for (size_t k=0; k<nl; k++) {
					size_t cell = j + k * off;
					if (std::isnan(v[cell])) {
						w[j]++;
					}
				}
			}
			if (!out.writeBlock(w, i)) return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::isfinite(bool falseNA, SpatOptions &opt) {
	SpatRaster out = geometry();
	out.setValueType(3);
    if (!hasValues()) return out;

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	if (falseNA) {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = std::isfinite(d) ? 1 : NAN;
			if (!out.writeBlock(v, i)) return out;
		}
	} else {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = std::isfinite(d);
			if (!out.writeBlock(v, i)) return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}


SpatRaster SpatRaster::isinfinite(bool falseNA, SpatOptions &opt) {
	SpatRaster out = geometry();
	out.setValueType(3);

    if (!hasValues()) return out;

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	if (falseNA) {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = std::isinf(d) ? 1 : NAN;
			if (!out.writeBlock(v, i)) return out;
		}
	} else {
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> v;
			readBlock(v, out.bs, i);
			for (double &d : v) d = std::isinf(d);
			if (!out.writeBlock(v, i)) return out;
		}
	}
	readStop();
	out.writeStop();
	return(out);
}


std::vector<std::vector<double>> SpatRaster::where(std::string what, bool values, SpatOptions &opt) {

	unsigned nl = nlyr();
	std::vector<std::vector<double>> out(nl);

	std::vector<std::string> f {"min", "max"};
	if (std::find(f.begin(), f.end(), what) == f.end()) {
		setError("unknown where function");
		return out;
	}

  	if (!hasValues()) {
		setError("SpatRaster has no values");
		return out;
	}

	if (!readStart()) {
		return(out);
	}
	BlockSize bs = getBlockSize(opt);

	std::vector<double> val;
	bool do_min = what == "min";

	std::vector<bool> hr = hasRange();
	bool hasR = true;
	for (size_t i=0; i<hr.size(); i++) {
		if (!hr[i]) {
			hasR = false;
			break;
		}
	}
	if (hasR) {
		if (do_min) {
			val = range_min();
			for (double &d : val) d *= 1.00001;
		} else {
			val = range_max();
			for (double &d : val) d *= 0.99999;
		}
	} else {
		if (do_min) {
			val.resize(nl, std::numeric_limits<double>::max());
		} else {
			val.resize(nl, std::numeric_limits<double>::lowest());
		}
	}
	for (size_t i=0; i<bs.n; i++) {
		std::vector<double> v;
		readBlock(v, bs, i);
		size_t lyrsize = bs.nrows[i] * ncol();
		size_t boff = i * lyrsize;
		for (size_t j=0; j<nl; j++) {
			size_t off = j * lyrsize;
			if (do_min) {
				for (size_t k=off; k<(off+lyrsize); k++) {
					if (!std::isnan(v[k])) {
						if (v[k] < val[j]) {
							val[j] = v[k];
							out[j].resize(0);
							double cell = k - off + boff;
							out[j].push_back(cell);
						} else if (v[k] == val[j]) {
							double cell = k - off + boff;
							out[j].push_back(cell);
						}
					}
				}
			} else {
				for (size_t k=off; k<(off+lyrsize); k++) {
					if (!std::isnan(v[k])) {
						if (v[k] > val[j]) {
							val[j] = v[k];
							out[j].resize(0);
							double cell = k - off + boff;
							out[j].push_back(cell);
						} else if (v[k] == val[j]) {
							double cell = k - off + boff;
							out[j].push_back(cell);
						}
					}
				}
			}
			if (values) {
				std::vector<double> wval(out[j].size(), val[j]);
				out[j].insert(out[j].end(), wval.begin(), wval.end());
			}
		}
	}
	readStop();
	return(out);
}

