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

//#include <type_traits>
#include "spatRaster.h"
#include "vecmath.h"


template <typename T>
void cumsum(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] += v[i-1];
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) | is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] += v[i-1];
            }
        }
    }
}

template <typename T>
void cumprod(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] *= v[i-1];
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) | is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] *= v[i-1];
            }
        }
    }
}


template <typename T>
void cummax(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] = std::max(v[i], v[i-1]);
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) | is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] = std::max(v[i], v[i-1]);
            }
        }
    }
}


template <typename T>
void cummin(std::vector<T>& v, bool narm) {
    if (narm) {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i])) {
                v[i] = v[i-1];
            } else if (!is_NA(v[i-1])){
                v[i] = std::min(v[i], v[i-1]);
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) | is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] = std::min(v[i], v[i-1]);
            }
        }
    }
}


SpatRaster SpatRaster::cum(std::string fun, bool narm, SpatOptions &opt) {

	SpatRaster out = geometry();

	std::vector<std::string> f {"sum", "prod", "min", "max"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown cum function");
		return out;
	}
	if (!hasValues()) {
		out.setError("raster has no values");
		return out;
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
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
		if (!out.writeValues(a, out.bs.row[i])) return out;

	}
	out.writeStop();
	readStop();
	return(out);
}



SpatRaster SpatRaster::summary_numb(std::string fun, std::vector<double> add, bool narm, SpatOptions &opt) {

	SpatRaster out = geometry(1);

	std::vector<std::string> f {"sum", "mean", "min", "max", "range", "prod", "any", "all"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.setError("unknown summary function");
		return out;
	}

	if (fun == "range") {
		out.source[0].nlyr = 2;
		out.source[0].names.resize(2);
		out.source[0].names[0] = "range_min" ;
		out.source[0].names[1] = "range_max" ;
	} else {
		out.source[0].names[0] = fun;
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	unsigned nl = nlyr();
	std::vector<double> v(nl);
	v.insert( v.end(), add.begin(), add.end() );

	//unsigned nc;
	unsigned nlout = out.nlyr();

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		unsigned nc = out.bs.nrows[i] * out.ncol();
		std::vector<double> b(nc * nlout);
		for (size_t j=0; j<nc; j++) {
			for (size_t k=0; k<nl; k++) {
				v[k] = a[j+k*nc];
			}
			if (fun == "sum") {
				b[j] = vsum(v, narm);
			} else if (fun == "mean") {
				b[j] = vmean(v, narm);
			} else if (fun == "prod") {
				b[j] = vprod(v, narm);
			} else if (fun == "min") {
				b[j] = vmin(v, narm);
			} else if (fun == "max") {
				b[j] = vmax(v, narm);
			} else if (fun == "any") {
				b[j] = vany(v, narm);
			} else if (fun == "all") {
				b[j] = vall(v, narm);
			} else if (fun == "range") {
                std::vector<double> rng = vrange(v, narm);
				b[j] = rng[0];
				b[j+nc] = rng[1];
			}
		}
		if (!out.writeValues(b, out.bs.row[i])) return out;

	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::summary(std::string fun, bool narm, SpatOptions &opt) {
	std::vector<double> add;
	return summary_numb(fun, add, narm, opt);
}

