using namespace std;
#include "spat.h"
#include <type_traits>
#include <vector>


template <typename T>
bool is_NA(const T v) {
    if (std::is_floating_point<T>::value) {
        return isnan(v);
    }
    return (v == std::numeric_limits<T>::min()) ;
}


template <class T> class NA {
public:
    static constexpr T value = std::is_floating_point<T>::value ? NAN : std::numeric_limits<T>::min();
};


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
                v[i] = max(v[i], v[i-1]);
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) | is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] = max(v[i], v[i-1]);
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
                v[i] = min(v[i], v[i-1]);
            }
        }
    } else {
        for (size_t i=1; i<v.size(); i++) {
            if (is_NA(v[i]) | is_NA(v[i-1])) {
                v[i] = NA<T>::value;
            } else {
                v[i] = min(v[i], v[i-1]);
            }
        }
    }
}



template <typename T>
T vsum(std::vector<T>& v, bool narm) {
	T x = v[0];
	T nodata = NA<T>::value;
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (x == nodata) {
				x = v[i];
			} else if (v[i] != nodata) {
				x += v[i];
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (x != nodata) {
				if (v[i] == nodata) {
					x = nodata;
				} else {
					x += v[i];
				}
			}
		}
	}
	return x;
}


template <typename T>
T vprod(std::vector<T>& v, bool narm) {
	T x = v[0];
	T nodata = NA<T>::value;
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (x == nodata) {
				x = v[i];
			} else if (v[i] != nodata) {
				x *= v[i];
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (x != nodata) {
				if (v[i] == nodata) {
					x = nodata;
				} else {
					x *= v[i];
				}
			}
		}
	}
	return x;
}



template <typename T>
T vmean(std::vector<T>& v, bool narm) {
	T x = 0;
	unsigned d = 0;
	T nodata = NA<T>::value;
	if (narm) {
		for (size_t i=0; i<v.size(); i++) {
			if (v[i] != nodata) {
				x += v[i];
				d++;
			}
		}
	} else {
		for (size_t i=0; i<v.size(); i++) {
			if (x != nodata) {
				if (v[i] == nodata) {
					x = nodata;
					d = 0;
				} else {
					x += v[i];
					d++;
				}
			}
		}
	}
	if (d > 0) {
		x /= d;
	} else {
		x = nodata;
	}
	return x;
}




template <typename T>
T vmin(std::vector<T>& v, bool narm) {
	T x = v[0];
	T nodata = NA<T>::value;
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (v[i] != nodata) {
				if (x == nodata) {
					x = v[i];
				} else {
					x = min(x, v[i]);
				}
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (x != nodata) {
				if (v[i] == nodata) {
					x = nodata;
				} else {
					x = min(x, v[i]);
				}
			}
		}
	}
	return x;
}


template <typename T>
T vmax(std::vector<T>& v, bool narm) {
	T x = v[0];
	T nodata = NA<T>::value;

	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (v[i] != nodata) {
				if (x == nodata) {
					x = v[i];
				} else {
					x = max(x, v[i]);
				}
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (x != nodata) {
				if (v[i] == nodata) {
					x = nodata;
				} else {
					x = max(x, v[i]);
				}
			}
		}
	}
	return x;
}


// problematic; should be ok for int and float but
// won't work with bool values (nodata == 0)
template <typename T>
T vany(std::vector<T>& v, bool narm) {
	T nodata = NA<T>::value;
	T x = nodata;
	bool hasnd = false;
	for (size_t i=1; i<v.size(); i++) {
		if (v[i] != nodata) {
			if (v[i] != 0) {
				x = 1;
				break;
			} else {
				x = 0;
			}
		} else {
			hasnd = true;
		}
	}
	if (hasnd & (x == 0) & (!narm)) {
		x = nodata;
	}

	return x;
}

template <typename T>
T vall(std::vector<T>& v, bool narm) {
	T nodata = NA<T>::value;
	T x;
	for (size_t i=1; i<v.size(); i++) {
		if (narm) {
			x = -1;
			if (v[i] != nodata) {
				if (v[i] == 0) {
					x = 0;
					break;
				} else {
					x = 1;
				}
			}
			x = x < 0 ? nodata : x;
		} else {
			x = 1;
			if ((v[i] == nodata) | (v[i] == 0)) {
				x = 0;
				break;
			}
		}
	}
	return x;
}


template <typename T>
std::vector<T> vrange(std::vector<T>& v, bool narm) {
	std::vector<T> x = { v[0], v[0] };
	T nodata = NA<T>::value;

	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (v[i] != nodata) {
				if (x[0] == nodata) {
					x[0] = v[i];
					x[1] = v[i];
				} else {
					x[0] = min(x[0], v[i]);
					x[1] = max(x[1], v[i]);
				}
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (x[0] != nodata) {
				if (v[i] == nodata) {
					x[0] = nodata;
					x[1] = nodata;
				} else {
					x[0] = min(x[0], v[i]);
					x[1] = max(x[1], v[i]);
				}
			}
		}
	}
	return x;
}




SpatRaster SpatRaster::cum(std::string fun, bool narm, std::string filename, bool overwrite) {

	SpatRaster out = geometry();

	std::vector<string> f {"sum", "prod", "min", "max"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.error = true;
		out.error_message = "unknown cum function";
		return out;
	}

  	out.writeStart(filename, overwrite);
	readStart();
	unsigned nl = out.nlyr();
	std::vector<double> v(nl);
	unsigned nc;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		nc = out.bs.nrows[i] * out.ncol;
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
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::summary(std::string fun, bool narm, std::string filename, bool overwrite) {

	SpatRaster out = geometry(1);

	std::vector<string> f {"sum", "mean", "min", "max", "range", "prod", "any", "all"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.error = true;
		out.error_message = "unknown summary function";
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

  	out.writeStart(filename, overwrite);
	readStart();
	unsigned nl = nlyr();
	std::vector<double> v(nl);
	unsigned nc;
	unsigned nlout = out.nlyr();
	
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readBlock(out.bs, i);
		nc = out.bs.nrows[i] * out.ncol;
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
                std:vector<double> rng = vrange(v, narm);
				b[j] = rng[0];
				b[j+nc] = rng[1];
			}
		}
		out.writeValues(b, out.bs.row[i]);
	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::summary_numb(std::string fun, std::vector<double> add, bool narm, std::string filename, bool overwrite) {

	SpatRaster out = geometry();
	out.source[0].nlyr=1;
	std::vector<string> f {"sum", "mean", "min", "max", "range", "prod", "any", "all"};
	if (std::find(f.begin(), f.end(), fun) == f.end()) {
		out.error = true;
		out.error_message = "unknown summary function";
		return out;
	}
	if (fun == "range") {
		out.source[0].nlyr = 2;
	}
	//std:vector<double> rng(2);

	return out;
}

