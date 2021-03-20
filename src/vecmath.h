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

#ifndef VECMATH_GUARD
#define VECMATH_GUARD

#include <type_traits>
#include <vector>
#include "NA.h"
#include <math.h>



template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}


static inline double interpolate(double x, double y1, double y2, unsigned x1, unsigned x2) {
	double denom = (x2-x1);
	return y1 + (x-x1) * (y2-y1)/denom;
}



static inline std::vector<double> vquantile(std::vector<double> v, const std::vector<double>& probs, bool narm) {
	size_t n = v.size();
    if (n==0) {
        return std::vector<double>(probs.size(), NAN);
    }
    if (n == 1) {
        return std::vector<double>(probs.size(), v[0]);
    }

	//na_omit(v);
	v.erase(std::remove_if(std::begin(v), std::end(v),
        [](const double& value) { return std::isnan(value); }),
        std::end(v));

	if ((!narm) & (v.size() < n)) {
        return std::vector<double>(probs.size(), NAN);
	}
	n = v.size();
    std::sort(v.begin(), v.end());

	size_t pn = probs.size();
	std::vector<double> q(pn);

    for (size_t i = 0; i < pn; ++i) {
		double x = probs[i] * (n-1);
		unsigned x1 = std::floor(x);
		unsigned x2 = std::ceil(x);
		if (x1 == x2) {
			q[i] = v[x1];
		} else {
			q[i] = interpolate(x, v[x1], v[x2], x1, x2);
		}
    }
    return q;
}



template <typename T>
std::vector<T> vunique(std::vector<T> d) {
	std::sort(d.begin(), d.end());
	d.erase(std::unique(d.begin(), d.end()), d.end());
	return d;
}

template <typename T>
std::vector<std::string> vtostring(std::vector<T>& v) {
	std::vector<std::string> s;
	std::transform(std::begin(v),
           std::end(v), std::back_inserter(s),
           [](double d) { return std::to_string(d); } 
        );
	return s;
}



template <typename T>
T vmedian(std::vector<T>& v, bool narm) {
	size_t n = v.size();
	std::vector<T> vv;
	vv.reserve(n);
	for (size_t i=0; i<n; i++) {
        if (!is_NA(v[i])) {
            vv.push_back(v[i]);
        } else if (!narm) {
            return NA<T>::value;
        }
	}
	n = vv.size();
	if (n == 0) {
		return(NA<T>::value);
	}
	size_t n2 = n / 2;
	std::nth_element(vv.begin(), vv.begin()+n2, vv.end());
	T med = vv[n2];
	if (n % 2 == 1) {
		return med;
	} else {
		std::nth_element(vv.begin(), vv.begin()+n2-1, vv.end());
		return 0.5 * (med + vv[n2-1] );
	}
}



template <typename T>
T vsum(std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {		
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(x)) {
				x = v[i];
			} else if (!is_NA(v[i])) {
				x += v[i];
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(x)) {
				if (is_NA(v[i])) {
					x = NA<T>::value;
					break;
				} else {
					x += v[i];
				}
			}
		}
	}
	return x;
}

template <typename T>
T vsum2(std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {		
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(x)) {
				x = v[i] * v[i];
			} else if (!is_NA(v[i])) {
				x += v[i] * v[i];
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(x)) {
				if (is_NA(v[i])) {
					x = NA<T>::value;
					break;
				} else {
					x += v[i] * v[i];
				}
			}
		}
	}
	return x;
}


template <typename T>
T vprod(std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(x)) {
				x = v[i];
			} else if (!is_NA(v[i])) {
				x *= v[i];
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(x)) {
				if (is_NA(v[i])) {
					x = NA<T>::value;
					break;
				} else {
					x *= v[i];
				}
			}
		}
	}
	return x;
}



template <typename T>
double vmean(std::vector<T>& v, bool narm) {
	double x = 0;
	unsigned d = 0;
	if (narm) {
		for (size_t i=0; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				x += v[i];
				d++;
			}
		}
	} else {
		for (size_t i=0; i<v.size(); i++) {
			if (!std::isnan(x)) {
				if (is_NA(v[i])) {
					x = NAN;
					d = 0;
					break;
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
		x = NAN;
	}
	return x;
}

template <typename T>
double vsd(std::vector<T>& v, bool narm) {
	double m = vmean(v, narm);
	if (std::isnan(m)) return m;
	double x = v[0];
	size_t n = 0;
	for (size_t i=0; i<v.size(); i++) {
		if (!is_NA(v[i])) {
			double d = (v[i] - m);
			x += d * d;
			n++;
		}
	}
	n--;
	if (n==0) return NAN;
	x = sqrt(x / n);
	return x;
}



template <typename T>
double vsdpop(std::vector<T>& v, bool narm) {
	double m = vmean(v, narm);
	if (std::isnan(m)) return m;
	double x = v[0];
	size_t n = 0;
	for (size_t i=0; i<v.size(); i++) {
		if (!is_NA(v[i])) {
			double d = (v[i] - m);
			x += d * d;
			n++;
		}
	}
	x = sqrt(x / n);
	return x;
}



template <typename T>
T vmin(std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(x)) {
					x = v[i];
				} else {
					x = std::min(x, v[i]);
				}
			}
		}
	} else {
		if (is_NA(x)) return x;
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				x = std::min(x, v[i]);
			}
		}
	}
	return x;
}


template <typename T>
T vfirst(std::vector<T>& v, bool narm) {
	if (narm) {
		for (size_t i=0; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				return v[i];
			}
		}
	} 
	return v[0];
}


template <typename T>
T vmax(std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(x)) {
					x = v[i];
				} else {
					x = std::max(x, v[i]);
				}
			}
		}
	} else {
		if (is_NA(x)) return x;
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				x = std::max(x, v[i]);
			}
		}
	}
	return x;
}


template <typename T>
double vwhich(std::vector<T>& v, bool narm) {
	double out;
	for (size_t i=0; i<v.size(); i++) {
		if ((!is_NA(v[i])) && v[i] != 0) {
			out = i+1;
			return out;
		}
	}
	out = NAN;
	return out;
}



template <typename T>
T vwhichmin(std::vector<T>& v, bool narm) {
	T x = v[0];
	T out;
	if (is_NA(x)) {
		out = NA<T>::value;
	} else {
		out = 0;		
	}
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(out)) {
					x = v[i];
					out = i;
				} else if (v[i] < x) {
					x = v[i];
					out = i;
				}
			}
		}
	} else {
		if (is_NA(x)) { return out; }
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				if (v[i] < x) {
					x = v[i];
					out = i;
				}
			}
		}
	}
	if (is_NA(out)) {
		return out;
	} else {
		return (out + 1);  // +1 for R
	}	
}


template <typename T>
T vwhichmax(std::vector<T>& v, bool narm) {

	T x = v[0];
	T out;
	if (is_NA(x)) {
		out = NA<T>::value;
	} else {
		out = 0;		
	}
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(out)) {
					x = v[i];
					out = i;
				} else if (v[i] > x) {
					x = v[i];
					out = i;
				}
			}
		}
	} else {
		if (is_NA(x)) { return out; }
		for (size_t i=0; i<v.size(); i++) {
			if (is_NA(v[i])) {
				return NA<T>::value;
			} else {
				if (v[i] > x) {
					x = v[i];
					out = i;
				}
			}
		}
	}
	if (is_NA(out)) {
		return out;
	} else {
		return (out + 1);  // +1 for R
	}	
}


// problematic; should be ok for int and float but
// won't work with bool values (nodata == 0)
template <typename T>
T vall(std::vector<T>& v, bool narm) {
	T x;
	if (narm) {
		x = 1;
        for (size_t i=0; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (v[i] == 0) {
					x = 0;
					break;
				} else {
					x = 1;
				}
			}
        }
		//x = x < 0 ? NA<T>::value : x;
    } else {
		x = 1;
        for (size_t i=0; i<v.size(); i++) {
            if (is_NA(v[i]) | (v[i] == 0)) {
                x = v[i];
                break;
			}
		}
	}
	return x;
}


template <typename T>
T vany(std::vector<T>& v, bool narm) {
	T x = 0;
	bool hasnd = false;
	for (size_t i=0; i<v.size(); i++) {
		if (!is_NA(v[i])) {
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
		x = NA<T>::value;;
	}
	return x;
}



template <typename T>
std::vector<T> vrange(std::vector<T>& v, bool narm) {
	std::vector<T> x = { v[0], v[0] };

	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(v[i])) {
				if (is_NA(x[0])) {
					x[0] = v[i];
					x[1] = v[i];
				} else {
					x[0] = std::min(x[0], v[i]);
					x[1] = std::max(x[1], v[i]);
				}
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(x[0])) {
				if (is_NA(v[i])) {
					x[0] = NA<T>::value;
					x[1] = NA<T>::value;
				} else {
					x[0] = std::min(x[0], v[i]);
					x[1] = std::max(x[1], v[i]);
				}
			}
		}
	}
	return x;
}



template <typename T>
T vmodal(std::vector<T>& v, bool narm) {

	size_t n = v.size();
    std::vector<unsigned> counts(n, 0);

	std::sort(v.begin(), v.end());

    for (size_t i = 0; i < n; ++i) {
        counts[i] = 0;
        size_t j = 0;
        while ((j < i) && (v[i] != v[j])) {
            ++j;
        }
        ++(counts[j]);
    }
	
    size_t maxCount = 0;
	for (size_t i = 1; i < n; ++i) {
		if (counts[i] > counts[maxCount]) {
			maxCount = i;
		}
	}
	
    return v[maxCount];
}



template <typename T>
std::vector<bool> visna(std::vector<T>& v) {
	std::vector<bool> x(v.size(), false);
	for (size_t i=0; i<v.size(); i++) {
		if (is_NA(v[i])) {
			x[i] = true;
		}
	}
	return x;
}


template <typename T>
std::vector<bool> visnotna(std::vector<T>& v) {
	std::vector<bool> x(v.size(), true);
	for (size_t i=0; i<v.size(); i++) {
		if (is_NA(v[i])) {
			x[i] = false;
		}
	}
	return x;
}




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


#endif

