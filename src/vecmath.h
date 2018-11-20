// Copyright (c) 2018  Robert J. Hijmans
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

#include <type_traits>
#include <vector>
#include "NA.h"


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
T vmean(std::vector<T>& v, bool narm) {
	T x = 0;
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
			if (!is_NA(x)) {
				if (is_NA(v[i])) {
					x = NA<T>::value;
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
		x = NA<T>::value;
	}
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
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(x)) {
				if (is_NA(v[i])) {
					x = NA<T>::value;
					break;
				} else {
					x = std::min(x, v[i]);
				}
			}
		}
	}
	return x;
}


template <typename T>
T vmax(std::vector<T>& v, bool narm) {
	T x = v[0];
	if (narm) {
		for (size_t i=1; i<v.size(); i++) {
			if (is_NA(v[i])) {
				if (is_NA(x)) {
					x = v[i];
				} else {
					x = std::max(x, v[i]);
				}
			}
		}
	} else {
		for (size_t i=1; i<v.size(); i++) {
			if (!is_NA(x)) {
				if (is_NA(v[i])) {
					x = NA<T>::value;
					break;
				} else {
					x = std::max(x, v[i]);
				}
			}
		}
	}
	return x;
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

