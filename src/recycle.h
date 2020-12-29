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

#include <stddef.h>

template <typename T>
void recycle(std::vector<T> &v, unsigned n) {
	size_t s = v.size();
	v.resize(n);
	for (size_t i=s; i<n; i++) {
		v[i] = v[i % s];
	}
}


/*
template <typename T>
void recycle(std::vector<T> &v, unsigned n) {
	size_t s = v.size();
	if (s > n) {
		v.resize(n);
	} else if (s < n) {
		v.reserve(n);
		for (size_t i=s; i<n; i++) {
			v.push_back(v[i % s]);
		}
	}
}
*/

template <typename T>
void recycle(std::vector<T> &x, std::vector<T> &y) {
	size_t xsize = x.size();
	size_t ysize = y.size();
	if (xsize != ysize) {
		size_t n = std::max(xsize, ysize);
		if (xsize > ysize) {
			y.resize(n);
			for (size_t i=ysize; i<n; i++) {
				y[i] = y[i % ysize];
			} 				
		} else {
			x.resize(n);
			for (size_t i=xsize; i<n; i++) {
				x[i] = x[i % xsize];
			} 				
		}
	}
}


template <typename T>
void rep(std::vector<T> &v, unsigned n) {
	size_t s = v.size();
	v.reserve(n * s);
	for (size_t i=1; i<n; i++) {
		for (size_t j=0; j<s; j++) {
			v.push_back(v[j]);
		}
	}
}

template <typename T>
void rep_each(std::vector<T> &v, unsigned n) {
	std::vector<T> vv = v;
	size_t s = v.size();
	v.resize(0);
	v.reserve(n * s);
	for (size_t j=0; j<s; j++) {
		for (size_t i=0; i<n; i++) {
			v.push_back(vv[j]);
		}
	}
}

template <typename T>
void rep_each_vect(std::vector<T> &v, std::vector<size_t> n) {
	std::vector<T> vv = v;
	v.resize(0);
	size_t nsum = 0;
	for (size_t i=0; i<n.size(); i++) nsum += n[i];
	v.reserve(nsum);
	for (size_t i=0; i<vv.size(); i++) {
		for (size_t j=0; j<n[i]; j++) {
			v.push_back(vv[i]);
		}
	}
}

template <typename T>
std::vector<T> seq(T start, T end, T increment) {
	size_t s = floor((end - start) / increment);
	std::vector<T> out;
	out.reserve(s);
	for (size_t i=0; i<=s; i++) {
		out.push_back(start + i * increment);
	}
	return out;
}

