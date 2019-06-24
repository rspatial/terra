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

#include <stddef.h>

template <typename T>
void recycle(std::vector<T> &v, unsigned n) {
	size_t s = v.size();
	v.resize(n);
	for (size_t i=s; i<n; i++) {
		v[i] = v[i % s];
	}
}


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

