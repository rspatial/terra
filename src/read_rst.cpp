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

/*
  C++ for reading gridfiles
  Robert Hijmans
  January 2008
  r.hijmans@gmail.com
*/

#include <vector>
#include <fstream>
//#include <iostream>
#include "spatRaster.h"


/*
https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
#include <climits>
template <typename T>
T swap_endian(T u) {
    static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");
    union {
        T u;
        unsigned char u8[sizeof(T)];
    } source, dest;

    source.u = u;
    for (size_t k = 0; k < sizeof(T); k++) {
        dest.u8[k] = source.u8[sizeof(T) - k - 1];
    }
    return dest.u;
}
swap_endian<double>(42).
*/



std::vector<double> readINT2(std::string file, unsigned long cell, unsigned n) {
	const int dsize = 2;
	std::vector<short> v(n);
	short* value = &v[0];

	std::ifstream f (file, std::ios::in | std::ios::binary);
	f.seekg ( cell * dsize, std::ios::beg);
	f.read ((char*)value, dsize*n);
	f.close();

	std::vector<double> vv(v.begin(), v.end());
	return vv;
}



std::vector<double> readINT4(std::string file, unsigned long cell, unsigned n) {
	const int dsize = 4;
	std::vector<long> v(n);
	long* value = &v[0];

	std::ifstream f (file, std::ios::in | std::ios::binary);
	f.seekg ( cell * dsize, std::ios::beg);
	f.read ((char*)value, dsize*n);
	f.close();

	std::vector<double> vv(v.begin(), v.end());
	return vv;
}



std::vector<double> readFLT4(std::string file, std::string order, unsigned long start, unsigned n) {

	const int dsize = 4;
	size_t nlyr = 1;
//	unsigned lyr = 0;
	std::vector<float> v(n * nlyr);
	float* value = &v[0];

	start = start * dsize;
	n = dsize * n;
	std::ifstream f (file, std::ios::in | std::ios::binary);
	if (order == "BIL") {
		f.seekg (start, std::ios::beg);
		f.read ((char*)value, n*nlyr);
	} 
	
	
//	if (order == "BSQ") {
//		for (size_t i = 0; i < nlyr; i++) {
//			f.seekg (start * lyr, std::ios::beg);
//			f.read ((char*)value, n);
//		}


	f.close();
	std::vector<double> vv(v.begin(), v.end());
	return vv;
}


std::vector<double> readFLT8(std::string file, std::string order, unsigned long start, unsigned n) {

	const int dsize = 8;
	size_t nlyr = 1;
//	unsigned lyr = 0;
	std::vector<double> v(n * nlyr);
	double* value = &v[0];

	start = start * dsize;
	n = dsize * n;
	std::ifstream f (file, std::ios::in | std::ios::binary);
	if (order == "BIL") {
		f.seekg (start, std::ios::beg);
		f.read ((char*)value, n*nlyr);
	} 
	f.close();
	std::vector<double> vv(v.begin(), v.end());
	return vv;
}

/*std::vector<double> readFLT8(std::string file, unsigned long cell, unsigned n, unsigned nlyr, string order) {
	const int dsize = 8;
	std::vector<double> v(n);
	double* value = &v[0];

	ifstream f (file, std::ios::in | std::ios::binary);
	f.seekg ( cell * dsize, std::ios/::beg);
	f.read ((char*)value, dsize*n);
	f.close();

	return v;
}
*/

