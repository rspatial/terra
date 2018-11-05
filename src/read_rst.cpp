/*
  C++ for reading gridfiles
  Robert Hijmans
  January 2008
  r.hijmans@gmail.com
*/

#include <vector>
#include <fstream>
#include <iostream>
#include "spatraster.h"


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
	unsigned lyr = 0;
	std::vector<float> v(n * nlyr);
	float* value = &v[0];

	start = start * dsize;
	n = dsize * n;
	std::ifstream f (file, std::ios::in | std::ios::binary);
//	if (order == "BSQ") {
		for (size_t i = 0; i < nlyr; i++) { 
			f.seekg (start * lyr, std::ios::beg);
			f.read ((char*)value, n); 
		}

	
	f.close();
	std::vector<double> vv(v.begin(), v.end());
	return vv;
}


std::vector<double> readFLT8(std::string file, std::string order, unsigned long start, unsigned n) {

	const int dsize = 8;
	size_t nlyr = 1;
	unsigned lyr = 0;
	std::vector<double> v(n * nlyr);
	double* value = &v[0];

	start = start * dsize;
	n = dsize * n;
	std::ifstream f (file, std::ios::in | std::ios::binary);
//	if (order == "BSQ") {
		for (size_t i = 0; i < nlyr; i++) { 
			f.seekg (start * lyr, std::ios::beg);
			f.read ((char*)value, n); 
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

