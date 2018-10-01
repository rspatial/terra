/*
  C++ for reading gridfiles
  Robert Hijmans
  January 2008
  r.hijmans@gmail.com
*/


using namespace std;

#include "spat.h"
#include <vector>
#include <fstream>
#include <iostream>


std::vector<double> readINT2(string file, unsigned long cell, unsigned n) {
	const int dsize = 2;
	std::vector<short> v(n);
	short* value = &v[0];

	ifstream f (file, ios::in | ios::binary);
	f.seekg ( cell * dsize, ios::beg);
	f.read ((char*)value, dsize*n); 
	f.close();

	std::vector<double> vv(v.begin(), v.end());
	return vv;
}



std::vector<double> readINT4(string file, unsigned long cell, unsigned n) {
	const int dsize = 4;
	std::vector<long> v(n);
	long* value = &v[0];

	ifstream f (file, ios::in | ios::binary);
	f.seekg ( cell * dsize, ios::beg);
	f.read ((char*)value, dsize*n); 
	f.close();

	std::vector<double> vv(v.begin(), v.end());
	return vv;
}


 
std::vector<double> readFLT4(string file, string order, unsigned long start, unsigned n, std::vector<unsigned> lyrs) {

	const int dsize = 4;
	size_t nlyr = lyrs.size();
	std::vector<float> v(n * nlyr);
	float* value = &v[0];

	start = start * dsize;
	n = dsize * n;
	ifstream f (file, ios::in | ios::binary);
//	if (order == "BSQ") {
		for (size_t i = 0; i < nlyr; i++) { 
			f.seekg (start * lyrs[i], ios::beg);
			f.read ((char*)value, n); 
		}

	
	f.close();
	std::vector<double> vv(v.begin(), v.end());
	return vv;
}

std::vector<double> readFLT8(string file, unsigned long cell, unsigned n, unsigned nlyr, string order) {
	const int dsize = 8;
	std::vector<double> v(n);
	double* value = &v[0];

	ifstream f (file, ios::in | ios::binary);
	f.seekg ( cell * dsize, ios::beg);
	f.read ((char*)value, dsize*n); 
	f.close();

	return v;
}


