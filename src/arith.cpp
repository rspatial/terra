using namespace std;
#include <vector>
#include "spat.h"

#include <algorithm>
#include <functional>
#include <cassert>


template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<T>());
    return result;
}


template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<T>());
    return result;
}


template <typename T>
std::vector<T> operator/(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::divides<T>());
    return result;
}

template <typename T>
std::vector<T> operator*(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::multiplies<T>());
    return result;
}


/*
template <typename T>
std::vector<T> operator%(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::modulus<T>());
    return result;
}
*/


// the comparisons do not return what I want
// for NAN == NAN (or x > NAN, etc). Should be NAN, not false
template <typename T>
std::vector<T> operator==(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::equal_to<T>());
    return result;
}

template <typename T>
std::vector<T> operator!=(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::not_equal_to<T>());
    return result;
}

template <typename T>
std::vector<T> operator>=(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::greater<T>());
    return result;
}

template <typename T>
std::vector<T> operator<=(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::less<T>());
    return result;
}


template <typename T>
std::vector<T> operator>(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::greater<T>());
    return result;
}

template <typename T>
std::vector<T> operator<(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::less<T>());
    return result;
}


SpatRaster SpatRaster::arith(SpatRaster x, std::string oper, std::string filename, bool overwrite) {
	SpatRaster out = *this;
	out.values.resize(0);
  	out.writeStart(filename, overwrite);
	readStart();
	x.readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol);
		std::vector<double> b = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol);
		if (oper == "+") {
			a = a + b; 
		} else if (oper == "-") {
			a = a - b; 
		} else if (oper == "*") {
			a = a * b; 
		} else if (oper == "/") {
			a = a / b; 
		} else if (oper == "%") {
			// a = a % b; 
		} else if (oper == "==") {
			a = a == b; 
		} else if (oper == "!=") {
			a = a == b; 
		} else if (oper == ">=") {
			a = a >= b; 
		} else if (oper == "<=") {
			a = a <= b; 
		} else if (oper == ">") {
			a = a > b; 
		} else if (oper == "<") {
			a = a < b; 
		} else {
			// stop
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();	
	x.readStop();	
	return(out);
}



SpatRaster SpatRaster::arith(double x, std::string oper, std::string filename, bool overwrite) {

	SpatRaster out = *this;
	out.values.resize(0);
  	out.writeStart(filename, overwrite);
	readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol);
		if (oper == "+") {
			for(double& d : a)  d += x;
		} else if (oper == "-") {
			for(double& d : a)  d -= x;
		} else if (oper == "*") {
			for(double& d : a)  d *= x;
		} else if (oper == "/") {
			for(double& d : a)  d /= x;
		} else if (oper == "%") {
			for(double& d : a) std::fmod(d,x);
		} else if (oper == "==") {
			for(double& d : a) d = d == x;
		} else if (oper == "!=") {
			for(double& d : a) d = d != x;
		} else if (oper == ">=") {
			for(double& d : a) d = d >= x;
		} else if (oper == "<=") {
			for(double& d : a) d = d <= x;
		} else if (oper == ">") {
			for(double& d : a) d = d > x;
		} else if (oper == "<") {
			for(double& d : a) d = d < x;
		} else {
			// stop
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();		
	return(out);
}


SpatRaster SpatRaster::arith_rev(double x, std::string oper, std::string filename, bool overwrite) {

	SpatRaster out = *this;
	out.values.resize(0);
  	out.writeStart(filename, overwrite);
	readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> a = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol);
		if (oper == "+") {
			for(double& d : a)  d = x + d;
		} else if (oper == "-") {
			for(double& d : a)  d = x - d;
		} else if (oper == "*") {
			for(double& d : a)  d = x * d;
		} else if (oper == "/") {
			for(double& d : a)  d = x / d;
		} else if (oper == "%") {
			for(double& d : a)  std::fmod(x,d);
		} else if (oper == "==") {
			for(double& d : a) d = d == x;
		} else if (oper == "!=") {
			for(double& d : a) d = d != x;
		} else if (oper == ">") {
			for(double& d : a)  d = x > d;
		} else if (oper == "<") {
			for(double& d : a)  d = x < d;
		} else {
			// stop
		}
		out.writeValues(a, out.bs.row[i]);
	}
	out.writeStop();
	readStop();		
	return(out);
}


/*

std::vector<std::vector<double> > matrix(int nrow, int ncol) {
	std::vector<std::vector<double> > m (nrow, std::vector<double>(ncol));
	return(m);
}


int main() {
	std::vector<vector<double> > d = matrix(10, 2);
	std::vector<double> m (10);
	m[1] = 1;
	m[5] = 1;
	d = mask(d, m, 1, 9, false);
	for (int i=0; i < d.size(); i++) {
		for (int j=0; j < d[0].size(); j++) {
			std::cout << ' ' << d[i][j];
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}
*/