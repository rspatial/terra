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



template <typename T>
std::vector<T> operator%(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
//    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::modulus<T>());
	for (size_t i=0; i<a.size(); i++) {
		result[i] = std::fmod(a[i], b[i]);
	}
    return result;
}




template <typename T>
std::vector<T> operator==(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::equal_to<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
    return result;
}

template <typename T>
std::vector<T> operator!=(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::not_equal_to<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
    return result;
}

template <typename T>
std::vector<T> operator>=(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::greater<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
    return result;
}

template <typename T>
std::vector<T> operator<=(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::less<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
    return result;
}


template <typename T>
std::vector<T> operator>(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::greater<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
    return result;
}

template <typename T>
std::vector<T> operator<(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::less<T>());
	for (size_t i=0; i<a.size(); i++) {
		if (std::isnan(a[i]) || std::isnan(b[i])) {
			result[i] = NAN;
		} 
	}
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
			 a = a % b; 
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
		if (std::isnan(x)) {
			for(double& d : a)  d = NAN;
		} else if (oper == "+") {
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
			for(double& d : a) if (!std::isnan(d)) d = d == x;
		} else if (oper == "!=") {
			for(double& d : a) if (!std::isnan(d)) d = d != x;
		} else if (oper == ">=") {
			for(double& d : a) if (!std::isnan(d)) d = d >= x;
		} else if (oper == "<=") {
			for(double& d : a) if (!std::isnan(d)) d = d <= x;
		} else if (oper == ">") {
			for(double& d : a) if (!std::isnan(d)) d = d > x;
		} else if (oper == "<") {
			for(double& d : a) if (!std::isnan(d)) d = d < x;
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
		} else if (oper == ">=") {
			for(double& d : a) d = x >= d;
		} else if (oper == "<=") {
			for(double& d : a) d = x <= d;			
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


