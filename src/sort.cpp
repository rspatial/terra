
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <limits>
#include <cmath>

std::vector<std::size_t> sort_order_nan_a(const std::vector<double> &x){
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (
                std::isnan(x[i]) ? false :
                std::isnan(x[j]) ? true :
                x[i] < x[j]); });
	return p;
}

std::vector<std::size_t> sort_order_nan_d(const std::vector<double> &x){
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (
                std::isnan(x[i]) ? false :
                std::isnan(x[j]) ? true :
                x[i] > x[j]); });
	return p;
}

std::vector<std::size_t> sort_order_nal_a(const std::vector<long> &x){
	long NAL = std::numeric_limits<long>::min();
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (
                x[i] == NAL ? false :
                x[j] == NAL ? true :
                x[i] < x[j]); });
	return p;
}

std::vector<std::size_t> sort_order_nal_d(const std::vector<long> &x){
	long NAL = std::numeric_limits<long>::min();
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (
                x[i] == NAL ? false :
                x[j] == NAL ? true :
                x[i] > x[j]); });
	return p;
}

std::vector<std::size_t> sort_order_nas_a(const std::vector<std::string> &x){
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (
                x[i] == "____NA_+" ? false :
                x[j] == "____NA_+" ? true :
                x[i] < x[j]); });
	return p;
}

std::vector<std::size_t> sort_order_nas_d(const std::vector<std::string> &x){
	std::vector<std::size_t> p(x.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return (
                x[i] == "____NA_+" ? false :
                x[j] == "____NA_+" ? true :
                x[i] > x[j]); });
	return p;
}


