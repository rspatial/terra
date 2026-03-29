// Copyright (c) 2018-2026  Robert J. Hijmans
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


