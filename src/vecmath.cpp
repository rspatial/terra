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
#include <algorithm>
#include <string>
#include <functional>
#include <cmath>
#include <numeric>
#include "vecmath.h"


bool haveFun(std::string fun) {
	std::vector<std::string> f {"sum", "mean", "median", "modal", "which", "which.min", "which.max", "min", "max", "prod", "any", "all", "none", "sd", "std", "first", "expH"};
	auto it = std::find(f.begin(), f.end(), fun);
	if (it == f.end()) {
		return false;
	}
	return true;
}

std::function<double(std::vector<double>&, bool)> getFun(std::string fun) {
	std::function<double(std::vector<double>&, bool)> theFun;
	if (fun == "mean") {
		theFun = vmean<double>;
	} else if (fun == "sum") {
		theFun = vsum<double>;
	} else if (fun == "sum2") {
		theFun = vsum2<double>;
	} else if (fun == "min") {
		theFun = vmin<double>;
	} else if (fun == "max") {
		theFun = vmax<double>;
	} else if (fun == "median") {
		theFun = vmedian<double>;
	} else if (fun == "modal") {
		theFun = vmodal<double>;
	} else if (fun == "prod") {
		theFun = vprod<double>;
	} else if (fun == "which") {
		theFun = vwhich<double>;
	} else if (fun == "which.min") {
		theFun = vwhichmin<double>;
	} else if (fun == "which.max") {
		theFun = vwhichmax<double>;
	} else if (fun == "any") {
		theFun = vany<double>;
	} else if (fun == "none") {
		theFun = vnone<double>;
	} else if (fun == "all") {
		theFun = vall<double>;
	} else if (fun == "sd") {
		theFun = vsd<double>;
	} else if (fun == "std") {
		theFun = vsdpop<double>;
	} else if (fun == "first") {
		theFun = vfirst<double>;
	} else {
		theFun = vmean<double>;
	}
	return theFun;
}


static bool vselect_better(double a, size_t ia, double b, size_t ib, bool low) {
	if (a != b) {
		return low ? (a < b) : (a > b);
	}
	return ia < ib;
}


std::vector<size_t> vselect_index(const std::vector<double>& v, size_t n, bool low) {
	std::vector<size_t> idx;
	idx.reserve(v.size());
	for (size_t i=0; i<v.size(); i++) {
		if (!std::isnan(v[i])) {
			idx.push_back(i);
		}
	}
	if ((n == 0) || idx.empty()) {
		return std::vector<size_t>();
	}
	if (n > idx.size()) {
		n = idx.size();
	}
	std::nth_element(idx.begin(), idx.begin()+n, idx.end(),
		[&](size_t i, size_t j) { return vselect_better(v[i], i, v[j], j, low); });
	idx.resize(n);
	return idx;
}


std::vector<double> vselect(const std::vector<double>& v, size_t n, bool low) {
	std::vector<size_t> idx = vselect_index(v, n, low);
	std::vector<double> out;
	out.reserve(idx.size());
	for (size_t i=0; i<idx.size(); i++) {
		out.push_back(v[idx[i]]);
	}
	return out;
}


void vselect_keep(std::vector<double>& val, std::vector<size_t>& cell, size_t n, bool low) {
	size_t m = val.size();
	if ((n >= m) || (m == 0) || (n == 0)) {
		if (n == 0) {
			val.clear();
			cell.clear();
		}
		return;
	}
	std::vector<size_t> ord(m);
	std::iota(ord.begin(), ord.end(), 0);
	std::nth_element(ord.begin(), ord.begin()+n, ord.end(),
		[&](size_t i, size_t j) { return vselect_better(val[i], cell[i], val[j], cell[j], low); });
	std::vector<double> nv;
	std::vector<size_t> nc;
	nv.reserve(n);
	nc.reserve(n);
	for (size_t i=0; i<n; i++) {
		nv.push_back(val[ord[i]]);
		nc.push_back(cell[ord[i]]);
	}
	val.swap(nv);
	cell.swap(nc);
}


bool ball(const std::vector<bool>& v) {
    for (size_t i=0; i<v.size(); i++) {
		if (!v[i]) return false;
	}
	return true;
}

bool bany(const std::vector<bool>& v) {
    for (size_t i=0; i<v.size(); i++) {
		if (v[i]) return true;
	}
	return false;
}

bool bnone(const std::vector<bool>& v) {
    for (size_t i=0; i<v.size(); i++) {
		if (v[i]) return false;
	}
	return true;
}
