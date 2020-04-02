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

#include <vector>
#include <algorithm>
#include <string>
#include <functional>
#include "vecmath.h"


std::function<double(std::vector<double>&, bool)> getFun(std::string fun) {
	std::function<double(std::vector<double>&, bool)> theFun;
	if (fun == "mean") {
		theFun = vmean<double>;
	} else if (fun == "sum") {
		theFun = vsum<double>;
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
	} else if (fun == "which.min") {
		theFun = vwhichmin<double>;
	} else if (fun == "which.max") {
		theFun = vwhichmax<double>;
	} else if (fun == "any") {
		theFun = vany<double>;
	} else if (fun == "all") {
		theFun = vall<double>;
	} else {
		theFun = vmean<double>;
	}
	return theFun;
}

