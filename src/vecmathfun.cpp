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
#include <functional>

#include "vecmath.h"


std::function<double(std::vector<double>&, bool)> getFun(std::string fun) {
	std::function<double(std::vector<double>&, bool)> agFun;
	if (fun == "mean") {
		agFun = vmean<double>;
	} else if (fun == "sum") {
		agFun = vsum<double>;
	} else if (fun == "min") {
		agFun = vmin<double>;
	} else if (fun == "max") {
		agFun = vmax<double>;
	} else if (fun == "median") {
		agFun = vmedian<double>;
	} else if (fun == "modal") {
		agFun = vmodal<double>;
	} else {
		agFun = vmean<double>;
	}
	return agFun;
}

