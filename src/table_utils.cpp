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

#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include "spatLookup.h"


std::map<double, size_t> table(std::vector<double> &v) {
	SpatFrequencyTable<double> count;
	for (double val : v) {
		if(!std::isnan(val)) count[val]++;
	}
	return std::map<double, size_t>(count.begin(), count.end());
}


std::map<double, size_t> combine_tables(std::map<double, size_t> &x, std::map<double, size_t> &y) {
	for(auto p : y) {
		x[p.first] += p.second;
	}
	return(x);
}


std::vector<double> table2vector(std::map<double, size_t> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}
	out[0].insert(out[0].end(), out[1].begin(), out[1].end());
	return out[0];
}

std::vector<std::vector<double>> table2vector2(std::map<double, size_t> &x) {
	std::vector<std::vector<double>> out(2);
	for( auto p : x ) {
		out[0].push_back(p.first);
		out[1].push_back(p.second);
	}
	return out;
}
