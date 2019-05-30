// Copyright (c) 2018-2019  Robert J. Hijmans
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

#include "spatRaster.h"
#include "string_utils.h"


template <typename T>
std::vector<long unsigned> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<long unsigned> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](long unsigned i1, long unsigned i2) {return v[i1] < v[i2];});

  return idx;
}

void make_valid(std::vector<std::string> &s) {
    for (size_t i=0; i<s.size(); i++) {
        lrtrim(s[i]);
        if (s[i] == "") s[i] = "X";
        if (isdigit(s[i][0])) s[i] = "X" + s[i];
//        if ((s[i][0] == ".") & (s[i].size() > 1)) {
//			if (isdigit(s[i][1])) s[i] = "X" + s[i];
//		}
        std::replace(s[i].begin(), s[i].end(), ' ', '.');
    }
}

void make_unique(std::vector<std::string> &s) {
    std::vector<long unsigned> x = sort_indexes(s);
    std::sort(s.begin(), s.end());
    std::vector<std::string> ss = s;
    unsigned j = 1;
    for (size_t i=1; i<s.size(); i++) {
        if (s[i] == s[i-1]) {
            ss[i-1] = s[i-1] + "_" + std::to_string(j);
            ss[i] = s[i] + "_" + std::to_string(j + 1);
            j++;
        } else {
            j = 1;
        }
    }
    for (size_t i=0; i<s.size(); i++) {
        s[x[i]] = ss[i];
    }
}


std::vector<std::string> SpatRaster::getNames() {
	std::vector<std::string> x;
	for (size_t i=0; i<source.size(); i++) {
		x.insert(x.end(), source[i].names.begin(), source[i].names.end());
	}
	return(x);
}


bool SpatRaster::setNames(std::vector<std::string> names) {
	if (names.size() != nlyr()) {
		return false;
	} else {
        make_valid(names);
        make_unique(names);
        size_t begin=0;
        size_t end;
        for (size_t i=0; i<source.size(); i++)	{
            end = begin + source[i].nlyr;
            source[i].names = std::vector<std::string> (names.begin() + begin, names.begin() + end);
            begin = end;
        }
        return true;
	}
}

