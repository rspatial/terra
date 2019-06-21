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
        make_valid_names(names);
        make_unique_names(names);
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

