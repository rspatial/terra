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


#include <vector>
#include "spatRaster.h"


void SpatRaster::createCategories(unsigned layer) {
	// subset to layer
	SpatOptions opt;
	std::vector<unsigned> lyrs(1, layer);
	SpatRaster r = subset(lyrs, opt);
	std::vector<std::vector<double>> u = r.unique(false);
    std::vector<unsigned> sl = findLyr(layer);
	source[sl[0]].cats[sl[1]].id = u[0];
	source[sl[0]].hasCategories[sl[1]] = true;
}


std::vector<bool> SpatRaster::hasCategories() {
	std::vector<bool> b(nlyr());
	std::vector<unsigned> ns = nlyrBySource();
	unsigned k = 0;
	for (size_t i=0; i<source.size(); i++) {
		for (size_t j=0; j<ns[i]; j++) {
			b[k] = source[i].hasCategories[j];
			k++;
		}
	}
	return b;
}



void SpatRaster::setCategoriesDF(unsigned layer, SpatDataFrame df) {
    std::vector<unsigned> sl = findLyr(layer);
	if (df.nrow() == source[sl[0]].cats[sl[1]].id.size()) {
		source[sl[0]].cats[sl[1]].df = df;
	} else {
		setError("nrow(df) does not match number of categories");
	}
}


SpatCategories SpatRaster::getLayerCategories(unsigned layer) {
    std::vector<unsigned> sl = findLyr(layer);
	SpatCategories cat = source[sl[0]].cats[sl[1]];
	return cat;
}

std::vector<SpatCategories> SpatRaster::getCategories() {
	std::vector<SpatCategories> cats;
	for (size_t i=0; i<source.size(); i++) {
		cats.insert(cats.end(), source[i].cats.begin(), source[i].cats.end());
	}
	return cats;
}

