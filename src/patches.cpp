// Copyright (c) 2018-2025  Robert J. Hijmans
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

void patch_search(const std::vector<double> &m, std::vector<double> &patches, const int &i, const long &ncol, const int &patch, const size_t &dirs) {
// DFS

	std::vector<long> directions;
	if (dirs==4) {
		directions = {-ncol, ncol, -1, 1};		
	} else {
		directions = {-ncol, ncol, -1, 1, -ncol-1, -ncol+1, ncol-1, ncol+1};
	}
	
	size_t ncell = m.size();	
    patches[i] = patch; 
    for (size_t d=0; d<dirs; d++) {
        size_t j = i + directions[d];
        if (j >= 0 && j < ncell && (!std::isnan(m[j])) && std::isnan(patches[j]) && m[j] == m[i]) {
            patch_search(m, patches, j, ncol, patch, dirs); 
        }
    }
}


SpatRaster SpatRaster::patches(size_t dirs, SpatOptions &opt) {

	SpatRaster out = geometry(1, false);
	if (!hasValues()) {
		out.setError("cannot compute surfaceArea for a raster with no values");
		return out;
	}
	if (nlyr() != 1) {
		out.setError("can only compute surfaceArea for a single raster layer");
		return out;		
	}
	if (!((dirs == 4) || (dirs == 8))) {
		out.setError("directions should be 4 or 8");
		return out;		
	}

	if (!canProcessInMemory(opt)) {
		out.setError("cannot do this for large rasters");
		return out;		
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

    size_t patch = 1; 
	size_t nc = ncol();
	std::vector<double> v;

/*
	std::vector<double> patches(out.bs.nrows[0] * nc, NAN);
	for (size_t i = 0; i < out.bs.n; i++) {
		if (i > 0) {
			readValues(v, out.bs.row[i]-1, out.bs.nrows[i]+1, 0, nc);

			std::vector<double> old_p(patches.end()-nc, patches.end());
			patches = std::vector<double>(v.size(), NAN);

			for (size_t j=0; j<nc; j++) {
				if (!std::isnan(v[j])) {
					patch_search(v, patches, j, nc, old_p[j], dirs); 
				}
			}
			patches.erase(patches.begin(), patches.begin()+nc);
			v.erase(v.begin(), v.begin()+nc);
		} else {
			readBlock(v, out.bs, i);
		}
		for (size_t j=0; j<v.size(); j++) {
			if ((!std::isnan(v[j])) && std::isnan(patches[j])) {
				patch_search(v, patches, j, nc, patch, dirs); 
				patch++; 
			}
		}
		
		if (!out.writeBlock(patches, i)) return out;
	}
*/

	std::vector<double> patches(nrow() * nc, NAN);
	readValues(v, 0, nrow(), 0, nc);
	for (size_t j=0; j<v.size(); j++) {
		if ((!std::isnan(v[j])) && std::isnan(patches[j])) {
			patch_search(v, patches, j, nc, patch, dirs); 
			patch++; 
		}
	}
	if (!out.writeValues(patches, 0, nrow())) return out;
	
	readStop();
	out.writeStop();
	return(out);
}

