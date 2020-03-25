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

#include "spatRaster.h"
#include "string_utils.h"
#include "file_utils.h"
//#include "hdr.h"

bool SpatRaster::constructFromFile(std::string fname) {

	if (!file_exists(fname)) {
		setError("file does not exist");
		return false;
	}

	std::string ext = getFileExt(fname);
    #ifdef useGDAL
	return constructFromFileGDAL(fname);
	#else 		
	return false;
    #endif // useGDAL
}



bool SpatRaster::constructFromFiles(std::vector<std::string> fnames) {

	SpatRaster r = SpatRaster(fnames[0]);
	setSource(r.source[0]);
	for (size_t i=1; i<fnames.size(); i++) {
		r = SpatRaster(fnames[i]);
		if (!compare_geom(r, false, true, true)) {
			setError("geometry of " + fnames[i] + " does not match previous sources");
			return false;
		} else {
			addSource(r);
		}
	}
	return true;
}

