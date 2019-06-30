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
#include "math_utils.h"
#include "file_utils.h"
#include "hdr.h"


template <typename T>
bool write_bin(std::string filename, std::vector<double> v, unsigned offset,
		std::vector<double> &range_min, std::vector<double> &range_max, 
		unsigned nlyr, std::string bandorder) {

	std::vector<T> values(v.begin(), v.end());
	size_t step = values.size() / nlyr;
	double vmin, vmax;
	
	for (size_t i=0; i<nlyr; i++) {
		size_t start = step * i;
		minmax(values.begin()+start, values.begin()+start+step, vmin, vmax);
		range_min[i] = std::min(range_min[i], vmin);
		range_max[i] = std::max(range_max[i], vmax);
	}

	std::ofstream fs(filename, std::ios::out | std::ios::binary);
	size_t size = sizeof(T);
	//make sure we write at the right place if not in order.
	//also deal with BIL/BIP/BSQ
	//long cursize = fs.tellp() / sizeof(double);
    //long needpos = row * ncol();
    //std::cout << cursize << " " << needpos << "\n";

	//fs.write((char*)&values[0], values.size() * sizeof(T));
	fs.seekp(offset * size);
	fs.write(reinterpret_cast<const char*>(&values[0]),
								values.size()*size);
	fs.close();
	return(true);
}



bool SpatRaster::writeStartBinary(std::string filename, std::string datatype, std::string bandorder, bool overwrite) {

	SpatMessages m = can_write(filename, overwrite);
	if (m.has_error) {
		msg = m;
		return(false);
	}
	source[0].driver = "raster";
	for (size_t i =0; i<nlyr(); i++) {
		source[0].range_min[i] = std::numeric_limits<double>::max();
		source[0].range_max[i] = std::numeric_limits<double>::lowest();
	}
	source[0].filename = setFileExt(filename, ".gri");
//	source[0].fsopen(filename);
//  create file and close it
//  std::string griname = setFileExt(source[0].filename, ".gri");
//  std::ofstream fs(griname, std::ios::out | std::ios::binary);
//  fs.close();
	return true;
}


bool SpatRaster::writeValuesBinary(std::vector<double>& vals, unsigned startrow, unsigned nrows, unsigned startcol, unsigned ncols){

	bool success = true;
    std::string datatype = source[0].datatype;
	std::string fname = source[0].filename;
	unsigned nl = source[0].nlyr;
	std::string bdo = source[0].bandorder;
// Ojo -- writing by row, ignoring startcol and ncols
	size_t offset = startrow * ncol();
	if (datatype == "INT2S") {
		success = write_bin<short>(fname, vals, offset, source[0].range_min, source[0].range_max, nl, bdo);
	} else if (datatype == "INT4S") {
		success = write_bin<long>(fname, vals, offset, source[0].range_min, source[0].range_max, nl, bdo);
	} else if (datatype == "INT8S") {
		success = write_bin<long long int>(fname, vals, offset, source[0].range_min, source[0].range_max, nl, bdo);
	} else if (datatype == "FLT4S") {
		success = write_bin<float>(fname, vals, offset, source[0].range_min, source[0].range_max, nl, bdo);
	} else if (datatype == "FLT8S") {
		success = write_bin<double>(fname, vals, offset, source[0].range_min, source[0].range_max, nl, bdo);
	} else {
		setError("datatype not yet implemented: " + source[0].datatype);
		success = false;
	}
	return success;
}

//bool SpatRaster::writeStopBinary() {
//	return true;
//}


bool SpatRaster::writeRasterBinary(std::string filename, std::string datatype, std::string bandorder, bool overwrite) {
	bool success = true;
    if (!writeStartBinary(filename, datatype, bandorder, overwrite)) return(false); 
	std::vector<double> v = getValues();
		//source[0].fsopen(filename);
		//source[0].fswrite(v);
		//source[0].fsclose();
//	    std::string grifile = setFileExt(filename, ".gri");
//		std::ofstream fs(grifile, std::ios::out | std::ios::binary);
//		fs.write((char*)&v[0], v.size() * sizeof(double));
//		fs.close();
    if (!writeValues(v, 0, nrow(), 0, ncol()))  { return(false); }
    if (!writeStop())  { return(false); }
    if (!writeHDR(filename))  { return(false); }
	return success;
}




bool SpatRaster::writeHDR(std::string filename) {
	std::vector<std::string> s(14);
	s[0] = std::to_string(extent.xmin);
	s[1] = std::to_string(extent.xmax);
	s[2] = std::to_string(extent.ymin);
	s[3] = std::to_string(extent.ymax);
	s[4] = crs;
	s[5] = std::to_string(nrow());
	s[6] = std::to_string(ncol());
	s[7] = std::to_string(nlyr());
	s[8] = concatenate(getNames(), std::string(":|:"));
	s[9] = "FLT8S"; // double
	s[10] = std::to_string(-1 * std::numeric_limits<double>::max());
	s[11] = concatenate(dbl2str(range_min()), std::string(":|:"));
	s[12] = concatenate(dbl2str(range_max()), std::string(":|:"));
	s[13] = filename;
	return hdr_write(s);
}

