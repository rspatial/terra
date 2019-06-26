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

#include <random>
#include <chrono>
#include "spatRaster.h"
#include "file_utils.h"
#include "string_utils.h"
#include "math_utils.h"
#include "hdr.h"

std::string tempFile(SpatOptions &opt, double seed) {
    std::vector<char> characters = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K',
    'L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m',
    'n','o','p','q','r','s','t','u','v','w','x','y','z' };
    std::default_random_engine generator(std::random_device{}());
    generator.seed(seed);
    std::uniform_int_distribution<> distrib(0, characters.size()-1);
    auto draw = [ characters, &distrib, &generator ]() {
		return characters[ distrib(generator) ];
	};
    std::string filename(15, 0);
    std::generate_n(filename.begin(), 15, draw);
	std::string tmpdir = opt.get_tempdir();
	#ifdef useGDAL
	std::string extension = ".tif";
	#else
	std::string extension = ".grd";	
	#endif
	filename = tmpdir + "/spat_" + filename + extension;
	return filename;
}


bool SpatRaster::isSource(std::string filename) {
	std::vector<std::string> ff = filenames();
	for (size_t i=0; i<ff.size(); i++) {
		if (ff[i] == filename) {
			return true;
		}
	}
	return false;
}


bool SpatRaster::writeRaster(SpatOptions &opt) {

	bool success = true;
	std::string filename = opt.get_filename();
	bool overwrite = opt.get_overwrite();

	
	if (filename == "") {
		double seed = std::chrono::system_clock::now().time_since_epoch().count();
		filename = tempFile(opt, seed);
	} else if (file_exists(filename)) {
		SpatMessages m = can_write(filename, overwrite);
		if (m.has_error) {
			msg = m;
			return false;
		}
	}

	std::string ext = getFileExt(filename);
	lowercase(ext);
	std::string datatype = opt.get_datatype();

	SpatRaster out = geometry(nlyr());
	if (ext == ".grd") {
		std::string bandorder = opt.get_bandorder();
	    if (!out.writeRasterBinary(filename, datatype, bandorder, overwrite)) {
			msg = out.msg; 
			return(false);
		}
		setSources(out.source);

	} else {
		std::string format = opt.get_filetype();
        #ifdef useGDAL
        return writeRasterGDAL(filename, format, datatype, overwrite);
		#else
		setError("GDAL is not available");
	    return false;
        #endif
	}
	return success;
}


bool SpatRaster::writeStart(SpatOptions &opt) {

	std::string filename = opt.get_filename();

	if (filename == "") {
		if (!canProcessInMemory(4)) {
			double seed = std::chrono::system_clock::now().time_since_epoch().count();
			filename = tempFile(opt, seed);
		}
	}

	if (filename != "") {

		std::string ext = getFileExt(filename);
		std::string dtype = opt.get_datatype();
		source[0].datatype = dtype;
		bool overwrite = opt.get_overwrite();

		lowercase(ext);
		if (ext == ".grd") {
			std::string bandorder = opt.get_bandorder();
			if (! writeStartBinary(filename, dtype, bandorder, overwrite) ) {
				return false; 
			}
		} else {
			// open GDAL filestream
			#ifdef useGDAL
			if (! writeStartGDAL(filename, opt.get_filetype(), dtype, overwrite) ) {
				return false; 
			}
			#else
			setError("GDAL is not available");
			return false;
			#endif
		}
	}
	if (source[0].open_write) {
		addWarning("file was already open");
	}
	source[0].open_write = true;
	source[0].filename = filename;
	bs = getBlockSize(opt.get_blocksizemp());

    #ifdef useRcpp
	pbar = new Progress(bs.n, opt.do_progress(bs.n));
	#endif
	return true;
}


bool SpatRaster::writeValues(std::vector<double> &vals, unsigned row) {
	bool success = true;
	
	if (!source[0].open_write) {
		setError("cannot write (no open file)");
		return false;
	}
	if (source[0].driver == "raster") {

		success = writeValuesBinary(vals, row);

	} else if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeValuesGDAL(vals, row);
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
		size_t sz = source[0].values.size();
		size_t start = row * ncol() * nlyr();
		if (sz == 0) { // first or all
			source[0].values = vals;
		} else if (sz == start) { // in chunks
			source[0].values.insert(source[0].values.end(), vals.begin(), vals.end());
		} else { // async
			if (start+vals.size() > sz) {
				source[0].values.resize(start+vals.size());
			}
			for (size_t i=0; i<vals.size(); i++) {
				source[0].values[start+i] = vals[i];
			}
		}
	}

#ifdef useRcpp
	if (Progress::check_abort()) {
		pbar->cleanup();
		setError("aborted");
		return(false);
	}
	pbar->increment();
#endif

	return success;
}

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

bool SpatRaster::writeValues2(std::vector<std::vector<double>> &vals, unsigned row){
    std::vector<double> vv = flatten(vals);
    return writeValues(vv, row);
}

bool SpatRaster::writeStop(){
	if (!source[0].open_write) {
		setError("cannot close a file that is not open");
		return false;
	}
	source[0].open_write = false;
	bool success = true;
	source[0].memory = false;
	if (source[0].driver == "raster") {
		//source[0].fsclose();
		writeHDR(source[0].filename);
		setRange();
	} else if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeStopGDAL();
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
		source[0].setRange();
		source[0].driver = "memory";
		source[0].memory = true;
	}
	source[0].hasValues = true;

#ifdef useRcpp
	//pbar->cleanup();
	delete pbar;
#endif

	return success;
}


bool SpatRaster::setValues(std::vector<double> _values) {
	bool result = false;
	SpatRaster g = geometry();

	if (_values.size() == g.size()) {

		RasterSource s = g.source[0];
		s.values = _values;
		s.hasValues = true;
		s.memory = true;
		s.names = getNames();
		s.driver = "memory";

		s.setRange();
		setSource(s);
		result = true;
	} else {
		setError("incorrect number of values");
	}
	return (result);
}

void SpatRaster::setRange() {
	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].memory) { // for now. should read from files as needed
			source[i].setRange();
		}
	}
}

void RasterSource::setRange() {
	double vmin, vmax;
	size_t nc = ncol * nrow;
	size_t start;
	range_min.resize(nlyr);
	range_max.resize(nlyr);
	hasRange.resize(nlyr);
	for (size_t i=0; i<nlyr; i++) {
		start = nc * i;
		minmax(values.begin()+start, values.begin()+start+nc, vmin, vmax);
		range_min[i] = vmin;
		range_max[i] = vmax;
		hasRange[i] = true;
	}
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





/*
bool SpatRaster::writeStartFs(std::string filename, std::string format, std::string datatype, bool overwrite,  fstream& f) {

	lrtrim(filename);
	if (filename == "") {
		if (!canProcessInMemory(4)) {
			filename = "random_file_name.grd";
		}
	}

	if (filename == "") {
		source.driver = {"memory"};

	} else {
		// if (!overwrite) check if file exists
		string ext = getFileExt(filename);
		lowercase(ext);
		if (ext == ".grd") {
			source.driver = {"raster"};
			string fname = setFileExt(filename, ".gri");
			f.open(fname, ios::out | ios::binary);
			(*fs).open(fname, ios::out | ios::binary);
			fs = &f;
		} else {
			source.driver = {"gdal"} ;
			// open GDAL filestream
		}
	}

	source.filename = {filename};
	bs = getBlockSize(4);
	return true;
}
*/
