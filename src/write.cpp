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
#include "file_utils.h"
#include "string_utils.h"
#include "math_utils.h"


void SpatRaster::fill(double x) {
 	if (source[0].driver == "gdal") {
		#ifdef useGDAL
		fillValuesGDAL(x);
		#endif
	} else if (source[0].driver == "memory") {
		source[0].values.resize(size(), x);
	}
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


SpatRaster SpatRaster::writeRaster(SpatOptions &opt) {

	std::string filename = opt.get_filename();
	bool overwrite = opt.get_overwrite();
	SpatRaster out = geometry();

	if (filename == "") {
		#ifdef useGDAL
		std::string extension = ".tif";
		#else
		std::string extension = ".grd";
		#endif
		filename = tempFile(opt.get_tempdir(), extension);
//	} else if (file_exists(filename)) {
	}
	SpatMessages m = can_write(filename, overwrite);
	if (m.has_error) {
		out.msg = m;
		return out;
	}

	std::string ext = getFileExt(filename);
	lowercase(ext);
	std::string datatype = opt.get_datatype();

	if (opt.names.size() == nlyr()) {
		setNames(opt.names);
	}
	std::string format = opt.get_filetype();
    #ifdef useGDAL
    out = writeRasterGDAL(filename, format, datatype, true, opt);
	#else
	out.setError("GDAL is not available");
    return out;
    #endif
	return out;
}


bool SpatRaster::writeStart(SpatOptions &opt) {

	std::string filename = opt.get_filename();

	if (opt.names.size() == nlyr()) {
		setNames(opt.names);
	}

	if (filename == "") {
		if (!canProcessInMemory(4) || opt.get_todisk()) {
			#ifdef useGDAL
			std::string extension = ".tif";
			#else
			std::string extension = ".grd";
			#endif
			filename = tempFile(opt.get_tempdir(), extension);
		}
	}

	if (filename != "") {

		std::string ext = getFileExt(filename);
		std::string dtype = opt.get_datatype();
		source[0].datatype = dtype;
		bool overwrite = opt.get_overwrite();

		lowercase(ext);
		// open GDAL filestream
		#ifdef useGDAL
		if (! writeStartGDAL(filename, opt.get_filetype(), dtype, overwrite, opt) ) {
			return false;
		}
		#else
		setError("GDAL is not available");
		return false;
		#endif
	}
	if (source[0].open_write) {
		addWarning("file was already open");
	}
	source[0].open_write = true;
	source[0].filename = filename;
	bs = getBlockSize(opt.get_blocksizemp(), opt.get_steps());

    #ifdef useRcpp
	pbar = new Progress(bs.n+2, opt.do_progress(bs.n));
	pbar->increment();
	#endif
	return true;
}



bool SpatRaster::writeValues(std::vector<double> &vals, unsigned startrow, unsigned nrows, unsigned startcol, unsigned ncols) {
	bool success = true;

	if (!source[0].open_write) {
		setError("cannot write (no open file)");
		return false;
	}

	if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeValuesGDAL(vals, startrow, nrows, startcol, ncols);
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
		success = writeValuesMem(vals, startrow, nrows, startcol, ncols);
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


bool SpatRaster::writeValues2(std::vector<std::vector<double>> &vals, unsigned startrow, unsigned nrows, unsigned startcol, unsigned ncols) {
    std::vector<double> vv = flatten(vals);
    return writeValues(vv, startrow, nrows, startcol, ncols);
}

bool SpatRaster::writeStop(){
	if (!source[0].open_write) {
		setError("cannot close a file that is not open");
		return false;
	}
	source[0].open_write = false;
	bool success = true;
	source[0].memory = false;
	if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeStopGDAL();
		//source[0].hasValues = true;
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
   		source[0].setRange();
		source[0].driver = "memory";
		source[0].memory = true;
		if (source[0].values.size() > 0) {
			source[0].hasValues = true;
		}
	}

#ifdef useRcpp
	if (Progress::check_abort()) {
		pbar->cleanup();
		setError("aborted");
		return(false);
	}
	pbar->increment();
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
	if (values.size() == nc * nlyr) {
		for (size_t i=0; i<nlyr; i++) {
			start = nc * i;
			minmax(values.begin()+start, values.begin()+start+nc, vmin, vmax);
			range_min[i] = vmin;
			range_max[i] = vmax;
			hasRange[i] = true;
		}
	}
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
