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
#include "string_utils.h"
#include "math_utils.h"
#include "hdr.h"
#include "spatRaster.h"

bool canWrite(std::string filename) {
	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == NULL) {
		return false;				
	}
	fclose(fp);
	return true;
}


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
	std::string extension = ".tif";
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

	std::string filename = opt.get_filename();
	std::string datatype = opt.get_datatype();
	bool overwrite = opt.get_overwrite();

	if (filename == "") {
		double seed = std::chrono::system_clock::now().time_since_epoch().count();
		filename = tempFile(opt, seed);
	} else if (file_exists(filename)) {
		if (overwrite) {
			if (isSource(filename)) {
				setError("cannot overwrite object to its own file source");
				return false;
			}
			remove(filename.c_str());
		} else {
			setError("file exists");
			return false;
		}
	}
	if (!canWrite(filename)) {
		setError("cannot write file");
		return false;				
	}

	std::string ext = getFileExt(filename);
	lowercase(ext);

	if (ext == ".grd") {
		std::vector<double> v = getValues();
		//source[0].fsopen(filename);
		//source[0].fswrite(v);
		//source[0].fsclose();

//	    std::string grifile = setFileExt(filename, ".gri");
//		std::ofstream fs(grifile, std::ios::out | std::ios::binary);
//		fs.write((char*)&v[0], v.size() * sizeof(double));
//		fs.close();

        writeStart(opt);
        writeValues(v, 0);
        writeStop();
        return writeHDR(filename);
	} else {
		std::string format = opt.get_filetype();
        #ifdef useGDAL
        return writeRasterGDAL(filename, format, datatype, overwrite);
		#else
		setError("GDAL is not available");
	    return false;
        #endif // useGDAL
	}
}


bool SpatRaster::writeStart(SpatOptions &opt) {

//	double inf = std::numeric_limits<double>::infinity();
//	s.min_range = inf;
//	s.max_range = -inf;
	bool success = true;
	std::string filename = opt.get_filename();
	std::string datatype = opt.get_datatype();

	if (filename == "") {
		if (!canProcessInMemory(4)) {
			double seed = std::chrono::system_clock::now().time_since_epoch().count();
			filename = tempFile(opt, seed);
		}
	}

	if (filename == "") {
		source[0].driver = "memory";
		source[0].memory = true;

	} else {
		source[0].memory = false;
		if (file_exists(filename)) {
			if (opt.get_overwrite()) {
				remove(filename.c_str());
			} else {
				setError("file exists");
				return false;
			}
		}
		if (!canWrite(filename)) {
			setError("cannot write file");
			return false;				
		}
		
		std::string ext = getFileExt(filename);
		lowercase(ext);
		if (ext == ".grd") {
			source[0].driver = "raster";
//			source[0].fsopen(filename);
            // create file and close it
            std::string griname = setFileExt(source[0].filename, ".gri");
            std::ofstream fs(griname, std::ios::out | std::ios::binary);
            fs.close();
		} else {
			// open GDAL filestream
			#ifdef useGDAL
			source[0].driver = "gdal" ;
			std::string format = opt.get_filetype();
			success = writeStartGDAL(filename, format, datatype);
			if (!success) {
				setError("file exists");
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
	//bs = getBlockSize(4);


    #ifdef useRcpp
    //SpatProgressBar pb;
	//pbar = new Progress(bs.n, opt.do_progress(bs.n), pb);
	pbar = new Progress(bs.n, opt.do_progress(bs.n));
	#endif
	return success;
}



bool SpatRaster::writeValues(std::vector<double> &vals, unsigned row) {
	if (!source[0].open_write) {
		setError("cannot write (no open file)");
		return false;
	}
	bool success = true;
	if (source[0].driver == "raster") {
		//source[0].fswrite(vals);

		unsigned sz = vals.size();
		std::string fname = setFileExt(source[0].filename, ".gri");

		// make sure we write at the right place if not in order.
		// also deal with BIL/BIP/BSQ
		std::ofstream fs(fname, std::ios::ate | std::ios::binary);
		long cursize = fs.tellp() / sizeof(double);
        long needpos = row * ncol();
        //std::cout << cursize << " " << needpos << "\n";

		fs.write(reinterpret_cast<const char*>(&vals[0]), sz*sizeof(double));
		fs.close();

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
	if (source[0].driver == "raster") {
		//source[0].fsclose();
		writeHDR(source[0].filename);
	} else if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeStopGDAL();
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
		source[0].setRange();
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
