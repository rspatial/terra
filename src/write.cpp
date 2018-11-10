#include "spatraster.h"
#include "SimpleIni/SimpleIni.h"
#include "string_utils.h"
#include "math_utils.h"

bool SpatRaster::isSource(std::string filename) {
	std::vector<std::string> ff = filenames();
	for (size_t i=0; i<ff.size(); i++) {
		if (ff[i] == filename) {
			return true;
		}
	}
	return false;
}


bool SpatRaster::writeRaster(std::string filename, std::string format, std::string datatype, bool overwrite) {
	lrtrim(filename);
	if (filename == "") {
		filename = "random_file_name.grd";
	} else if (file_exists(filename)) {
		if (overwrite) {
			if (isSource(filename)) {
				setError("cannot overwrite object to itself");
			}
			remove(filename.c_str());
		} else {
			setError("file exists");
		}
	}

	std::string ext = getFileExt(filename);
	lowercase(ext);

	if (ext == ".grd") {
		std::ofstream fs(filename, std::ios::ate | std::ios::binary);
		std::vector<double> v = getValues();
		fs.write((char*)&v[0], v.size() * sizeof(double));
		fs.close();
        return writeHDR(filename);
	} else {
        #ifdef useGDAL
        return writeRasterGDAL(filename, format, datatype, overwrite);
		#else
		setError("GDAL is not available");
	    return false;
        #endif // useGDAL
	}
}


bool SpatRaster::writeStart(std::string filename, std::string format, std::string datatype, bool overwrite) {

//	double inf = std::numeric_limits<double>::infinity();
//	s.min_range = inf;
//	s.max_range = -inf;
	bool success = true;
	lrtrim(filename);
	if (filename == "") {
		if (!canProcessInMemory(4)) {
			filename = "random_file_name.grd";
		}
	}

	if (filename == "") {
		source[0].driver = "memory";
		source[0].memory = true;

	} else {
		source[0].memory = false;
		bool exists = file_exists(filename);
		std::string ext = getFileExt(filename);
		lowercase(ext);
		if (ext == ".grd") {
			source[0].driver = "raster";
			if (exists) {
				if (overwrite) {
					remove(filename.c_str());
				} else {
					// stop()
				}
			}
			//(*fs).open(fname, ios::out | ios::binary);
		} else {
			// open GDAL filestream
			#ifdef useGDAL
			source[0].driver = "gdal" ;
			success = writeStartGDAL(filename, format, datatype, overwrite);
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
	source[0].filename = {filename};
	bs = getBlockSize(4);
	return success;
}



bool SpatRaster::writeValues(std::vector<double> vals, unsigned row){
	if (!source[0].open_write) {
		setError("cannot write (no open file)");
		return false;
	}
	bool success = true;
	if (source[0].driver == "raster") {
		unsigned size = vals.size();
		//(*fs).write(reinterpret_cast<const char*>(&vals[0]), size*sizeof(double));
		std::string fname = setFileExt(source[0].filename, ".gri");
		std::ofstream fs(fname, std::ios::ate | std::ios::binary);
		fs.write(reinterpret_cast<const char*>(&vals[0]), size*sizeof(double));
		fs.close();

	} else if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeValuesGDAL(vals, row);
		#else
		setError("GDAL is not available");
		return false;
		#endif
	} else {
        source[0].values = vals;
        source[0].hasValues = true;
        source[0].memory = true;
        source[0].setRange();
	}
	return success;
}


bool SpatRaster::writeStop(){
	if (!source[0].open_write) {
		setError("cannot close a file that is not open");
		return false;
	}
	source[0].open_write = false;
	bool success = true;

	if (source[0].driver == "raster") {
		//(*fs).close();
		writeHDR(source[0].filename);
	} else if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeStopGDAL();
		#else
		setError("GDAL is not available");
		return false;
		#endif
	}
	source[0].hasValues = true;
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
	unsigned nc = ncol * nrow;
	unsigned start;
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
	CSimpleIniA ini;
	ini.SetValue("version", NULL, NULL);
	ini.SetValue("version", "version", "2");
	ini.SetValue("georeference", NULL, NULL);
	ini.SetValue("georeference", "xmin", std::to_string(extent.xmin).c_str());
	ini.SetValue("georeference", "xmax", std::to_string(extent.xmax).c_str());
	ini.SetValue("georeference", "ymin", std::to_string(extent.ymin).c_str());
	ini.SetValue("georeference", "ymax", std::to_string(extent.ymax).c_str());
	ini.SetValue("georeference", "crs", crs.c_str());
	ini.SetValue("dimensions", "nrow", std::to_string(nrow).c_str());
	ini.SetValue("dimensions", "ncol", std::to_string(ncol).c_str());
	ini.SetValue("dimensions", "nlyr", std::to_string(nlyr()).c_str());
	ini.SetValue("dimensions", "names", concatenate(getNames(), std::string(":|:")).c_str());
	ini.SetValue("data", NULL, NULL);
	ini.SetValue("data", "datatype", "FLT8S"); // double
	ini.SetValue("data", "nodata", std::to_string(-1 * std::numeric_limits<double>::max()).c_str());
	ini.SetValue("data", "range_min", concatenate(dbl2str(range_min()), std::string(":|:")).c_str());
	ini.SetValue("data", "range_max", concatenate(dbl2str(range_max()), std::string(":|:")).c_str());

	//string f = setFileExt(filename, ".grd");
	SI_Error rc = ini.SaveFile(filename.c_str());
	return rc >= 0 ;
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
