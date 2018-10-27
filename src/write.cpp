#include "spat.h"
#include "SimpleIni.h"
#include "util.h"
using namespace std;

bool SpatRaster::isSource(std::string filename) {
	std::vector<string> ff = filenames();
	for (size_t i=0; i<ff.size(); i++) {
		if (ff[i] == filename) {
			return true;
		}
	}
	return false;
}


bool SpatRaster::writeRaster(std::string filename, bool overwrite) {
	lrtrim(filename);
	if (filename == "") {
		filename = "random_file_name.grd";
	} else if (file_exists(filename)) {
		if (overwrite) {
			if (isSource(filename)) {
				error = true;
				error_message = "cannot overwrite object to itself";				
			}
			remove(filename.c_str());
		} else {
			error = true;
			error_message = "file exists";
		    return false;
		}
	}

	string ext = getFileExt(filename);
	lowercase(ext);
	
	if (ext == ".grd") {
		ofstream fs(filename, ios::ate | ios::binary);
		std::vector<double> v = getValues();
		fs.write((char*)&v[0], v.size() * sizeof(double));
		fs.close();
        return writeHDR(filename);
	} else {
        #ifdef useGDAL
        return writeRasterGDAL(filename, overwrite);
		#else
		error = true;
		error_message = "gdal is not available";
	    return false;			
        #endif // useGDAL
	}
}


bool SpatRaster::writeStart(std::string filename, bool overwrite) {

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

	} else {
		bool exists = file_exists(filename);
		string ext = getFileExt(filename);
		lowercase(ext);
		if (ext == ".grd") {
			source[0].driver = {"raster"};
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
			source[0].driver = {"gdal"} ;
			success = writeStartGDAL(filename, overwrite);
			#else 
			error = true;
			error_message = "gdal is not available";	
			return false;
			#endif
		}
		
		if (open_write) {
			warning = true;
			warning_message.push_back("file was already open");
		}
		open_write = true;		
	}
	source[0].filename = {filename};
	bs = getBlockSize(4);
	return success;
}



bool SpatRaster::writeValues(std::vector<double> vals, unsigned row){
	if (!open_write) {
		error = true;
		error_message = "cannot write (no open file)";
		return false;
	}
	bool success = true;
	if (source[0].driver == "raster") {
		unsigned size = vals.size();
		//(*fs).write(reinterpret_cast<const char*>(&vals[0]), size*sizeof(double));
		string fname = setFileExt(source[0].filename, ".gri");
		ofstream fs(fname, ios::ate | ios::binary);
		fs.write(reinterpret_cast<const char*>(&vals[0]), size*sizeof(double));
		fs.close();

	} else if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeValuesGDAL(vals, row);
		#else
		error = true;
		error_message = "gdal is not available";	
		return false;
		#endif
	} else {
		setValues(vals);
	}
	return success;
}

bool SpatRaster::writeStop(){
	if (!open_write) {
		error = true;
		error_message = "cannot close a file that is not open";
		return false;
	}
	open_write = false;
	bool success = true;

	if (source[0].driver == "raster") {
		//(*fs).close();
		writeHDR(source[0].filename);
	} else if (source[0].driver == "gdal") {
		#ifdef useGDAL
		success = writeStopGDAL();
		#else
		error = true;
		error_message = "gdal is not available";	
		return false;
		#endif
	}
	source[0].hasValues = true;
	// return this new raster instead?
	return success;
}


void SpatRaster::setValues(std::vector<double> _values) {
	//bool result = false;
	//if (_values.size() == size()) {
		source[0].values = _values;
		source[0].hasValues = true;
		source[0].memory = true;

		// todo clear source...
		setRange();

		source[0].names = std::vector<string> {"layer"};
		//result = true;
	//}
	//return (result);
}



void SpatRaster::setRange() {
	double vmin, vmax;
	int imin, imax;
	// for each layer {
		vector_minmax(source[0].values, vmin, imin, vmax, imax);
		source[0].range_min = std::vector<double> {vmin};
		source[0].range_max = std::vector<double> {vmax};
		source[0].hasRange = std::vector<bool> {true};
	//}
}



bool SpatRaster::writeHDR(string filename) {
	CSimpleIniA ini;
	ini.SetValue("version", NULL, NULL);
	ini.SetValue("version", "version", "2");
	ini.SetValue("georeference", NULL, NULL);
	ini.SetValue("georeference", "xmin", to_string(extent.xmin).c_str());
	ini.SetValue("georeference", "xmax", to_string(extent.xmax).c_str());
	ini.SetValue("georeference", "ymin", to_string(extent.ymin).c_str());
	ini.SetValue("georeference", "ymax", to_string(extent.ymax).c_str());
	ini.SetValue("georeference", "crs", crs.c_str());
	ini.SetValue("dimensions", "nrow", to_string(nrow).c_str());
	ini.SetValue("dimensions", "ncol", to_string(ncol).c_str());
	ini.SetValue("dimensions", "nlyr", to_string(nlyr()).c_str());
	ini.SetValue("dimensions", "names", concatenate(getNames(), std::string(":|:")).c_str());
	ini.SetValue("data", NULL, NULL);
	ini.SetValue("data", "datatype", "FLT8S"); // double
	ini.SetValue("data", "nodata", to_string(-1 * numeric_limits<double>::max()).c_str());
	ini.SetValue("data", "range_min", concatenate(dbl2str(range_min()), std::string(":|:")).c_str());
	ini.SetValue("data", "range_max", concatenate(dbl2str(range_max()), std::string(":|:")).c_str());

	//string f = setFileExt(filename, ".grd");
	SI_Error rc = ini.SaveFile(filename.c_str());
	return rc >= 0 ;
}





/*
bool SpatRaster::writeStartFs(std::string filename, bool overwrite,  fstream& f) {

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
