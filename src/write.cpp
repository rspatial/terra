#include "spat.h"
#include "SimpleIni.h"
#include "util.h"
//#include <fstream>
using namespace std;


bool canProcessInMemory() {
	return true;
}

BlockSize SpatRaster::getBlockSize() {
	BlockSize bs;
		
	if (source.filename[0] == "") {
	// in memory
		bs.row = {0};
		bs.nrows = {nrow};
		bs.n = 1;
		
	} else {
	
	// to be improved, see raster::blockSize
		bs.row = {0, unsigned(floor(nrow/2))};
		bs.nrows = {bs.row[1], nrow-bs.row[1]};
		bs.n = 2;
	}
	return bs;
}


SpatRaster SpatRaster::writeRaster(std::string filename, bool overwrite) {
//	if ((filename == "") { stop}
//	if ((!hasValues) {stop }
	writeStart(filename, overwrite);
	writeValues(getValues(), 0);
	writeStop();
	return SpatRaster(filename);
}

bool SpatRaster::writeStart(std::string filename, bool overwrite) {
	
	lrtrim(filename);
	if (filename == "") {
		if (!canProcessInMemory()) {
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
			source.driver = {"native"};
			string fname = setFileExt(filename, ".gri");
		} else {
			source.driver = {"gdal"} ;			
			// open GDAL filestream		
		}
	}
	
	source.filename = {filename};
	bs = getBlockSize();
	return true;
}

bool SpatRaster::writeStartFs(std::string filename, bool overwrite,  fstream& f) {
	
	lrtrim(filename);
	if (filename == "") {
		if (!canProcessInMemory()) {
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
			source.driver = {"native"};
			string fname = setFileExt(filename, ".gri");
			f.open(fname, ios::out | ios::binary);
			fs = &f;
		} else {
			source.driver = {"gdal"} ;			
			// open GDAL filestream		
		}
	}
	
	source.filename = {filename};
	bs = getBlockSize();
	return true;
}


bool SpatRaster::writeStop(){

	if (source.driver[0] == "native") {
		(*fs).close();
		writeHDR();
	} else if (source.driver[0] == "gdal") {
	
	}
	return true;
}

bool SpatRaster::writeValues(std::vector<double> vals, unsigned row){
	
	if (source.driver[0] == "native") {
		size_t size = vals.size();
		(*fs).write(reinterpret_cast<const char*>(&vals[0]), size*sizeof(double));
	} else if (source.driver[0] == "gdal") {
		// write with gdal

	} else {
		setValues(vals);
	}
	return true;
}


void SpatRaster::setValues(std::vector<double> _values) {
	//bool result = false;
	//if (_values.size() == size()) {
		values = _values;
		hasValues = true;
		source.memory = {true};
		
		// todo clear source...
		setRange();

		names = std::vector<string> {"layer"};
		//result = true;
	//} 
	//return (result);
}



void vector_minmax(std::vector<double> v, double &min, int &imin, double &max, int &imax) {
    std::vector<double>::size_type p=0;
    imax = -1; imin=-1;
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    for (auto &val : v) {
		if (!std::isnan(val)) {
			if (val > max) {
				imax = p;
				max = val;
			}
			if (val < min) {
				imin = p;
				min = val;
			}
		}
        p++;
    }
}


void SpatRaster::setRange() {
	double vmin, vmax;
	int imin, imax;
	// for each layer {
		vector_minmax(values, vmin, imin, vmax, imax); 
		range_min = std::vector<double> {vmin};
		range_max = std::vector<double> {vmax};
		hasRange = std::vector<bool> {true};
	//}
}



bool SpatRaster::writeHDR() { 
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
	ini.SetValue("dimensions", "nlyr", to_string(nlyr).c_str());
	ini.SetValue("dimensions", "names", concatenate(names, std::string(":|:")).c_str());		
	ini.SetValue("data", NULL, NULL);
	ini.SetValue("data", "datatype", "FLT8S"); // double
	ini.SetValue("data", "nodata", to_string(-1 * numeric_limits<double>::max()).c_str());
	ini.SetValue("data", "range_min", concatenate(dbl2str(range_min), std::string(":|:")).c_str());
	ini.SetValue("data", "range_max", concatenate(dbl2str(range_max), std::string(":|:")).c_str());

	string f = setFileExt(source.filename[0], ".grd");
	SI_Error rc = ini.SaveFile(f.c_str());
	if (rc < 0) {
		return false;
	} else {
		return true;
	}	
}

