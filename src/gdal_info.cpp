#include <vector>
#include <string>
#include "string_utils.h"
#include "file_utils.h"
#include "spatRaster.h"

#include "cpl_port.h"
#include "cpl_conv.h" // CPLFree()
#include "gdal_version.h"

#include "gdal_errors.h" 

// code adapted from the 'sf' package by Edzer Pebesma et al

#if (!(GDAL_VERSION_MAJOR == 2 && GDAL_VERSION_MINOR < 1))
# include "gdal_utils.h" // requires >= 2.1

std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> openopts) {
	
	std::string ret = "";
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(openopts); // open options
	GDALInfoOptions* opt = GDALInfoOptionsNew(options_char.data(), NULL);
	GDALDatasetH ds = GDALOpenEx((const char *) filename.c_str(), GA_ReadOnly, NULL, oo_char.data(), NULL);
	if (ds == NULL) return ret; // #nocov
	char *ret_val = GDALInfo(ds, opt);
	ret = ret_val;
	CPLFree(ret_val);
	GDALInfoOptionsFree(opt);
	GDALClose(ds);
	return ret;
}
#else
std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> oo) {
	std::string out = "GDAL version >= 2.1 required for gdalinfo");
	return out;
}

#endif



std::vector<std::vector<std::string>> sdinfo(std::string fname) {
	std::vector<std::vector<std::string>> out(5);
	GDALDataset *poDataset;
    poDataset = (GDALDataset *) GDALOpen( fname.c_str(), GA_ReadOnly );
    if( poDataset == NULL ) {
		if (!file_exists(fname)) {
			out[0] = std::vector<std::string> {"no such file"};
		} else {
			out[0] = std::vector<std::string> {"cannot open file"};
		}
		return out;
	}
	char **metadata = poDataset->GetMetadata("SUBDATASETS");
	if (metadata == NULL) {
		out[0] = std::vector<std::string> {"no subdatasets"};
		GDALClose( (GDALDatasetH) poDataset );	
		return out;		
	}
	std::vector<std::string> meta;
	for (size_t i=0; metadata[i] != NULL; i++) {
		meta.push_back(metadata[i]);
	}
	if (meta.size() == 0) {
		GDALClose( (GDALDatasetH) poDataset );	
		out[0] = std::vector<std::string> {"no subdatasets"};
		return out;
	}	
	SpatRaster sub;
	std::vector<std::string> name, desc, nr, nc, nl;
	std::string ndelim = "NAME=";
	std::string ddelim = "DESC=";

    for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(ndelim);
		if (pos != std::string::npos) {
			s.erase(0, pos + ndelim.length());
			name.push_back(s);
			if (sub.constructFromFile(s, -1, "")) {
				nr.push_back( std::to_string(sub.nrow()));
				nc.push_back(std::to_string(sub.ncol()));
				nl.push_back(std::to_string(sub.nlyr()));
			}
		} else {
			size_t pos = s.find(ddelim);
			if (pos != std::string::npos) {
				s.erase(0, pos + ddelim.length());
				desc.push_back(s);
			} else {
				desc.push_back("");				
			}
		}
	}
	GDALClose( (GDALDatasetH) poDataset );
	out[0] = name;
	out[1] = desc;
	out[2] = nr;
	out[3] = nc;
	out[4] = nl;
	return out;
}




