#include <vector>
#include <string>
#include "string_utils.h"
#include "file_utils.h"
#include "spatRaster.h"

#include "cpl_port.h"
#include "cpl_conv.h" // CPLFree()
//#include "gdal_version.h"


std::string sectostr(int x) {
	char buffer[20];
	time_t now = x;
	tm *utc = gmtime(&now);
	strftime (buffer, 20, "%Y-%m-%d %H:%M:%S", utc); 
	std::string s = buffer;
	return s;
}


std::vector<std::string> get_metadata(std::string filename) {
	
	std::vector<std::string> out;
	
    GDALDataset *poDataset;
    poDataset = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly );
    if( poDataset == NULL )  {
		return out;
	}
	
	char **m = poDataset->GetMetadata();
	while (*m != nullptr) {
		out.push_back(*m++);
	}
	
	GDALClose( (GDALDatasetH) poDataset );
	return out;	
}


std::vector<std::string> get_metadata_sds(std::string filename) {
	std::vector<std::string> meta;
    GDALDataset *poDataset;
    poDataset = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly );
    if( poDataset == NULL )  {
		return meta;
	}
	char **metadata = poDataset->GetMetadata("SUBDATASETS");
	if (metadata != NULL) {
		for (size_t i=0; metadata[i] != NULL; i++) {
			meta.push_back(metadata[i]);
		}
	}
	GDALClose( (GDALDatasetH) poDataset );	
	return meta;		
}

std::vector<std::vector<std::string>> parse_metadata_sds(std::vector<std::string> meta) {
	
	std::vector<std::string> name, var, desc, nr, nc, nl;
	std::string ndelim = "NAME=";
	std::string ddelim = "DESC=";

    for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(ndelim);
		if (pos != std::string::npos) {
			s.erase(0, pos + ndelim.length());
			name.push_back(s);
			std::string vdelim = "\":";
			size_t pos = s.find(vdelim);

			if (pos != std::string::npos) {
				s.erase(0, pos + vdelim.length());
				var.push_back(s);
			}
		} else {
			size_t pos = s.find(ddelim);
			if (pos != std::string::npos) {
				s.erase(0, pos + ddelim.length());
				pos = s.find("]");
				std::string dims = s.substr(1, pos-1);
				
				std::vector<std::string> d = strsplit(dims, "x");
				if (d.size() < 2) {
					nl.push_back("0");	
					nr.push_back("0");	
					nc.push_back("0");		
				} else if (d.size() == 2) {
					nl.push_back("1");	
					nr.push_back(d[0]);	
					nc.push_back(d[1]);	
				} else {
					size_t ds = d.size()-1;
					size_t nls = stoi(d[ds-2]);
					for (size_t i=0; i<(ds-2); i++) {
						nls *= stoi(d[i]);
					}
					nl.push_back(std::to_string(nls));	
					nr.push_back(d[ds-1]);	
					nc.push_back(d[ds]);	
				}
				//desc.push_back(std::string(pos, s.size()));
				s = s.substr(pos+2, s.size());
				pos = s.find(" ");
				s = s.substr(0, pos);
				desc.push_back(s);
				

			//	nr.push_back( std::to_string(sub.nrow()));
			//	nc.push_back(std::to_string(sub.ncol()));
			//	nl.push_back(std::to_string(sub.nlyr()));

			} else {
				desc.push_back("");				
			}
		}
	}
	std::vector<std::vector<std::string>> out(6);
	out[0] = name;
	out[1] = var;
	out[2] = desc;
	out[3] = nr;
	out[4] = nc;
	out[5] = nl;
	return out;
}



std::vector<std::vector<std::string>> sdinfo(std::string fname) {

	std::vector<std::vector<std::string>> out(6);
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
	std::vector<std::string> name, var, desc, nr, nc, nl;
	std::string ndelim = "NAME=";
	std::string ddelim = "DESC=";

    for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(ndelim);
		if (pos != std::string::npos) {
			s.erase(0, pos + ndelim.length());
			name.push_back(s);
			std::string vdelim = "\":";
			size_t pos = s.find(vdelim);
			if (sub.constructFromFile(s, {-1}, {""})) {
				nr.push_back( std::to_string(sub.nrow()));
				nc.push_back(std::to_string(sub.ncol()));
				nl.push_back(std::to_string(sub.nlyr()));
			}
			if (pos != std::string::npos) {
				s.erase(0, pos + vdelim.length());
				var.push_back(s);
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
	//var.resize(name.size());
	out[0] = name;
	out[1] = var;
	out[2] = desc;
	out[3] = nr;
	out[4] = nc;
	out[5] = nl;
	return out;
}





#if (!(GDAL_VERSION_MAJOR == 2 && GDAL_VERSION_MINOR < 1))
# include "gdal_utils.h" // requires >= 2.1

std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> openopts) {
// adapted from the 'sf' package by Edzer Pebesma et al
	
	std::string out = "";
	std::vector <char *> options_char = string_to_charpnt(options);
	std::vector <char *> oo_char = string_to_charpnt(openopts); // open options
	GDALInfoOptions* opt = GDALInfoOptionsNew(options_char.data(), NULL);
	GDALDatasetH ds = GDALOpenEx((const char *) filename.c_str(), GA_ReadOnly, NULL, oo_char.data(), NULL);
	if (ds == NULL) return out;
	char *val = GDALInfo(ds, opt);
	out = val;
	CPLFree(val);
	GDALInfoOptionsFree(opt);
	GDALClose(ds);
	return out;
}
#else
std::string gdalinfo(std::string filename, std::vector<std::string> options, std::vector<std::string> oo) {
	std::string out = "GDAL version >= 2.1 required for gdalinfo");
	return out;
}

#endif

