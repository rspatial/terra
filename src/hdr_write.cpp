#include <vector>
#include <string>
#include "SimpleIni.h"

bool hdr_write(std::vector<std::string> s) {
	CSimpleIniA ini;
	ini.SetValue("version", NULL, NULL);
	ini.SetValue("version", "version", "2");
	ini.SetValue("georeference", NULL, NULL);
	ini.SetValue("georeference", "xmin", s[0].c_str());
	ini.SetValue("georeference", "xmax", s[1].c_str());
	ini.SetValue("georeference", "ymin", s[2].c_str());
	ini.SetValue("georeference", "ymax", s[3].c_str());
	ini.SetValue("georeference", "crs", s[4].c_str());
	ini.SetValue("dimensions", "nrow", s[5].c_str());
	ini.SetValue("dimensions", "ncol", s[6].c_str());
	ini.SetValue("dimensions", "nlyr", s[7].c_str());
	ini.SetValue("dimensions", "names", s[8].c_str());
	ini.SetValue("data", NULL, NULL);
	ini.SetValue("data", "datatype", s[9].c_str()); 
	ini.SetValue("data", "nodata", s[10].c_str());
	ini.SetValue("data", "range_min", s[11].c_str());
	ini.SetValue("data", "range_max", s[12].c_str());

	SI_Error rc = ini.SaveFile(s[13].c_str());
	return rc >= 0 ;
}


