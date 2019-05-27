#include <vector>
#include <string>
#include "string_utils.h"
#include "SimpleIni.h"

std::vector<std::string> hdr_read(std::string filename) {
	CSimpleIniA ini(true, false, false);
	char ss[filename.length()];
	strcpy(ss, filename.c_str());
	SI_Error rc = ini.LoadFile(ss);
	std::vector<std::string> s(16);
	s[15] = "false";
	if (rc < 0) {
		return s;
	} 
	s[0] = ini.GetValue("georeference", "xmin");
	s[1] = ini.GetValue("georeference", "xmax");
	s[2] = ini.GetValue("georeference", "ymin");
	s[3] = ini.GetValue("georeference", "ymax");
	s[4] = ini.GetValue("data", "datatype");
	s[5] = ini.GetValue("data", "bandorder", "");
	s[6] = ini.GetValue("data", "byteorder", "");
	s[7] = ini.GetValue("version", "version", "1");
	if (s[7] == "1") {
		s[8] = ini.GetValue("georeference", "nrows");
		s[9] =	ini.GetValue("georeference", "ncols");
		s[10] = ini.GetValue("data", "nbands");
		s[11] = ini.GetValue("georeference", "projection");
		s[12] = ini.GetValue("data", "nodatavalue");
		s[13] = ini.GetValue("data", "minvalue");
		s[14] = ini.GetValue("data", "maxvalue");
		s[15] = ini.GetValue("description", "layername");
	} else {  // version 2
		s[8] = ini.GetValue("dimensions", "nrow");
		s[9] = ini.GetValue("dimensions", "ncol");
		s[10] = ini.GetValue("dimensions", "nlyr");
		s[11] = ini.GetValue("georeference", "crs");
		s[12] = ini.GetValue("data", "nodata");
		s[13] = ini.GetValue("data", "range_min");
		s[14] = ini.GetValue("data", "range_max");
		s[15] = ini.GetValue("dimensions", "names");
	}
	return(s);
}



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


