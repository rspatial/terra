
#include <vector>
#include <string>
#include "string_utils.h"
#include "SimpleIni.h"

std::vector<std::string> readIni(std::string filename) {
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

