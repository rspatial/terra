using namespace std;
#include "spat.h"
#include "SimpleIni.h"
#include "util.h"


bool SpatRaster::constructFromFile(std::string fname) {
	
	string ext = getFileExt(fname);
	
	if (ext != ".grd") {
		
		return constructFromFileGDAL(fname);		
		
	} else {
		
		CSimpleIniA ini(true, false, false);
		char ss[fname.length()];
		strcpy(ss, fname.c_str());
		SI_Error rc = ini.LoadFile(ss);
		if (rc < 0) {
			return false;
			
		} else {

			RasterSource s;
			double xmin = atof(ini.GetValue("georeference", "xmin"));
			double xmax = atof(ini.GetValue("georeference", "xmax"));
			double ymin = atof(ini.GetValue("georeference", "ymin"));
			double ymax = atof(ini.GetValue("georeference", "ymax"));	
			SpatExtent e(xmin, xmax, ymin, ymax);
			s.extent = e;
			s.datatype = ini.GetValue("data", "datatype");
			s.bandorder = ini.GetValue("data", "bandorder");
			s.byteorder = ini.GetValue("data", "byteorder");
			
			string smin, smax, snames;
			string version = ini.GetValue("version", "version", "1");
			if (version == "1") {
				s.nrow = atoi(ini.GetValue("georeference", "nrows"));
				s.ncol = atoi(ini.GetValue("georeference", "ncols"));
				s.crs = ini.GetValue("georeference", "projection");
				s.nlyr = atoi(ini.GetValue("data", "nbands"));
				s.NAflag  = atof(ini.GetValue("data", "nodatavalue"));
				smin = ini.GetValue("data", "minvalue");
				smax = ini.GetValue("data", "maxvalue");
				snames = ini.GetValue("description", "layername");
			} else {  // version 2			
				s.nrow = atoi(ini.GetValue("dimensions", "nrow"));
				s.ncol = atoi(ini.GetValue("dimensions", "ncol"));
				s.nlyr = atoi(ini.GetValue("dimensions", "nlyr"));
				snames = ini.GetValue("dimensions", "names");
				s.crs = ini.GetValue("georeference", "crs");
				smin = ini.GetValue("data", "range_min");
				smax = ini.GetValue("data", "range_max");		
				s.NAflag = atof(ini.GetValue("data", "nodata"));
			}

			s.names = { strsplit(snames, ":|:") }; 
			s.range_min = str2dbl(strsplit(smin, ":"));
			s.range_max = str2dbl(strsplit(smax, ":"));	
			s.hasValues = true; 
			s.hasRange = { true };
			s.filename = setFileExt(fname, ".gri");
			s.memory = false;		
			s.driver = "raster";
			
			setSource(s);
			
			return true;
		}
   }
   return true;
}


