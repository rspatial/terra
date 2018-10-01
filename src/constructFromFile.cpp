using namespace std;
#include <string>
#include <vector>
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

			double xmin = atof(ini.GetValue("georeference", "xmin"));
			double xmax = atof(ini.GetValue("georeference", "xmax"));
			double ymin = atof(ini.GetValue("georeference", "ymin"));
			double ymax = atof(ini.GetValue("georeference", "ymax"));	
			SpatExtent e(xmin, xmax, ymin, ymax);

			unsigned nlyrs;
			double vna;
			string dtp, crs, smin, smax, bandorder, byteorder;
			std::vector<string> nms;
			std::vector<double> dmin, dmax; 			

			string version = ini.GetValue("version", "version", "1");
			if (version == "1") {
				unsigned nrows = atoi(ini.GetValue("georeference", "nrows"));
				unsigned ncols = atoi(ini.GetValue("georeference", "ncols"));
				SpatExtent e(xmin, xmax, ymin, ymax);
				setExtent(e, false);
				ncol = ncols;
				nrow = nrows;
				crs = ini.GetValue("georeference", "projection");
				dtp = ini.GetValue("data", "datatype");
				nlyrs = atoi(ini.GetValue("data", "nbands"));
				vna  = atof(ini.GetValue("data", "nodatavalue"));
				bandorder = ini.GetValue("data", "bandorder");
				byteorder = ini.GetValue("data", "byteorder");
				smin = ini.GetValue("data", "minvalue");
				smax = ini.GetValue("data", "maxvalue");
				dmin = str2dbl(strsplit(smin, ":"));
				dmax = str2dbl(strsplit(smax, ":"));	
				string snames = ini.GetValue("description", "layername");
				nms = { strsplit(snames, ":") }; 
				
			} else {  // version 2
			
				unsigned nrows = atoi(ini.GetValue("dimensions", "nrow"));
				unsigned ncols = atoi(ini.GetValue("dimensions", "ncol"));
				setExtent(e, false);
				ncol = ncols;
				nrow = nrows;
				nlyrs = atoi(ini.GetValue("dimensions", "nlyr"));
				string snames = ini.GetValue("dimensions", "names");
				nms = { strsplit(snames, ":|:") }; 
				crs = ini.GetValue("georeference", "crs");
				dtp = ini.GetValue("data", "datatype");
				bandorder = ini.GetValue("data", "bandorder");
				byteorder = ini.GetValue("data", "byteorder");
				smin = ini.GetValue("data", "range_min");
				smax = ini.GetValue("data", "range_max");		
				dmin = str2dbl(strsplit(smin, ":|:"));
				dmax = str2dbl(strsplit(smax, ":|:"));	
				vna = atof(ini.GetValue("data", "nodata"));
			}

			RasterSource s;
			s.datatype =  dtp ;
			s.nlyr =  nlyrs ;		
	    	s.bandorder =  bandorder ;
	    	s.endian =  byteorder ;

			s.hasValues = true; 
			s.NAflag =  vna ;
			s.hasRange = { true };
			s.range_min = dmin;
			s.range_max = dmax;		
			s.filename = setFileExt(fname, ".gri");
			s.memory = false;
			
			s.driver = "raster";
			s.crs = crs;
			s.names = nms;
			
			source = { s };
			
			//setnlyr();
			return true;
		}
   }
   return true;
}


