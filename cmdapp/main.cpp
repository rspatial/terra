#include <iostream>
#include <cstdlib>
#include "spatRaster.h"
#include "file_utils.h"
#include "show.h"
#include "fun.h"
#include "gdal_init.h"

int main(int argc, char *argv[]) {

	std::vector<std::string> arguments = std::vector<std::string>(argv, argv + argc);

 	arguments[0] = arguments[0].substr(arguments[0].find_last_of("/\\") + 1);
 	if (arguments.size() < 2) {
		std::string msg = "usage: " +  arguments[0] + " method input output parameters";
		std::cout << msg << std::endl;
		SpatRaster out;
		SpatExtent e(0,10,0,10);
		SpatOptions opt;
		out = out.crop(e, "near", false, opt);
		show(out);
		
        return 1;
	}

	const char *gdal_data = std::getenv("GDAL_DATA");
#ifdef TERRA_GDAL_DATA
	if (!gdal_data || !gdal_data[0])
		gdal_data = TERRA_GDAL_DATA;
#endif
	const char *proj_path = std::getenv("PROJ_LIB");
	if (!proj_path || !proj_path[0])
		proj_path = std::getenv("PROJ_DATA");
#ifdef TERRA_PROJ_DATA
	if (!proj_path || !proj_path[0])
		proj_path = TERRA_PROJ_DATA;
#endif
	const char *gdal_plugindir = std::getenv("GDAL_DRIVER_PATH");
#ifdef TERRA_GDAL_PLUGINDIR
	if (!gdal_plugindir || !gdal_plugindir[0])
		gdal_plugindir = TERRA_GDAL_PLUGINDIR;
#endif
	terra_gdal_app_init(gdal_data, proj_path, gdal_plugindir);

	SpatRaster out;
	std::string method = arguments[1];
	//SpatRaster input(arguments[2], {-1}, {""});
    //show(input);
	
	
    if (method == "show") out = SpatRaster(arguments[2], {-1}, {""}, {}, {}, false, true, {});
    if (method == "aggregate") out = aggregate(arguments);
    if (method == "") out = SpatRaster();
 
    if (out.hasError()) {
       	std::cout <<  out.msg.error << std::endl;
       	return 1;
    } else {
		show(out);
    }
	return 0;
}
