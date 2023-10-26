#include <iostream>
#include "spatRaster.h"
#include "file_utils.h"
#include "show.h"
#include "fun.h"

int main(int argc, char *argv[]) {

	std::vector<std::string> arguments = std::vector<std::string>(argv, argv + argc);

 	arguments[0] = arguments[0].substr(arguments[0].find_last_of("/\\") + 1);
 	if (arguments.size() < 2) {
		std::string msg = "usage: " +  arguments[0] + " method input output parameters";
		std::cout << msg << std::endl;
        return 1;
	}
	GDALAllRegister();
	SpatRaster out;
	std::string method = arguments[1];
	//SpatRaster input(arguments[2], {-1}, {""});
    //show(input);
    if (method == "show") out = SpatRaster(arguments[2], {-1}, {""}, {""}, {""});
    if (method == "aggregate") out = aggregate(arguments);
 
    if (out.hasError()) {
       	std::cout <<  out.msg.error << std::endl;
       	return 1;
    } else {
		show(out);
    }
	return 0;
}

