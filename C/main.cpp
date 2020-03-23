#include <iostream>
#include "spatRaster.h"
#include "gdal_frmts.h"

int main(int argc, char *argv[]) {

    SpatRaster test;
    std::cout <<  test.nrow() << std::endl;


    if(argc < 2) {
        std::string path = (std::string)argv[0];
        std::string basename = path.substr(path.find_last_of("/\\") + 1);
        std::string usage = "Usage: " + basename + " method  args ...";
        std::cout << usage << std::endl;
        return 1;
    } else {
	GDALAllRegister();
        std::string x = (std::string)argv[1];
        SpatRaster r(x);
        std::cout << r.ncol() << std::endl;
        return 0;
    }
	return 0;
}

