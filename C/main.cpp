#include <vector>
#include <fstream>
#include "spatRaster.h"


int main(int argc, char *argv[]) {

    if(argc < 2) {
        std::string path = (std::string)argv[0];
        std::string basename = path.substr(path.find_last_of("/\\") + 1);
        std::string usage = "Usage: " + basename + "method args";
        std::cout << usage << std::endl;
        return 1;
    } else {
        return 0;
    }
	return 0;
}

