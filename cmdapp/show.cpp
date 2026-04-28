#include <iostream>
#include <sstream>
#include <iterator>
#include <vector>
#include <cmath>
#include "spatRaster.h"

std::string join(const std::vector<std::string>& vec, const char* delim) {
    std::stringstream s;
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::string>(s, delim));
    return s.str();
}


void show(SpatRaster &r) {
    if (r.msg.has_error) {
        std::cout << "Error: " << r.msg.error << std::endl;
    } else {
        std::cout << r.show();
    }
}


void show(SpatVector &v) {
    if (v.msg.has_error) {
        std::cout << "Error: " << v.msg.error << std::endl;
    } else {
        std::cout << v.show();
    }
}


void showValues(SpatRaster d) {
	SpatOptions opt;
    std::vector<double> v = d.getValues(0, opt);
    for (size_t k=0; k<d.nlyr(); k++) {
        for (size_t i=0; i<d.nrow(); i++) {
            for (size_t j=0; j<d.ncol(); j++) {
                size_t cell = k*d.ncell() + i*d.ncol()+j;
                if (std::isnan(v[cell])) {
                    std::cout << "NA" << " ";
                } else {
                    std::cout << v[cell] << " ";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
      std::cout << "\n";
    }
      std::cout << "\n";
}


void showValues(std::vector<double> v, unsigned nr, unsigned nc) {
    for (size_t i=0; i<nr; i++) {
        for (size_t j=0; j<nc; j++) {
            size_t cell = i*nc+j;
            if (std::isnan(v[cell])) {
                std::cout << "NA" << " ";
            } else {
                std::cout << v[cell] << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void showValues(std::vector<std::vector<double>> v) {
    for (size_t i=0; i<v.size(); i++) {
        for (size_t j=0; j<v[i].size(); j++) {
            if (std::isnan(v[i][j])) {
                std::cout << "NA" << " ";
            } else {
                std::cout << v[i][j] << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}



void showValues(std::vector<std::vector<std::vector<double>>> v) {
    for (size_t i=0; i<v.size(); i++) {
        for (size_t j=0; j<v[i].size(); j++) {
            for (size_t k=0; k<v[i][j].size(); k++) {
                if (std::isnan(v[i][j][k])) {
                    std::cout << "NA" << " ";
                } else {
                    std::cout << v[i][j][k] << " ";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

